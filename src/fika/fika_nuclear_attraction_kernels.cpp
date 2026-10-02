//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT

#include "fika_nuclear_attraction_kernels.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <tuple>
#include <variant>

#include "fika_gaussian_normalization.hpp"
#include "fika_overlap_screening.hpp"
#include "fika_boys.hpp"
#include "fika_gaunt.hpp"

namespace fika::detail {

namespace {

auto squared_distance(const Point3D<double>& a, const Point3D<double>& b) -> double {
  const double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;
  return x * x + y * y + z * z;
}

/// Dense contraction matrix of a shell, K x N row-major (transpose = false) or N x K (true).
auto coefficient_matrix(const BasisShell& shell, bool transpose) -> std::vector<double> {
  return std::visit(
      [transpose](const auto& s) {
        const std::size_t primitives = s.primitive_count();
        const std::size_t contractions = s.contraction_count();
        std::vector<double> matrix(contractions * primitives, 0.0);
        for (std::size_t k = 0; k < contractions; ++k) {
          for_each_coefficient(s, k, [&](std::size_t i, double coefficient) {
            matrix[transpose ? k * primitives + i : i * contractions + k] = coefficient;
          });
        }
        return matrix;
      },
      shell);
}

auto coupling_count(int l, int l_prime) -> std::size_t {
  return static_cast<std::size_t>(std::min(l, l_prime) + 1);
}

}  // namespace

auto far_field_threshold(int n) -> double {
  // Q(s, T) <= T^(s-1) exp(-T) / Gamma(s) / (1 - (s-1)/T) for T > s - 1, s = n + 1/2.
  const double s = n + 0.5;
  const double target = std::ldexp(1.0, -53);
  const auto bound = [s](double t) {
    return std::exp((s - 1.0) * std::log(t) - t - std::lgamma(s)) / (1.0 - (s - 1.0) / t);
  };
  double low = std::max(1.0, s);  // bound decreases beyond its largest pole-free point
  double high = low + 1.0;
  while (bound(high) > target) {
    low = high;
    high *= 2.0;
  }
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double middle = 0.5 * (low + high);
    (bound(middle) > target ? low : high) = middle;
  }
  return high;
}

SourceFields::SourceFields(const Molecule<double>& molecule, const MolecularBasis& basis,
                           const PointSources& sources, const FarFieldExpansion* far_field)
    : sources_(sources), atoms_(molecule.coordinates().begin(), molecule.coordinates().end()) {
  assert(sources.charges.size() == sources.charge_coordinates.size());
  assert(sources.dipoles.size() == sources.dipole_coordinates.size());
  assert(sources.quadrupoles.size() == sources.quadrupole_coordinates.size());
  thetas_.reserve(sources.quadrupoles.size());
  for (const Quadrupole& quadrupole : sources.quadrupoles) {
    thetas_.push_back(quadrupole_components(quadrupole));
  }
  const auto charges = sources.charges;
  const auto coordinates = sources.charge_coordinates;
  const auto dipoles = sources.dipoles;
  const auto dipole_coordinates = sources.dipole_coordinates;
  const auto quadrupoles = sources.quadrupoles;
  const auto quadrupole_coordinates = sources.quadrupole_coordinates;
  const std::size_t atom_count = atoms_.size();
  std::vector<int> orders(atom_count, 0);
  std::vector<double> far_radius_squared(atom_count, 0.0);
  std::vector<double> dipole_far_radius_squared(atom_count, 0.0);
  std::vector<double> quadrupole_far_radius_squared(atom_count, 0.0);
  omega_offsets_.assign(atom_count + 1, 0);
  for (std::size_t atom = 0; atom < atom_count; ++atom) {
    const auto shells = basis.atom_basis(atom).shells();
    int l_max = 0;
    double smallest = std::numeric_limits<double>::infinity();
    for (const BasisShell& shell : shells) {
      l_max = std::max(l_max, angular_momentum(shell));
      smallest = std::min(smallest, min_exponent(shell));
    }
    orders[atom] = shells.empty() ? 0 : 2 * l_max;
    if (!shells.empty()) {
      far_radius_squared[atom] = far_field_threshold(orders[atom]) / (2.0 * smallest);
      dipole_far_radius_squared[atom] = far_field_threshold(orders[atom] + 1) / (2.0 * smallest);
      quadrupole_far_radius_squared[atom] =
          far_field_threshold(orders[atom] + 2) / (2.0 * smallest);
    }
    const auto order = static_cast<std::size_t>(orders[atom]);
    omega_offsets_[atom + 1] = omega_offsets_[atom] + (order + 1) * (order + 1);
  }
  omega_.assign(omega_offsets_.back(), 0.0);

  // Near lists (counted, then filled) and Omega of each atom, in parallel over atoms. With a
  // far-field expansion, only the near charges of the atom's leaf (which include every charge
  // near the atom) are visited, in index order, and the other charges enter Omega through the
  // leaf's local expansion.
  std::vector<std::vector<std::uint32_t>> near_lists(atom_count);
  std::vector<std::vector<std::uint32_t>> near_dipole_lists(atom_count);
  std::vector<std::vector<std::uint32_t>> near_quadrupole_lists(atom_count);
#pragma omp parallel
  {
    std::vector<Point3D<double>> far_separations;
    std::vector<double> far_charges;
    std::vector<Dipole> far_dipoles;
    std::vector<std::array<double, 5>> far_thetas;
    SolidHarmonics harmonics;
    std::vector<std::size_t> counts;
    std::vector<std::uint32_t> candidates;
    ExpansionWorkspace expansion;
    std::vector<double> expanded;
#pragma omp for schedule(dynamic, 4)
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
      const int order = orders[atom];
      double* omega = omega_.data() + omega_offsets_[atom];
      far_separations.clear();
      far_charges.clear();
      candidates.clear();
      if (far_field != nullptr) {
        assert(far_field->contains(atoms_[atom]));
        const std::size_t leaf = far_field->leaf(atoms_[atom]);
        const auto near = far_field->near(leaf);
        candidates.assign(near.begin(), near.end());
        std::sort(candidates.begin(), candidates.end());
        expanded.resize(static_cast<std::size_t>((order + 1) * (order + 1)));
        local_field_tensor(far_field->local(leaf), far_field->centre(leaf), far_field->order(),
                           atoms_[atom], order, expanded, expansion);
        std::copy(expanded.begin(), expanded.end(),
                  omega_.begin() + static_cast<std::ptrdiff_t>(omega_offsets_[atom]));
      }
      const std::size_t candidate_count =
          far_field != nullptr ? candidates.size() : coordinates.size();
      for (std::size_t i = 0; i < candidate_count; ++i) {
        const std::size_t c = far_field != nullptr ? candidates[i] : i;
        const Point3D<double>& position = coordinates[c];
        const double r2 = squared_distance(position, atoms_[atom]);
        if (r2 < far_radius_squared[atom]) {
          near_lists[atom].push_back(static_cast<std::uint32_t>(c));
        } else if (charges[c] != 0.0) {
          far_separations.push_back({position.x - atoms_[atom].x, position.y - atoms_[atom].y,
                                     position.z - atoms_[atom].z});
          far_charges.push_back(charges[c]);
        }
      }
      if (!far_separations.empty()) {
        counts.assign(static_cast<std::size_t>(order) + 1, far_separations.size());
        harmonics.compute(far_separations, counts);
        const auto r = harmonics.distances();
        for (std::size_t c = 0; c < far_separations.size(); ++c) {
          const double inverse = 1.0 / r[c];
          double scale = far_charges[c] * inverse;  // q / |W|^(2L + 1)
          for (int big_l = 0; big_l <= order; ++big_l) {
            for (int m = -big_l; m <= big_l; ++m) {
              omega[big_l * big_l + m + big_l] += scale * harmonics.values(big_l, m)[c];
            }
            scale *= inverse * inverse;
          }
        }
      }

      // Dipoles: near ones listed, far ones add -(2L + 1) A+_LM(W, mu) / W^(2L+3); with a
      // far-field expansion only the near dipoles of the atom's leaf (the rest is in the local
      // expansion).
      far_separations.clear();
      far_dipoles.clear();
      candidates.clear();
      if (far_field != nullptr) {
        const auto near = far_field->near_dipoles(far_field->leaf(atoms_[atom]));
        candidates.assign(near.begin(), near.end());
        std::sort(candidates.begin(), candidates.end());
      }
      const std::size_t dipole_candidates =
          far_field != nullptr ? candidates.size() : dipoles.size();
      for (std::size_t i = 0; i < dipole_candidates; ++i) {
        const std::size_t d = far_field != nullptr ? candidates[i] : i;
        const Point3D<double>& position = dipole_coordinates[d];
        const double r2 = squared_distance(position, atoms_[atom]);
        if (r2 < dipole_far_radius_squared[atom]) {
          near_dipole_lists[atom].push_back(static_cast<std::uint32_t>(d));
        } else if (dipoles[d] != Dipole{}) {
          far_separations.push_back({position.x - atoms_[atom].x, position.y - atoms_[atom].y,
                                     position.z - atoms_[atom].z});
          far_dipoles.push_back(dipoles[d]);
        }
      }
      if (!far_separations.empty()) {
        counts.assign(static_cast<std::size_t>(order) + 2, far_separations.size());
        harmonics.compute(far_separations, counts);
        const auto r = harmonics.distances();
        for (std::size_t d = 0; d < far_separations.size(); ++d) {
          const double inverse = 1.0 / r[d];
          double scale = -inverse * inverse * inverse;  // -1 / |W|^(2L + 3)
          for (int big_l = 0; big_l <= order; ++big_l) {
            const double factor = (2 * big_l + 1) * scale;
            for (const DipoleCouplingTerm& term : dipole_couplings(big_l).plus) {
              omega[big_l * big_l + term.big_m + big_l] +=
                  factor * term.value * dipole_component(far_dipoles[d], term.m) *
                  harmonics.values(big_l + 1, term.harmonic_m)[d];
            }
            scale *= inverse * inverse;
          }
        }
      }

      // Quadrupoles: near ones listed, far ones add
      // (2L + 1)(2L + 3) / 2 B^(L+2)_LM(W, theta) / W^(2L+5); with a far-field expansion only
      // the near quadrupoles of the atom's leaf.
      far_separations.clear();
      far_thetas.clear();
      candidates.clear();
      if (far_field != nullptr) {
        const auto near = far_field->near_quadrupoles(far_field->leaf(atoms_[atom]));
        candidates.assign(near.begin(), near.end());
        std::sort(candidates.begin(), candidates.end());
      }
      const std::size_t quadrupole_candidates =
          far_field != nullptr ? candidates.size() : quadrupoles.size();
      for (std::size_t i = 0; i < quadrupole_candidates; ++i) {
        const std::size_t e = far_field != nullptr ? candidates[i] : i;
        const Point3D<double>& position = quadrupole_coordinates[e];
        const double r2 = squared_distance(position, atoms_[atom]);
        if (r2 < quadrupole_far_radius_squared[atom]) {
          near_quadrupole_lists[atom].push_back(static_cast<std::uint32_t>(e));
        } else {
          far_separations.push_back({position.x - atoms_[atom].x, position.y - atoms_[atom].y,
                                     position.z - atoms_[atom].z});
          far_thetas.push_back(thetas_[e]);
        }
      }
      if (!far_separations.empty()) {
        counts.assign(static_cast<std::size_t>(order) + 3, far_separations.size());
        harmonics.compute(far_separations, counts);
        const auto r = harmonics.distances();
        for (std::size_t e = 0; e < far_separations.size(); ++e) {
          const double inverse = 1.0 / r[e];
          double scale = std::pow(inverse, 5);  // 1 / |W|^(2L + 5)
          for (int big_l = 0; big_l <= order; ++big_l) {
            const double factor = 0.5 * (2 * big_l + 1) * (2 * big_l + 3) * scale;
            for (const QuadrupoleCouplingTerm& term : quadrupole_couplings(big_l).raised) {
              omega[big_l * big_l + term.big_m + big_l] +=
                  factor * term.value * far_thetas[e][static_cast<std::size_t>(term.mu + 2)] *
                  harmonics.values(big_l + 2, term.harmonic_m)[e];
            }
            scale *= inverse * inverse;
          }
        }
      }
    }
  }
  const auto flatten = [&](const std::vector<std::vector<std::uint32_t>>& lists,
                           std::vector<std::size_t>& offsets, std::vector<std::uint32_t>& flat) {
    offsets.assign(atom_count + 1, 0);
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
      offsets[atom + 1] = offsets[atom] + lists[atom].size();
    }
    flat.reserve(offsets.back());
    for (const auto& list : lists) {
      flat.insert(flat.end(), list.begin(), list.end());
    }
  };
  flatten(near_lists, near_offsets_, near_);
  flatten(near_dipole_lists, near_dipole_offsets_, near_dipoles_);
  flatten(near_quadrupole_lists, near_quadrupole_offsets_, near_quadrupoles_);
}

auto make_same_atom_shell_pair(const BasisShell& bra, const BasisShell& ket) -> SameAtomShellPair {
  SameAtomShellPair pair;
  pair.l = angular_momentum(bra);
  pair.l_prime = angular_momentum(ket);
  const auto alphas = exponents(bra);
  const auto betas = exponents(ket);
  pair.bra_primitives = alphas.size();
  pair.ket_primitives = betas.size();
  pair.bra_contractions = contraction_count(bra);
  pair.ket_contractions = contraction_count(ket);
  for (const double alpha : alphas) {
    for (const double beta : betas) {
      pair.p.push_back(alpha + beta);
    }
  }
  const int n = pair.l + pair.l_prime;
  for (int big_l = std::abs(pair.l - pair.l_prime); big_l <= n; big_l += 2) {
    const int k = (n - big_l) / 2;
    const double gamma = std::tgamma(k + big_l + 1.5) / (2 * big_l + 1);
    for (const double p : pair.p) {
      pair.far_moments.push_back(gamma * std::pow(p, -(k + big_l + 1.5)));
    }
  }
  pair.bra_coefficients = coefficient_matrix(bra, true);
  pair.ket_coefficients = coefficient_matrix(ket, false);
  return pair;
}

void same_atom_potential_values(const SameAtomShellPair& pair, const SourceFields& fields,
                                std::span<const AtomPair> pairs, std::size_t n,
                                KernelWorkspace& workspace, std::span<double> values) {
  const int l = pair.l;
  const int l_prime = pair.l_prime;
  const int order = l + l_prime;
  const int lowest = std::abs(l - l_prime);
  const std::size_t couplings = coupling_count(l, l_prime);
  const std::size_t ka = pair.bra_primitives, kb = pair.ket_primitives;
  const std::size_t na = pair.bra_contractions, nb = pair.ket_contractions;
  const std::size_t primitive_pairs = ka * kb;
  const std::size_t contracted_pairs = na * nb;
  const auto ket_components = static_cast<std::size_t>(2 * l_prime + 1);
  const std::size_t components = static_cast<std::size_t>(2 * l + 1) * ket_components;
  assert(values.size() >= contracted_pairs * components * n);

  // Field rows (L, M): coupling c has L = lowest + 2c and rows first_row[c] .. first_row[c] + 2L.
  std::array<std::size_t, max_angular_momentum + 2> first_row{};
  for (std::size_t c = 0; c < couplings; ++c) {
    first_row[c + 1] =
        first_row[c] + static_cast<std::size_t>(2 * (lowest + 2 * static_cast<int>(c)) + 1);
  }
  const std::size_t field_rows = first_row[couplings];
  workspace.field.resize(field_rows * primitive_pairs);
  workspace.contracted.resize(field_rows * contracted_pairs);
  workspace.boys.resize(static_cast<std::size_t>(order) + 3);  // quadrupoles: two orders more
  std::fill_n(values.begin(), contracted_pairs * components * n, 0.0);
  const auto charges = fields.charges();
  const auto coordinates = fields.charge_coordinates();
  const auto dipoles = fields.dipoles();
  const auto dipole_coordinates = fields.dipole_coordinates();
  const auto quadrupoles = fields.quadrupoles();
  const auto quadrupole_coordinates = fields.quadrupole_coordinates();
  std::vector<std::size_t> counts;

  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t atom = pairs[i].bra;
    assert(pairs[i].ket == atom);
    double* field = workspace.field.data();
    std::fill_n(field, field_rows * primitive_pairs, 0.0);

    // Near charges: Y^(L)_M,ab += q S_LM(W) G_kL(p_ab, W^2), with T = p W^2 and
    // G_kL = k! p^-(k+1) sum_{j=0..k} T^j / j! F_(L+j)(T) (all terms positive).
    const auto near = fields.near(atom);
    if (!near.empty()) {
      const Point3D<double>& centre = fields.atom_position(atom);
      workspace.charge_separations.clear();
      for (const std::uint32_t c : near) {
        workspace.charge_separations.push_back({coordinates[c].x - centre.x,
                                                coordinates[c].y - centre.y,
                                                coordinates[c].z - centre.z});
      }
      counts.assign(static_cast<std::size_t>(order) + 1, near.size());
      workspace.charge_harmonics.compute(workspace.charge_separations, counts);
      const auto w2 = workspace.charge_harmonics.distances_squared();
      const double* f = workspace.boys.data();
      for (std::size_t c = 0; c < near.size(); ++c) {
        const double q = charges[near[c]];
        if (q == 0.0) {
          continue;
        }
        for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
          const double p = pair.p[ab];
          const double t = p * w2[c];
          boys(order, t, workspace.boys);
          const double inverse_p = 1.0 / p;
          for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
            const int big_l = lowest + 2 * static_cast<int>(coupling);
            const int k = (order - big_l) / 2;
            double term = 1.0;
            double sum = f[big_l];
            double scale = q * inverse_p;  // q k! p^-(k+1)
            for (int j = 1; j <= k; ++j) {
              term *= t / j;
              sum += term * f[big_l + j];
              scale *= j * inverse_p;
            }
            const double g = scale * sum;
            double* rows = field + first_row[coupling] * primitive_pairs + ab;
            for (int m = -big_l; m <= big_l; ++m) {
              rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                  g * workspace.charge_harmonics.values(big_l, m)[c];
            }
          }
        }
      }
    }

    // Near dipoles (dipole notes Eq. 17): Y^(L)_M,ab += H+_kL(p_ab, W^2) A+_LM(W, mu) +
    // H-_kL(p_ab, W^2) A-_LM(W, mu), with Boys orders up to n + 1.
    const auto near_dipoles = fields.near_dipoles(atom);
    if (!near_dipoles.empty()) {
      const Point3D<double>& centre = fields.atom_position(atom);
      workspace.charge_separations.clear();
      for (const std::uint32_t d : near_dipoles) {
        workspace.charge_separations.push_back({dipole_coordinates[d].x - centre.x,
                                                dipole_coordinates[d].y - centre.y,
                                                dipole_coordinates[d].z - centre.z});
      }
      counts.assign(static_cast<std::size_t>(order) + 2, near_dipoles.size());
      workspace.charge_harmonics.compute(workspace.charge_separations, counts);
      const auto w2 = workspace.charge_harmonics.distances_squared();
      // A+ and A- of each field row (L, M) for the current dipole.
      workspace.charge_weights.resize(2 * field_rows);
      double* a_plus = workspace.charge_weights.data();
      double* a_minus = a_plus + field_rows;
      for (std::size_t d = 0; d < near_dipoles.size(); ++d) {
        const Dipole& mu = dipoles[near_dipoles[d]];
        if (mu == Dipole{}) {
          continue;
        }
        std::fill_n(a_plus, 2 * field_rows, 0.0);
        for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
          const int big_l = lowest + 2 * static_cast<int>(coupling);
          const auto& terms = dipole_couplings(big_l);
          for (const DipoleCouplingTerm& term : terms.plus) {
            a_plus[first_row[coupling] + static_cast<std::size_t>(term.big_m + big_l)] +=
                term.value * dipole_component(mu, term.m) *
                workspace.charge_harmonics.values(big_l + 1, term.harmonic_m)[d];
          }
          for (const DipoleCouplingTerm& term : terms.minus) {
            a_minus[first_row[coupling] + static_cast<std::size_t>(term.big_m + big_l)] +=
                term.value * dipole_component(mu, term.m) *
                workspace.charge_harmonics.values(big_l - 1, term.harmonic_m)[d];
          }
        }
        for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
          const double p = pair.p[ab];
          boys(order + 1, p * w2[d], workspace.boys);
          for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
            const int big_l = lowest + 2 * static_cast<int>(coupling);
            const auto kernels =
                dipole_kernels((order - big_l) / 2, big_l, p, w2[d], workspace.boys);
            double* rows = field + first_row[coupling] * primitive_pairs + ab;
            for (std::size_t m = 0; m < static_cast<std::size_t>(2 * big_l + 1); ++m) {
              const std::size_t row = first_row[coupling] + m;
              rows[m * primitive_pairs] +=
                  kernels.plus * a_plus[row] + kernels.minus * a_minus[row];
            }
          }
        }
      }
    }

    // Near quadrupoles (quadrupole notes Eq. 17): Y^(L)_M,ab += 2 (K0 B^(L+2) + K1 B^(L) +
    // K2 B^(L-2))_M at (p_ab, W^2), with Boys orders up to n + 2.
    const auto near_quadrupoles = fields.near_quadrupoles(atom);
    if (!near_quadrupoles.empty()) {
      const Point3D<double>& centre = fields.atom_position(atom);
      workspace.charge_separations.clear();
      for (const std::uint32_t e : near_quadrupoles) {
        workspace.charge_separations.push_back({quadrupole_coordinates[e].x - centre.x,
                                                quadrupole_coordinates[e].y - centre.y,
                                                quadrupole_coordinates[e].z - centre.z});
      }
      counts.assign(static_cast<std::size_t>(order) + 3, near_quadrupoles.size());
      workspace.charge_harmonics.compute(workspace.charge_separations, counts);
      const auto w2 = workspace.charge_harmonics.distances_squared();
      // B^(L+2), B^(L) and B^(L-2) of each field row (L, M) for the current quadrupole.
      workspace.charge_weights.resize(3 * field_rows);
      double* b_raised = workspace.charge_weights.data();
      double* b_same = b_raised + field_rows;
      double* b_lowered = b_same + field_rows;
      for (std::size_t e = 0; e < near_quadrupoles.size(); ++e) {
        const auto theta = quadrupole_components(quadrupoles[near_quadrupoles[e]]);
        std::fill_n(b_raised, 3 * field_rows, 0.0);
        for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
          const int big_l = lowest + 2 * static_cast<int>(coupling);
          const auto& terms = quadrupole_couplings(big_l);
          const auto add = [&](const std::vector<QuadrupoleCouplingTerm>& list, int degree,
                               double* out) {
            for (const QuadrupoleCouplingTerm& term : list) {
              out[first_row[coupling] + static_cast<std::size_t>(term.big_m + big_l)] +=
                  term.value * theta[static_cast<std::size_t>(term.mu + 2)] *
                  workspace.charge_harmonics.values(degree, term.harmonic_m)[e];
            }
          };
          add(terms.raised, big_l + 2, b_raised);
          add(terms.same, big_l, b_same);
          add(terms.lowered, big_l - 2, b_lowered);
        }
        for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
          const double p = pair.p[ab];
          boys(order + 2, p * w2[e], workspace.boys);
          for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
            const int big_l = lowest + 2 * static_cast<int>(coupling);
            const auto kernels =
                quadrupole_kernels((order - big_l) / 2, big_l, p, w2[e], workspace.boys);
            double* rows = field + first_row[coupling] * primitive_pairs + ab;
            for (std::size_t m = 0; m < static_cast<std::size_t>(2 * big_l + 1); ++m) {
              const std::size_t row = first_row[coupling] + m;
              rows[m * primitive_pairs] +=
                  2.0 * (kernels.k0 * b_raised[row] + kernels.k1 * b_same[row] +
                         kernels.k2 * b_lowered[row]);
            }
          }
        }
      }
    }

    // Far sources: Y^(L)_M,ab += Gamma(k + L + 3/2) / (2L + 1) p_ab^-(k + L + 3/2) Omega_LM.
    const auto omega = fields.omega(atom);
    if (!omega.empty()) {
      for (std::size_t coupling = 0; coupling < couplings; ++coupling) {
        const int big_l = lowest + 2 * static_cast<int>(coupling);
        const double* moments = pair.far_moments.data() + coupling * primitive_pairs;
        for (int m = -big_l; m <= big_l; ++m) {
          const double value = omega[static_cast<std::size_t>(big_l * big_l + m + big_l)];
          if (value == 0.0) {
            continue;
          }
          double* row =
              field + (first_row[coupling] + static_cast<std::size_t>(m + big_l)) * primitive_pairs;
          for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
            row[ab] += moments[ab] * value;
          }
        }
      }
    }

    // Contraction c^T Y d per row, then the Gaunt assembly 2 pi sum_LM C^{LM}_{lm,l'm'} Y^(L)_M.
    double* contracted = workspace.contracted.data();
    for (std::size_t row = 0; row < field_rows; ++row) {
      const double* y = field + row * primitive_pairs;
      for (std::size_t bra = 0; bra < na; ++bra) {
        for (std::size_t ket = 0; ket < nb; ++ket) {
          double sum = 0.0;
          for (std::size_t a = 0; a < ka; ++a) {
            const double c = pair.bra_coefficients[bra * ka + a];
            if (c == 0.0) {
              continue;
            }
            double inner = 0.0;
            for (std::size_t b = 0; b < kb; ++b) {
              inner += y[a * kb + b] * pair.ket_coefficients[b * nb + ket];
            }
            sum += c * inner;
          }
          contracted[row * contracted_pairs + bra * nb + ket] = sum;
        }
      }
    }
    for (const GauntEntry& entry : gaunt_coefficients(l, l_prime)) {
      const std::size_t component = static_cast<std::size_t>(entry.m + l) * ket_components +
                                    static_cast<std::size_t>(entry.m_prime + l_prime);
      const std::size_t row = first_row[static_cast<std::size_t>((entry.big_l - lowest) / 2)] +
                              static_cast<std::size_t>(entry.big_m + entry.big_l);
      const double scale = 2.0 * std::numbers::pi * entry.value;
      for (std::size_t ij = 0; ij < contracted_pairs; ++ij) {
        values[(ij * components + component) * n + i] +=
            scale * contracted[row * contracted_pairs + ij];
      }
    }
  }
}

namespace {

/// sum_i x[i] y[i] with four partial sums (independent chains that vectorize without
/// reassociating a single sum).
auto dot(const double* x, const double* y, std::size_t count) -> double {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  std::size_t i = 0;
  for (; i + 4 <= count; i += 4) {
    s0 += x[i] * y[i];
    s1 += x[i + 1] * y[i + 1];
    s2 += x[i + 2] * y[i + 2];
    s3 += x[i + 3] * y[i + 3];
  }
  for (; i < count; ++i) {
    s0 += x[i] * y[i];
  }
  return (s0 + s1) + (s2 + s3);
}

/// (l over lambda)!! = (2l - 1)!! / ((2 lambda - 1)!! (2l - 2 lambda - 1)!!), (-1)!! = 1.
auto double_factorial_binomial(int l, int lambda) -> double {
  const auto odd = [](int k) {  // (2k - 1)!!
    double value = 1.0;
    for (int i = 1; i <= k; ++i) {
      value *= 2 * i - 1;
    }
    return value;
  };
  return odd(l) / (odd(lambda) * odd(l - lambda));
}

auto build_translation_terms(int l) -> std::vector<TranslationTerm> {
  std::vector<TranslationTerm> terms;
  for (int lambda = 0; lambda <= l; ++lambda) {
    const double binomial = double_factorial_binomial(l, lambda);
    for (const GauntEntry& entry : gaunt_coefficients(lambda, l - lambda)) {
      if (entry.big_l == l) {
        terms.push_back({entry.big_m, lambda, entry.m, entry.m_prime, binomial * entry.value});
      }
    }
  }
  return terms;
}

/// Scheme II structure of a pair (l, l') (notes Eqs. 22-23), independent of exponents and
/// geometry: the distinct (Lambda, kappa) of the charge field, the combinations
/// (lambda, lambda', Lambda) of the contraction, and the kernel terms
/// Gamma = sum T^lm_{lambda mu,nu} T^l'm'_{lambda' mu',nu'} C^{Lambda M}_{lambda mu,lambda' mu'}
/// S_(l-lambda),nu(R) S_(l'-lambda'),nu'(R), summed over mu, mu'.
struct SchemeTwoTable {
  using Distinct = FieldRowBlock;  // field rows first_row .. first_row + 2 Lambda
  struct Combination {
    int lambda;
    int lambda_prime;
    int big_l;
    std::size_t distinct;
    std::size_t first_row;  // contracted rows first_row .. first_row + 2 Lambda
  };
  struct Term {
    std::size_t component;  // (m + l)(2l' + 1) + m' + l'
    std::size_t row;        // contracted row of (lambda, lambda', Lambda, M)
    int bra_l, bra_m;       // S_(l-lambda),nu(R)
    int ket_l, ket_m;       // S_(l'-lambda'),nu'(R)
    double value;
  };
  std::vector<Distinct> distinct;
  std::size_t field_rows = 0;
  std::vector<Combination> combinations;
  std::size_t contracted_rows = 0;
  std::vector<Term> terms;
};

auto build_scheme_two_table(int l, int l_prime) -> SchemeTwoTable {
  SchemeTwoTable table;
  for (int lambda = 0; lambda <= l; ++lambda) {
    for (int lambda_prime = 0; lambda_prime <= l_prime; ++lambda_prime) {
      for (int big_l = std::abs(lambda - lambda_prime); big_l <= lambda + lambda_prime;
           big_l += 2) {
        const int kappa = (lambda + lambda_prime - big_l) / 2;
        auto found = std::ranges::find_if(
            table.distinct, [&](const auto& d) { return d.big_l == big_l && d.kappa == kappa; });
        if (found == table.distinct.end()) {
          table.distinct.push_back({big_l, kappa, table.field_rows});
          table.field_rows += static_cast<std::size_t>(2 * big_l + 1);
          found = table.distinct.end() - 1;
        }
        table.combinations.push_back({lambda, lambda_prime, big_l,
                                      static_cast<std::size_t>(found - table.distinct.begin()),
                                      table.contracted_rows});
        table.contracted_rows += static_cast<std::size_t>(2 * big_l + 1);
      }
    }
  }
  const auto bra_terms = translation_terms(l);
  const auto ket_terms = translation_terms(l_prime);
  std::map<std::tuple<std::size_t, std::size_t, int, int, int, int>, double> terms;
  for (const auto& combination : table.combinations) {
    for (const GauntEntry& gaunt :
         gaunt_coefficients(combination.lambda, combination.lambda_prime)) {
      if (gaunt.big_l != combination.big_l) {
        continue;
      }
      for (const TranslationTerm& bra : bra_terms) {
        if (bra.lambda != combination.lambda || bra.mu != gaunt.m) {
          continue;
        }
        for (const TranslationTerm& ket : ket_terms) {
          if (ket.lambda != combination.lambda_prime || ket.mu != gaunt.m_prime) {
            continue;
          }
          const auto component =
              static_cast<std::size_t>((bra.m + l) * (2 * l_prime + 1) + ket.m + l_prime);
          const std::size_t row =
              combination.first_row + static_cast<std::size_t>(gaunt.big_m + combination.big_l);
          terms[{component, row, l - combination.lambda, bra.nu, l_prime - combination.lambda_prime,
                 ket.nu}] += bra.value * ket.value * gaunt.value;
        }
      }
    }
  }
  for (const auto& [key, value] : terms) {
    if (std::abs(value) > 1e-14) {
      const auto& [component, row, bra_l, bra_m, ket_l, ket_m] = key;
      table.terms.push_back({component, row, bra_l, bra_m, ket_l, ket_m, value});
    }
  }
  return table;
}

template <typename Table, typename Build>
auto lazy_table(int l, int l_prime, Build&& build) -> const Table& {
  struct Slot {
    std::once_flag built;
    Table table;
  };
  static std::array<std::array<Slot, max_angular_momentum + 1>, max_angular_momentum + 1> slots;
  Slot& slot = slots[static_cast<std::size_t>(l)][static_cast<std::size_t>(l_prime)];
  std::call_once(slot.built, [&] { slot.table = build(l, l_prime); });
  return slot.table;
}

auto scheme_two_table(int l, int l_prime) -> const SchemeTwoTable& {
  return lazy_table<SchemeTwoTable>(l, l_prime, build_scheme_two_table);
}

/// Factorized Scheme II assembly of a pair (l, l'):
///   Z^(lambda lambda')_mu,mu' = sum_{Lambda M} C^{Lambda M}_{lambda mu,lambda' mu'} Y^(lambda
///   lambda' Lambda)_M, G^(lambda)_m,mu = sum_nu T^lm_{lambda mu,nu} S_(l-lambda),nu(R),
///   H^(lambda') likewise, V = 2 pi sum_{lambda lambda'} G^(lambda) Z^(lambda lambda')
///   H^(lambda')^T.
struct FactorizedTable {
  struct GauntTerm {
    std::size_t z;    // mu (2 lambda' + 1) + mu' (indices from 0)
    std::size_t row;  // contracted row of (lambda, lambda', Lambda, M)
    double value;
  };
  struct TranslationEntry {
    std::size_t cell;  // m (2 lambda + 1) + mu (indices from 0)
    int harmonic_l;
    int harmonic_m;
    double value;
  };
  std::vector<std::vector<GauntTerm>> gaunt;        // per lambda (l' + 1) + lambda'
  std::vector<std::vector<TranslationEntry>> bra;   // per lambda
  std::vector<std::vector<TranslationEntry>> ket;   // per lambda'
  std::vector<std::vector<std::size_t>> bra_cells;  // nonzero cells of G^(lambda)
  std::vector<std::vector<std::size_t>> ket_cells;  // nonzero cells of H^(lambda')
  std::size_t operations = 0;  // estimated multiply-adds per contracted pair and atom pair
};

auto build_factorized_table(int l, int l_prime) -> FactorizedTable {
  const SchemeTwoTable& table = scheme_two_table(l, l_prime);
  FactorizedTable factorized;
  const auto lambdas = static_cast<std::size_t>(l + 1);
  const auto lambda_primes = static_cast<std::size_t>(l_prime + 1);
  factorized.gaunt.resize(lambdas * lambda_primes);
  for (const auto& combination : table.combinations) {
    auto& terms = factorized.gaunt[static_cast<std::size_t>(combination.lambda) * lambda_primes +
                                   static_cast<std::size_t>(combination.lambda_prime)];
    for (const GauntEntry& entry :
         gaunt_coefficients(combination.lambda, combination.lambda_prime)) {
      if (entry.big_l == combination.big_l) {
        terms.push_back(
            {static_cast<std::size_t>((entry.m + combination.lambda) *
                                          (2 * combination.lambda_prime + 1) +
                                      entry.m_prime + combination.lambda_prime),
             combination.first_row + static_cast<std::size_t>(entry.big_m + combination.big_l),
             entry.value});
      }
    }
  }
  const auto translation = [](int shell_l,
                              std::vector<std::vector<FactorizedTable::TranslationEntry>>& entries,
                              std::vector<std::vector<std::size_t>>& cells) {
    entries.resize(static_cast<std::size_t>(shell_l + 1));
    cells.resize(static_cast<std::size_t>(shell_l + 1));
    for (const TranslationTerm& term : translation_terms(shell_l)) {
      const auto lambda = static_cast<std::size_t>(term.lambda);
      const auto cell = static_cast<std::size_t>((term.m + shell_l) * (2 * term.lambda + 1) +
                                                 term.mu + term.lambda);
      entries[lambda].push_back({cell, shell_l - term.lambda, term.nu, term.value});
      if (std::ranges::find(cells[lambda], cell) == cells[lambda].end()) {
        cells[lambda].push_back(cell);
      }
    }
  };
  translation(l, factorized.bra, factorized.bra_cells);
  translation(l_prime, factorized.ket, factorized.ket_cells);

  // Z, then W = G Z (nonzero G cells times the 2 lambda' + 1 columns), then V = W H^T.
  for (std::size_t lambda = 0; lambda < lambdas; ++lambda) {
    for (std::size_t lambda_prime = 0; lambda_prime < lambda_primes; ++lambda_prime) {
      factorized.operations += factorized.gaunt[lambda * lambda_primes + lambda_prime].size() +
                               factorized.bra_cells[lambda].size() * (2 * lambda_prime + 1);
    }
  }
  for (std::size_t lambda_prime = 0; lambda_prime < lambda_primes; ++lambda_prime) {
    factorized.operations +=
        static_cast<std::size_t>(2 * l + 1) * factorized.ket_cells[lambda_prime].size();
  }
  for (std::size_t lambda = 0; lambda < lambdas; ++lambda) {
    factorized.operations += factorized.bra[lambda].size();
  }
  for (std::size_t lambda_prime = 0; lambda_prime < lambda_primes; ++lambda_prime) {
    factorized.operations += factorized.ket[lambda_prime].size();
  }
  return factorized;
}

auto factorized_table(int l, int l_prime) -> const FactorizedTable& {
  return lazy_table<FactorizedTable>(l, l_prime, build_factorized_table);
}

}  // namespace

auto make_source_far_field(const Molecule<double>& molecule, const MolecularBasis& basis,
                           const PointSources& sources, double threshold) -> FarFieldExpansion {
  double charge_sum = 0.0;
  for (const double q : sources.charges) {
    charge_sum += std::abs(q);
  }
  double dipole_sum = 0.0;
  for (const Dipole& mu : sources.dipoles) {
    dipole_sum += std::hypot(mu(0), mu(1), mu(2));
  }
  double quadrupole_sum = 0.0;
  for (const Quadrupole& q : sources.quadrupoles) {
    quadrupole_sum += traceless_quadrupole_norm(q);
  }
  int l_max = 0;
  double smallest = std::numeric_limits<double>::infinity();
  for (const AtomBasis& atom_basis : basis.unique_bases()) {
    for (const BasisShell& shell : atom_basis.shells()) {
      l_max = std::max(l_max, angular_momentum(shell));
      smallest = std::min(smallest, min_exponent(shell));
    }
  }
  const int rank = 2 * l_max;
  const auto bases = basis.unique_bases();
  std::vector<double> moments(static_cast<std::size_t>(rank) + 1, 0.0);
  std::vector<double> cutoffs(bases.size() * bases.size(), -1.0);
  for (const AtomBasisPair& pair : basis.unique_basis_pairs()) {
    double cutoff = -1.0;
    for (const BasisShell& bra : bases[pair.bra].shells()) {
      for (const BasisShell& ket : bases[pair.ket].shells()) {
        cutoff = std::max(cutoff, nuclear_attraction_shell_pair_bound(bra, ket, charge_sum)
                                      .cutoff_distance_squared(threshold));
        if (dipole_sum > 0.0) {
          cutoff = std::max(cutoff, dipole_potential_shell_pair_bound(bra, ket, dipole_sum)
                                        .cutoff_distance_squared(threshold));
        }
        if (quadrupole_sum > 0.0) {
          cutoff = std::max(cutoff, quadrupole_potential_shell_pair_bound(bra, ket, quadrupole_sum)
                                        .cutoff_distance_squared(threshold));
        }
        for (int big_l = 0; big_l <= rank; ++big_l) {
          auto& moment = moments[static_cast<std::size_t>(big_l)];
          moment = std::max(moment, far_field_moment_bound(bra, ket, big_l));
        }
      }
    }
    cutoffs[pair.bra * bases.size() + pair.ket] = cutoff;
    cutoffs[pair.ket * bases.size() + pair.bra] = cutoff;
  }
  FarFieldOptions options;
  options.rank = rank;
  for (const double moment : moments) {
    options.tolerances.push_back(moment > 0.0 ? threshold / ((rank + 1) * moment) : threshold);
  }
  const int boys_order = !sources.quadrupoles.empty() ? rank + 2
                         : !sources.dipoles.empty()   ? rank + 1
                                                      : rank;
  options.penetration_radius = std::sqrt(far_field_threshold(boys_order) / (2.0 * smallest));

  std::vector<AtomPair> pairs;
  const auto atoms = molecule.coordinates();
  const auto indices = basis.atom_basis_indices();
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    for (std::size_t b = 0; b <= a; ++b) {
      const double x = atoms[a].x - atoms[b].x, y = atoms[a].y - atoms[b].y,
                   z = atoms[a].z - atoms[b].z;
      if (a == b || x * x + y * y + z * z <= cutoffs[indices[a] * bases.size() + indices[b]]) {
        pairs.push_back({a, b});
      }
    }
  }
  return FarFieldExpansion(sources, atoms, pairs, options);
}

auto coulomb_kernel(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> double {
  assert(boys.size() > static_cast<std::size_t>(big_l + kappa));
  const double t = p * u2;
  double term = 1.0;  // T^i / i!
  double sum = boys[static_cast<std::size_t>(big_l)];
  double scale = 1.0 / p;  // kappa! p^-(kappa + 1)
  for (int i = 1; i <= kappa; ++i) {
    term *= t / i;
    sum += term * boys[static_cast<std::size_t>(big_l + i)];
    scale *= i / p;
  }
  return scale * sum;
}

auto dipole_kernels(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> DipoleKernels {
  assert(boys.size() > static_cast<std::size_t>(big_l + kappa + 1));
  const double top = boys[static_cast<std::size_t>(big_l + kappa + 1)];
  const double u_power = std::pow(u2, kappa);  // U^(2 kappa)
  DipoleKernels kernels;
  kernels.plus = -2.0 * u_power * top;
  if (big_l > 0) {
    kernels.minus = (2 * big_l - 1) * (coulomb_kernel(kappa, big_l, p, u2, boys) -
                                       2.0 * u_power * u2 * top / (2 * big_l + 1));
  }
  return kernels;
}

auto dipole_couplings(int big_l) -> const DipoleCouplings& {
  constexpr int max_rank = 2 * max_angular_momentum;
  assert(big_l >= 0 && big_l <= max_rank);
  struct Slot {
    std::once_flag built;
    DipoleCouplings couplings;
  };
  static std::array<Slot, max_rank + 1> slots;
  Slot& slot = slots[static_cast<std::size_t>(big_l)];
  std::call_once(slot.built, [&] {
    // A+: S_1m S_Lambda,M = sum C^{Lambda+1,M'}_{1m,Lambda M} S_(Lambda+1),M' + ... .
    for (const GauntEntry& entry : multipole_gaunt_coefficients(1, big_l)) {
      if (entry.big_l == big_l + 1) {
        slot.couplings.plus.push_back({entry.m_prime, entry.m, entry.big_m, entry.value});
      }
    }
    // A-: C^{Lambda M}_{1m,(Lambda-1)M'} from S_1m S_(Lambda-1),M'.
    if (big_l > 0) {
      for (const GauntEntry& entry : multipole_gaunt_coefficients(1, big_l - 1)) {
        if (entry.big_l == big_l) {
          slot.couplings.minus.push_back({entry.big_m, entry.m, entry.m_prime, entry.value});
        }
      }
    }
  });
  return slot.couplings;
}

auto quadrupole_components(const Quadrupole& q) -> std::array<double, 5> {
  const double two_over_sqrt3 = 2.0 / std::sqrt(3.0);
  return {two_over_sqrt3 * q(0, 1),                   // mu = -2: xy
          two_over_sqrt3 * q(1, 2),                   // mu = -1: yz
          (2.0 * q(2, 2) - q(0, 0) - q(1, 1)) / 3.0,  // mu = 0
          two_over_sqrt3 * q(0, 2),                   // mu = 1: xz
          (q(0, 0) - q(1, 1)) / std::sqrt(3.0)};      // mu = 2: xx - yy
}

auto quadrupole_kernels(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> QuadrupoleKernels {
  assert(boys.size() > static_cast<std::size_t>(big_l + kappa + 2));
  const double g = coulomb_kernel(kappa, big_l, p, u2, boys);
  const double first = boys[static_cast<std::size_t>(big_l + kappa + 1)];
  const double second = boys[static_cast<std::size_t>(big_l + kappa + 2)];
  const double u_power = std::pow(u2, kappa);  // x^kappa
  const double g1 = -u_power * first;          // G'
  // G'' = -kappa x^(kappa-1) F_(L+k+1) + p x^kappa F_(L+k+2); x G'' without dividing by x.
  const double x_g2 = -kappa * u_power * first + p * u_power * u2 * second;
  const double g2 =
      (kappa > 0 ? -kappa * std::pow(u2, kappa - 1) * first : 0.0) + p * u_power * second;
  QuadrupoleKernels kernels;
  kernels.k0 = g2;
  kernels.k1 = x_g2 + (big_l + 1.5) * g1;
  kernels.k2 = u2 * x_g2 + (2 * big_l + 1) * u2 * g1 + (big_l + 0.5) * (big_l - 0.5) * g;
  return kernels;
}

auto quadrupole_couplings(int big_l) -> const QuadrupoleCouplings& {
  constexpr int max_rank = 2 * max_angular_momentum;
  assert(big_l >= 0 && big_l <= max_rank);
  struct Slot {
    std::once_flag built;
    QuadrupoleCouplings couplings;
  };
  static std::array<Slot, max_rank + 1> slots;
  Slot& slot = slots[static_cast<std::size_t>(big_l)];
  std::call_once(slot.built, [&] {
    // S_2mu S_Lambda,M = sum C^{L'M'}_{2mu,Lambda M} r^(2k') S_L'M', L' = Lambda + 2 - 2k'.
    for (const GauntEntry& entry : multipole_gaunt_coefficients(2, big_l)) {
      const QuadrupoleCouplingTerm term{entry.m_prime, entry.m, entry.big_m, entry.value};
      if (entry.big_l == big_l + 2) {
        slot.couplings.raised.push_back(term);
      } else if (entry.big_l == big_l) {
        slot.couplings.same.push_back(term);
      } else {
        slot.couplings.lowered.push_back(term);
      }
    }
  });
  return slot.couplings;
}

auto translation_terms(int l) -> std::span<const TranslationTerm> {
  assert(l >= 0 && l <= max_angular_momentum);
  struct Slot {
    std::once_flag built;
    std::vector<TranslationTerm> terms;
  };
  static std::array<Slot, max_angular_momentum + 1> slots;
  Slot& slot = slots[static_cast<std::size_t>(l)];
  std::call_once(slot.built, [&] { slot.terms = build_translation_terms(l); });
  return slot.terms;
}

auto make_two_atom_shell_pair(const BasisShell& bra, const BasisShell& ket, const SourceSums& sums,
                              double threshold, AssemblyForm form) -> TwoAtomShellPair {
  TwoAtomShellPair pair;
  pair.l = angular_momentum(bra);
  pair.l_prime = angular_momentum(ket);
  const auto alphas = exponents(bra);
  const auto betas = exponents(ket);
  pair.bra_exponents.assign(alphas.begin(), alphas.end());
  pair.ket_exponents.assign(betas.begin(), betas.end());
  pair.bra_contractions = contraction_count(bra);
  pair.ket_contractions = contraction_count(ket);
  pair.bra_coefficients = coefficient_matrix(bra, true);
  pair.ket_coefficients = coefficient_matrix(ket, false);
  pair.far_threshold = far_field_threshold(pair.l + pair.l_prime);
  pair.dipole_far_threshold = far_field_threshold(pair.l + pair.l_prime + 1);
  pair.quadrupole_far_threshold = far_field_threshold(pair.l + pair.l_prime + 2);
  const auto largest = [](const std::vector<double>& matrix, std::size_t rows, std::size_t columns,
                          bool per_column) {
    std::vector<double> result(per_column ? columns : rows, 0.0);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < columns; ++j) {
        double& value = result[per_column ? j : i];
        value = std::max(value, std::abs(matrix[i * columns + j]));
      }
    }
    return result;
  };
  // c^T is N_A x K_A (primitives are columns); d is K_B x N_B (primitives are rows).
  pair.bra_largest = largest(pair.bra_coefficients, pair.bra_contractions, alphas.size(), true);
  pair.ket_largest = largest(pair.ket_coefficients, betas.size(), pair.ket_contractions, false);
  pair.bound_prefactor = (pair.l + pair.l_prime + 1) * 2.0 * std::numbers::pi * sums.charges;
  pair.dipole_prefactor =
      (pair.l + pair.l_prime + 1) * 2.0 * std::numbers::pi * 2.0 * dipole_boys_bound * sums.dipoles;
  pair.quadrupole_prefactor = (pair.l + pair.l_prime + 1) * 4.0 * std::numbers::pi *
                              quadrupole_boys_bound * sums.quadrupoles;
  pair.primitive_threshold = threshold / static_cast<double>(alphas.size() * betas.size());
  // Build the shared tables outside the kernels and choose the assembly form.
  const std::size_t flat_operations = scheme_two_table(pair.l, pair.l_prime).terms.size();
  const std::size_t factorized_operations = factorized_table(pair.l, pair.l_prime).operations;
  pair.factorized = form == AssemblyForm::factorized ||
                    (form == AssemblyForm::automatic && factorized_operations < flat_operations);
  return pair;
}

auto far_field_moment_bound(const BasisShell& bra, const BasisShell& ket, int rank) -> double {
  const int l = angular_momentum(bra);
  const int l_prime = angular_momentum(ket);
  if (rank > l + l_prime) {
    return 0.0;
  }
  const auto alphas = exponents(bra);
  const auto betas = exponents(ket);
  const auto largest = [](const BasisShell& shell) {
    const std::size_t primitives = exponents(shell).size();
    const auto matrix = coefficient_matrix(shell, false);  // K x N
    const std::size_t contractions = matrix.size() / primitives;
    std::vector<double> result(primitives, 0.0);
    for (std::size_t i = 0; i < primitives; ++i) {
      for (std::size_t k = 0; k < contractions; ++k) {
        result[i] = std::max(result[i], std::abs(matrix[i * contractions + k]));
      }
    }
    return result;
  };
  const auto bra_largest = largest(bra);
  const auto ket_largest = largest(ket);
  const auto binomial = [](int n, int k) {
    double value = 1.0;
    for (int i = 1; i <= k; ++i) {
      value = value * (n + 1 - i) / i;
    }
    return value;
  };
  double total = 0.0;
  for (std::size_t a = 0; a < alphas.size(); ++a) {
    for (std::size_t b = 0; b < betas.size(); ++b) {
      const double p = alphas[a] + betas[b];
      const double mu = alphas[a] * betas[b] / p;
      double sum = 0.0;
      for (int i = 0; i <= l; ++i) {
        for (int j = 0; j <= l_prime; ++j) {
          // (b/p)^(l-i) (a/p)^(l'-j) sup_R R^k exp(-mu R^2), k = l - i + l' - j.
          const int k = l - i + l_prime - j;
          const double peak = k == 0 ? 1.0 : std::pow(k / (2.0 * mu * std::numbers::e), 0.5 * k);
          const double s = 0.5 * (i + j + rank + 3);
          sum += binomial(l, i) * binomial(l_prime, j) * std::pow(betas[b] / p, l - i) *
                 std::pow(alphas[a] / p, l_prime - j) * peak * 2.0 * std::numbers::pi *
                 std::tgamma(s) * std::pow(p, -s);
        }
      }
      total += bra_largest[a] * ket_largest[b] * sum;
    }
  }
  return total;
}

auto field_row_blocks(int l, int l_prime) -> std::span<const FieldRowBlock> {
  return scheme_two_table(l, l_prime).distinct;
}

auto field_row_count(int l, int l_prime) -> std::size_t {
  return scheme_two_table(l, l_prime).field_rows;
}

void primitive_pair_field_weights(const TwoAtomShellPair& pair, const SolidHarmonics& harmonics,
                                  std::size_t a, std::size_t b, std::span<double> weights) {
  const int l = pair.l;
  const int l_prime = pair.l_prime;
  const SchemeTwoTable& table = scheme_two_table(l, l_prime);
  const std::size_t rows = table.field_rows;
  const auto components = static_cast<std::size_t>((2 * l + 1) * (2 * l_prime + 1));
  assert(weights.size() >= components * rows);
  std::fill_n(weights.begin(), components * rows, 0.0);
  const double alpha = pair.bra_exponents[a];
  const double beta = pair.ket_exponents[b];
  const double p = alpha + beta;
  const double t = beta / p;
  const double r2 = harmonics.distances_squared()[0];
  const double exponential = 2.0 * std::numbers::pi * std::exp(-alpha * beta / p * r2);
  // Per combination (lambda, lambda', Lambda): prefactor; per contracted row: its combination.
  thread_local std::vector<double> prefactors;
  thread_local std::vector<std::size_t> row_combination;
  prefactors.resize(table.combinations.size());
  row_combination.resize(table.contracted_rows);
  for (std::size_t k = 0; k < table.combinations.size(); ++k) {
    const auto& combination = table.combinations[k];
    prefactors[k] = exponential * std::pow(-t, l - combination.lambda) *
                    std::pow(1.0 - t, l_prime - combination.lambda_prime);
    for (int m = 0; m <= 2 * combination.big_l; ++m) {
      row_combination[combination.first_row + static_cast<std::size_t>(m)] = k;
    }
  }
  for (const auto& term : table.terms) {
    const std::size_t k = row_combination[term.row];
    const auto& combination = table.combinations[k];
    const std::size_t field_row =
        table.distinct[combination.distinct].first_row + (term.row - combination.first_row);
    weights[term.component * rows + field_row] += prefactors[k] * term.value *
                                                  harmonics.values(term.bra_l, term.bra_m)[0] *
                                                  harmonics.values(term.ket_l, term.ket_m)[0];
  }
}

void two_atom_potential_values(const TwoAtomShellPair& pair, const SourceFields& fields,
                               const SolidHarmonics& harmonics, std::span<const AtomPair> pairs,
                               std::size_t n, KernelWorkspace& workspace, std::span<double> values,
                               const FarFieldExpansion* far_field) {
  const int l = pair.l;
  const int l_prime = pair.l_prime;
  const int order = l + l_prime;
  const SchemeTwoTable& table = scheme_two_table(l, l_prime);
  const std::size_t ka = pair.bra_exponents.size(), kb = pair.ket_exponents.size();
  const std::size_t na = pair.bra_contractions, nb = pair.ket_contractions;
  const std::size_t primitive_pairs = ka * kb;
  const std::size_t contracted_pairs = na * nb;
  const std::size_t components = static_cast<std::size_t>((2 * l + 1) * (2 * l_prime + 1));
  const bool segmented = contracted_pairs == 1;
  assert(values.size() >= contracted_pairs * components * n);
  std::fill_n(values.begin(), contracted_pairs * components * n, 0.0);
  workspace.field.resize(table.field_rows * primitive_pairs);
  workspace.contracted.resize(table.contracted_rows * contracted_pairs * n);  // [row][IJ][i]
  workspace.boys.resize(static_cast<std::size_t>(order) + 3);  // quadrupoles: two orders more
  // Per primitive pair: far-field factors Gamma(kappa + Lambda + 3/2)/(2 Lambda + 1)
  // p^-(kappa + Lambda + 3/2) of each distinct (Lambda, kappa), and the contraction prefactors.
  workspace.powers.resize(table.distinct.size() * primitive_pairs);
  workspace.leading.resize(table.combinations.size() * primitive_pairs);
  for (std::size_t d = 0; d < table.distinct.size(); ++d) {
    const auto& distinct = table.distinct[d];
    const double gamma =
        std::tgamma(distinct.kappa + distinct.big_l + 1.5) / (2 * distinct.big_l + 1);
    for (std::size_t a = 0; a < ka; ++a) {
      for (std::size_t b = 0; b < kb; ++b) {
        const double p = pair.bra_exponents[a] + pair.ket_exponents[b];
        workspace.powers[d * primitive_pairs + a * kb + b] =
            gamma * std::pow(p, -(distinct.kappa + distinct.big_l + 1.5));
      }
    }
  }
  // Charges summed per primitive pair: all of them, or the near charges of the product centre's
  // leaf (gathered) with the rest through the leaf's local expansion.
  std::span<const double> charges = fields.charges();
  std::span<const Point3D<double>> coordinates = fields.charge_coordinates();
  std::vector<std::size_t> counts(static_cast<std::size_t>(order) + 1, coordinates.size());
  const double* f = workspace.boys.data();
  const std::size_t tensor_size = static_cast<std::size_t>((order + 1) * (order + 1));
  workspace.far_field.resize(2 * tensor_size);
  double* phi = workspace.far_field.data();
  double* expanded = phi + tensor_size;
  // Dipoles likewise: all of them, or the near dipoles of the leaf (gathered).
  std::span<const Dipole> dipoles = fields.dipoles();
  std::span<const Point3D<double>> dipole_coordinates = fields.dipole_coordinates();
  std::vector<std::size_t> dipole_counts(static_cast<std::size_t>(order) + 2, dipoles.size());
  // Quadrupoles likewise (theta_mu and positions).
  std::span<const std::array<double, 5>> thetas = fields.thetas();
  std::span<const Point3D<double>> quadrupole_coordinates = fields.quadrupole_coordinates();
  std::vector<std::size_t> quadrupole_counts(static_cast<std::size_t>(order) + 3, thetas.size());

  for (std::size_t i = 0; i < n; ++i) {
    const Point3D<double>& a_position = fields.atom_position(pairs[i].bra);
    const Point3D<double>& b_position = fields.atom_position(pairs[i].ket);
    const Point3D<double> r{a_position.x - b_position.x, a_position.y - b_position.y,
                            a_position.z - b_position.z};
    const double r2 = r.x * r.x + r.y * r.y + r.z * r.z;
    const double r_power = std::max(1.0, std::pow(std::sqrt(r2), order));  // max(1, R^L)
    double* field = workspace.field.data();
    std::fill_n(field, table.field_rows * primitive_pairs, 0.0);
    std::fill_n(workspace.leading.begin(), table.combinations.size() * primitive_pairs, 0.0);

    // Charge field at each product centre P = A - t R: Y^(Lambda kappa)_M,ab =
    // sum_C q S_Lambda,M(U) G_kappa,Lambda(p, U^2), U = C - P.
    for (std::size_t a = 0; a < ka; ++a) {
      for (std::size_t b = 0; b < kb; ++b) {
        const std::size_t ab = a * kb + b;
        const double p = pair.bra_exponents[a] + pair.ket_exponents[b];
        const double t = pair.ket_exponents[b] / p;
        const double mu_r2 = pair.bra_exponents[a] * pair.ket_exponents[b] / p * r2;
        if (pair.bra_largest[a] * pair.ket_largest[b] * pair.bound_prefactor / p * r_power *
                    std::exp(-mu_r2) +
                pair.bra_largest[a] * pair.ket_largest[b] * pair.dipole_prefactor / std::sqrt(p) *
                    r_power * std::exp(-mu_r2) +
                pair.bra_largest[a] * pair.ket_largest[b] * pair.quadrupole_prefactor * r_power *
                    std::exp(-mu_r2) <
            pair.primitive_threshold) {
          continue;  // negligible primitive pair: its field and prefactors stay zero
        }
        const Point3D<double> centre{a_position.x - t * r.x, a_position.y - t * r.y,
                                     a_position.z - t * r.z};
        if (far_field != nullptr) {
          assert(far_field->contains(centre));
          const std::size_t leaf = far_field->leaf(centre);
          local_field_tensor(far_field->local(leaf), far_field->centre(leaf), far_field->order(),
                             centre, order, std::span(expanded, tensor_size), workspace.expansion);
          const auto near = far_field->near(leaf);
          workspace.gathered_charges.resize(near.size());
          workspace.gathered_coordinates.resize(near.size());
          for (std::size_t g = 0; g < near.size(); ++g) {
            workspace.gathered_charges[g] = fields.charges()[near[g]];
            workspace.gathered_coordinates[g] = fields.charge_coordinates()[near[g]];
          }
          charges = workspace.gathered_charges;
          coordinates = workspace.gathered_coordinates;
          std::fill(counts.begin(), counts.end(), coordinates.size());
          const auto near_dipoles = far_field->near_dipoles(leaf);
          workspace.gathered_dipoles.resize(near_dipoles.size());
          workspace.gathered_dipole_coordinates.resize(near_dipoles.size());
          for (std::size_t g = 0; g < near_dipoles.size(); ++g) {
            workspace.gathered_dipoles[g] = fields.dipoles()[near_dipoles[g]];
            workspace.gathered_dipole_coordinates[g] = fields.dipole_coordinates()[near_dipoles[g]];
          }
          dipoles = workspace.gathered_dipoles;
          dipole_coordinates = workspace.gathered_dipole_coordinates;
          std::fill(dipole_counts.begin(), dipole_counts.end(), dipoles.size());
          const auto near_quadrupoles = far_field->near_quadrupoles(leaf);
          workspace.thetas.resize(near_quadrupoles.size());
          workspace.gathered_quadrupole_coordinates.resize(near_quadrupoles.size());
          for (std::size_t g = 0; g < near_quadrupoles.size(); ++g) {
            workspace.thetas[g] = fields.thetas()[near_quadrupoles[g]];
            workspace.gathered_quadrupole_coordinates[g] =
                fields.quadrupole_coordinates()[near_quadrupoles[g]];
          }
          thetas = workspace.thetas;
          quadrupole_coordinates = workspace.gathered_quadrupole_coordinates;
          std::fill(quadrupole_counts.begin(), quadrupole_counts.end(),
                    quadrupole_coordinates.size());
        }
        const std::size_t dipole_count = dipoles.size();
        const std::size_t quadrupole_count = quadrupole_coordinates.size();
        const std::size_t charge_count = coordinates.size();
        workspace.charge_separations.resize(charge_count);
        workspace.charge_weights.resize(2 * charge_count);
        double* weight = workspace.charge_weights.data();
        double* inverse_square = weight + charge_count;
        for (std::size_t c = 0; c < charge_count; ++c) {
          workspace.charge_separations[c] = {coordinates[c].x - centre.x,
                                             coordinates[c].y - centre.y,
                                             coordinates[c].z - centre.z};
        }
        workspace.charge_harmonics.compute(workspace.charge_separations, counts);
        const auto u2 = workspace.charge_harmonics.distances_squared();
        const auto u = workspace.charge_harmonics.distances();
        const double inverse_p = 1.0 / p;

        // Far charges (T = p U^2 >= T_far): the field tensor Phi_Lambda,M = sum q S_Lambda,M(U) /
        // U^(2 Lambda + 1), Lambda <= l + l', as dense dot products over branch-free weights
        // (zero for near charges); Y^(Lambda kappa)_M += Gamma(kappa + Lambda + 3/2) / (2 Lambda
        // + 1) p^-(kappa + Lambda + 3/2) Phi_Lambda,M. Near charges are listed for the exact loop.
        workspace.near_charges.clear();
        for (std::size_t c = 0; c < charge_count; ++c) {
          const bool far = p * u2[c] >= pair.far_threshold;
          weight[c] = far ? charges[c] / u[c] : 0.0;
          inverse_square[c] = far ? 1.0 / u2[c] : 0.0;  // near ones may sit on P
          if (!far && charges[c] != 0.0) {
            workspace.near_charges.push_back(static_cast<std::uint32_t>(c));
          }
        }
        for (int big_l = 0; big_l <= order; ++big_l) {
          for (int m = -big_l; m <= big_l; ++m) {
            const auto row = static_cast<std::size_t>(big_l * big_l + m + big_l);
            phi[row] =
                dot(weight, workspace.charge_harmonics.values(big_l, m).data(), charge_count);
            if (far_field != nullptr) {
              phi[row] += expanded[row];
            }
          }
          if (big_l < order) {
            for (std::size_t c = 0; c < charge_count; ++c) {
              weight[c] *= inverse_square[c];
            }
          }
        }
        for (std::size_t d = 0; d < table.distinct.size(); ++d) {
          const int big_l = table.distinct[d].big_l;
          const double factor = workspace.powers[d * primitive_pairs + ab];
          double* rows = field + table.distinct[d].first_row * primitive_pairs + ab;
          for (int m = -big_l; m <= big_l; ++m) {
            rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                factor * phi[static_cast<std::size_t>(big_l * big_l + m + big_l)];
          }
        }

        // Near charges: exact kernel G_kappa,Lambda = kappa! p^-(kappa + 1)
        // sum_{j=0..kappa} T^j / j! F_(Lambda+j)(T).
        for (const std::uint32_t c : workspace.near_charges) {
          const double q = charges[c];
          const double big_t = p * u2[c];
          boys(order, big_t, workspace.boys);
          for (std::size_t d = 0; d < table.distinct.size(); ++d) {
            const int big_l = table.distinct[d].big_l;
            const int kappa = table.distinct[d].kappa;
            double term = 1.0;
            double sum = f[big_l];
            double scale = q * inverse_p;  // q kappa! p^-(kappa + 1)
            for (int j = 1; j <= kappa; ++j) {
              term *= big_t / j;
              sum += term * f[big_l + j];
              scale *= j * inverse_p;
            }
            const double g = scale * sum;
            double* rows = field + table.distinct[d].first_row * primitive_pairs + ab;
            for (int m = -big_l; m <= big_l; ++m) {
              rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                  g * workspace.charge_harmonics.values(big_l, m)[c];
            }
          }
        }
        // Dipoles (dipole notes Eqs. 22 and 24): far ones (p U^2 >= T_far(n + 1)) form the field
        // tensor Psi_Lambda,M = -(2 Lambda + 1) sum A+_Lambda,M(U, mu) / U^(2 Lambda + 3) (dense
        // dot products over the Racah components mu_(m) / U^(2 Lambda + 3), zero for near
        // dipoles), entering Y like the charge tensor; near ones add H+ A+ + H- A- with Boys
        // orders up to n + 1.
        if (dipole_count > 0) {
          workspace.dipole_separations.resize(dipole_count);
          for (std::size_t d = 0; d < dipole_count; ++d) {
            workspace.dipole_separations[d] = {dipole_coordinates[d].x - centre.x,
                                               dipole_coordinates[d].y - centre.y,
                                               dipole_coordinates[d].z - centre.z};
          }
          workspace.dipole_harmonics.compute(workspace.dipole_separations, dipole_counts);
          const auto du2 = workspace.dipole_harmonics.distances_squared();
          const auto du = workspace.dipole_harmonics.distances();
          workspace.dipole_weights.resize(4 * dipole_count);
          double* dipole_weight = workspace.dipole_weights.data();  // m = -1, 0, 1 rows
          double* dipole_inverse_square = dipole_weight + 3 * dipole_count;
          workspace.near_dipoles.clear();
          for (std::size_t d = 0; d < dipole_count; ++d) {
            const bool far = p * du2[d] >= pair.dipole_far_threshold;
            const double scale = far ? 1.0 / (du2[d] * du[d]) : 0.0;  // 1 / U^3
            for (int m = -1; m <= 1; ++m) {
              dipole_weight[static_cast<std::size_t>(m + 1) * dipole_count + d] =
                  scale * dipole_component(dipoles[d], m);
            }
            dipole_inverse_square[d] = far ? 1.0 / du2[d] : 0.0;  // near ones may sit on P
            if (!far && dipoles[d] != Dipole{}) {
              workspace.near_dipoles.push_back(static_cast<std::uint32_t>(d));
            }
          }
          workspace.dipole_couplings.resize(3 * tensor_size);
          double* psi = workspace.dipole_couplings.data() + 2 * tensor_size;
          std::fill_n(psi, tensor_size, 0.0);
          for (int big_l = 0; big_l <= order; ++big_l) {
            const double factor = -(2.0 * big_l + 1.0);
            for (const DipoleCouplingTerm& term : dipole_couplings(big_l).plus) {
              psi[static_cast<std::size_t>(big_l * big_l + term.big_m + big_l)] +=
                  factor * term.value *
                  dot(dipole_weight + static_cast<std::size_t>(term.m + 1) * dipole_count,
                      workspace.dipole_harmonics.values(big_l + 1, term.harmonic_m).data(),
                      dipole_count);
            }
            if (big_l < order) {
              for (std::size_t w = 0; w < 3 * dipole_count; ++w) {
                dipole_weight[w] *= dipole_inverse_square[w % dipole_count];
              }
            }
          }
          for (std::size_t k = 0; k < table.distinct.size(); ++k) {
            const int big_l = table.distinct[k].big_l;
            const double factor = workspace.powers[k * primitive_pairs + ab];
            double* rows = field + table.distinct[k].first_row * primitive_pairs + ab;
            for (int m = -big_l; m <= big_l; ++m) {
              rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                  factor * psi[static_cast<std::size_t>(big_l * big_l + m + big_l)];
            }
          }
          double* a_plus = workspace.dipole_couplings.data();
          double* a_minus = a_plus + tensor_size;
          for (const std::uint32_t d : workspace.near_dipoles) {
            const Dipole& mu = dipoles[d];
            std::fill_n(a_plus, 2 * tensor_size, 0.0);
            for (int big_l = 0; big_l <= order; ++big_l) {
              const auto& terms = dipole_couplings(big_l);
              for (const DipoleCouplingTerm& term : terms.plus) {
                a_plus[big_l * big_l + term.big_m + big_l] +=
                    term.value * dipole_component(mu, term.m) *
                    workspace.dipole_harmonics.values(big_l + 1, term.harmonic_m)[d];
              }
              for (const DipoleCouplingTerm& term : terms.minus) {
                a_minus[big_l * big_l + term.big_m + big_l] +=
                    term.value * dipole_component(mu, term.m) *
                    workspace.dipole_harmonics.values(big_l - 1, term.harmonic_m)[d];
              }
            }
            boys(order + 1, p * du2[d], workspace.boys);
            for (std::size_t k = 0; k < table.distinct.size(); ++k) {
              const int big_l = table.distinct[k].big_l;
              const auto kernels =
                  dipole_kernels(table.distinct[k].kappa, big_l, p, du2[d], workspace.boys);
              double* rows = field + table.distinct[k].first_row * primitive_pairs + ab;
              for (int m = -big_l; m <= big_l; ++m) {
                const auto row = static_cast<std::size_t>(big_l * big_l + m + big_l);
                rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                    kernels.plus * a_plus[row] + kernels.minus * a_minus[row];
              }
            }
          }
        }

        // Quadrupoles (quadrupole notes Eqs. 23 and 24): far ones (p U^2 >= T_far(n + 2)) form
        // the field tensor Psi_Lambda,M = (2 Lambda + 1)(2 Lambda + 3) / 2 sum B^(Lambda+2)_M(U,
        // theta) / U^(2 Lambda + 5) (dense dot products over theta_mu / U^(2 Lambda + 5), zero
        // for near ones), entering Y like the charge tensor; near ones add 2 (K0 B^(Lambda+2) +
        // K1 B^(Lambda) + K2 B^(Lambda-2)) with Boys orders up to n + 2.
        if (quadrupole_count > 0) {
          workspace.quadrupole_separations.resize(quadrupole_count);
          for (std::size_t e = 0; e < quadrupole_count; ++e) {
            workspace.quadrupole_separations[e] = {quadrupole_coordinates[e].x - centre.x,
                                                   quadrupole_coordinates[e].y - centre.y,
                                                   quadrupole_coordinates[e].z - centre.z};
          }
          workspace.quadrupole_harmonics.compute(workspace.quadrupole_separations,
                                                 quadrupole_counts);
          const auto qu2 = workspace.quadrupole_harmonics.distances_squared();
          const auto qu = workspace.quadrupole_harmonics.distances();
          workspace.quadrupole_weights.resize(6 * quadrupole_count);
          double* theta_weight = workspace.quadrupole_weights.data();  // mu = -2..2 rows
          double* quadrupole_inverse_square = theta_weight + 5 * quadrupole_count;
          workspace.near_quadrupoles.clear();
          for (std::size_t e = 0; e < quadrupole_count; ++e) {
            const bool far = p * qu2[e] >= pair.quadrupole_far_threshold;
            const double scale = far ? 1.0 / (qu2[e] * qu2[e] * qu[e]) : 0.0;  // 1 / U^5
            for (std::size_t mu = 0; mu < 5; ++mu) {
              theta_weight[mu * quadrupole_count + e] = scale * thetas[e][mu];
            }
            quadrupole_inverse_square[e] = far ? 1.0 / qu2[e] : 0.0;  // near ones may sit on P
            if (!far) {
              workspace.near_quadrupoles.push_back(static_cast<std::uint32_t>(e));
            }
          }
          workspace.quadrupole_couplings.resize(4 * tensor_size);
          double* psi_quadrupole = workspace.quadrupole_couplings.data() + 3 * tensor_size;
          std::fill_n(psi_quadrupole, tensor_size, 0.0);
          for (int big_l = 0; big_l <= order; ++big_l) {
            const double factor = 0.5 * (2 * big_l + 1) * (2 * big_l + 3);
            for (const QuadrupoleCouplingTerm& term : quadrupole_couplings(big_l).raised) {
              psi_quadrupole[static_cast<std::size_t>(big_l * big_l + term.big_m + big_l)] +=
                  factor * term.value *
                  dot(theta_weight + static_cast<std::size_t>(term.mu + 2) * quadrupole_count,
                      workspace.quadrupole_harmonics.values(big_l + 2, term.harmonic_m).data(),
                      quadrupole_count);
            }
            if (big_l < order) {
              for (std::size_t w = 0; w < 5 * quadrupole_count; ++w) {
                theta_weight[w] *= quadrupole_inverse_square[w % quadrupole_count];
              }
            }
          }
          for (std::size_t k = 0; k < table.distinct.size(); ++k) {
            const int big_l = table.distinct[k].big_l;
            const double factor = workspace.powers[k * primitive_pairs + ab];
            double* rows = field + table.distinct[k].first_row * primitive_pairs + ab;
            for (int m = -big_l; m <= big_l; ++m) {
              rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                  factor * psi_quadrupole[static_cast<std::size_t>(big_l * big_l + m + big_l)];
            }
          }
          double* b_raised = workspace.quadrupole_couplings.data();
          double* b_same = b_raised + tensor_size;
          double* b_lowered = b_same + tensor_size;
          for (const std::uint32_t e : workspace.near_quadrupoles) {
            const auto& theta = thetas[e];
            std::fill_n(b_raised, 3 * tensor_size, 0.0);
            for (int big_l = 0; big_l <= order; ++big_l) {
              const auto& terms = quadrupole_couplings(big_l);
              const auto add = [&](const std::vector<QuadrupoleCouplingTerm>& list, int degree,
                                   double* out) {
                for (const QuadrupoleCouplingTerm& term : list) {
                  out[big_l * big_l + term.big_m + big_l] +=
                      term.value * theta[static_cast<std::size_t>(term.mu + 2)] *
                      workspace.quadrupole_harmonics.values(degree, term.harmonic_m)[e];
                }
              };
              add(terms.raised, big_l + 2, b_raised);
              add(terms.same, big_l, b_same);
              add(terms.lowered, big_l - 2, b_lowered);
            }
            boys(order + 2, p * qu2[e], workspace.boys);
            for (std::size_t k = 0; k < table.distinct.size(); ++k) {
              const int big_l = table.distinct[k].big_l;
              const auto kernels =
                  quadrupole_kernels(table.distinct[k].kappa, big_l, p, qu2[e], workspace.boys);
              double* rows = field + table.distinct[k].first_row * primitive_pairs + ab;
              for (int m = -big_l; m <= big_l; ++m) {
                const auto row = static_cast<std::size_t>(big_l * big_l + m + big_l);
                rows[static_cast<std::size_t>(m + big_l) * primitive_pairs] +=
                    2.0 * (kernels.k0 * b_raised[row] + kernels.k1 * b_same[row] +
                           kernels.k2 * b_lowered[row]);
              }
            }
          }
        }

        // Contraction prefactors e^(-mu R^2) (-t)^(l - lambda) (1 - t)^(l' - lambda'), with the
        // coefficients folded in for segmented pairs.
        const double exponential = std::exp(-mu_r2);
        const double coefficients =
            segmented ? pair.bra_coefficients[a] * pair.ket_coefficients[b] : 1.0;
        for (std::size_t k = 0; k < table.combinations.size(); ++k) {
          const auto& combination = table.combinations[k];
          workspace.leading[k * primitive_pairs + ab] =
              coefficients * exponential * std::pow(-t, l - combination.lambda) *
              std::pow(1.0 - t, l_prime - combination.lambda_prime);
        }
      }
    }

    // Contraction: Y^(lambda lambda' Lambda)_M = c^T (prefactor * Y^(Lambda kappa)_M) d, kept for
    // all atom pairs of the batch (atom pair innermost) for the assembly below.
    double* contracted = workspace.contracted.data();
    for (std::size_t k = 0; k < table.combinations.size(); ++k) {
      const auto& combination = table.combinations[k];
      const double* prefactor = workspace.leading.data() + k * primitive_pairs;
      const std::size_t rows = static_cast<std::size_t>(2 * combination.big_l + 1);
      for (std::size_t m = 0; m < rows; ++m) {
        const double* y =
            field + (table.distinct[combination.distinct].first_row + m) * primitive_pairs;
        double* out = contracted + (combination.first_row + m) * contracted_pairs * n + i;
        if (segmented) {
          double sum = 0.0;
          for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
            sum += prefactor[ab] * y[ab];
          }
          out[0] = sum;
          continue;
        }
        for (std::size_t bra = 0; bra < na; ++bra) {
          for (std::size_t ket = 0; ket < nb; ++ket) {
            double sum = 0.0;
            for (std::size_t a = 0; a < ka; ++a) {
              const double c = pair.bra_coefficients[bra * ka + a];
              if (c == 0.0) {
                continue;
              }
              double inner = 0.0;
              for (std::size_t b = 0; b < kb; ++b) {
                inner +=
                    prefactor[a * kb + b] * y[a * kb + b] * pair.ket_coefficients[b * nb + ket];
              }
              sum += c * inner;
            }
            out[(bra * nb + ket) * n] = sum;
          }
        }
      }
    }
  }

  // Assembly over the whole batch, streaming contiguously over the n atom pairs.
  const double* contracted = workspace.contracted.data();
  if (!pair.factorized) {
    // Flat: V_(IJ),(mm') = 2 pi sum Gamma Y^(lambda lambda' Lambda)_M, one term at a time.
    for (const auto& term : table.terms) {
      const double* bra_harmonic = harmonics.values(term.bra_l, term.bra_m).data();
      const double* ket_harmonic = harmonics.values(term.ket_l, term.ket_m).data();
      const double scale = 2.0 * std::numbers::pi * term.value;
      for (std::size_t ij = 0; ij < contracted_pairs; ++ij) {
        const double* y = contracted + (term.row * contracted_pairs + ij) * n;
        double* out = values.data() + (ij * components + term.component) * n;
        for (std::size_t i = 0; i < n; ++i) {
          out[i] += scale * bra_harmonic[i] * ket_harmonic[i] * y[i];
        }
      }
    }
    return;
  }

  // Factorized: geometric translation matrices G^(lambda), H^(lambda') (shared by all contracted
  // pairs), then per contracted pair Z = Gaunt contraction, W = G Z, V += 2 pi W H^T.
  const FactorizedTable& factorized = factorized_table(l, l_prime);
  const auto lambdas = static_cast<std::size_t>(l + 1);
  const auto lambda_primes = static_cast<std::size_t>(l_prime + 1);
  const auto bra_components = static_cast<std::size_t>(2 * l + 1);
  const auto ket_components = static_cast<std::size_t>(2 * l_prime + 1);
  const auto translation_matrices =
      [&](const std::vector<std::vector<FactorizedTable::TranslationEntry>>& entries,
          std::size_t rows, std::vector<double>& matrices, std::vector<std::size_t>& offsets) {
        offsets.assign(entries.size() + 1, 0);
        for (std::size_t lambda = 0; lambda < entries.size(); ++lambda) {
          offsets[lambda + 1] = offsets[lambda] + rows * (2 * lambda + 1);
        }
        matrices.assign(offsets.back() * n, 0.0);
        for (std::size_t lambda = 0; lambda < entries.size(); ++lambda) {
          for (const auto& entry : entries[lambda]) {
            const double* s_values = harmonics.values(entry.harmonic_l, entry.harmonic_m).data();
            double* out = matrices.data() + (offsets[lambda] + entry.cell) * n;
            for (std::size_t i = 0; i < n; ++i) {
              out[i] += entry.value * s_values[i];
            }
          }
        }
      };
  std::vector<std::size_t> bra_offsets, ket_offsets;
  translation_matrices(factorized.bra, bra_components, workspace.translation_bra, bra_offsets);
  translation_matrices(factorized.ket, ket_components, workspace.translation_ket, ket_offsets);
  const double* g = workspace.translation_bra.data();
  const double* h = workspace.translation_ket.data();
  const double two_pi = 2.0 * std::numbers::pi;
  for (std::size_t ij = 0; ij < contracted_pairs; ++ij) {
    for (std::size_t lambda_prime = 0; lambda_prime < lambda_primes; ++lambda_prime) {
      const std::size_t columns = 2 * lambda_prime + 1;
      workspace.bra_transformed.assign(bra_components * columns * n, 0.0);
      double* w = workspace.bra_transformed.data();
      for (std::size_t lambda = 0; lambda < lambdas; ++lambda) {
        const auto& gaunt = factorized.gaunt[lambda * lambda_primes + lambda_prime];
        if (gaunt.empty()) {
          continue;
        }
        const std::size_t inner = 2 * lambda + 1;
        workspace.gaunt_contracted.assign(inner * columns * n, 0.0);
        double* z = workspace.gaunt_contracted.data();
        for (const auto& term : gaunt) {
          const double* y = contracted + (term.row * contracted_pairs + ij) * n;
          double* out = z + term.z * n;
          for (std::size_t i = 0; i < n; ++i) {
            out[i] += term.value * y[i];
          }
        }
        for (const std::size_t cell : factorized.bra_cells[lambda]) {
          const std::size_t m = cell / inner;
          const std::size_t mu = cell % inner;
          const double* g_row = g + (bra_offsets[lambda] + cell) * n;
          for (std::size_t mu_prime = 0; mu_prime < columns; ++mu_prime) {
            const double* z_row = z + (mu * columns + mu_prime) * n;
            double* w_row = w + (m * columns + mu_prime) * n;
            for (std::size_t i = 0; i < n; ++i) {
              w_row[i] += g_row[i] * z_row[i];
            }
          }
        }
      }
      for (const std::size_t cell : factorized.ket_cells[lambda_prime]) {
        const std::size_t m_prime = cell / columns;
        const std::size_t mu_prime = cell % columns;
        const double* h_row = h + (ket_offsets[lambda_prime] + cell) * n;
        for (std::size_t m = 0; m < bra_components; ++m) {
          const double* w_row = w + (m * columns + mu_prime) * n;
          double* out = values.data() + (ij * components + m * ket_components + m_prime) * n;
          for (std::size_t i = 0; i < n; ++i) {
            out[i] += two_pi * w_row[i] * h_row[i];
          }
        }
      }
    }
  }
}

}  // namespace fika::detail
