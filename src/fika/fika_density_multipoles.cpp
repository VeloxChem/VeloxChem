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

#include "fika_density_multipoles.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

#include "fika_basis_shell.hpp"
#include "fika_nuclear_attraction_kernels.hpp"
#include "fika_boys.hpp"
#include "fika_solid_harmonics.hpp"

namespace fika::detail {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::density_multipoles: " + reason);
}

/// Shell-pair blocks of the density: atoms A <= B, and shells s <= s' when A == B.
struct Block {
  std::size_t bra_atom;
  std::size_t ket_atom;
  std::size_t bra_shell;
  std::size_t ket_shell;
  double weight;  // 2 for a mirrored block, else 1
};

auto density_blocks(const MolecularBasis& basis) -> std::vector<Block> {
  std::vector<Block> blocks;
  const std::size_t atoms = basis.atom_count();
  for (std::size_t a = 0; a < atoms; ++a) {
    const std::size_t bra_shells = basis.atom_basis(a).shells().size();
    for (std::size_t b = a; b < atoms; ++b) {
      const std::size_t ket_shells = basis.atom_basis(b).shells().size();
      for (std::size_t s = 0; s < bra_shells; ++s) {
        for (std::size_t t = a == b ? s : 0; t < ket_shells; ++t) {
          blocks.push_back({a, b, s, t, a == b && s == t ? 1.0 : 2.0});
        }
      }
    }
  }
  return blocks;
}

}  // namespace

namespace {

/// Scratch storage of append_block.
struct BuildScratch {
  SolidHarmonics harmonics;
  std::vector<double> block_density;
  std::vector<double> effective;
  std::vector<double> primitive;
  std::vector<double> weights;
  std::vector<double> moments;
};

/// Appends the sources of one shell-pair block to `sources` (`pair`: its shell pair, built with
/// no screening; `threshold`: the block's potential-element threshold).
void append_block(const Block& block, const MolecularBasis& basis,
                  std::span<const Point3D<double>> coordinates, const DenseMatrix& density,
                  const TwoAtomShellPair& pair, std::size_t block_count, double accuracy,
                  DensityMultipoles& sources, BuildScratch& scratch) {
  const int l = pair.l;
  const int l_prime = pair.l_prime;
  const std::size_t na = pair.bra_contractions;
  const std::size_t nb = pair.ket_contractions;
  const auto ket_components = static_cast<std::size_t>(2 * l_prime + 1);
  const std::size_t components = static_cast<std::size_t>(2 * l + 1) * ket_components;
  const std::span<const double> alphas = pair.bra_exponents;
  const std::span<const double> betas = pair.ket_exponents;
  sources.primitive_pairs += alphas.size() * betas.size();
  // D_(Im),(Jm') at [(I nb + J) components + (m + l)(2l' + 1) + m' + l'].
  auto& block_density = scratch.block_density;
  block_density.assign(na * nb * components, 0.0);
  double norm = 0.0;
  for (std::size_t i = 0; i < na; ++i) {
    for (std::size_t j = 0; j < nb; ++j) {
      for (int m = -l; m <= l; ++m) {
        for (int m_prime = -l_prime; m_prime <= l_prime; ++m_prime) {
          const double value =
              density(basis.function_index(block.bra_atom, block.bra_shell, i, m),
                      basis.function_index(block.ket_atom, block.ket_shell, j, m_prime));
          block_density[(i * nb + j) * components +
                        static_cast<std::size_t>(m + l) * ket_components +
                        static_cast<std::size_t>(m_prime + l_prime)] = value;
          norm += std::abs(value);
        }
      }
    }
  }
  if (norm == 0.0) {
    return;
  }
  // The block's threshold shared by its primitive pairs (as make_two_atom_shell_pair does).
  const double primitive_threshold = accuracy /
                                     (static_cast<double>(block_count) * block.weight * norm) /
                                     static_cast<double>(alphas.size() * betas.size());
  const Point3D<double>& a_position = coordinates[block.bra_atom];
  const Point3D<double>& b_position = coordinates[block.ket_atom];
  const Point3D<double> r{a_position.x - b_position.x, a_position.y - b_position.y,
                          a_position.z - b_position.z};
  const int order = l + l_prime;
  const std::vector<std::size_t> counts(static_cast<std::size_t>(order) + 1, 1);
  scratch.harmonics.compute(std::span(&r, 1), counts);
  const double r2 = scratch.harmonics.distances_squared()[0];
  const double r_power = std::max(1.0, std::pow(std::sqrt(r2), order));
  const auto moment_rows = static_cast<std::size_t>((order + 1) * (order + 1));
  const std::size_t rows = field_row_count(l, l_prime);
  const auto row_blocks = field_row_blocks(l, l_prime);
  auto& effective = scratch.effective;
  auto& primitive = scratch.primitive;
  auto& weights = scratch.weights;
  auto& moments = scratch.moments;
  primitive.resize(components * rows);
  effective.resize(components);
  for (std::size_t a = 0; a < alphas.size(); ++a) {
    for (std::size_t b = 0; b < betas.size(); ++b) {
      const double p = alphas[a] + betas[b];
      const double mu_r2 = alphas[a] * betas[b] / p * r2;
      if (pair.bra_largest[a] * pair.ket_largest[b] * pair.dipole_prefactor / std::sqrt(p) *
              r_power * std::exp(-mu_r2) <
          primitive_threshold) {
        continue;  // the dipole-potential kernels skip it too
      }
      // D contracted with the coefficients: sum_IJ c_aI d_bJ D_(Im),(Jm').
      std::fill(effective.begin(), effective.end(), 0.0);
      for (std::size_t i = 0; i < na; ++i) {
        const double c = pair.bra_coefficients[i * alphas.size() + a];
        for (std::size_t j = 0; j < nb; ++j) {
          const double cd = block.weight * c * pair.ket_coefficients[b * nb + j];
          const double* d = block_density.data() + (i * nb + j) * components;
          for (std::size_t k = 0; k < components; ++k) {
            effective[k] += cd * d[k];
          }
        }
      }
      primitive_pair_field_weights(pair, scratch.harmonics, a, b, primitive);
      weights.assign(rows, 0.0);
      for (std::size_t k = 0; k < components; ++k) {
        for (std::size_t row = 0; row < rows; ++row) {
          weights[row] += effective[k] * primitive[k * rows + row];
        }
      }
      if (std::ranges::all_of(weights, [](double w) { return w == 0.0; })) {
        continue;
      }
      moments.assign(moment_rows, 0.0);
      for (const FieldRowBlock& row_block : row_blocks) {
        const int big_l = row_block.big_l;
        const double s = row_block.kappa + big_l + 1.5;
        const double factor = std::tgamma(s) / (2 * big_l + 1) * std::pow(p, -s);
        for (int m = -big_l; m <= big_l; ++m) {
          moments[static_cast<std::size_t>(big_l * big_l + m + big_l)] +=
              factor * weights[row_block.first_row + static_cast<std::size_t>(m + big_l)];
        }
      }
      const double t = betas[b] / p;
      sources.centres.push_back(
          {a_position.x - t * r.x, a_position.y - t * r.y, a_position.z - t * r.z});
      sources.exponents.push_back(p);
      sources.penetration_squared.push_back(pair.dipole_far_threshold / p);
      sources.bra_l.push_back(l);
      sources.ket_l.push_back(l_prime);
      sources.weights.insert(sources.weights.end(), weights.begin(), weights.end());
      sources.weight_offsets.push_back(sources.weights.size());
      sources.moments.insert(sources.moments.end(), moments.begin(), moments.end());
      sources.offsets.push_back(sources.moments.size());
    }
  }
}

/// Appends `part` (offsets from zero) to `sources`.
void append_sources(DensityMultipoles& sources, const DensityMultipoles& part) {
  const std::size_t weight_base = sources.weights.size();
  const std::size_t moment_base = sources.moments.size();
  sources.centres.insert(sources.centres.end(), part.centres.begin(), part.centres.end());
  sources.exponents.insert(sources.exponents.end(), part.exponents.begin(), part.exponents.end());
  sources.penetration_squared.insert(sources.penetration_squared.end(),
                                     part.penetration_squared.begin(),
                                     part.penetration_squared.end());
  sources.bra_l.insert(sources.bra_l.end(), part.bra_l.begin(), part.bra_l.end());
  sources.ket_l.insert(sources.ket_l.end(), part.ket_l.begin(), part.ket_l.end());
  sources.weights.insert(sources.weights.end(), part.weights.begin(), part.weights.end());
  sources.moments.insert(sources.moments.end(), part.moments.begin(), part.moments.end());
  for (std::size_t s = 1; s < part.weight_offsets.size(); ++s) {
    sources.weight_offsets.push_back(weight_base + part.weight_offsets[s]);
    sources.offsets.push_back(moment_base + part.offsets[s]);
  }
  sources.primitive_pairs += part.primitive_pairs;
}

}  // namespace

auto density_multipoles(const Molecule<double>& molecule, const MolecularBasis& basis,
                        const DenseMatrix& density, double accuracy) -> DensityMultipoles {
  const std::size_t n = basis.function_count();
  if (molecule.size() != basis.atom_count()) {
    fail("molecule has " + std::to_string(molecule.size()) + " atoms, basis " +
         std::to_string(basis.atom_count()));
  }
  if (density.rows() != n || density.columns() != n ||
      density.symmetry() != MatrixSymmetry::symmetric) {
    fail("the density must be a symmetric " + std::to_string(n) + " x " + std::to_string(n) +
         " matrix");
  }
  if (!(accuracy > 0.0) || !std::isfinite(accuracy)) {
    fail("the accuracy must be positive and finite");
  }
  const auto blocks = density_blocks(basis);
  // Shell pairs (exponents, coefficients, bounds) per pair of shells of the unique bases.
  const auto bases = basis.atom_basis_indices();
  std::map<std::array<std::size_t, 4>, std::size_t> pair_index;
  std::vector<TwoAtomShellPair> pairs;
  std::vector<std::size_t> block_pairs(blocks.size());
  for (std::size_t k = 0; k < blocks.size(); ++k) {
    const Block& block = blocks[k];
    const std::array<std::size_t, 4> key{bases[block.bra_atom], block.bra_shell,
                                         bases[block.ket_atom], block.ket_shell};
    const auto [found, inserted] = pair_index.try_emplace(key, pairs.size());
    if (inserted) {
      pairs.push_back(make_two_atom_shell_pair(
          basis.atom_basis(block.bra_atom).shells()[block.bra_shell],
          basis.atom_basis(block.ket_atom).shells()[block.ket_shell], {.dipoles = 1.0}, 0.0));
    }
    block_pairs[k] = found->second;
  }
  // Fixed chunks of blocks in parallel, concatenated in order (independent of the thread count).
  constexpr std::size_t chunk = 64;
  const std::size_t chunks = (blocks.size() + chunk - 1) / chunk;
  std::vector<DensityMultipoles> parts(chunks);
  const auto coordinates = molecule.coordinates();
#pragma omp parallel
  {
    BuildScratch scratch;
#pragma omp for schedule(dynamic, 1)
    for (std::size_t c = 0; c < chunks; ++c) {
      DensityMultipoles& part = parts[c];
      part.offsets.push_back(0);
      part.weight_offsets.push_back(0);
      for (std::size_t k = c * chunk; k < std::min(blocks.size(), (c + 1) * chunk); ++k) {
        append_block(blocks[k], basis, coordinates, density, pairs[block_pairs[k]], blocks.size(),
                     accuracy, part, scratch);
      }
    }
  }
  DensityMultipoles sources;
  sources.offsets.push_back(0);
  sources.weight_offsets.push_back(0);
  for (const DensityMultipoles& part : parts) {
    append_sources(sources, part);
  }
  return sources;
}

namespace {

/// Adds to e (Racah components m = -1, 0, 1) the field at separation U = site - P (|U| = u,
/// harmonics of `harmonics` column i up to order rank + 1) of source s: far through its moments,
/// E_(m) = -sum_L (2L + 1) / U^(2L + 3) sum_M Q_LM A+_LM(U, e_(m)); near through its weights,
/// E_(m) = sum_rows W_row (H+ A+_row(U, e_(m)) + H- A-_row(U, e_(m))).
void add_source_field(const DensityMultipoles& sources, std::size_t s,
                      const SolidHarmonics& harmonics, std::size_t i, bool far,
                      std::vector<double>& boys_values, double* e) {
  const int rank = sources.rank(s);
  const double u2 = harmonics.distances_squared()[i];
  if (far) {
    const auto moments = sources.moments_of(s);
    double inverse = 1.0 / (u2 * harmonics.distances()[i]);  // U^-(2L + 3) for L = 0
    for (int big_l = 0; big_l <= rank; ++big_l) {
      const double scale = -(2 * big_l + 1) * inverse;
      for (const DipoleCouplingTerm& term : dipole_couplings(big_l).plus) {
        e[term.m + 1] += scale *
                         moments[static_cast<std::size_t>(big_l * big_l + term.big_m + big_l)] *
                         term.value * harmonics.values(big_l + 1, term.harmonic_m)[i];
      }
      inverse /= u2;
    }
    return;
  }
  const double p = sources.exponents[s];
  boys_values.resize(static_cast<std::size_t>(rank) + 2);
  boys(rank + 1, p * u2, boys_values);
  const auto weights = sources.weights_of(s);
  for (const FieldRowBlock& row_block : field_row_blocks(sources.bra_l[s], sources.ket_l[s])) {
    const int big_l = row_block.big_l;
    const DipoleKernels kernels = dipole_kernels(row_block.kappa, big_l, p, u2, boys_values);
    const double* w = weights.data() + row_block.first_row + big_l;  // w[M]
    const auto& couplings = dipole_couplings(big_l);
    for (const DipoleCouplingTerm& term : couplings.plus) {
      e[term.m + 1] += kernels.plus * w[term.big_m] * term.value *
                       harmonics.values(big_l + 1, term.harmonic_m)[i];
    }
    for (const DipoleCouplingTerm& term : couplings.minus) {
      e[term.m + 1] += kernels.minus * w[term.big_m] * term.value *
                       harmonics.values(big_l - 1, term.harmonic_m)[i];
    }
  }
}

}  // namespace

void add_density_field_block(const DensityMultipoles& sources, std::span<const std::uint32_t> which,
                             std::span<const Point3D<double>> sites, std::span<double> sums,
                             DensityFieldWorkspace& workspace) {
  const std::size_t count = sites.size();
  assert(sums.size() >= 3 * count);
  workspace.separations.resize(count);
  for (const std::uint32_t s : which) {
    const Point3D<double>& centre = sources.centres[s];
    for (std::size_t i = 0; i < count; ++i) {
      workspace.separations[i] = {sites[i].x - centre.x, sites[i].y - centre.y,
                                  sites[i].z - centre.z};
    }
    const std::vector<std::size_t> counts(static_cast<std::size_t>(sources.rank(s)) + 2, count);
    workspace.harmonics.compute(workspace.separations, counts);
    const auto u2 = workspace.harmonics.distances_squared();
    for (std::size_t i = 0; i < count; ++i) {
      add_source_field(sources, s, workspace.harmonics, i, u2[i] >= sources.penetration_squared[s],
                       workspace.boys, sums.data() + 3 * i);
    }
  }
}

void add_density_field(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                       std::span<Point3D<double>> field) {
  if (field.size() != sites.size()) {
    throw std::invalid_argument("fika::add_density_field: field has " +
                                std::to_string(field.size()) + " entries for " +
                                std::to_string(sites.size()) + " sites");
  }
  std::vector<std::uint32_t> all(sources.size());
  for (std::size_t s = 0; s < all.size(); ++s) {
    all[s] = static_cast<std::uint32_t>(s);
  }
  constexpr std::size_t chunk = 64;
  const std::size_t chunks = (sites.size() + chunk - 1) / chunk;
#pragma omp parallel
  {
    std::vector<double> sums;
    DensityFieldWorkspace workspace;
#pragma omp for schedule(dynamic, 1)
    for (std::size_t c = 0; c < chunks; ++c) {
      const std::size_t begin = c * chunk;
      const std::size_t count = std::min(chunk, sites.size() - begin);
      sums.assign(3 * count, 0.0);
      add_density_field_block(sources, all, sites.subspan(begin, count), sums, workspace);
      for (std::size_t i = 0; i < count; ++i) {
        field[begin + i].x += sums[3 * i + 2];
        field[begin + i].y += sums[3 * i];
        field[begin + i].z += sums[3 * i + 1];
      }
    }
  }
}

void add_far_density_field(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                           std::span<Point3D<double>> field) {
  if (field.size() != sites.size()) {
    throw std::invalid_argument("fika::add_far_density_field: field has " +
                                std::to_string(field.size()) + " entries for " +
                                std::to_string(sites.size()) + " sites");
  }
  std::vector<Point3D<double>> separations(sites.size());
  SolidHarmonics harmonics;
  for (std::size_t s = 0; s < sources.size(); ++s) {
    const Point3D<double>& centre = sources.centres[s];
    const int rank = sources.rank(s);
    for (std::size_t i = 0; i < sites.size(); ++i) {
      separations[i] = {sites[i].x - centre.x, sites[i].y - centre.y, sites[i].z - centre.z};
      const double u2 = separations[i].x * separations[i].x + separations[i].y * separations[i].y +
                        separations[i].z * separations[i].z;
      if (u2 < sources.penetration_squared[s]) {
        throw std::domain_error("fika::add_far_density_field: site " + std::to_string(i) +
                                " within the penetration radius of source " + std::to_string(s));
      }
    }
    const std::vector<std::size_t> counts(static_cast<std::size_t>(rank) + 2, sites.size());
    harmonics.compute(separations, counts);
    std::vector<double> boys_values;  // unused by far sources
    for (std::size_t i = 0; i < sites.size(); ++i) {
      double e[3] = {0.0, 0.0, 0.0};  // m = -1, 0, 1
      add_source_field(sources, s, harmonics, i, true, boys_values, e);
      field[i].x += e[2];
      field[i].y += e[0];
      field[i].z += e[1];
    }
  }
}

}  // namespace fika::detail
