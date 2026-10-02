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

#include "fika_two_centre_kernels.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <mutex>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <variant>

#include "fika_gaussian_normalization.hpp"
#include "fika_gaunt.hpp"
#include "fika_vector_math.hpp"

namespace fika::detail {

namespace {

/// Exponents and effective coefficients of an uncontracted or segmented shell.
struct Primitives {
  std::span<const double> exponents;
  std::span<const double> coefficients;
};

auto primitives_of(const BasisShell& shell) -> Primitives {
  return std::visit(
      [](const auto& s) -> Primitives {
        using Shell = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<Shell, GeneralShell>) {
          throw std::logic_error("fika: overlap kernel does not handle general contractions");
        } else {
          return {s.exponents(), s.coefficients()};
        }
      },
      shell);
}

/// Angular-assembly term: row `component` += value * T_coupling * S_{L,M}.
struct AssemblyTerm {
  std::size_t component;
  std::size_t coupling;  // index of L in |l - l'|, |l - l'| + 2, ...
  int big_l;
  int big_m;
  double value;
};

/// Gaunt terms of the correction (L < l + l'), for l == l' only components with m <= m'.
auto build_assembly_terms(int l, int l_prime) -> std::vector<AssemblyTerm> {
  std::vector<AssemblyTerm> terms;
  const int lowest = std::abs(l - l_prime);
  for (const GauntEntry& entry : gaunt_coefficients(l, l_prime)) {
    if (entry.big_l >= l + l_prime || (l == l_prime && entry.m_prime < entry.m)) {
      continue;
    }
    terms.push_back(
        {static_cast<std::size_t>((entry.m + l) * (2 * l_prime + 1) + entry.m_prime + l_prime),
         static_cast<std::size_t>((entry.big_l - lowest) / 2), entry.big_l, entry.big_m,
         entry.value});
  }
  return terms;
}

/// One pair's assembly terms, built on first request.
struct AssemblySlot {
  std::once_flag built;
  std::vector<AssemblyTerm> terms;
};

auto assembly_terms(int l, int l_prime) -> const std::vector<AssemblyTerm>& {
  static std::array<std::array<AssemblySlot, max_angular_momentum + 1>, max_angular_momentum + 1>
      table;
  AssemblySlot& slot = table[static_cast<std::size_t>(l)][static_cast<std::size_t>(l_prime)];
  std::call_once(slot.built, [&] { slot.terms = build_assembly_terms(l, l_prime); });
  return slot.terms;
}

/// values[(m, m')] = W0 S_lm S_l'm' + sum_L T_L sum_M C^{LM}_{lm,l'm'} S_LM over n columns.
template <int L, int LPrime>
void assemble(const SolidHarmonics& harmonics, std::size_t n, const double* w0,
              const double* radial, std::span<double> values) {
  constexpr int ket_components = 2 * LPrime + 1;
  for (int m = -L; m <= L; ++m) {
    const double* bra = harmonics.values(L, m).data();
    for (int m_prime = (L == LPrime ? m : -LPrime); m_prime <= LPrime; ++m_prime) {
      const double* ket = harmonics.values(LPrime, m_prime).data();
      double* row =
          values.data() + static_cast<std::size_t>((m + L) * ket_components + m_prime + LPrime) * n;
      for (std::size_t i = 0; i < n; ++i) {
        row[i] = w0[i] * bra[i] * ket[i];
      }
    }
  }
  for (const AssemblyTerm& term : assembly_terms(L, LPrime)) {
    double* row = values.data() + term.component * n;
    const double* t = radial + term.coupling * n;
    const double* s = harmonics.values(term.big_l, term.big_m).data();
    const double c = term.value;
    for (std::size_t i = 0; i < n; ++i) {
      row[i] += c * t[i] * s[i];
    }
  }
}

void assemble_dispatch(int l, int l_prime, const SolidHarmonics& harmonics, std::size_t n,
                       const double* w0, const double* radial, std::span<double> values) {
  dispatch_angular_momentum(l, [&]<int A>(std::integral_constant<int, A>) {
    dispatch_angular_momentum(l_prime, [&]<int B>(std::integral_constant<int, B>) {
      assemble<A, B>(harmonics, n, w0, radial, values);
    });
  });
}

/// Primitive factor (-1)^l b^l a^l' p^-(l+l') (pi/p)^(3/2) without contraction coefficients.
auto primitive_factor(int l, int l_prime, double alpha, double beta) -> double {
  const double p = alpha + beta;
  const double q = std::numbers::pi / p;
  const double sign = l % 2 == 0 ? 1.0 : -1.0;
  return sign * std::pow(beta, l) * std::pow(alpha, l_prime) / std::pow(p, l + l_prime) * q *
         std::sqrt(q);
}

/// D^(t)_{kappa L} = binom(kappa, t) (kappa + L + 1/2)(kappa + L - 1/2)...(kappa - t + L + 3/2)
/// for t = 1..kappa, appended to `d`.
void append_d_coefficients(int kappa, int big_l, std::vector<double>& d) {
  double binomial = 1.0;
  double product = 1.0;
  for (int t = 1; t <= kappa; ++t) {
    binomial = binomial * (kappa - t + 1) / t;
    product *= kappa + big_l + 1.5 - t;
    d.push_back(binomial * product);
  }
}

auto sign(int t) -> double {
  return t % 2 == 0 ? 1.0 : -1.0;
}

/// Operator-dependent part of the common form (see two_centre_kernels.hpp).
auto make_form(int l, int l_prime, TwoCentreIntegral integral) -> IntegralForm {
  const int n = l + l_prime;
  const int k_max = std::min(l, l_prime);
  IntegralForm form;
  form.integral = integral;
  const bool kinetic = integral == TwoCentreIntegral::kinetic;
  form.first_power = kinetic ? -2 : 0;
  form.power_count = static_cast<std::size_t>(kinetic ? k_max + 2 : k_max + 1);
  form.w0_quadratic = kinetic ? -2.0 : 0.0;
  form.w0_constant = kinetic ? 2.0 * (n + 1.5) : 1.0;
  form.w0_constant_row = kinetic ? 1 : 0;
  form.radial_offsets.push_back(0);
  std::vector<double> d;
  for (int big_l = std::abs(l - l_prime); big_l <= n - 2; big_l += 2) {
    const int k = (n - big_l) / 2;
    d.clear();
    if (kinetic) {
      append_d_coefficients(k, big_l, d);  // D^(1)_kL
      form.radial_coefficients.push_back(2.0 * d[0]);
      d.clear();
      append_d_coefficients(k + 1, big_l, d);  // D^(t)_(k+1)L
      for (int t = 2; t <= k + 1; ++t) {
        form.radial_coefficients.push_back(-2.0 * (sign(t) * d[static_cast<std::size_t>(t - 1)]));
      }
    } else {
      append_d_coefficients(k, big_l, d);
      for (int t = 1; t <= k; ++t) {
        form.radial_coefficients.push_back(sign(t) * d[static_cast<std::size_t>(t - 1)]);
      }
    }
    form.radial_offsets.push_back(static_cast<int>(form.radial_coefficients.size()));
  }
  return form;
}

auto coupling_count(const IntegralForm& form) -> std::size_t {
  return form.radial_offsets.size() - 1;
}

/// Primitive factors of power rows j = first_power.. (factor times mu^-j), row-major over
/// `primitive_pairs`, at column `index`.
void set_power_factors(const IntegralForm& form, double factor, double mu, std::size_t index,
                       std::size_t primitive_pairs, std::vector<double>& factors) {
  for (int power = form.first_power; power < 0; ++power) {
    factor *= mu;
  }
  for (std::size_t row = 0; row < form.power_count; ++row) {
    factors[row * primitive_pairs + index] = factor;
    factor /= mu;
  }
}

/// w0 over n: w0_quadratic R^2 W(0) + w0_constant W(w0_constant_row), or W(0) itself for the
/// overlap (returned without copying).
auto leading_term(const IntegralForm& form, std::span<const double> r2, const double* w,
                  std::size_t n, std::vector<double>& leading) -> const double* {
  if (form.w0_quadratic == 0.0 && form.w0_constant == 1.0 && form.w0_constant_row == 0) {
    return w;
  }
  leading.resize(n);
  const double* constant_row = w + form.w0_constant_row * n;
  for (std::size_t i = 0; i < n; ++i) {
    leading[i] = form.w0_quadratic * r2[i] * w[i] + form.w0_constant * constant_row[i];
  }
  return leading.data();
}

/// Dense transposed contraction matrix of a shell: N x K, row-major (effective coefficients).
auto transposed_coefficients(const BasisShell& shell) -> std::vector<double> {
  return std::visit(
      [](const auto& s) {
        const std::size_t primitives = s.primitive_count();
        std::vector<double> matrix(s.contraction_count() * primitives, 0.0);
        for (std::size_t k = 0; k < s.contraction_count(); ++k) {
          for_each_coefficient(s, k, [&](std::size_t i, double coefficient) {
            matrix[k * primitives + i] = coefficient;
          });
        }
        return matrix;
      },
      shell);
}

}  // namespace

auto make_segmented_shell_pair(const BasisShell& bra, const BasisShell& ket,
                               TwoCentreIntegral integral) -> SegmentedShellPair {
  const Primitives a = primitives_of(bra);
  const Primitives b = primitives_of(ket);
  SegmentedShellPair pair;
  pair.l = angular_momentum(bra);
  pair.l_prime = angular_momentum(ket);
  pair.form = make_form(pair.l, pair.l_prime, integral);
  pair.primitive_pairs = a.exponents.size() * b.exponents.size();
  pair.mu.reserve(pair.primitive_pairs);
  pair.factors.resize(pair.form.power_count * pair.primitive_pairs);
  std::size_t index = 0;
  for (std::size_t i = 0; i < a.exponents.size(); ++i) {
    const double alpha = a.exponents[i];
    for (std::size_t j = 0; j < b.exponents.size(); ++j, ++index) {
      const double beta = b.exponents[j];
      const double mu = alpha * beta / (alpha + beta);
      pair.mu.push_back(mu);
      set_power_factors(pair.form,
                        a.coefficients[i] * b.coefficients[j] *
                            primitive_factor(pair.l, pair.l_prime, alpha, beta),
                        mu, index, pair.primitive_pairs, pair.factors);
    }
  }
  return pair;
}

void segmented_shell_pair_values(const SegmentedShellPair& pair, const SolidHarmonics& harmonics,
                                 std::size_t n, KernelWorkspace& workspace,
                                 std::span<double> values) {
  assert(values.size() >= static_cast<std::size_t>((2 * pair.l + 1) * (2 * pair.l_prime + 1)) * n);
  if (n == 0) {
    return;
  }
  const IntegralForm& form = pair.form;
  const auto r2 = harmonics.distances_squared();
  const std::size_t couplings = coupling_count(form);
  const std::size_t rows = form.power_count;
  workspace.accumulators.resize(rows * n);
  workspace.radial.resize(std::max<std::size_t>(couplings, 1) * n);
  double* w = workspace.accumulators.data();
  double* t = workspace.radial.data();
  const double* w0 = w;

  if (pair.primitive_pairs == 1) {
    // Simplified path: W(t) = q_t exp(-mu R^2), so w0 and T_L are polynomials in R^2 times one
    // exponential. For the overlap (q_t = q_0 mu^-t) T_L = W(0) sum_t r_t mu^-t R^(2(deg - t)).
    const bool overlap = form.integral == TwoCentreIntegral::overlap;
    exp_scaled_negative(r2, pair.mu[0], std::span(w, n));  // exp(-mu R^2)
    const double* multiplier = w;
    if (overlap) {
      const double q0 = pair.factors[0];
      for (std::size_t i = 0; i < n; ++i) {
        w[i] *= q0;
      }
    } else {
      workspace.leading.resize(n);
      const double quadratic = form.w0_quadratic * pair.factors[0];
      const double constant = form.w0_constant * pair.factors[form.w0_constant_row];
      for (std::size_t i = 0; i < n; ++i) {
        workspace.leading[i] = (quadratic * r2[i] + constant) * w[i];
      }
      w0 = workspace.leading.data();
    }
    const double inverse_mu = 1.0 / pair.mu[0];
    for (std::size_t c = 0; c < couplings; ++c) {
      const auto first = static_cast<std::size_t>(form.radial_offsets[c]);
      const auto deg = static_cast<std::size_t>(form.radial_offsets[c + 1]) - first;
      std::array<double, max_angular_momentum + 2> coefficient{};  // r_t times its row factor
      double mu_power = 1.0;
      for (std::size_t k = 1; k <= deg; ++k) {
        mu_power *= inverse_mu;
        coefficient[k] =
            form.radial_coefficients[first + k - 1] * (overlap ? mu_power : pair.factors[k]);
      }
      double* row = t + c * n;
      for (std::size_t i = 0; i < n; ++i) {
        double polynomial = coefficient[1];  // Horner in R^2, highest power first
        for (std::size_t k = 2; k <= deg; ++k) {
          polynomial = polynomial * r2[i] + coefficient[k];
        }
        row[i] = multiplier[i] * polynomial;
      }
    }
  } else {
    // Stage 1: exponentials and accumulators W(t) = sum_ab q_ab,t exp(-mu_ab R^2).
    workspace.exponentials.resize(pair.primitive_pairs * n);
    double* e = workspace.exponentials.data();
    for (std::size_t ab = 0; ab < pair.primitive_pairs; ++ab) {
      exp_scaled_negative(r2, pair.mu[ab], std::span(e + ab * n, n));
    }
    std::fill(w, w + rows * n, 0.0);
    for (std::size_t k = 0; k < rows; ++k) {
      double* row = w + k * n;
      const double* q = pair.factors.data() + k * pair.primitive_pairs;
      for (std::size_t ab = 0; ab < pair.primitive_pairs; ++ab) {
        const double factor = q[ab];
        const double* exponential = e + ab * n;
        for (std::size_t i = 0; i < n; ++i) {
          row[i] += factor * exponential[i];
        }
      }
    }
    // Stage 2: radial factors T_L = sum_{t=1..deg} r_t R^(2(deg - t)) W(t) (Horner in R^2).
    for (std::size_t c = 0; c < couplings; ++c) {
      const auto first = static_cast<std::size_t>(form.radial_offsets[c]);
      const auto deg = static_cast<std::size_t>(form.radial_offsets[c + 1]) - first;
      const double* d = form.radial_coefficients.data() + first;  // d[t - 1] for t = 1..deg
      double* row = t + c * n;
      for (std::size_t i = 0; i < n; ++i) {
        double sum = d[0] * w[n + i];
        for (std::size_t k = 2; k <= deg; ++k) {
          sum = sum * r2[i] + d[k - 1] * w[k * n + i];
        }
        row[i] = sum;
      }
    }
    w0 = leading_term(form, r2, w, n, workspace.leading);
  }

  // Stage 3: angular assembly.
  assemble_dispatch(pair.l, pair.l_prime, harmonics, n, w0, t, values);
}

auto make_general_shell_pair(const BasisShell& bra, const BasisShell& ket,
                             TwoCentreIntegral integral) -> GeneralShellPair {
  const auto alphas = exponents(bra);
  const auto betas = exponents(ket);
  GeneralShellPair pair;
  pair.l = angular_momentum(bra);
  pair.l_prime = angular_momentum(ket);
  pair.form = make_form(pair.l, pair.l_prime, integral);
  pair.bra_primitives = alphas.size();
  pair.ket_primitives = betas.size();
  pair.bra_contractions = contraction_count(bra);
  pair.ket_contractions = contraction_count(ket);

  const std::size_t primitive_pairs = alphas.size() * betas.size();
  pair.mu.reserve(primitive_pairs);
  pair.factors.resize(pair.form.power_count * primitive_pairs);
  std::size_t index = 0;
  for (const double alpha : alphas) {
    for (const double beta : betas) {
      const double mu = alpha * beta / (alpha + beta);
      pair.mu.push_back(mu);
      set_power_factors(pair.form, primitive_factor(pair.l, pair.l_prime, alpha, beta), mu, index,
                        primitive_pairs, pair.factors);
      ++index;
    }
  }
  pair.bra_coefficients = transposed_coefficients(bra);
  pair.ket_coefficients = transposed_coefficients(ket);

  // Multiply-adds per power and separation of the two orders of the half-transformations.
  const std::size_t bra_first_cost =
      pair.bra_contractions * pair.ket_primitives * (pair.bra_primitives + pair.ket_contractions);
  const std::size_t ket_first_cost =
      pair.bra_primitives * pair.ket_contractions * (pair.ket_primitives + pair.bra_contractions);
  pair.bra_first = bra_first_cost <= ket_first_cost;
  return pair;
}

void general_shell_pair_values(const GeneralShellPair& pair, const SolidHarmonics& harmonics,
                               std::size_t n, KernelWorkspace& workspace,
                               std::span<double> values) {
  const IntegralForm& form = pair.form;
  const auto bra_components = static_cast<std::size_t>(2 * pair.l + 1);
  const auto ket_components = static_cast<std::size_t>(2 * pair.l_prime + 1);
  const std::size_t components = bra_components * ket_components;
  const std::size_t na = pair.bra_contractions;
  const std::size_t nb = pair.ket_contractions;
  const std::size_t ka = pair.bra_primitives;
  const std::size_t kb = pair.ket_primitives;
  const std::size_t contracted_pairs = na * nb;
  assert(values.size() >= contracted_pairs * components * n);
  if (n == 0) {
    return;
  }
  const auto r2 = harmonics.distances_squared();
  const std::size_t k_count = form.power_count;
  const std::size_t primitive_pairs = ka * kb;

  // Stage 1: exponentials and primitive factors w_ab,j exp(-mu_ab R^2), stored [a][j][b][i] when
  // the bra primitives are contracted first and [b][j][a][i] otherwise.
  workspace.exponentials.resize(primitive_pairs * n);
  workspace.weighted.resize(k_count * primitive_pairs * n);
  double* e = workspace.exponentials.data();
  double* weighted = workspace.weighted.data();
  for (std::size_t ab = 0; ab < primitive_pairs; ++ab) {
    exp_scaled_negative(r2, pair.mu[ab], std::span(e + ab * n, n));
  }
  for (std::size_t a = 0; a < ka; ++a) {
    for (std::size_t b = 0; b < kb; ++b) {
      const std::size_t ab = a * kb + b;
      const double* exponential = e + ab * n;
      for (std::size_t j = 0; j < k_count; ++j) {
        const std::size_t row =
            pair.bra_first ? (a * k_count + j) * kb + b : (b * k_count + j) * ka + a;
        const double factor = pair.factors[j * primitive_pairs + ab];
        double* out = weighted + row * n;
        for (std::size_t i = 0; i < n; ++i) {
          out[i] = factor * exponential[i];
        }
      }
    }
  }

  // Stage 2: two-index transformation W_j = c^T w_j d, stored [j][I][J][i]. The first
  // half-transformation covers all j in one product.
  workspace.accumulators.resize(k_count * contracted_pairs * n);
  double* w = workspace.accumulators.data();
  const double* c = pair.bra_coefficients.data();  // N_A x K_A
  const double* d = pair.ket_coefficients.data();  // N_B x K_B
  if (pair.bra_first) {
    const std::size_t width = k_count * kb * n;
    workspace.half.resize(na * width);
    double* half = workspace.half.data();  // [I][j][b][i]
    gemm(na, width, ka, c, ka, weighted, width, half, width);
    for (std::size_t bra = 0; bra < na; ++bra) {
      for (std::size_t j = 0; j < k_count; ++j) {
        gemm(nb, n, kb, d, kb, half + (bra * k_count + j) * kb * n, n, w + (j * na + bra) * nb * n,
             n);
      }
    }
  } else {
    const std::size_t width = k_count * ka * n;
    workspace.half.resize(nb * width);
    double* half = workspace.half.data();  // [J][j][a][i]
    gemm(nb, width, kb, d, kb, weighted, width, half, width);
    for (std::size_t ket = 0; ket < nb; ++ket) {
      for (std::size_t j = 0; j < k_count; ++j) {
        gemm(na, n, ka, c, ka, half + (ket * k_count + j) * ka * n, n,
             w + (j * contracted_pairs + ket) * n, nb * n);
      }
    }
  }

  // Stage 3: kernels Psi_t over n (one per accumulator row), only m <= m' for l == l':
  // the product P = S_lm S_l'm' enters row 0 (times w0_quadratic R^2) and row w0_constant_row
  // (times w0_constant), and each Gaunt term of T_L enters rows t = 1..deg.
  const bool mirrored = pair.l == pair.l_prime;
  workspace.kernels.resize(k_count * components * n);
  double* psi = workspace.kernels.data();
  std::fill(psi, psi + k_count * components * n, 0.0);
  for (int m = -pair.l; m <= pair.l; ++m) {
    const double* bra = harmonics.values(pair.l, m).data();
    for (int m_prime = mirrored ? m : -pair.l_prime; m_prime <= pair.l_prime; ++m_prime) {
      const double* ket = harmonics.values(pair.l_prime, m_prime).data();
      const std::size_t component = static_cast<std::size_t>(m + pair.l) * ket_components +
                                    static_cast<std::size_t>(m_prime + pair.l_prime);
      double* constant_row = psi + (form.w0_constant_row * components + component) * n;
      for (std::size_t i = 0; i < n; ++i) {
        constant_row[i] += form.w0_constant * (bra[i] * ket[i]);
      }
      if (form.w0_quadratic != 0.0) {
        double* row = psi + component * n;
        for (std::size_t i = 0; i < n; ++i) {
          row[i] += form.w0_quadratic * r2[i] * (bra[i] * ket[i]);
        }
      }
    }
  }
  // R^(2p) for p = 1..power_count - 2 (the largest power in a T_L is deg - 1 <= power_count - 2).
  const std::size_t power_count = k_count > 2 ? k_count - 2 : 0;
  workspace.powers.resize(power_count * n);
  double* powers = workspace.powers.data();
  for (std::size_t p = 0; p < power_count; ++p) {
    for (std::size_t i = 0; i < n; ++i) {
      powers[p * n + i] = p == 0 ? r2[i] : powers[(p - 1) * n + i] * r2[i];
    }
  }
  for (const AssemblyTerm& term : assembly_terms(pair.l, pair.l_prime)) {
    const auto first = static_cast<std::size_t>(form.radial_offsets[term.coupling]);
    const auto deg = static_cast<std::size_t>(form.radial_offsets[term.coupling + 1]) - first;
    const double* s = harmonics.values(term.big_l, term.big_m).data();
    for (std::size_t k = 1; k <= deg; ++k) {
      const double scale = term.value * form.radial_coefficients[first + k - 1];
      double* row = psi + (k * components + term.component) * n;
      if (k == deg) {
        for (std::size_t i = 0; i < n; ++i) {
          row[i] += scale * s[i];
        }
      } else {
        const double* power = powers + (deg - k - 1) * n;  // R^(2(deg - k))
        for (std::size_t i = 0; i < n; ++i) {
          row[i] += scale * power[i] * s[i];
        }
      }
    }
  }

  // Stage 4: values_(IJ),(mm') = sum_t W(t)_IJ Psi_t,mm'.
  for (std::size_t m = 0; m < bra_components; ++m) {
    for (std::size_t m_prime = mirrored ? m : 0; m_prime < ket_components; ++m_prime) {
      const std::size_t component = m * ket_components + m_prime;
      for (std::size_t ij = 0; ij < contracted_pairs; ++ij) {
        double* out = values.data() + (ij * components + component) * n;
        const double* w0 = w + ij * n;
        const double* psi0 = psi + component * n;
        for (std::size_t i = 0; i < n; ++i) {
          out[i] = w0[i] * psi0[i];
        }
        for (std::size_t j = 1; j < k_count; ++j) {
          const double* wj = w + (j * contracted_pairs + ij) * n;
          const double* psij = psi + (j * components + component) * n;
          for (std::size_t i = 0; i < n; ++i) {
            out[i] += wj[i] * psij[i];
          }
        }
      }
    }
  }
}

}  // namespace fika::detail
