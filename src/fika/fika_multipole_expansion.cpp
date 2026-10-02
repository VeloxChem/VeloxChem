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

#include "fika_multipole_expansion.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace fika::detail {

namespace {

using Complex = std::complex<double>;

// Complex product without the C99 Annex G NaN/infinity recovery of operator*.
inline auto times(Complex a, Complex b) noexcept -> Complex {
  return {a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real()};
}

// n! for n <= 2 max_expansion_order.
const std::array<double, 2 * max_expansion_order + 1> factorials = [] {
  std::array<double, 2 * max_expansion_order + 1> result{};
  result[0] = 1.0;
  for (std::size_t n = 1; n < result.size(); ++n) {
    result[n] = result[n - 1] * static_cast<double>(n);
  }
  return result;
}();

auto difference(const Point3D<double>& a, const Point3D<double>& b) noexcept -> Point3D<double> {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// All components m = -l..l of an expansion (order p), real and imaginary parts apart, at index
// l^2 + l + m (reversed: l^2 + l - m): the operators sum over contiguous, branch-free ranges of m.
void expand_components(std::span<const Complex> values, int order, std::vector<double>& real,
                       std::vector<double>& imaginary, bool reversed = false) {
  const auto size = static_cast<std::size_t>((order + 1) * (order + 1));
  real.resize(size);
  imaginary.resize(size);
  for (int l = 0; l <= order; ++l) {
    const auto centre = static_cast<std::size_t>(l * l + l);
    for (int m = 0; m <= l; ++m) {
      const Complex value = values[expansion_index(l, m)];
      const std::size_t plus =
          reversed ? centre - static_cast<std::size_t>(m) : centre + static_cast<std::size_t>(m);
      const std::size_t minus =
          reversed ? centre + static_cast<std::size_t>(m) : centre - static_cast<std::size_t>(m);
      real[plus] = value.real();
      imaginary[plus] = value.imag();
      const double sign = (m & 1) != 0 ? -1.0 : 1.0;  // X_(l,-m) = (-1)^m conj(X_lm)
      real[minus] = sign * value.real();
      imaginary[minus] = -sign * value.imag();
    }
  }
}

}  // namespace

void regular_harmonics(const Point3D<double>& r, int order, std::span<Complex> values) {
  assert(order >= 0 && order <= max_expansion_order);
  assert(values.size() >= expansion_size(order));
  const double r2 = r.x * r.x + r.y * r.y + r.z * r.z;
  const Complex xy{r.x, r.y};
  values[0] = 1.0;
  for (int l = 0; l < order; ++l) {
    // R_(l+1),(l+1) = -(x + iy) / (2l + 2) R_ll;
    // R_(l+1),m = ((2l + 1) z R_lm - r^2 R_(l-1),m) / ((l + m + 1)(l - m + 1)).
    values[expansion_index(l + 1, l + 1)] =
        times(xy, values[expansion_index(l, l)]) * (-1.0 / (2 * l + 2));
    for (int m = 0; m <= l; ++m) {
      Complex value = values[expansion_index(l, m)] * ((2 * l + 1) * r.z);
      if (m < l) {
        value -= values[expansion_index(l - 1, m)] * r2;
      }
      values[expansion_index(l + 1, m)] = value / static_cast<double>((l + m + 1) * (l - m + 1));
    }
  }
}

void irregular_harmonics(const Point3D<double>& r, int order, std::span<Complex> values) {
  assert(order >= 0 && order <= max_expansion_order);
  assert(values.size() >= expansion_size(order));
  const double r2 = r.x * r.x + r.y * r.y + r.z * r.z;
  assert(r2 > 0.0);
  const double inverse_r2 = 1.0 / r2;
  const Complex xy{r.x * inverse_r2, r.y * inverse_r2};
  const double z = r.z * inverse_r2;
  values[0] = 1.0 / std::sqrt(r2);
  for (int l = 0; l < order; ++l) {
    // I_(l+1),(l+1) = -(2l + 1)(x + iy) / r^2 I_ll;
    // I_(l+1),m = ((2l + 1) z I_lm - (l^2 - m^2) I_(l-1),m) / r^2.
    values[expansion_index(l + 1, l + 1)] =
        times(xy, values[expansion_index(l, l)]) * static_cast<double>(-(2 * l + 1));
    for (int m = 0; m <= l; ++m) {
      Complex value = values[expansion_index(l, m)] * ((2 * l + 1) * z);
      if (m < l) {
        value -=
            values[expansion_index(l - 1, m)] * (static_cast<double>(l * l - m * m) * inverse_r2);
      }
      values[expansion_index(l + 1, m)] = value;
    }
  }
}

namespace {

// Lanes of the batched source operators.
constexpr std::size_t batch = 8;

// Batched recursion of regular_harmonics over `batch` lanes (layout: coefficient major, lanes
// innermost), from the values R_00 already in re[0..batch), im[0..batch). The arrays do not
// overlap (__restrict lets the lane loops vectorize).
[[gnu::always_inline]] inline void regular_recursion(int order, const double* __restrict x,
                                                     const double* __restrict y,
                                                     const double* __restrict z,
                                                     const double* __restrict r2,
                                                     double* __restrict re, double* __restrict im) {
  for (int l = 0; l < order; ++l) {
    const double* top_re = re + expansion_index(l, l) * batch;
    const double* top_im = im + expansion_index(l, l) * batch;
    double* next_re = re + expansion_index(l + 1, l + 1) * batch;
    double* next_im = im + expansion_index(l + 1, l + 1) * batch;
    const double factor = -1.0 / (2 * l + 2);
#pragma omp simd
    for (std::size_t b = 0; b < batch; ++b) {  // R_(l+1),(l+1) = -(x + iy)/(2l + 2) R_ll
      next_re[b] = factor * (x[b] * top_re[b] - y[b] * top_im[b]);
      next_im[b] = factor * (x[b] * top_im[b] + y[b] * top_re[b]);
    }
    for (int m = 0; m <= l; ++m) {
      const double* current_re = re + expansion_index(l, m) * batch;
      const double* current_im = im + expansion_index(l, m) * batch;
      double* out_re = re + expansion_index(l + 1, m) * batch;
      double* out_im = im + expansion_index(l + 1, m) * batch;
      const auto a = static_cast<double>(2 * l + 1);
      const double inverse = 1.0 / static_cast<double>((l + m + 1) * (l - m + 1));
      if (m < l) {
        const double* previous_re = re + expansion_index(l - 1, m) * batch;
        const double* previous_im = im + expansion_index(l - 1, m) * batch;
#pragma omp simd
        for (std::size_t b = 0; b < batch; ++b) {
          out_re[b] = (a * z[b] * current_re[b] - r2[b] * previous_re[b]) * inverse;
          out_im[b] = (a * z[b] * current_im[b] - r2[b] * previous_im[b]) * inverse;
        }
      } else {
#pragma omp simd
        for (std::size_t b = 0; b < batch; ++b) {
          out_re[b] = a * z[b] * current_re[b] * inverse;
          out_im[b] = a * z[b] * current_im[b] * inverse;
        }
      }
    }
  }
}

// Batched recursion of irregular_harmonics (x, y, z already divided by r^2), from I_00 in
// re[0..batch), im[0..batch).
[[gnu::always_inline]] inline void irregular_recursion(
    int order, const double* __restrict x, const double* __restrict y, const double* __restrict z,
    const double* __restrict inverse_r2, double* __restrict re, double* __restrict im) {
  for (int l = 0; l < order; ++l) {
    double* top_re = re + expansion_index(l, l) * batch;
    double* top_im = im + expansion_index(l, l) * batch;
    double* next_re = re + expansion_index(l + 1, l + 1) * batch;
    double* next_im = im + expansion_index(l + 1, l + 1) * batch;
    const auto factor = static_cast<double>(-(2 * l + 1));
#pragma omp simd
    for (std::size_t b = 0; b < batch; ++b) {  // I_(l+1),(l+1) = -(2l + 1)(x + iy)/r^2 I_ll
      next_re[b] = factor * (x[b] * top_re[b] - y[b] * top_im[b]);
      next_im[b] = factor * (x[b] * top_im[b] + y[b] * top_re[b]);
    }
    for (int m = 0; m <= l; ++m) {
      const double* current_re = re + expansion_index(l, m) * batch;
      const double* current_im = im + expansion_index(l, m) * batch;
      double* out_re = re + expansion_index(l + 1, m) * batch;
      double* out_im = im + expansion_index(l + 1, m) * batch;
      const auto a = static_cast<double>(2 * l + 1);
      if (m < l) {
        const double* previous_re = re + expansion_index(l - 1, m) * batch;
        const double* previous_im = im + expansion_index(l - 1, m) * batch;
        const auto c = static_cast<double>(l * l - m * m);
#pragma omp simd
        for (std::size_t b = 0; b < batch; ++b) {
          out_re[b] = a * z[b] * current_re[b] - c * inverse_r2[b] * previous_re[b];
          out_im[b] = a * z[b] * current_im[b] - c * inverse_r2[b] * previous_im[b];
        }
      } else {
#pragma omp simd
        for (std::size_t b = 0; b < batch; ++b) {
          out_re[b] = a * z[b] * current_re[b];
          out_im[b] = a * z[b] * current_im[b];
        }
      }
    }
  }
}

}  // namespace

void add_charges_to_multipole(std::span<const double> charges,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order,
                              std::span<Complex> multipole, ExpansionWorkspace& workspace) {
  assert(charges.size() == coordinates.size());
  const std::size_t size = expansion_size(order);
  assert(multipole.size() >= size);
  // Batches of charges as in add_charges_to_local: q R_lm of every charge by the recursion of
  // regular_harmonics, the charges innermost (vectorized), then conjugated and summed per
  // coefficient. Lanes past the last charge hold a zero charge; fewer charges than a batch take
  // the scalar path.
  if (charges.size() < batch) {
    workspace.harmonics.resize(size);
    for (std::size_t c = 0; c < charges.size(); ++c) {
      regular_harmonics(difference(coordinates[c], centre), order, workspace.harmonics);
      for (std::size_t i = 0; i < size; ++i) {
        multipole[i] += charges[c] * std::conj(workspace.harmonics[i]);
      }
    }
    return;
  }
  workspace.real.resize(size * batch);
  workspace.imaginary.resize(size * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], r2[batch];
  for (std::size_t first = 0; first < charges.size(); first += batch) {
    const std::size_t count = std::min(batch, charges.size() - first);
    for (std::size_t b = 0; b < batch; ++b) {
      const bool used = b < count;
      const Point3D<double> r =
          used ? difference(coordinates[first + b], centre) : Point3D<double>{0.0, 0.0, 0.0};
      x[b] = r.x;
      y[b] = r.y;
      z[b] = r.z;
      r2[b] = r.x * r.x + r.y * r.y + r.z * r.z;
      re[b] = used ? charges[first + b] : 0.0;
      im[b] = 0.0;
    }
    regular_recursion(order, x, y, z, r2, re, im);
    for (std::size_t i = 0; i < size; ++i) {
      double sum_re = 0.0;
      double sum_im = 0.0;
#pragma omp simd reduction(+ : sum_re, sum_im)
      for (std::size_t b = 0; b < batch; ++b) {
        sum_re += re[i * batch + b];
        sum_im += im[i * batch + b];
      }
      multipole[i] += Complex{sum_re, -sum_im};  // conj
    }
  }
}

namespace {

// M'_lm += sum_jk conj(R_jk(O - O')) M_(l-j),(m-k) with the shift harmonics R expanded to all
// components; the source is expanded in reversed m order, so both operands of the sum over k
// are contiguous forward ranges (SIMD reduction).
void translate_multipole(std::span<const Complex> source, const double* r_re, const double* r_im,
                         int order, std::span<Complex> target, ExpansionWorkspace& workspace) {
  expand_components(source, order, workspace.real, workspace.imaginary, true);
  const double* m_re = workspace.real.data();
  const double* m_im = workspace.imaginary.data();
  for (int l = 0; l <= order; ++l) {
    for (int m = 0; m <= l; ++m) {
      double sum_re = 0.0;
      double sum_im = 0.0;
      for (int j = 0; j <= l; ++j) {
        const int k_low = std::max(-j, m - (l - j));
        const int k_high = std::min(j, m + (l - j));
        const double* b_re = r_re + j * j + j;  // R_jk at k
        const double* b_im = r_im + j * j + j;
        const int n = l - j;
        const double* a_re = m_re + n * n + n - m;  // M_n,(m-k) at k (reversed)
        const double* a_im = m_im + n * n + n - m;
#pragma omp simd reduction(+ : sum_re, sum_im)
        for (int k = k_low; k <= k_high; ++k) {  // a conj(b)
          sum_re += a_re[k] * b_re[k] + a_im[k] * b_im[k];
          sum_im += a_im[k] * b_re[k] - a_re[k] * b_im[k];
        }
      }
      target[expansion_index(l, m)] += Complex{sum_re, sum_im};
    }
  }
}

}  // namespace

auto make_multipole_shift(const Point3D<double>& from, const Point3D<double>& to, int order)
    -> MultipoleShift {
  MultipoleShift shift;
  shift.order = order;
  std::vector<Complex> harmonics(expansion_size(order));
  regular_harmonics(difference(from, to), order, harmonics);
  expand_components(harmonics, order, shift.real, shift.imaginary);
  return shift;
}

void translate_multipole(std::span<const Complex> source, const MultipoleShift& shift,
                         std::span<Complex> target, ExpansionWorkspace& workspace) {
  translate_multipole(source, shift.real.data(), shift.imaginary.data(), shift.order, target,
                      workspace);
}

void translate_multipole(std::span<const Complex> source, const Point3D<double>& from,
                         const Point3D<double>& to, int order, std::span<Complex> target,
                         ExpansionWorkspace& workspace) {
  workspace.harmonics.resize(expansion_size(order));
  regular_harmonics(difference(from, to), order, workspace.harmonics);
  expand_components(workspace.harmonics, order, workspace.kernel_real, workspace.kernel_imaginary);
  translate_multipole(source, workspace.kernel_real.data(), workspace.kernel_imaginary.data(),
                      order, target, workspace);
}

void multipole_to_local(std::span<const Complex> multipole, const Point3D<double>& from,
                        const Point3D<double>& to, int order, std::span<Complex> local,
                        ExpansionWorkspace& workspace) {
  // L_lm = (-1)^l sum_(j + l <= p) sum_k M_jk I_(j+l),(k+m)(Q - O).
  workspace.harmonics.resize(expansion_size(order));
  irregular_harmonics(difference(to, from), order, workspace.harmonics);
  expand_components(multipole, order, workspace.real, workspace.imaginary);
  expand_components(workspace.harmonics, order, workspace.kernel_real, workspace.kernel_imaginary);
  const double* m_re = workspace.real.data();
  const double* m_im = workspace.imaginary.data();
  const double* i_re = workspace.kernel_real.data();
  const double* i_im = workspace.kernel_imaginary.data();
  // The sums over k vectorize as OpenMP SIMD reductions (fixed order for a given build).
  for (int l = 0; l <= order; ++l) {
    for (int m = 0; m <= l; ++m) {
      double sum_re = 0.0;
      double sum_im = 0.0;
      for (int j = 0; j <= order - l; ++j) {
        const double* a_re = m_re + j * j;  // k = -j..j
        const double* a_im = m_im + j * j;
        const int kernel_start = (j + l) * (j + l) + (j + l) + m - j;
        const double* b_re = i_re + kernel_start;
        const double* b_im = i_im + kernel_start;
#pragma omp simd reduction(+ : sum_re, sum_im)
        for (int k = 0; k <= 2 * j; ++k) {
          sum_re += a_re[k] * b_re[k] - a_im[k] * b_im[k];
          sum_im += a_re[k] * b_im[k] + a_im[k] * b_re[k];
        }
      }
      const double sign = (l & 1) != 0 ? -1.0 : 1.0;
      local[expansion_index(l, m)] += Complex{sign * sum_re, sign * sum_im};
    }
  }
}

void add_real_multipole_to_multipole(std::span<const double> moments, int rank,
                                     const Point3D<double>& position, const Point3D<double>& centre,
                                     int order, std::span<Complex> multipole,
                                     ExpansionWorkspace& workspace) {
  assert(rank >= 0 && order >= 0 && order <= max_expansion_order);
  assert(moments.size() >= static_cast<std::size_t>((rank + 1) * (rank + 1)));
  assert(multipole.size() >= expansion_size(order));
  const int top = std::min(rank, order);
  // Complex coefficients about the position, all components m = -n..n (index n^2 + n + m).
  workspace.shifted.assign(expansion_size(top), 0.0);
  for (int n = 0; n <= top; ++n) {
    const auto row = static_cast<std::size_t>(n * n + n);
    workspace.shifted[expansion_index(n, 0)] =
        moments[row] / factorials[static_cast<std::size_t>(n)];
    for (int m = 1; m <= n; ++m) {
      const double scale =
          ((m & 1) != 0 ? -1.0 : 1.0) /
          (std::sqrt(2.0) * std::sqrt(factorials[static_cast<std::size_t>(n - m)] *
                                      factorials[static_cast<std::size_t>(n + m)]));
      workspace.shifted[expansion_index(n, m)] =
          Complex{scale * moments[row + static_cast<std::size_t>(m)],
                  -scale * moments[row - static_cast<std::size_t>(m)]};
    }
  }
  expand_components(workspace.shifted, top, workspace.real, workspace.imaginary);
  // M_lm(centre) += sum_(n <= k) sum_k' M_nk' conj(R_(l-n),(m-k')(position - centre)).
  workspace.harmonics.resize(expansion_size(order));
  regular_harmonics(difference(position, centre), order, workspace.harmonics);
  expand_components(workspace.harmonics, order, workspace.kernel_real, workspace.kernel_imaginary);
  const double* m_re = workspace.real.data();
  const double* m_im = workspace.imaginary.data();
  const double* r_re = workspace.kernel_real.data();
  const double* r_im = workspace.kernel_imaginary.data();
  for (int l = 0; l <= order; ++l) {
    for (int m = 0; m <= l; ++m) {
      double sum_re = 0.0;
      double sum_im = 0.0;
      for (int n = 0; n <= std::min(top, l); ++n) {
        const int j = l - n;
        for (int k = std::max(-n, m - j); k <= std::min(n, m + j); ++k) {
          const double a_re = m_re[n * n + n + k];
          const double a_im = m_im[n * n + n + k];
          const double b_re = r_re[j * j + j + m - k];
          const double b_im = r_im[j * j + j + m - k];
          sum_re += a_re * b_re + a_im * b_im;  // a conj(b)
          sum_im += a_im * b_re - a_re * b_im;
        }
      }
      multipole[expansion_index(l, m)] += Complex{sum_re, sum_im};
    }
  }
}

void add_charges_to_local(std::span<const double> charges,
                          std::span<const Point3D<double>> coordinates,
                          const Point3D<double>& centre, int order, std::span<Complex> local,
                          ExpansionWorkspace& workspace) {
  assert(charges.size() == coordinates.size());
  const std::size_t size = expansion_size(order);
  assert(local.size() >= size);
  // Batches of charges: q I_lm of every charge by the recursion of irregular_harmonics, the
  // charges innermost (vectorized), then summed per coefficient. Lanes past the last charge
  // hold a zero charge at unit distance; fewer charges than a batch take the scalar path.
  if (charges.size() < batch) {
    workspace.harmonics.resize(size);
    for (std::size_t c = 0; c < charges.size(); ++c) {
      irregular_harmonics(difference(coordinates[c], centre), order, workspace.harmonics);
      for (std::size_t i = 0; i < size; ++i) {
        local[i] += charges[c] * workspace.harmonics[i];
      }
    }
    return;
  }
  workspace.real.resize(size * batch);
  workspace.imaginary.resize(size * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], inverse_r2[batch];
  for (std::size_t first = 0; first < charges.size(); first += batch) {
    const std::size_t count = std::min(batch, charges.size() - first);
    for (std::size_t b = 0; b < batch; ++b) {
      const bool used = b < count;
      const Point3D<double> r =
          used ? difference(coordinates[first + b], centre) : Point3D<double>{1.0, 0.0, 0.0};
      const double r2 = r.x * r.x + r.y * r.y + r.z * r.z;
      assert(r2 > 0.0);
      inverse_r2[b] = 1.0 / r2;
      x[b] = r.x * inverse_r2[b];
      y[b] = r.y * inverse_r2[b];
      z[b] = r.z * inverse_r2[b];
      re[b] = used ? charges[first + b] / std::sqrt(r2) : 0.0;
      im[b] = 0.0;
    }
    irregular_recursion(order, x, y, z, inverse_r2, re, im);
    for (std::size_t i = 0; i < size; ++i) {
      double sum_re = 0.0;
      double sum_im = 0.0;
#pragma omp simd reduction(+ : sum_re, sum_im)
      for (std::size_t b = 0; b < batch; ++b) {
        sum_re += re[i * batch + b];
        sum_im += im[i * batch + b];
      }
      local[i] += Complex{sum_re, sum_im};
    }
  }
}

namespace {

// L'_jk = sum_in L_(j+i),(k+n) conj(R_in(Q' - Q)) for j <= rank: the local expansion about Q'
// (rows j <= rank of `target`, accumulated).
void shift_local(std::span<const Complex> source, const Point3D<double>& from,
                 const Point3D<double>& to, int order, int rank, std::span<Complex> target,
                 ExpansionWorkspace& workspace) {
  workspace.harmonics.resize(expansion_size(order));
  regular_harmonics(difference(to, from), order, workspace.harmonics);
  expand_components(source, order, workspace.real, workspace.imaginary);
  expand_components(workspace.harmonics, order, workspace.kernel_real, workspace.kernel_imaginary);
  const double* l_re = workspace.real.data();
  const double* l_im = workspace.imaginary.data();
  const double* r_re = workspace.kernel_real.data();
  const double* r_im = workspace.kernel_imaginary.data();
  for (int j = 0; j <= rank; ++j) {
    for (int k = 0; k <= j; ++k) {
      double sum_re = 0.0;
      double sum_im = 0.0;
      for (int i = 0; i <= order - j; ++i) {
        const int source_start = (j + i) * (j + i) + (j + i) + k - i;  // n = -i..i
        const double* a_re = l_re + source_start;
        const double* a_im = l_im + source_start;
        const double* b_re = r_re + i * i;
        const double* b_im = r_im + i * i;
#pragma omp simd reduction(+ : sum_re, sum_im)
        for (int n = 0; n <= 2 * i; ++n) {  // a conj(b)
          sum_re += a_re[n] * b_re[n] + a_im[n] * b_im[n];
          sum_im += a_im[n] * b_re[n] - a_re[n] * b_im[n];
        }
      }
      target[expansion_index(j, k)] += Complex{sum_re, sum_im};
    }
  }
}

}  // namespace

void translate_local(std::span<const Complex> source, const Point3D<double>& from,
                     const Point3D<double>& to, int order, std::span<Complex> target,
                     ExpansionWorkspace& workspace) {
  shift_local(source, from, to, order, order, target, workspace);
}

void local_field_tensor(std::span<const Complex> local, const Point3D<double>& centre, int order,
                        const Point3D<double>& point, int rank, std::span<double> phi,
                        ExpansionWorkspace& workspace) {
  assert(rank >= 0 && rank <= order);
  assert(phi.size() >= static_cast<std::size_t>((rank + 1) * (rank + 1)));
  workspace.shifted.assign(expansion_size(rank), 0.0);
  shift_local(local, centre, point, order, rank, workspace.shifted, workspace);
  // The local coefficients at x are sum q I_Lambda,M(s - x); in real form
  // Phi_Lambda,0 = Re / Lambda!, Phi_Lambda,(+-M) = (-1)^M sqrt(2) (Re, Im) / sqrt((L-M)! (L+M)!).
  for (int big_l = 0; big_l <= rank; ++big_l) {
    const std::size_t row = static_cast<std::size_t>(big_l * big_l + big_l);
    phi[row] = workspace.shifted[expansion_index(big_l, 0)].real() /
               factorials[static_cast<std::size_t>(big_l)];
    for (int m = 1; m <= big_l; ++m) {
      const double scale = ((m & 1) != 0 ? -std::sqrt(2.0) : std::sqrt(2.0)) /
                           std::sqrt(factorials[static_cast<std::size_t>(big_l - m)] *
                                     factorials[static_cast<std::size_t>(big_l + m)]);
      const Complex value = workspace.shifted[expansion_index(big_l, m)];
      phi[row + static_cast<std::size_t>(m)] = scale * value.real();
      phi[row - static_cast<std::size_t>(m)] = scale * value.imag();
    }
  }
}

namespace {

// X_(l, m) of a batched expansion (lanes innermost) for |m| <= l, from the stored m >= 0:
// X_(l,-m) = (-1)^m conj(X_lm).
struct BatchedValue {
  double re;
  double im;
};

inline auto batched(const double* re, const double* im, int l, int m, std::size_t b) noexcept
    -> BatchedValue {
  if (m >= 0) {
    const std::size_t i = expansion_index(l, m) * batch + b;
    return {re[i], im[i]};
  }
  const std::size_t i = expansion_index(l, -m) * batch + b;
  const double sign = (m & 1) != 0 ? -1.0 : 1.0;
  return {sign * re[i], -sign * im[i]};
}

// Loads a batch of dipoles: separations r = s - centre (unit x for unused lanes, zero dipole).
void load_dipoles(std::span<const Dipole> dipoles, std::span<const Point3D<double>> coordinates,
                  const Point3D<double>& centre, std::size_t first, double* x, double* y, double* z,
                  double* r2, double* mu_x, double* mu_y, double* mu_z) {
  const std::size_t count = std::min(batch, dipoles.size() - first);
  for (std::size_t b = 0; b < batch; ++b) {
    const bool used = b < count;
    const Point3D<double> r =
        used ? difference(coordinates[first + b], centre) : Point3D<double>{1.0, 0.0, 0.0};
    x[b] = r.x;
    y[b] = r.y;
    z[b] = r.z;
    r2[b] = r.x * r.x + r.y * r.y + r.z * r.z;
    mu_x[b] = used ? dipoles[first + b](0) : 0.0;
    mu_y[b] = used ? dipoles[first + b](1) : 0.0;
    mu_z[b] = used ? dipoles[first + b](2) : 0.0;
  }
}

}  // namespace

void add_dipoles_to_multipole(std::span<const Dipole> dipoles,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order,
                              std::span<Complex> multipole, ExpansionWorkspace& workspace) {
  assert(dipoles.size() == coordinates.size());
  assert(multipole.size() >= expansion_size(order));
  if (dipoles.empty() || order == 0) {
    return;  // a dipole has no monopole
  }
  // M_lm += mu . grad_s conj(R_lm(s - O)) = conj(mu_z R_(l-1),m + (mu_x - i mu_y) / 2
  // R_(l-1),(m+1) - (mu_x + i mu_y) / 2 R_(l-1),(m-1)): regular harmonics of order p - 1.
  const int lower = order - 1;
  workspace.real.resize(expansion_size(lower) * batch);
  workspace.imaginary.resize(expansion_size(lower) * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], r2[batch], mu_x[batch], mu_y[batch], mu_z[batch];
  for (std::size_t first = 0; first < dipoles.size(); first += batch) {
    load_dipoles(dipoles, coordinates, centre, first, x, y, z, r2, mu_x, mu_y, mu_z);
    for (std::size_t b = 0; b < batch; ++b) {
      re[b] = 1.0;
      im[b] = 0.0;
    }
    regular_recursion(lower, x, y, z, r2, re, im);
    for (int l = 1; l <= order; ++l) {
      for (int m = 0; m <= l; ++m) {
        double sum_re = 0.0;
        double sum_im = 0.0;
        for (std::size_t b = 0; b < batch; ++b) {
          double value_re = 0.0;
          double value_im = 0.0;
          if (m <= l - 1) {
            const auto r = batched(re, im, l - 1, m, b);
            value_re += mu_z[b] * r.re;
            value_im += mu_z[b] * r.im;
          }
          if (m + 1 <= l - 1) {  // (mu_x - i mu_y) / 2 R_(l-1),(m+1)
            const auto r = batched(re, im, l - 1, m + 1, b);
            value_re += 0.5 * (mu_x[b] * r.re + mu_y[b] * r.im);
            value_im += 0.5 * (mu_x[b] * r.im - mu_y[b] * r.re);
          }
          if (std::abs(m - 1) <= l - 1) {  // -(mu_x + i mu_y) / 2 R_(l-1),(m-1)
            const auto r = batched(re, im, l - 1, m - 1, b);
            value_re -= 0.5 * (mu_x[b] * r.re - mu_y[b] * r.im);
            value_im -= 0.5 * (mu_x[b] * r.im + mu_y[b] * r.re);
          }
          sum_re += value_re;
          sum_im += value_im;
        }
        multipole[expansion_index(l, m)] += Complex{sum_re, -sum_im};  // conj
      }
    }
  }
}

void add_dipoles_to_local(std::span<const Dipole> dipoles,
                          std::span<const Point3D<double>> coordinates,
                          const Point3D<double>& centre, int order, std::span<Complex> local,
                          ExpansionWorkspace& workspace) {
  assert(dipoles.size() == coordinates.size());
  assert(local.size() >= expansion_size(order));
  assert(order + 1 <= max_expansion_order);
  if (dipoles.empty()) {
    return;
  }
  // L_lm += mu . grad_s I_lm(s - Q) = -mu_z I_(l+1),m + (mu_x - i mu_y) / 2 I_(l+1),(m+1)
  // - (mu_x + i mu_y) / 2 I_(l+1),(m-1): irregular harmonics of order p + 1.
  const int upper = order + 1;
  workspace.real.resize(expansion_size(upper) * batch);
  workspace.imaginary.resize(expansion_size(upper) * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], r2[batch], inverse_r2[batch], mu_x[batch],
      mu_y[batch], mu_z[batch];
  for (std::size_t first = 0; first < dipoles.size(); first += batch) {
    load_dipoles(dipoles, coordinates, centre, first, x, y, z, r2, mu_x, mu_y, mu_z);
    for (std::size_t b = 0; b < batch; ++b) {
      assert(r2[b] > 0.0);
      inverse_r2[b] = 1.0 / r2[b];
      re[b] = 1.0 / std::sqrt(r2[b]);
      im[b] = 0.0;
      x[b] *= inverse_r2[b];
      y[b] *= inverse_r2[b];
      z[b] *= inverse_r2[b];
    }
    irregular_recursion(upper, x, y, z, inverse_r2, re, im);
    for (int l = 0; l <= order; ++l) {
      for (int m = 0; m <= l; ++m) {
        double sum_re = 0.0;
        double sum_im = 0.0;
        for (std::size_t b = 0; b < batch; ++b) {
          const auto centre_term = batched(re, im, l + 1, m, b);
          const auto raised = batched(re, im, l + 1, m + 1, b);
          const auto lowered = batched(re, im, l + 1, m - 1, b);
          sum_re += -mu_z[b] * centre_term.re + 0.5 * (mu_x[b] * raised.re + mu_y[b] * raised.im) -
                    0.5 * (mu_x[b] * lowered.re - mu_y[b] * lowered.im);
          sum_im += -mu_z[b] * centre_term.im + 0.5 * (mu_x[b] * raised.im - mu_y[b] * raised.re) -
                    0.5 * (mu_x[b] * lowered.im + mu_y[b] * lowered.re);
        }
        local[expansion_index(l, m)] += Complex{sum_re, sum_im};
      }
    }
  }
}

namespace {

// Loads a batch of quadrupoles: separations r = s - centre (unit x for unused lanes, zero
// quadrupole) and the coefficients of 1/2 Theta : grad grad = c++ D+^2 + conj(c++) D-^2 +
// c+z D+ d_z + conj(c+z) D- d_z + czz d_z^2 on harmonic functions (D+- = d_x +- i d_y):
// c++ = (Q_xx - Q_yy) / 8 - i Q_xy / 4, c+z = (Q_xz - i Q_yz) / 2,
// czz = (Q_zz - (Q_xx + Q_yy) / 2) / 2 (the trace of Q cancels).
void load_quadrupoles(std::span<const Quadrupole> quadrupoles,
                      std::span<const Point3D<double>> coordinates, const Point3D<double>& centre,
                      std::size_t first, double* x, double* y, double* z, double* r2,
                      double* cpp_re, double* cpp_im, double* cpz_re, double* cpz_im, double* czz) {
  const std::size_t count = std::min(batch, quadrupoles.size() - first);
  for (std::size_t b = 0; b < batch; ++b) {
    const bool used = b < count;
    const Point3D<double> r =
        used ? difference(coordinates[first + b], centre) : Point3D<double>{1.0, 0.0, 0.0};
    x[b] = r.x;
    y[b] = r.y;
    z[b] = r.z;
    r2[b] = r.x * r.x + r.y * r.y + r.z * r.z;
    const Quadrupole q = used ? quadrupoles[first + b] : Quadrupole{};
    cpp_re[b] = (q(0, 0) - q(1, 1)) / 8.0;
    cpp_im[b] = -q(0, 1) / 4.0;
    cpz_re[b] = q(0, 2) / 2.0;
    cpz_im[b] = -q(1, 2) / 2.0;
    czz[b] = (q(2, 2) - 0.5 * (q(0, 0) + q(1, 1))) / 2.0;
  }
}

// a b for complex numbers given as (re, im).
inline auto multiply(double a_re, double a_im, BatchedValue b) noexcept -> BatchedValue {
  return {a_re * b.re - a_im * b.im, a_re * b.im + a_im * b.re};
}

}  // namespace

void add_quadrupoles_to_multipole(std::span<const Quadrupole> quadrupoles,
                                  std::span<const Point3D<double>> coordinates,
                                  const Point3D<double>& centre, int order,
                                  std::span<Complex> multipole, ExpansionWorkspace& workspace) {
  assert(quadrupoles.size() == coordinates.size());
  assert(multipole.size() >= expansion_size(order));
  if (quadrupoles.empty() || order < 2) {
    return;  // a quadrupole has no monopole or dipole moment
  }
  // M_lm += conj(1/2 Theta : grad grad R_lm(s - O)) = conj(c++ R_(l-2),(m+2) + conj(c++)
  // R_(l-2),(m-2) + c+z R_(l-2),(m+1) - conj(c+z) R_(l-2),(m-1) + czz R_(l-2),m): regular
  // harmonics of order p - 2.
  const int lower = order - 2;
  workspace.real.resize(expansion_size(lower) * batch);
  workspace.imaginary.resize(expansion_size(lower) * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], r2[batch], cpp_re[batch], cpp_im[batch],
      cpz_re[batch], cpz_im[batch], czz[batch];
  for (std::size_t first = 0; first < quadrupoles.size(); first += batch) {
    load_quadrupoles(quadrupoles, coordinates, centre, first, x, y, z, r2, cpp_re, cpp_im, cpz_re,
                     cpz_im, czz);
    for (std::size_t b = 0; b < batch; ++b) {
      re[b] = 1.0;
      im[b] = 0.0;
    }
    regular_recursion(lower, x, y, z, r2, re, im);
    for (int l = 2; l <= order; ++l) {
      const int n = l - 2;
      for (int m = 0; m <= l; ++m) {
        double sum_re = 0.0;
        double sum_im = 0.0;
        for (std::size_t b = 0; b < batch; ++b) {
          const auto add = [&](BatchedValue v, double sign) {
            sum_re += sign * v.re;
            sum_im += sign * v.im;
          };
          if (std::abs(m + 2) <= n) {
            add(multiply(cpp_re[b], cpp_im[b], batched(re, im, n, m + 2, b)), 1.0);
          }
          if (std::abs(m - 2) <= n) {
            add(multiply(cpp_re[b], -cpp_im[b], batched(re, im, n, m - 2, b)), 1.0);
          }
          if (std::abs(m + 1) <= n) {
            add(multiply(cpz_re[b], cpz_im[b], batched(re, im, n, m + 1, b)), 1.0);
          }
          if (std::abs(m - 1) <= n) {
            add(multiply(cpz_re[b], -cpz_im[b], batched(re, im, n, m - 1, b)), -1.0);
          }
          if (m <= n) {
            add(multiply(czz[b], 0.0, batched(re, im, n, m, b)), 1.0);
          }
        }
        multipole[expansion_index(l, m)] += Complex{sum_re, -sum_im};  // conj
      }
    }
  }
}

void add_quadrupoles_to_local(std::span<const Quadrupole> quadrupoles,
                              std::span<const Point3D<double>> coordinates,
                              const Point3D<double>& centre, int order, std::span<Complex> local,
                              ExpansionWorkspace& workspace) {
  assert(quadrupoles.size() == coordinates.size());
  assert(local.size() >= expansion_size(order));
  assert(order + 2 <= max_expansion_order);
  if (quadrupoles.empty()) {
    return;
  }
  // L_lm += 1/2 Theta : grad grad I_lm(s - Q) = c++ I_(l+2),(m+2) + conj(c++) I_(l+2),(m-2)
  // - c+z I_(l+2),(m+1) + conj(c+z) I_(l+2),(m-1) + czz I_(l+2),m: irregular harmonics of
  // order p + 2.
  const int upper = order + 2;
  workspace.real.resize(expansion_size(upper) * batch);
  workspace.imaginary.resize(expansion_size(upper) * batch);
  double* re = workspace.real.data();
  double* im = workspace.imaginary.data();
  alignas(64) double x[batch], y[batch], z[batch], r2[batch], inverse_r2[batch], cpp_re[batch],
      cpp_im[batch], cpz_re[batch], cpz_im[batch], czz[batch];
  for (std::size_t first = 0; first < quadrupoles.size(); first += batch) {
    load_quadrupoles(quadrupoles, coordinates, centre, first, x, y, z, r2, cpp_re, cpp_im, cpz_re,
                     cpz_im, czz);
    for (std::size_t b = 0; b < batch; ++b) {
      assert(r2[b] > 0.0);
      inverse_r2[b] = 1.0 / r2[b];
      re[b] = 1.0 / std::sqrt(r2[b]);
      im[b] = 0.0;
      x[b] *= inverse_r2[b];
      y[b] *= inverse_r2[b];
      z[b] *= inverse_r2[b];
    }
    irregular_recursion(upper, x, y, z, inverse_r2, re, im);
    for (int l = 0; l <= order; ++l) {
      const int n = l + 2;
      for (int m = 0; m <= l; ++m) {
        double sum_re = 0.0;
        double sum_im = 0.0;
        for (std::size_t b = 0; b < batch; ++b) {
          const auto add = [&](BatchedValue v, double sign) {
            sum_re += sign * v.re;
            sum_im += sign * v.im;
          };
          add(multiply(cpp_re[b], cpp_im[b], batched(re, im, n, m + 2, b)), 1.0);
          add(multiply(cpp_re[b], -cpp_im[b], batched(re, im, n, m - 2, b)), 1.0);
          add(multiply(cpz_re[b], cpz_im[b], batched(re, im, n, m + 1, b)), -1.0);
          add(multiply(cpz_re[b], -cpz_im[b], batched(re, im, n, m - 1, b)), 1.0);
          add(multiply(czz[b], 0.0, batched(re, im, n, m, b)), 1.0);
        }
        local[expansion_index(l, m)] += Complex{sum_re, sum_im};
      }
    }
  }
}

auto field_tensor_error_bound(int order, int rank, double charge_sum, double source_radius,
                              double target_radius, double distance) -> double {
  assert(rank >= 0 && rank <= order);
  const double theta = (source_radius + target_radius) / distance;
  assert(theta >= 0.0 && theta < 1.0);
  double binomial = 1.0;  // C(p + 1, Lambda)
  for (int i = 1; i <= rank; ++i) {
    binomial = binomial * (order + 2 - i) / i;
  }
  return charge_sum * binomial * std::pow(theta, order - rank + 1) /
         std::pow((1.0 - theta) * distance, rank + 1);
}

auto dipole_field_tensor_error_bound(int order, int rank, double dipole_sum, double source_radius,
                                     double target_radius, double distance) -> double {
  // The expansion error e(s) of a unit charge at s is harmonic in s, so |mu . grad_s e| <=
  // |mu| (3 / delta) sup_{|s' - s| <= delta} |e(s')| <= |mu| (3 / delta) B(r_S + delta); the
  // smallest over delta = gap / 2^k, gap = d - r_S - r_T.
  const double gap = distance - source_radius - target_radius;
  assert(gap > 0.0);
  double best = std::numeric_limits<double>::infinity();
  for (int k = 1; k <= 10; ++k) {
    const double delta = gap / static_cast<double>(1 << k);
    best =
        std::min(best, 3.0 / delta *
                           field_tensor_error_bound(order, rank, dipole_sum, source_radius + delta,
                                                    target_radius, distance));
  }
  return best;
}

auto quadrupole_field_tensor_error_bound(int order, int rank, double quadrupole_sum,
                                         double source_radius, double target_radius,
                                         double distance) -> double {
  // The expansion error e(s) of a unit charge at s is harmonic in s. Its degree-2 part about s,
  // h_2(y) = y^T A y, obeys |h_2| <= (10 / (3 sqrt 3)) sup_{|s' - s| <= delta} |e(s')| on the
  // sphere |y| = delta (projection with P_2), so ||A||_2 <= (10 / (3 sqrt 3)) sup / delta^2 and
  // |1/2 Theta : grad grad e| = |Theta : A| <= 2 ||Theta||_2 ||A||_2 (traceless Theta); the
  // smallest over delta = gap / 2^k, gap = d - r_S - r_T.
  const double gap = distance - source_radius - target_radius;
  assert(gap > 0.0);
  const double constant = 20.0 / (3.0 * std::sqrt(3.0));
  double best = std::numeric_limits<double>::infinity();
  for (int k = 1; k <= 10; ++k) {
    const double delta = gap / static_cast<double>(1 << k);
    best = std::min(best,
                    constant / (delta * delta) *
                        field_tensor_error_bound(order, rank, quadrupole_sum, source_radius + delta,
                                                 target_radius, distance));
  }
  return best;
}

namespace {

// L_k = (2k + 1) / 2 integral_-1^1 |P_k(t)| dt, k = 0..16 (rounded up in the last digit).
constexpr std::array<double, 17> projection_constants = {1.0,
                                                         1.5,
                                                         1.9245008972987526,
                                                         2.2750000000000001,
                                                         2.5793099446611808,
                                                         2.8516999510581614,
                                                         3.1004147449826600,
                                                         3.3306903958680230,
                                                         3.5460879283452746,
                                                         3.7491558840341224,
                                                         3.9417912002007628,
                                                         4.1254503743358061,
                                                         4.3012802561879392,
                                                         4.4702028784552523,
                                                         4.6329726104096648,
                                                         4.7902159130285983,
                                                         4.9424597529666089};

}  // namespace

auto real_multipole_field_tensor_error_bound(int order, int rank,
                                             std::span<const double> moment_sums,
                                             double source_radius, double target_radius,
                                             double distance) -> double {
  assert(moment_sums.size() <= projection_constants.size());
  const double gap = distance - source_radius - target_radius;
  assert(gap > 0.0);
  double total = 0.0;
  for (std::size_t k = 0; k < moment_sums.size(); ++k) {
    if (moment_sums[k] == 0.0) {
      continue;
    }
    if (k == 0) {
      total += field_tensor_error_bound(order, rank, moment_sums[0], source_radius, target_radius,
                                        distance);
      continue;
    }
    const double constant = std::sqrt(2.0 * static_cast<double>(k) + 1.0) * projection_constants[k];
    double best = std::numeric_limits<double>::infinity();
    for (int j = 1; j <= 10; ++j) {
      const double delta = gap / static_cast<double>(1 << j);
      best = std::min(
          best, constant / std::pow(delta, static_cast<double>(k)) *
                    field_tensor_error_bound(order, rank, moment_sums[k], source_radius + delta,
                                             target_radius, distance));
    }
    total += best;
  }
  return total;
}

auto real_multipole_order(double accuracy, int rank, std::span<const double> moment_sums,
                          double source_radius, double target_radius, double distance,
                          int max_order) -> int {
  for (int order = rank; order <= max_order; ++order) {
    if (real_multipole_field_tensor_error_bound(order, rank, moment_sums, source_radius,
                                                target_radius, distance) <= accuracy) {
      return order;
    }
  }
  return -1;
}

}  // namespace fika::detail
