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

#include "fika_uniform_m2l.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "fika_point3d.hpp"
#include "fika_multipole_expansion.hpp"
#include "SerialDenseLinearAlgebra.hpp"

namespace fika::detail {

namespace {

using Complex = std::complex<double>;

auto build_offsets(int separation) -> std::vector<CellOffset> {
  const int reach = 2 * separation + 1;
  std::vector<CellOffset> offsets;
  for (int x = -reach; x <= reach; ++x) {
    for (int y = -reach; y <= reach; ++y) {
      for (int z = -reach; z <= reach; ++z) {
        if (std::max({std::abs(x), std::abs(y), std::abs(z)}) > separation) {
          offsets.push_back({x, y, z});
        }
      }
    }
  }
  return offsets;
}

/// Sign and conjugation of coefficient (l, m) under the axis reflections of an offset.
struct Reflection {
  bool x = false, y = false, z = false;

  auto conjugates() const noexcept -> bool { return x != y; }
  auto sign(int l, int m) const noexcept -> double {
    int parity = 0;
    if (x) {
      parity += m;
    }
    if (z) {
      parity += l + m;
    }
    return (parity & 1) != 0 ? -1.0 : 1.0;
  }
};

}  // namespace

UniformM2L::UniformM2L(int order, int separation) : order_(order), separation_(separation) {
  if (order < 1 || order > max_expansion_order) {
    throw std::invalid_argument("fika::UniformM2L: order " + std::to_string(order) +
                                " outside 1.." + std::to_string(max_expansion_order));
  }
  if (separation < 1 || separation > 2) {
    throw std::invalid_argument("fika::UniformM2L: separation must be 1 or 2");
  }
  size_ = expansion_size(order);
  offsets_ = build_offsets(separation);
  const std::size_t dim = 2 * size_;
  const int reach = 2 * separation + 1;
  kernel_index_.assign(static_cast<std::size_t>((reach + 1) * (reach + 1) * (reach + 1)), -1);
  const auto canonical_key = [&](int x, int y, int z) {
    return static_cast<std::size_t>((x * (reach + 1) + y) * (reach + 1) + z);
  };
  std::vector<Complex> harmonics(size_);
  for (int x = 0; x <= reach; ++x) {
    for (int y = 0; y <= reach; ++y) {
      for (int z = 0; z <= reach; ++z) {
        if (std::max({x, y, z}) <= separation) {
          continue;
        }
        kernel_index_[canonical_key(x, y, z)] = static_cast<int>(kernels_.size());
        irregular_harmonics(Point3D<double>{double(x), double(y), double(z)}, order, harmonics);
        const auto irregular = [&](int n, int mu) -> Complex {
          if (mu >= 0) {
            return harmonics[expansion_index(n, mu)];
          }
          const Complex value = std::conj(harmonics[expansion_index(n, -mu)]);
          return (mu & 1) != 0 ? -value : value;
        };
        std::vector<double> kernel(dim * dim, 0.0);
        // L_lm = (-1)^l sum_(j <= p - l) sum_k M_jk I_(j+l),(k+m)(v), with
        // M_(j,-k) = (-1)^k conj(M_jk): real and imaginary parts of the m >= 0 coefficients.
        for (int l = 0; l <= order; ++l) {
          const double sign_l = (l & 1) != 0 ? -1.0 : 1.0;
          for (int m = 0; m <= l; ++m) {
            const std::size_t re_out = expansion_index(l, m);
            const std::size_t im_out = size_ + re_out;
            for (int j = 0; j <= order - l; ++j) {
              const int n = j + l;
              for (int k = 0; k <= j; ++k) {
                const std::size_t re_in = expansion_index(j, k);
                const std::size_t im_in = size_ + re_in;
                const Complex first = irregular(n, m + k);
                double rr = first.real(), ri = -first.imag(), ir = first.imag(), ii = first.real();
                if (k > 0) {
                  const double s = (k & 1) != 0 ? -1.0 : 1.0;
                  const Complex second = irregular(n, m - k);
                  rr += s * second.real();
                  ri += s * second.imag();
                  ir += s * second.imag();
                  ii -= s * second.real();
                }
                kernel[re_out * dim + re_in] = sign_l * rr;
                kernel[re_out * dim + im_in] = sign_l * ri;
                kernel[im_out * dim + re_in] = sign_l * ir;
                kernel[im_out * dim + im_in] = sign_l * ii;
              }
            }
          }
        }
        kernels_.push_back(std::move(kernel));
      }
    }
  }
}

auto UniformM2L::interaction_offsets(int separation) -> std::vector<CellOffset> {
  return build_offsets(separation);
}

auto UniformM2L::interacts(const CellOffset& offset) const noexcept -> bool {
  const int largest = std::max({std::abs(offset[0]), std::abs(offset[1]), std::abs(offset[2])});
  return largest > separation_ && largest <= 2 * separation_ + 1;
}

void UniformM2L::apply(const CellOffset& offset, double edge,
                       std::span<const std::complex<double>> multipoles,
                       std::span<std::complex<double>> locals, std::size_t n,
                       UniformM2LWorkspace& workspace) const {
  if (!interacts(offset)) {
    throw std::invalid_argument("fika::UniformM2L: offset outside the interaction list");
  }
  const int ax = std::abs(offset[0]), ay = std::abs(offset[1]), az = std::abs(offset[2]);
  const int reach = 2 * separation_ + 1;
  if (multipoles.size() < n * size_ || locals.size() < n * size_ || !(edge > 0.0)) {
    throw std::invalid_argument("fika::UniformM2L: wrong sizes or edge");
  }
  if (n == 0) {
    return;
  }
  const Reflection reflection{offset[0] < 0, offset[1] < 0, offset[2] < 0};
  const std::vector<double>& kernel = kernels_[static_cast<std::size_t>(
      kernel_index_[static_cast<std::size_t>((ax * (reach + 1) + ay) * (reach + 1) + az)])];
  const std::size_t dim = 2 * size_;
  const double conjugate = reflection.conjugates() ? -1.0 : 1.0;
  // Per coefficient: input factor sign * h^-j (reflection, then scaling), output factor
  // sign * h^-(l+1).
  std::vector<double>& in_factor = workspace.in_factor;
  std::vector<double>& out_factor = workspace.out_factor;
  in_factor.resize(size_);
  out_factor.resize(size_);
  const double inverse_edge = 1.0 / edge;
  double power = 1.0;  // h^-l
  for (int l = 0; l <= order_; ++l) {
    for (int m = 0; m <= l; ++m) {
      const double sign = reflection.sign(l, m);
      in_factor[expansion_index(l, m)] = sign * power;
      out_factor[expansion_index(l, m)] = sign * power * inverse_edge;
    }
    power *= inverse_edge;
  }
  workspace.packed.resize(n * dim);
  workspace.product.resize(n * dim);
  for (std::size_t i = 0; i < n; ++i) {
    const Complex* source = multipoles.data() + i * size_;
    double* row = workspace.packed.data() + i * dim;
    for (std::size_t c = 0; c < size_; ++c) {
      row[c] = in_factor[c] * source[c].real();
      row[size_ + c] = in_factor[c] * conjugate * source[c].imag();
    }
  }
  sdenblas::serialMultABt(n, dim, dim, 1.0, workspace.packed.data(), dim, kernel.data(), dim, 0.0,
                          workspace.product.data(), dim);
  for (std::size_t i = 0; i < n; ++i) {
    Complex* target = locals.data() + i * size_;
    const double* row = workspace.product.data() + i * dim;
    for (std::size_t c = 0; c < size_; ++c) {
      target[c] += Complex{out_factor[c] * row[c], out_factor[c] * conjugate * row[size_ + c]};
    }
  }
}

}  // namespace fika::detail
