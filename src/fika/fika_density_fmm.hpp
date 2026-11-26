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

#ifndef fika_density_fmm_hpp
#define fika_density_fmm_hpp

// Internal: fast multipole evaluation of the field of a QM electron density (DensityMultipoles) at
// many sites.

#include <cstddef>
#include <span>

#include "fika_point3d.hpp"
#include "fika_density_multipoles.hpp"

namespace fika::detail {

/// Options of add_density_field_fmm.
struct DensityFmmOptions {
  double accuracy = 1e-10;       // field error target of the expansions (a.u., Euclidean norm)
  double theta = 0.35;           // cells interact through expansions when r_S + r_T <= theta d
  std::size_t source_leaf = 32;  // most sources per source leaf
  std::size_t target_leaf = 64;  // most sites per target leaf
  int max_order = 40;            // highest expansion order
  int order = 0;                 // starting order; 0: from the error bound (guaranteed)
  std::size_t sample_size = 64;  // sites checked against the direct sum
};

/// What add_density_field_fmm did.
struct DensityFmmReport {
  int order = 0;    // expansion order of the final evaluation
  int retries = 0;  // order increases after the sampled check
  std::size_t source_cells = 0;
  std::size_t target_cells = 0;
  std::size_t expansion_pairs = 0;  // (target cell, source cell) pairs through M2L
  std::size_t direct_pairs = 0;     // (site, source) pairs summed exactly
  double sampled_error = 0.0;       // largest field error on the sample
  double bound = 0.0;               // error bound at the final order (largest over target leaves)
};

/// Adds the field of the electrons of `sources` at `sites` (as add_density_field) through a
/// dual-tree fast multipole method: adaptive octrees over the source centres and over the sites;
/// a target cell T and a source cell S interact through P2M (add_real_multipole_to_multipole),
/// M2M, M2L, L2L and L2P when r_S + r_T <= theta |c_T - c_S| and no site of T lies within the
/// penetration radius of a source of S (|c_T - c_S| - r_T >= max_S(|P - c_S| + r_pen)); other
/// pairs split the larger cell, and leaf pairs are summed exactly (add_density_field_block). The
/// order is the smallest whose error bound (real_multipole_field_tensor_error_bound, summed over
/// the expansions reaching each target leaf) is at most `accuracy`; the field at `sample_size`
/// sites is then compared with the direct sum and the order raised by 2 while the error exceeds
/// the accuracy. Each site sums its contributions in a fixed order (independent of the thread
/// count). Throws std::invalid_argument for a wrong field size or invalid options and
/// std::runtime_error if max_order does not reach the accuracy.
auto add_density_field_fmm(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                           std::span<Point3D<double>> field, const DensityFmmOptions& options = {})
    -> DensityFmmReport;

}  // namespace fika::detail

#endif  // fika_density_fmm_hpp
