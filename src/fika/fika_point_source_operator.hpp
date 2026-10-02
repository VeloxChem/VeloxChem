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

#ifndef fika_point_source_operator_hpp
#define fika_point_source_operator_hpp

// Internal: the matrix of point sources of one kind (charges, dipoles or quadrupoles), shared by
// NuclearAttractionDriver, DipolePotentialDriver and QuadrupolePotentialDriver.

#include <cstddef>
#include <string_view>

#include "fika_molecular_basis.hpp"
#include "fika_summation.hpp"
#include "fika_block_sparse_matrix.hpp"
#include "fika_molecule.hpp"
#include "fika_far_field_expansion.hpp"

namespace fika::detail {

enum class PointSourceKind { charges, dipoles, quadrupoles };

struct PointSourceSettings {
  std::string_view caller;     // error-message prefix, e.g. "fika::DipolePotentialDriver"
  std::size_t block_size = 0;  // atom pairs per task; 0: automatic
  ChargeSummation summation = ChargeSummation::automatic;
  std::size_t multipole_count = 0;  // sources from which automatic summation uses multipoles
};

/// Symmetric matrix of the sources of `kind` in `sources` (the other kinds empty) over `basis`:
/// (A|V|A) blocks with ordering A, (A|V|B) blocks with Scheme II, far sources through a
/// FarFieldExpansion with multipole summation. Throws std::invalid_argument (prefixed with
/// settings.caller) for a negative or non-finite threshold, a basis of another molecule, source
/// counts that differ from the coordinate counts, non-finite sources or coordinates, or multipole
/// summation with threshold 0.
auto point_source_matrix(const Molecule<double>& molecule, const MolecularBasis& basis,
                         const PointSources& sources, PointSourceKind kind, double threshold,
                         const PointSourceSettings& settings) -> BlockSparseMatrix;

}  // namespace fika::detail

#endif  // fika_point_source_operator_hpp
