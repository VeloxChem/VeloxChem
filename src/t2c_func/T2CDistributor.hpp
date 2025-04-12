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
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef T2CDistributor_hpp
#define T2CDistributor_hpp

#include <cstdint>
#include <vector>

#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace t2cfunc {  // t2cfunc namespace

/**
 Distributes buffer of integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param indexes the compressed contracted GTOs indexes.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const std::vector<int64_t>& indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last) -> void;

/**
 Distributes buffer of sclaed integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param indexes the compressed contracted GTOs indexes.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const double                factor,
                const std::vector<int64_t>& indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last) -> void;

/**
 Distributes buffer of integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param mat_type the matrix type.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const std::vector<int64_t>& bra_indexes,
                const std::vector<int64_t>& ket_indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last,
                const mat_t                 mat_type) -> void;

/**
 Distributes buffer of scaled integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param mat_type the matrix type.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const double                factor,
                const std::vector<int64_t>& bra_indexes,
                const std::vector<int64_t>& ket_indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last,
                const mat_t                 mat_type) -> void;

/**
 Distributes buffer of integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const std::vector<int64_t>& bra_indexes,
                const std::vector<int64_t>& ket_indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last,
                const bool                  ang_order) -> void;

/**
 Distributes buffer of sclaed integrals into given matrix.

 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 */
auto distribute(CSubMatrix*                 matrix,
                const TDoubleArray&         buffer,
                const double                factor,
                const std::vector<int64_t>& bra_indexes,
                const std::vector<int64_t>& ket_indexes,
                const int64_t               bra_comp,
                const int64_t               ket_comp,
                const int64_t               bra_igto,
                const int64_t               ket_first,
                const int64_t               ket_last,
                const bool                  ang_order) -> void;

}  // namespace t2cfunc

#endif /* T2CDistributor_hpp */
