//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
