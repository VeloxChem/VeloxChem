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

#ifndef T2CUtils_hpp
#define T2CUtils_hpp

#include <cstddef>
#include <utility>

#include "Matrix.hpp"
#include "Point.hpp"
#include "SimdArray.hpp"
#include "SubMatrix.hpp"

namespace t2cfunc {  // t2cfunc namespace

/// @brief Computes R(AB) = A -B distances.
/// @param buffer The SIMD array containing R(AB) distances and Cartesian A, B coordinates.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_ab(CSimdArray<double>& buffer, const size_t index_ab, const size_t index_b, const TPoint<double>& r_a) -> void;

/// @brief Computes P center coordinates.
/// @param buffer The SIMD array containing R(AB) distances and Cartesian A, B coordinates.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The exponent on A center.
auto comp_coordinates_p(CSimdArray<double>& buffer, const size_t index_p, const size_t index_b, const TPoint<double>& r_a, const double a_exp)
    -> void;

/// @brief Computes R center coordinates.
/// @param buffer The SIMD array containing R(AB) distances and Cartesian A, B coordinates.
/// @param index_r The primary row index of  Cartesian R points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The exponent on A center.
/// @param c_exp The exponent on A center.
auto comp_coordinates_r(CSimdArray<double>& buffer, const size_t index_r, const size_t index_b, const TPoint<double>& r_a, const double a_exp, const double c_exp)
    -> void;

/// @brief Computes R(PB) = P - B distances.
/// @param buffer The SIMD array containing R(PB) distances.
/// @param index_pb The primary row index of R(PB) distances in SIMD array.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
auto comp_distances_pb(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_ab, const double a_exp) -> void;

/// @brief Computes R(PA) = P - A distances.
/// @param buffer The SIMD array containing R(PA) distances.
/// @param index_pa The primary row index of R(PA) distances in SIMD array.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
auto comp_distances_pa(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_ab, const double a_exp) -> void;

/// @brief Computes R(PB) = P - B distances.
/// @param buffer The SIMD array containing R(PB) distances.
/// @param index_pb The primary row index of R(PB) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
auto comp_distances_pb_from_p(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_p, const size_t index_b) -> void;

/// @brief Computes R(PA) = P - A distances.
/// @param buffer The SIMD array containing R(PA) distances.
/// @param index_pa The primary row index of R(PA) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_pa_from_p(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_p, const TPoint<double>& r_a) -> void;

/// @brief Computes R(RA) = R - A distances.
/// @param buffer The SIMD array containing R(RA) distances.
/// @param index_ra The primary row index of R(RA) distances in SIMD array.
/// @param index_r The primary row index of  Cartesian R points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_pa_from_p(CSimdArray<double>& buffer, const size_t index_ra, const size_t index_r, const TPoint<double>& r_a) -> void;

/// @brief Computes R(PC) = P - C distances.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_c The Cartesian C point coordinates.
auto comp_distances_pc(CSimdArray<double>& buffer, const size_t index_pc, const size_t index_p, const TPoint<double>& r_c) -> void;

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param a_exp The exponent on A center.
void comp_boys_args(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_pc, const double a_exp);

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param a_exp The exponent on A center.
/// @param omega The range separation factor.
void comp_boys_args(CSimdArray<double>& bf_data,
                    const size_t        index_args,
                    CSimdArray<double>& buffer,
                    const size_t        index_pc,
                    const double        a_exp,
                    const double        omega);

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(AB) distances.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
void comp_boys_args_with_rho(CSimdArray<double>& bf_data,
                             const size_t index_args,
                             CSimdArray<double>& buffer,
                             const size_t index_ab,
                             const double a_exp);

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const size_t ndims, const size_t nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param cposition The contracted array position.
/// @param pbuffer The primitive array.
/// @param pposition The starting position in primitive array.
/// @param nrows The number of rows to contract. 
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer,
            const size_t        cposition,
            CSimdArray<double>& pbuffer,
            const size_t pposition,
            const size_t nrows,
            const size_t ndims,
            const size_t nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param factor The scaling factor of primitive array.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer,
            CSimdArray<double>& pbuffer,
            const size_t        position,
            const double        factor,
            const size_t        ndims,
            const size_t        nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param facts The vector of scaling factors.
/// @param nfacts The number of scaling factors for single external center.
/// @param index The index of current scalling factor.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>&        cbuffer,
            CSimdArray<double>&        pbuffer,
            const size_t               position,
            const std::vector<double>& facts,
            const size_t               nfacts,
            const size_t               index,
            const size_t               ndims,
            const size_t               nblocks) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param indices The compressed contracted basis functions indexes.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param mat_type The matrix type.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const mat_t                      mat_type) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param ang_order The angular order of submatrix.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const bool                       ang_order) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param indices The compressed contracted basis functions indexes.
/// @param bra_angmom The angular momentum of integrals buffer on bra and ket sides.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       indices,
                const int                        bra_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_angmom The angular momentum of integrals buffer on bra side.
/// @param ket_angmom The angular momentum of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param mat_type The matrix type.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_angmom,
                const int                        ket_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const mat_t                      mat_type) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_angmom The angular momentum of integrals buffer on bra side.
/// @param ket_angmom The angular momentum of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param ang_order The angular order of submatrix.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_angmom,
                const int                        ket_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const bool                       ang_order) -> void;

}  // namespace t2cfunc

#endif /* T2CUtils_hpp */
