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

#ifndef T4CLocalDistributor_hpp
#define T4CLocalDistributor_hpp

#include <array>
#include <cstddef>

#include "Matrices.hpp"
#include "Matrix.hpp"
#include "SimdArray.hpp"

namespace t4cfunc {  // t2cfunc namespace

/// Distributes buffer of integrals into Fock matrix: restricted  2J - K.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_rest_jk(CMatrices&                       focks,
                              const std::string&               suffix,
                              const CMatrix*                   density,
                              const CSimdArray<double>&        buffer,
                              const size_t                     offset,
                              const std::vector<size_t>&       a_indices,
                              const std::vector<size_t>&       b_indices,
                              const std::vector<size_t>&       c_indices,
                              const std::vector<size_t>&       d_indices,
                              const std::vector<size_t>&       a_loc_indices,
                              const std::vector<size_t>&       b_loc_indices,
                              const std::vector<size_t>&       c_loc_indices,
                              const std::vector<size_t>&       d_loc_indices,
                              const int                        a_angmom,
                              const int                        b_angmom,
                              const int                        c_angmom,
                              const int                        d_angmom,
                              const size_t                     bra_igto,
                              const std::pair<size_t, size_t>& ket_range,
                              const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted scaled 2J - xK.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param factor The exchange contribution scaling factor.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_rest_jkx(CMatrices&                       focks,
                               const std::string&               suffix,
                               const CMatrix*                   density,
                               const CSimdArray<double>&        buffer,
                               const size_t                     offset,
                               const double                     factor,
                               const std::vector<size_t>&       a_indices,
                               const std::vector<size_t>&       b_indices,
                               const std::vector<size_t>&       c_indices,
                               const std::vector<size_t>&       d_indices,
                               const std::vector<size_t>&       a_loc_indices,
                               const std::vector<size_t>&       b_loc_indices,
                               const std::vector<size_t>&       c_loc_indices,
                               const std::vector<size_t>&       d_loc_indices,
                               const int                        a_angmom,
                               const int                        b_angmom,
                               const int                        c_angmom,
                               const int                        d_angmom,
                               const size_t                     bra_igto,
                               const std::pair<size_t, size_t>& ket_range,
                               const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  J.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_rest_j(CMatrices&                       focks,
                             const std::string&               suffix,
                             const CMatrix*                   density,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const std::vector<size_t>&       a_loc_indices,
                             const std::vector<size_t>&       b_loc_indices,
                             const std::vector<size_t>&       c_loc_indices,
                             const std::vector<size_t>&       d_loc_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range,
                             const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  K.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_rest_k(CMatrices&                       focks,
                             const std::string&               suffix,
                             const CMatrix*                   density,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const std::vector<size_t>&       a_loc_indices,
                             const std::vector<size_t>&       b_loc_indices,
                             const std::vector<size_t>&       c_loc_indices,
                             const std::vector<size_t>&       d_loc_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range,
                             const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  scaled xK.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param factor The exchange contribution scaling factor.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_rest_kx(CMatrices&                       focks,
                              const std::string&               suffix,
                              const CMatrix*                   density,
                              const CSimdArray<double>&        buffer,
                              const size_t                     offset,
                              const double                     factor,
                              const std::vector<size_t>&       a_indices,
                              const std::vector<size_t>&       b_indices,
                              const std::vector<size_t>&       c_indices,
                              const std::vector<size_t>&       d_indices,
                              const std::vector<size_t>&       a_loc_indices,
                              const std::vector<size_t>&       b_loc_indices,
                              const std::vector<size_t>&       c_loc_indices,
                              const std::vector<size_t>&       d_loc_indices,
                              const int                        a_angmom,
                              const int                        b_angmom,
                              const int                        c_angmom,
                              const int                        d_angmom,
                              const size_t                     bra_igto,
                              const std::pair<size_t, size_t>& ket_range,
                              const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: general  2J - K.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_gen_jk(CMatrices&                       focks,
                             const std::string&               suffix,
                             const CMatrix*                   density,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const std::vector<size_t>&       a_loc_indices,
                             const std::vector<size_t>&       b_loc_indices,
                             const std::vector<size_t>&       c_loc_indices,
                             const std::vector<size_t>&       d_loc_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range,
                             const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: general scaled  2J - xK.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param factor The exchange contribution scaling factor.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_gen_jkx(CMatrices&                       focks,
                              const std::string&               suffix,
                              const CMatrix*                   density,
                              const CSimdArray<double>&        buffer,
                              const size_t                     offset,
                              const double                     factor,
                              const std::vector<size_t>&       a_indices,
                              const std::vector<size_t>&       b_indices,
                              const std::vector<size_t>&       c_indices,
                              const std::vector<size_t>&       d_indices,
                              const std::vector<size_t>&       a_loc_indices,
                              const std::vector<size_t>&       b_loc_indices,
                              const std::vector<size_t>&       c_loc_indices,
                              const std::vector<size_t>&       d_loc_indices,
                              const int                        a_angmom,
                              const int                        b_angmom,
                              const int                        c_angmom,
                              const int                        d_angmom,
                              const size_t                     bra_igto,
                              const std::pair<size_t, size_t>& ket_range,
                              const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: general  J.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_gen_j(CMatrices&                       focks,
                            const std::string&               suffix,
                            const CMatrix*                   density,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const std::vector<size_t>&       a_loc_indices,
                            const std::vector<size_t>&       b_loc_indices,
                            const std::vector<size_t>&       c_loc_indices,
                            const std::vector<size_t>&       d_loc_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range,
                            const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: general  K.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_gen_k(CMatrices&                       focks,
                            const std::string&               suffix,
                            const CMatrix*                   density,
                            const CSimdArray<double>&        buffer,
                            const size_t                     offset,
                            const std::vector<size_t>&       a_indices,
                            const std::vector<size_t>&       b_indices,
                            const std::vector<size_t>&       c_indices,
                            const std::vector<size_t>&       d_indices,
                            const std::vector<size_t>&       a_loc_indices,
                            const std::vector<size_t>&       b_loc_indices,
                            const std::vector<size_t>&       c_loc_indices,
                            const std::vector<size_t>&       d_loc_indices,
                            const int                        a_angmom,
                            const int                        b_angmom,
                            const int                        c_angmom,
                            const int                        d_angmom,
                            const size_t                     bra_igto,
                            const std::pair<size_t, size_t>& ket_range,
                            const bool                       diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: general  scaled xK.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param buffer  The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param factor The exchange contribution scaling factor.
/// @param a_indices The compressed contracted GTOs indexes on center A.
/// @param b_indices The compressed contracted GTOs indexes on center B.
/// @param c_indices The compressed contracted GTOs indexes on center C.
/// @param d_indices The compressed contracted GTOs indexes on center D.
/// @param a_loc_indices The compressed local contracted GTOs indexes on center A.
/// @param b_loc_indices The compressed local contracted GTOs indexes on center B.
/// @param c_loc_indices The compressed local contracted GTOs indexes on center C.
/// @param d_loc_indices The compressed local contracted GTOs indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom The angular momentum of integrals buffer on center D.
/// @param bra_igto The index of GTO on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of GTOs on ket side.
/// @param diagonal The flag indicating to use block diagonal distribution pattern of  integrals.
auto local_distribute_gen_kx(CMatrices&                       focks,
                             const std::string&               suffix,
                             const CMatrix*                   density,
                             const CSimdArray<double>&        buffer,
                             const size_t                     offset,
                             const double                     factor,
                             const std::vector<size_t>&       a_indices,
                             const std::vector<size_t>&       b_indices,
                             const std::vector<size_t>&       c_indices,
                             const std::vector<size_t>&       d_indices,
                             const std::vector<size_t>&       a_loc_indices,
                             const std::vector<size_t>&       b_loc_indices,
                             const std::vector<size_t>&       c_loc_indices,
                             const std::vector<size_t>&       d_loc_indices,
                             const int                        a_angmom,
                             const int                        b_angmom,
                             const int                        c_angmom,
                             const int                        d_angmom,
                             const size_t                     bra_igto,
                             const std::pair<size_t, size_t>& ket_range,
                             const bool                       diagonal) -> void;

}  // namespace t4cfunc

#endif /* T4CLocalDistributor_hpp */
