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

/// Distributes buffer of integrals into Fock matrix: restricted  2J - xK.
/// @param focks  The local Fock matrices.
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

/// Distributes buffer of integrals into Fock matrix: restricted  xK.
/// @param focks  The local Fock matrices.
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
