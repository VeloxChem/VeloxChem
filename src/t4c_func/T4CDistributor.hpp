#ifndef T4CDistributor_hpp
#define T4CDistributor_hpp

#include <array>

#include "SimdArray.hpp"
#include "Matrix.hpp"

namespace t4cfunc {  // t2cfunc namespace


/// Distributes buffer of integrals into Fock matrix: restricted  2J - K.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_jk(      CMatrix*            fock,
                   const CMatrix*            density,
                   const CSimdArray<double>& buffer,
                   const int                 offset,
                   const std::vector<int>&   a_indices,
                   const std::vector<int>&   b_indices,
                   const std::vector<int>&   c_indices,
                   const std::vector<int>&   d_indices,
                   const int                 a_angmom,
                   const int                 b_angmom,
                   const int                 c_angmom,
                   const int                 d_angmom,
                   const int                 bra_igto,
                   const std::array<int, 2>& ket_range,
                   const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  2J - xK.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter factor : the scaling fator for exchange contribution to Fock matrix.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_jkx(      CMatrix*            fock,
                    const CMatrix*            density,
                    const CSimdArray<double>& buffer,
                    const int                 offset,
                    const double              factor,
                    const std::vector<int>&   a_indices,
                    const std::vector<int>&   b_indices,
                    const std::vector<int>&   c_indices,
                    const std::vector<int>&   d_indices,
                    const int                 a_angmom,
                    const int                 b_angmom,
                    const int                 c_angmom,
                    const int                 d_angmom,
                    const int                 bra_igto,
                    const std::array<int, 2>& ket_range,
                    const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  J.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_j(      CMatrix*            fock,
                  const CMatrix*            density,
                  const CSimdArray<double>& buffer,
                  const int                 offset,
                  const std::vector<int>&   a_indices,
                  const std::vector<int>&   b_indices,
                  const std::vector<int>&   c_indices,
                  const std::vector<int>&   d_indices,
                  const int                 a_angmom,
                  const int                 b_angmom,
                  const int                 c_angmom,
                  const int                 d_angmom,
                  const int                 bra_igto,
                  const std::array<int, 2>& ket_range,
                  const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  K.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_k(      CMatrix*            fock,
                  const CMatrix*            density,
                  const CSimdArray<double>& buffer,
                  const int                 offset,
                  const std::vector<int>&   a_indices,
                  const std::vector<int>&   b_indices,
                  const std::vector<int>&   c_indices,
                  const std::vector<int>&   d_indices,
                  const int                 a_angmom,
                  const int                 b_angmom,
                  const int                 c_angmom,
                  const int                 d_angmom,
                  const int                 bra_igto,
                  const std::array<int, 2>& ket_range,
                  const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  xK.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter factor : the scaling fator for exchange contribution to Fock matrix.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_kx(      CMatrix*            fock,
                   const CMatrix*            density,
                   const CSimdArray<double>& buffer,
                   const int                 offset,
                   const double              factor,
                   const std::vector<int>&   a_indices,
                   const std::vector<int>&   b_indices,
                   const std::vector<int>&   c_indices,
                   const std::vector<int>&   d_indices,
                   const int                 a_angmom,
                   const int                 b_angmom,
                   const int                 c_angmom,
                   const int                 d_angmom,
                   const int                 bra_igto,
                   const std::array<int, 2>& ket_range,
                   const bool                diagonal) -> void;


/// Distributes buffer of integrals into Fock matrix: restricted  general J.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_gen_j(      CMatrix*            fock,
                      const CMatrix*            density,
                      const CSimdArray<double>& buffer,
                      const int                 offset,
                      const std::vector<int>&   a_indices,
                      const std::vector<int>&   b_indices,
                      const std::vector<int>&   c_indices,
                      const std::vector<int>&   d_indices,
                      const int                 a_angmom,
                      const int                 b_angmom,
                      const int                 c_angmom,
                      const int                 d_angmom,
                      const int                 bra_igto,
                      const std::array<int, 2>& ket_range,
                      const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted general K.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_gen_k(      CMatrix*            fock,
                      const CMatrix*            density,
                      const CSimdArray<double>& buffer,
                      const int                 offset,
                      const std::vector<int>&   a_indices,
                      const std::vector<int>&   b_indices,
                      const std::vector<int>&   c_indices,
                      const std::vector<int>&   d_indices,
                      const int                 a_angmom,
                      const int                 b_angmom,
                      const int                 c_angmom,
                      const int                 d_angmom,
                      const int                 bra_igto,
                      const std::array<int, 2>& ket_range,
                      const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted general xK.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter factor : the scaling fator for exchange contribution to Fock matrix.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_gen_kx(      CMatrix*            fock,
                       const CMatrix*            density,
                       const CSimdArray<double>& buffer,
                       const int                 offset,
                       const double              factor,
                       const std::vector<int>&   a_indices,
                       const std::vector<int>&   b_indices,
                       const std::vector<int>&   c_indices,
                       const std::vector<int>&   d_indices,
                       const int                 a_angmom,
                       const int                 b_angmom,
                       const int                 c_angmom,
                       const int                 d_angmom,
                       const int                 bra_igto,
                       const std::array<int, 2>& ket_range,
                       const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  general 2J-K.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_gen_jk(      CMatrix*            fock,
                       const CMatrix*            density,
                       const CSimdArray<double>& buffer,
                       const int                 offset,
                       const std::vector<int>&   a_indices,
                       const std::vector<int>&   b_indices,
                       const std::vector<int>&   c_indices,
                       const std::vector<int>&   d_indices,
                       const int                 a_angmom,
                       const int                 b_angmom,
                       const int                 c_angmom,
                       const int                 d_angmom,
                       const int                 bra_igto,
                       const std::array<int, 2>& ket_range,
                       const bool                diagonal) -> void;

/// Distributes buffer of integrals into Fock matrix: restricted  general 2J-xK.
/// - Parameter fock : the pointer to Fock matrix.
/// - Parameter density : the pointer to AO density matrix.
/// - Parameter buffer : the integrals buffer.
/// - Parameter offset : the intgeral buffer offset.
/// - Parameter factor : the scaling fator for exchange contribution to Fock matrix.
/// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
/// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
/// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
/// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
/// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
/// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
/// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
/// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
/// - Parameter bra_igto: the index of GTO on bra side.
/// - Parameter ket_indices: the index of the range [ket_first, ket_last) of GTOs on ket side.
/// - Parameter diagonal: the flag indicating to use block diagonal distribution pattern of  integrals.
auto
distribute_rest_gen_jkx(      CMatrix*            fock,
                        const CMatrix*            density,
                        const CSimdArray<double>& buffer,
                        const int                 offset,
                        const double              factor, 
                        const std::vector<int>&   a_indices,
                        const std::vector<int>&   b_indices,
                        const std::vector<int>&   c_indices,
                        const std::vector<int>&   d_indices,
                        const int                 a_angmom,
                        const int                 b_angmom,
                        const int                 c_angmom,
                        const int                 d_angmom,
                        const int                 bra_igto,
                        const std::array<int, 2>& ket_range,
                        const bool                diagonal) -> void;


}  // namespace t4cfunc

#endif /* T4CDistributor_hpp */
