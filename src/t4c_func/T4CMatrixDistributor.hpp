#ifndef T4CMatrixDistributor_hpp
#define T4CMatrixDistributor_hpp

#include <vector>
#include <string>
#include <array>

#include "Matrix.hpp"
#include "SimdArray.hpp"

/// Class CT4CMatrixDistributor provides methods for distributing single Fock matrix associated
/// with single density matrix.
class CT4CMatrixDistributor
{
    /// The Fock matrix associated with distributor.
    CMatrix* _fock;
    
    /// The density matrix associated with distributor.
    const CMatrix* _density;
    
    /// The standard label of Fock matrix.
    std::string _label;
    
    /// The scalling factor for scaling exchange contribution.
    double _exchange_factor;
    
   public:
    
    /// Creates a Fock matrix distributor.
    /// - Parameter fock : the Fock matrix.
    /// - Parameter density : the density matrix.
    /// - Parameter label : the standard label of Fock matrix.
    /// - Parameter exchange_factor : the scaling factor of exchange contribution.
    CT4CMatrixDistributor(      CMatrix&     fock,
                          const CMatrix&     density,
                          const std::string& label,
                          const double       exchange_factor);
    
    /// Distributes buffer of integrals into Fock matrix.
    /// - Parameter buffer: the integrals buffer.
    /// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
    /// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
    /// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
    /// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
    /// - Parameter bra_range: the index of the range [bra_first, bra_last) of GTOs on bra side.
    /// - Parameter ket_range: the index of the range [ket_first, ket_last) of GTOs on ket side.
    auto distribute(const CSimdArray<double>& buffer,
                    const std::vector<int>&   a_indices,
                    const std::vector<int>&   b_indices,
                    const int                 a_angmom,
                    const int                 b_angmom,
                    const std::array<int, 2>& bra_range,
                    const std::array<int, 2>& ket_range) -> void;
    
    /// Distributes buffer of integrals into Fock matrix.
    /// - Parameter buffer: the integrals buffer.
    /// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
    /// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
    /// - Parameter c_indices: the compressed contracted GTOs indexes on center C.
    /// - Parameter d_indices: the compressed contracted GTOs indexes on center D.
    /// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
    /// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
    /// - Parameter c_angmom: the angular momentum of integrals buffer on center C.
    /// - Parameter d_angmom: the angular momentum of integrals buffer on center D.
    /// - Parameter bra_range: the index of the range [bra_first, bra_last) of GTOs on bra side.
    /// - Parameter ket_range: the index of the range [ket_first, ket_last) of GTOs on ket side.
    auto distribute(const CSimdArray<double>& buffer,
                    const std::vector<int>&   a_indices,
                    const std::vector<int>&   b_indices,
                    const std::vector<int>&   c_indices,
                    const std::vector<int>&   d_indices,
                    const int                 a_angmom,
                    const int                 b_angmom,
                    const int                 c_angmom,
                    const int                 d_angmom,
                    const std::array<int, 2>& bra_range,
                    const std::array<int, 2>& ket_range) -> void;
};

#endif /* T4CMatrixDistributor_hpp */
