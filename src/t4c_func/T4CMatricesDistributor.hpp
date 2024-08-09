#ifndef T4CMatricesDistributor_hpp
#define T4CMatricesDistributor_hpp

#include <vector>
#include <string>
#include <array>

#include "Matrices.hpp"
#include "SimdArray.hpp"

/// Class CT4CMatricesDistributor provides methods for distributing multiple Fock matrix associated
/// with multiple density matrix.
class CT4CMatricesDistributor
{
    /// The Fock matrix associated with distributor.
    CMatrices* _focks;
    
    /// The density matrix associated with distributor.
    const CMatrices* _densities;
    
    /// The vector of standard labels of Fock matrices.
    std::vector<std::string> _labels;
    
    /// The vector of scalling factors for scaling exchange contributions.
    std::vector<double> _exchange_factors;
    
   public:
    
    /// Creates a Fock matrix distributor.
    /// - Parameter focks : the Fock matrices.
    /// - Parameter densities : the density matrices.
    /// - Parameter label : the vector of standard labels of Fock matrices.
    /// - Parameter exchange_factor : the vector of scaling factors of exchange contributions.
    CT4CMatricesDistributor(      CMatrices&                focks,
                            const CMatrices&                densities,
                            const std::vector<std::string>& labels,
                            const std::vector<double>&      exchange_factors);
    
    /// Distributes buffer of integrals into Fock matrix.
    /// - Parameter buffer: the integrals buffer.
    /// - Parameter a_indices: the compressed contracted GTOs indexes on center A.
    /// - Parameter b_indices: the compressed contracted GTOs indexes on center B.
    /// - Parameter a_angmom: the angular momentum of integrals buffer on center A.
    /// - Parameter b_angmom: the angular momentum of integrals buffer on center B.
    /// - Parameter bra_igto: the index of GTO on bra side.
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
    /// - Parameter bra_igto: the index of GTO on bra side.
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


#endif /* T4CMatricesDistributor_hpp */
