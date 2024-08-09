#ifndef T4COrderedMatrixDistributor_hpp
#define T4COrderedMatrixDistributor_hpp

#include <vector>
#include <string>
#include <array>

#include "Matrix.hpp"
#include "Matrices.hpp"
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"

/// Class CT4COrderedMatrixDistributor provides methods for distributing single Fock matrix associated
/// with single density matrix.
class CT4COrderedMatrixDistributor
{
    /// The Fock matrix associated with distributor.
    CMatrix* _fock;
    
    /// The density matrix associated with distributor.
    const CMatrix* _density;
    
    /// The standard label of Fock matrix.
    std::string _label;
    
    /// The scalling factor for scaling exchange contribution.
    double _exchange_factor;
    
    /// The local storage matrices.
    CMatrices _matrices;
    
    /// The local indices for center A.
    std::vector<int> _a_loc_indices;
    
    /// The local indices for center B.
    std::vector<int> _b_loc_indices;
    
    /// The local indices for center C.
    std::vector<int> _c_loc_indices;
    
    /// The local indices for center D.
    std::vector<int> _d_loc_indices;
    
    /// Sets local matrices and their indices.
    /// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
    /// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
    auto _set_local_matrices(const CGtoPairBlock& bra_gto_pair_block,
                             const CGtoPairBlock& ket_gto_pair_block) -> void;
    
    /// Generates global indices map for given indices vector.
    /// - Parameter indices: the vector of indices.
    auto _get_global_indices(const std::vector<int>& indices) -> std::vector<int>;
    
   public:
    
    /// Creates a Fock matrix distributor.
    /// - Parameter fock : the Fock matrix.
    /// - Parameter density : the density matrix.
    /// - Parameter label : the standard label of Fock matrix.
    /// - Parameter gto_pair_block : the GTOs pair block.
    /// - Parameter exchange_factor : the scaling factor of exchange contribution.
    CT4COrderedMatrixDistributor(      CMatrix*       fock,
                                 const CMatrix*       density,
                                 const CGtoPairBlock& gto_pair_block,
                                 const std::string&   label,
                                 const double         exchange_factor);
    
    /// Creates a Fock matrix distributor.
    /// - Parameter fock : the Fock matrix.
    /// - Parameter density : the density matrix.
    /// - Parameter label : the standard label of Fock matrix.
    /// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
    /// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
    /// - Parameter exchange_factor : the scaling factor of exchange contribution.
    CT4COrderedMatrixDistributor(      CMatrix*       fock,
                                 const CMatrix*       density,
                                 const CGtoPairBlock& bra_gto_pair_block,
                                 const CGtoPairBlock& ket_gto_pair_block,
                                 const std::string&   label,
                                 const double         exchange_factor);
    
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
    
    /// Accumulates local fock matrices into Fock matrix.
    /// - Parameter gto_pair_block : the GTOs pair block.
    auto accumulate(const CGtoPairBlock& gto_pair_block) -> void;
    
    /// Accumulates local fock matrices into Fock matrix.
    /// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
    /// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
    auto accumulate(const CGtoPairBlock& bra_gto_pair_block,
                    const CGtoPairBlock& ket_gto_pair_block) -> void;
};

#endif /* T4COrderedMatrixDistributor_hpp */
