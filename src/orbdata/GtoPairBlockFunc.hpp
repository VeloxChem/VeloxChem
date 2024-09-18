#ifndef GtoPairBlockFunc_hpp
#define GtoPairBlockFunc_hpp

#include <vector>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis function pairs blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoPairBlock>;

/// @brief Creates vector of basis function pairs blocks for vector of basis functions blocks.
/// @param gto_blocks The vector of basis functions blocks on bra and ket sides.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<CGtoPairBlock>;

/// @brief Creates vector of basis function pairs blocks for pair of vectors of basis functions blocks.
/// @param bra_gto_blocks The vector of basis functions blocks on bra side.
/// @param ket_gto_blocks The vector of basis functions blocks on ket side.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const std::vector<CGtoBlock>& bra_gto_blocks, const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<CGtoPairBlock>; 

}  // namespace gtofunc

#endif /* GtoPairBlockFunc_hpp */
