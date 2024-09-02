#ifndef GtoPairBlockFunc_hpp
#define GtoPairBlockFunc_hpp

#include <vector>

#include "GtoPairBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis function pairs blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis function pairs blocks.
auto make_gto_pair_blocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoPairBlock>;

}  // namespace gtofunc

#endif /* GtoPairBlockFunc_hpp */
