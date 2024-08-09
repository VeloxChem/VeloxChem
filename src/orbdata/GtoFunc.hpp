#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <vector>

#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/// @brief Creates vector of basis functions blocks.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule) -> std::vector<CGtoBlock>;

/// @brief Creates vector of basis functions blocks for selected atoms.
/// @param basis The molecular basis.
/// @param molecule The molecule.
/// @param atoms The vector of atoms to select.
/// @return The vector of basis functions blocks.
auto make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms)
    -> std::vector<CGtoBlock>;

}  // namespace gtofunc

#endif /* GtoFunc_hpp */
