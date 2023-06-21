#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <cstdint>
#include <vector>

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "GtoBlock.hpp"

namespace gtofunc {  // gtofunc namespace

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @return the vector of contracted GTOs blocks.
 */
auto
makeGtoBlocks(const CMolecularBasis& basis,
              const CMolecule&       molecule) -> std::vector<CGtoBlock>;

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @param atoms the vector of atoms to select contracted GTOs.
 @return the vector of contracted GTOs blocks.
 */
auto
makeGtoBlocks(const CMolecularBasis&      basis,
              const CMolecule&            molecule,
              const std::vector<int64_t>& atoms) -> std::vector<CGtoBlock>;

}  // namespace gtofunc

#endif /* GtoFunc_hpp */
