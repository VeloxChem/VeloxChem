#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @return the vector of contracted GTOs blocks.
 */
auto makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoBlock>;

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @param atoms the vector of atoms to select contracted GTOs.
 @return the vector of contracted GTOs blocks.
 */
auto makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule, const std::vector<int64_t>& atoms) -> std::vector<CGtoBlock>;

/**
 Gets number of atomic orbitals from vector of contracted GTOs blocks.

 @param gto_blocks the vector of contracted GTOs blocks.
 @return the number of atomic orbitals.
 */
auto getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int64_t;

}  // namespace gtofunc

#endif /* GtoFunc_hpp */
