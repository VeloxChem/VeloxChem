#ifndef KineticEnergyDriver_hpp
#define KineticEnergyDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class CKineticEnergyDriver provides methods for computing two-center
 kinetic energy integrals.

 @author Z. Rinkevicius
 */
class CKineticEnergyDriver
{
   public:
    /**
     Creates an kinetic energy integrals driver.
     */
    CKineticEnergyDriver() = default;

    /**
     Computes kinetic energy matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @return the kinetic energy matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix;
};

#endif /* KineticEnergyDriver_hpp */
