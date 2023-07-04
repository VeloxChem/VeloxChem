#ifndef DipoleDriver_hpp
#define DipoleDriver_hpp

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class CDipoleDriver provides methods for computing two-center
 dipole integrals.

 @author Z. Rinkevicius
 */
class CDipoleDriver
{
   public:
    /**
     Creates a dipole integrals driver.
     */
    CDipoleDriver() = default;

    /**
     Computes dipole matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param point the external point of dipole reference.
     @return the dipole matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint3D& point) const -> CMatrices;
};

#endif /* DipoleDriver_hpp */
