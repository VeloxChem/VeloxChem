#ifndef OctupoleDriver_hpp
#define OctupoleDriver_hpp

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class COctupoleDriver provides methods for computing two-center
 octupole integrals.

 @author Z. Rinkevicius
 */
class COctupoleDriver
{
   public:
    /**
     Creates a octupole integrals driver.
     */
    COctupoleDriver() = default;

    /**
     Computes octupole matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param point the external point of dipole reference.
     @return the quadrupole matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint3D& point) const -> CMatrices;
};

#endif /* OctupoleDriver_hpp */
