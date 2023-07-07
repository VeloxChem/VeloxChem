#ifndef QuadrupoleDriver_hpp
#define QuadrupoleDriver_hpp

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class CQuadrupoleDriver provides methods for computing two-center
 quadrupole integrals.

 @author Z. Rinkevicius
 */
class CQuadrupoleDriver
{
   public:
    /**
     Creates a quadrupole integrals driver.
     */
    CQuadrupoleDriver() = default;

    /**
     Computes quadrupole matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param point the external point of dipole reference.
     @return the quadrupole matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint3D& point) const -> CMatrices;
};

#endif /* QuadrupoleDriver_hpp */
