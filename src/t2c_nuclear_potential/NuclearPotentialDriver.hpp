#ifndef NuclearPotentialDriver_hpp
#define NuclearPotentialDriver_hpp

#include "Matrices.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/**
 Class CNuclearPotentialDriver provides methods for computing two-center
 nuclear potential integrals.

 @author Z. Rinkevicius
 */
class CNuclearPotentialDriver
{
   public:
    /**
     Creates a nuclear potential integrals driver.
     */
    CNuclearPotentialDriver() = default;

    /**
     Computes nuclear potential matrices for given molecule, molecular basis, and vector of external points.

     @param basis the molecular basis.
     @param molecule the molecule.
     @param charges the vector of external charges.
     @param points the vector of external points.
     @return the nuclear potential matrices.
     */
    auto compute(const CMolecularBasis&       basis,
                 const CMolecule&             molecule,
                 const std::vector<double>&   charges,
                 const std::vector<TPoint3D>& points) const -> CMatrices;
};

#endif /* NuclearPotentialDriver_hpp */
