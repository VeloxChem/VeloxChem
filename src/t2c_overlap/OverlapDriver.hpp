#ifndef OverlapDriver_hpp
#define OverlapDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Matrix.hpp"

/**
 Class COverlapDriver provides methods for computing two-center
 overlap integrals.

 @author Z. Rinkevicius
 */
class COverlapDriver
{
   public:
    /**
     Creates an overlap integrals driver.
     */
    COverlapDriver() = default;
    
    /**
     Computes overlap matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto
    compute(const CMolecularBasis& basis,
            const CMolecule&       molecule) const -> CMatrix;
};

#endif /* OverlapDriver_hpp */
