#ifndef OverlapDriver_hpp
#define OverlapDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

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
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix;
    
    /**
     Computes overlap matrix for given molecule and pair of molecular bases.

     @param bra_basis the molecular basis on bra side.
     @param ket_basis the molecular basis on ket side.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto compute(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis, const CMolecule& molecule) const -> CMatrix;
};

#endif /* OverlapDriver_hpp */
