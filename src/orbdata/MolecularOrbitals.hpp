//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularOrbitals_hpp
#define MolecularOrbitals_hpp

#include <cstdint>
#include <vector>

#include "MolecularOrbitalsType.hpp"
#include "DenseMatrix.hpp"

/**
 Class CMolecularOrbitals stores data about molecular orbitals and provides set
 of methods for handling of molecular orbitals data.
 
 @author Z. Rinkevicius
 */
class CMolecularOrbitals
{
    /**
     The type of molecular orbitals.
     */
    morb _orbitalsType;
    
    /**
     The vector of dense matrices for storing molecular orbitals.
     */
    std::vector<CDenseMatrix> _orbitals;
    
public:
    
    /**
     Creates an empty molecular orbitals object.
     */
    CMolecularOrbitals();
    
    /**
     Creates a molecular orbitals object.
     
     @param orbitals the vector of dense matrices with molecular orbitals.
     @param orbitalsType the type of molecular orbitals.
     */
    CMolecularOrbitals(const std::vector<CDenseMatrix>& orbitals,
                       const morb                       orbitalsType);
    
    /**
     Creates a molecular orbitals object by copying other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals(const CMolecularOrbitals& source);
    
    /**
     Creates a molecular orbitals object by moving other molecular orbitals object.
     
     @param source the molecular basis object.
     */
    CMolecularOrbitals(CMolecularOrbitals&& source) noexcept;
    
    /**
     Destroys a molecular orbitals object.
     */
    ~CMolecularOrbitals();
    
    /**
     Assigns a molecular orbitals object by copying other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals& operator=(const CMolecularOrbitals& source);
    
    /**
     Assigns a molecular orbitals object by moving other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals& operator=(CMolecularOrbitals&& source) noexcept;
    
    /**
     Compares molecular orbitals object with other molecular orbitals object.
     
     @param other the molecular orbitals object.
     @return true if molecular orbitals objects are equal, false otherwise.
     */
    bool operator==(const CMolecularOrbitals& other) const;
    
    /**
     Compares molecular orbitals object with other molecular orbitals object.
     
     @param other the molecular orbitals object.
     @return true if molecular orbitals objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularOrbitals& other) const;
    
    /**
     Converts molecular orbitals object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the molecular orbitals object.
     */
    friend std::ostream& operator<<(      std::ostream&       output,
                                    const CMolecularOrbitals& source);
};

#endif /* MolecularOrbitals_hpp */
