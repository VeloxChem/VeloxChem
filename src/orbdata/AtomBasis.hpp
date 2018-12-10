//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AtomBasis_hpp
#define AtomBasis_hpp

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "mpi.h"

#include "BasisFunction.hpp"

/**
 Class CAtomBasis stores data about atomic basis set and provides set of methods
 for handling of atomic basis set data.
 
 @author Z. Rinkevicius
 */
class CAtomBasis
{
    /**
     The vector of basis function objects.
     */
    std::vector<CBasisFunction> _basisFunctions;

    /**
     The identifier of chemical element.
     */
    int32_t _idElemental;

    /**
     The maximum angular momentum.
     */
    int32_t _maxAngularMomentum;

public:

    /**
     Creates an empty atom basis object.
     */
    CAtomBasis();
   
    /**
     Creates a atom basis object by copying other atom basis object.
     
     @param source the atom basis object.
     */
    CAtomBasis(const CAtomBasis& source);

    /**
     Creates a atom basis object by moving other atom basis object.
     
     @param source the atom basis object.
     */
    CAtomBasis(CAtomBasis&& source) noexcept;

    /**
     Destroys a atom basis object.
     */
    ~CAtomBasis();
    
    /**
     Assigns a atom basis object by copying other atom basis object.
     
     @param source the atom basis object.
     */
    CAtomBasis& operator=(const CAtomBasis& source);
    
    /**
     Assigns a atom basis object by moving other atom basis object.
     
     @param source the atom basis object.
     */
    CAtomBasis& operator=(CAtomBasis&& source) noexcept;

    /**
     Compares atom basis object with other atom basis object.
     
     @param other the atom basis object.
     @return true if atom basis objects are equal, false otherwise.
     */
    bool operator==(const CAtomBasis& other) const;

    /**
     Compares atom basis object with other atom basis object.
     
     @param other the atom basis object.
     @return true if atom basis objects are not equal, false otherwise.
     */
    bool operator!=(const CAtomBasis& other) const;

    /**
     Set identifier of chemical element in atom basis.

     @param idElemental the identifier of chemical element.
     */
    void setIdElemental(const int32_t idElemental);

    /**
     Sets maximum angular momentum in atom basis.

     @param maxAngularMomentum the maximum angular momentum.
     */
    void setMaxAngularMomentum(const int32_t maxAngularMomentum);

    /**
     Adds basis function object to atom basis.

     @param basisFunction the basis function object.
     */
    void addBasisFunction(const CBasisFunction& basisFunction);

    /**
     Gets identifier of chemical element.

     @return the identifier of chemical element.
     */
    int32_t getIdElemental() const;

    /**
     Gets maximum angular momentum.

     @return the maximum angular momentum.
     */
    int32_t getMaxAngularMomentum() const;

    /**
     Gets number of basis functions with specific angular momentum.

     @param angularMomentum the angular momentum.
     @return the number of basis functions.
     */
    int32_t getNumberOfBasisFunctions(const int32_t angularMomentum) const;
    
    /**
     Gets reduced number of basis functions with specific angular momentum.
     
     @param angularMomentum the angular momentum.
     @return the reduced number of basis functions.
     */
    int32_t getNumberOfReducedBasisFunctions(const int32_t angularMomentum) const;

    /**
     Gets number of primitive Gaussain functions with requested angular momentum.

     @param angularMomentum the angular momentum.
     @return the number of primitive Gaussian functions.
     */
    int32_t getNumberOfPrimitiveFunctions(const int32_t angularMomentum) const;
    
    
    /**
     Gets number of normalization factors with requested angular momentum.

     @param angularMomentum the angular momentum.
     @return the number of normalization factors.
     */
    int32_t getNumberOfNormalizationFactors(const int32_t angularMomentum) const;

    /**
     Gets contraction string for atom basis.

     @return the contraction string.
     */
    std::string getContractionString() const;

    /**
     Gets primitives strig for atom basis.

     @return the primitives string.
     */
    std::string getPrimitivesString() const;

    /**
      Gets vector of basis function objects with specific angular momentum.

     @param angularMomentum the angular momentum.
     @return the vector of basis function objects.
     */
    std::vector<CBasisFunction> getBasisFunctions(const int32_t angularMomentum) const;

    /**
     Reduces atom basis set to valence atom basis

     @return the valence atom basis.
     */
    CAtomBasis reduceToValenceBasis() const;
    
    /**
     Broadcasts atom basis object within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t  rank,
                   MPI_Comm comm);

    /**
     Converts atom basis object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the atom basis object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CAtomBasis&   source);
};

#endif /* AtomBasis_hpp */
