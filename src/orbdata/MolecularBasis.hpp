//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularBasis_hpp
#define MolecularBasis_hpp

#include <cstdint>
#include <string>
#include <vector>

#include "mpi.h"

#include "AtomBasis.hpp"
#include "Molecule.hpp"

/**
Class CMolecularBasis stores data about molecular basis set and provides set of
methods for handling of molecular basis set data.

@author Z. Rinkevicius
*/
class CMolecularBasis
{
    /**
     The vector of atom basis objects.
     */
    std::vector<CAtomBasis> _atomicBasisSets;

    /**
     The maximum angular momentum.
     */
    int32_t _maxAngularMomentum;

    /**
     The name of molecular basis.
     */
    std::string _label;

   public:
    /**
     Creates an empty molecular basis object.
     */
    CMolecularBasis();

    /**
     Creates a molecular basis object by copying other molecular basis object.

     @param source the molecular basis object.
     */
    CMolecularBasis(const CMolecularBasis& source);

    /**
     Creates a molecular basis object by moving other molecular basis object.

     @param source the molecular basis object.
     */
    CMolecularBasis(CMolecularBasis&& source) noexcept;

    /**
     Destroys a molecular basis object.
     */
    ~CMolecularBasis();

    /**
     Assigns a molecular basis object by copying other molecular basis object.

     @param source the molecular basis object.
     */
    CMolecularBasis& operator=(const CMolecularBasis& source);

    /**
     Assigns a molecular basis object by moving other molecular basis object.

     @param source the molecular basis object.
     */
    CMolecularBasis& operator=(CMolecularBasis&& source) noexcept;

    /**
     Compares molecular basis object with other molecular basis object.

     @param other the molecular basis object.
     @return true if molecular basis objects are equal, false otherwise.
     */
    bool operator==(const CMolecularBasis& other) const;

    /**
     Compares molecular basis object with other molecular basis object.

     @param other the molecular basis object.
     @return true if molecular basis objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularBasis& other) const;

    /**
     Sets maximum angular momentum of molecular basis.

     @param maxAngularMomentum the maximum angular momentum.
     */
    void setMaxAngularMomentum(const int32_t maxAngularMomentum);

    /**
     Sets name of molecular basis.

     @param label the name of molecular basis.
     */
    void setLabel(const std::string& label);

    /**
     Adds atom basis object to molecular basis.

     @param atomBasis the atom basis object.
     */
    void addAtomBasis(const CAtomBasis& atomBasis);

    /**
     Reduces molecular basis to valence molecular basis.

     @return the valence molecular basis.
     */
    CMolecularBasis reduceToValenceBasis() const;

    /**
     Gets maximum angular momentum of molecular basis.

     @return the maximum angular momentum.
     */
    int32_t getMaxAngularMomentum() const;

    /**
     Gets maximum angular momentum of molecular basis for specific chemical
     element.

     @param idElemental the identifier of chemical element.
     @return the maximum angular momentum.
     */
    int32_t getMaxAngularMomentum(const int32_t idElemental) const;

    /**
     Gets maximum angular momentum of a molecule in this basis set.

     @param molecule the molecule.
     @return the maximum angular momentum.
     */
    int32_t getMolecularMaxAngularMomentum(const CMolecule& molecule) const;

    /**
     Gets name of molecular basis.

     @return the name of molecular basis.
     */
    std::string getLabel() const;

    /**
     Determines number of basis functions with specific angular momentum
     for selected chemical element in molecular basis.

     @param idElemental the identifier of chemical element.
     @param angularMomentum the angular momentum.
     @return number of basis functions.
     */
    int32_t getNumberOfBasisFunctions(const int32_t idElemental, const int32_t angularMomentum) const;

    /**
     Determines number of basis functions with specific angular momentum
     in molecular basis of selected molecule.

     @param molecule the molecule.
     @param angularMomentum the angular momentum.
     @return the number of basis functions.
     */
    int32_t getNumberOfBasisFunctions(const CMolecule& molecule, const int32_t angularMomentum) const;

    /**
     Determines number of basis functions with specific angular momentum
     in molecular basis of list of atoms in selected molecule.

     @param molecule the molecule.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     @param angularMomentum the angular momentum.
     @return the number of basis functions.
     */
    int32_t getNumberOfBasisFunctions(const CMolecule& molecule, const int32_t iAtom, const int32_t nAtoms, const int32_t angularMomentum) const;

    /**
     Determines number of primitive Gaussian functions with specific angular
     momentum in molecular basis of selected molecule.

     @param molecule the molecule.
     @param angularMomentum the angular momentum.
     @return the number of Gaussian functions.
     */
    int32_t getNumberOfPrimitiveBasisFunctions(const CMolecule& molecule, const int32_t angularMomentum) const;

    /**
     Determines number of primitive Gaussian functions with specific angular
     momentum in molecular basis of list of atoms in selected molecule.

     @param molecule the molecule.
     @param iAtom the index of first atom in list of atoms.
     @param nAtoms the number of atoms in list of atoms.
     @param angularMomentum the angular momentum.
     @return the number of Gaussian functions.
     */
    int32_t getNumberOfPrimitiveBasisFunctions(const CMolecule& molecule,
                                               const int32_t    iAtom,
                                               const int32_t    nAtoms,
                                               const int32_t    angularMomentum) const;

    /**
     Determines size of contracted AO basis for selected molecule.

     @param molecule the molecule.
     @return the size of contracted AO basis.
     */
    int32_t getDimensionsOfBasis(const CMolecule& molecule) const;

    /**
     Determines partial size up to specific angular momentum of contracted AO
     basis for selected molecule.

     @param molecule the molecule.
     @param angularMomentum the angular momentum.
     @return the partial size of contracted AO basis.
     */
    int32_t getPartialDimensionsOfBasis(const CMolecule& molecule, const int32_t angularMomentum) const;

    /**
     Determines size of primitive AO basis for selected molecule.

     @param molecule the molecule.
     @return the size of primitive AO basis.
     */
    int32_t getDimensionsOfPrimitiveBasis(const CMolecule& molecule) const;

    /**
     Gets atom basis object for specific chemical element from molecular basis.

     @param idElemental the identifier of chemical element.
     @return the atom basis object.
     */
    CAtomBasis getAtomBasis(const int32_t idElemental) const;

    /**
     Gets vector of basis function objects with specific angular momentum for
     specific chemical element from molecular basis.

     @param idElemental the identifier of chemical element.
     @param angularMomentum the angular momentum.
     @return  the vector of basis function objects.
     */
    std::vector<CBasisFunction> getBasisFunctions(const int32_t idElemental, const int32_t angularMomentum) const;

    /**
     Creates string representation map of basis functions.

     @param molecule the molecule.
     @return the string map of basis functions.
     */
    std::vector<std::string> getAOBasisMap(const CMolecule& molecule) const;

    /**
     Prints AO basis information to output stream for selected molecule.

     @param title the header line of AO basis output.
     @param molecule the molecule.
     */
    std::string printBasis(const char* title, const CMolecule& molecule) const;

    /**
     Broadcasts molecular basis object within domain of MPI communicator.

     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);

    /**
     Converts molecular basis object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the molecular basis object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CMolecularBasis& source);
};

#endif /* MolecularBasis_hpp */
