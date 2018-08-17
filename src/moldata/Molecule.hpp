//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef Molecule_hpp
#define Molecule_hpp

#include <cstdint>
#include <ostream>
#include <vector>
#include <set>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "OutputStream.hpp"

/**
 Class CMolecule stores data about single molecule and provides set of methods
 for handling of molecular data.
 
 @author Z. Rinkevicius
 */
class CMolecule
{
    /**
     The charge of molecule.
     */
    double _charge;

    /**
     The multiplicity of electronic ground state.
     */
    int32_t _multiplicity;
    
    /**
     The vector of atomic coordinates (x, y, z).
     */
    CMemBlock2D<double> _atomCoordinates;
    
    /**
     The vector of atomic charges.
     */
    CMemBlock<double> _atomCharges;
    
    /**
     The vector of atomic masses.
     */
    CMemBlock<double> _atomMasses;
    
    /**
     The vector of chemical element names.
     */
    std::vector<std::string> _atomLabels;
    
    /**
     The vector of global indexes.
     */
    CMemBlock<int32_t> _idsAtomic;
    
    /**
     The vector of chemical element identifiers.
     */
    CMemBlock<int32_t> _idsElemental;
    
    /**
     Prints atoms proximity error message to output stream, if distance between
     two atoms is shorter than given minimal distance.

     @param iAtomA the index of first atom.
     @param iAtomB the index of second atom.
     @param minDistance the minimal distance between atoms.
     @param oStream the output stream.
     */
    void _errorAtomsToClose(const int32_t        iAtomA,
                            const int32_t        iAtomB,
                            const double         minDistance,
                                  COutputStream& oStream) const;

public:

    /**
     Creates an empty molecule object.
     */
    CMolecule();
    
    /**
     Creates a molecule object.
     
     @param atomCoordinates the vector (all x,all y, all z) of atomic coordinates.
     @param atomCharges the vector of atomic charges.
     @param atomMasses the vector of atomic masses.
     @param atomLabels the vector of atomic names.
     @param idsElemental the vector of atomic identifiers.
     */
    CMolecule(const std::vector<double>&      atomCoordinates,
              const std::vector<double>&      atomCharges,
              const std::vector<double>&      atomMasses,
              const std::vector<std::string>& atomLabels,
              const std::vector<int32_t>&     idsElemental);

    /**
     Creates a molecule object by copying other molecule object.
     
     @param source the molecule object.
     */
    CMolecule(const CMolecule& source);

    /**
     Creates a molecule object by moving other molecule object.
     
     @param source the molecule object.
     */
    CMolecule(CMolecule&& source) noexcept;

    /**
     Destroys a molecule object.
     */
    ~CMolecule();

    /**
     Assigns a molecule object by copying other molecule object.
     
     @param source the molecule object.
     */
    CMolecule& operator=(const CMolecule& source);

    /**
     Assigns a molecule object by moving other molecule object.
     
     @param source the molecule object.
     */
    CMolecule& operator=(CMolecule&& source) noexcept;

    /**
     Compares molecule object with other molecule object.
     
     @param other the molecule object.
     @return true if molecule objects are equal, false otherwise.
     */
    bool operator==(const CMolecule& other) const;

    /**
     Compares molecule object with other molecule object.
     
     @param other the molecule object.
     @return true if molecule objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecule& other) const;
    
    /**
     Sets atomic indexes of all atoms in molecule object.

     @param startIndex the starting index of atomic indexes.
     */
    void setAtomicIndexes(const int32_t startIndex);
    
    /**
     Sets charge of molecule object.

     @param charge the charge of molecule.
     */
    void setCharge(const double charge);
    
    /**
     Sets multiplicity of molecule object.

     @param multiplicity the multiplicity (2S+1) of molecule.
     */
    void setMultiplicity(const int32_t multiplicity);
    
    /**
     Gets a charge of molecule.

     @return the charge of molecule.
     */
    double getCharge() const;
    
    /**
     Gets a multiplicity of molecule.

     @return the multiplicity of molecules.
     */
    int32_t getMultiplicity() const;
    
    /**
     Gets a total number of atoms in molecule.

     @return the total number of atoms.
     */
    int32_t getNumberOfAtoms() const;
    
    /**
     Gets number of atoms belonging to specific chemical element in molecule.

     @param idElemental the chemical element number.
     @return the number of atoms.
     */
    int32_t getNumberOfAtoms(const int32_t idElemental) const;
    
    /**
    Gets number of atoms belonging to specific chemical element in list of atoms
    in molecule.
    
    @param iAtom the index of first atom in list of atoms.
    @param nAtoms the number of atoms in list of atoms.
    @param idElemental the chemical element number.
    @return the number of atoms.
    */
    int32_t getNumberOfAtoms(const int32_t iAtom,
                             const int32_t nAtoms,
                             const int32_t idElemental) const;
    
    /**
     Gets set of unique chemical elements in molecule.

     @return the set of unique chemical elements.
     */
    std::set<int32_t> getElementalComposition() const;
    
    /**
     Gets a number of electrons in molecule.

     @return the number of electrons.
     */
    int32_t getNumberOfElectrons() const;
    
    /**
     Gets constant pointer to vector of chemical element identifiers.

     @return the constant pointer to chemical element indentifiers.
     */
    const int32_t* getIdsElemental() const;
    
    /**
     Gets vector Cartesian X coordinates of atoms in molecule.

     @return the constant pointer to vector of coordinates.
     */
    const double* getCoordinatesX() const;
    
    /**
     Gets vector Cartesian Y coordinates of atoms in molecule.
     
     @return the constant pointer to vector of coordinates.
     */
    const double* getCoordinatesY() const;
    
    /**
     Gets vector Cartesian Z coordinates of atoms in molecule.
     
     @return the constant pointer to vector of coordinates.
     */
    const double* getCoordinatesZ() const;
    
    /**
     Gets coordinates of all atoms in molecule.

     @return the 2D memory block object (all x, all y, all z).
     */
    CMemBlock2D<double> getCoordinates() const;
    
    /**
     Gets charges of all atoms in molecule.
     
     @return the memory block object (all charges).
     */
    CMemBlock<double> getCharges() const;
    
    /**
     Computes vector of distances to closest neighbouring atom for each atom.

     @return the vector of distances.
     */
    CMemBlock<double> getMinDistances() const;
    
    /**
     Prints geometry of molecule as table to output stream.

     @param oStream the output stream.
     */
    void printGeometry(COutputStream& oStream) const;
    
    /**
     Checks if any pair of atoms in molecule is closer than given minimal
     distance. Prints error message to output stream for first pair of atoms,
     which are to close.

     @param minDistance the minimal distance.
     @param oStream the output stream.
     @return true if proximity condition is not violated, false otherwise.
     */
    bool checkProximity(const double         minDistance,
                              COutputStream& oStream) const;
    
    /**
     Broadcasts molecule object within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);
    
    /**
     Converts molecule object to text output and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the molecule object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CMolecule&    source);
};

#endif /* Molecule_hpp */
