//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef Molecule_hpp
#define Molecule_hpp

#include <cstdint>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <mpi.h>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "AtomicRadii.hpp"

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
    double _charge{0.0};

    /**
     The multiplicity of electronic ground state.
     */
    int32_t _multiplicity{1};

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

   public:
    /**
     Creates an empty molecule object.
     */
    CMolecule() = default;

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
     Creates a molecule object by combining two molecule object.

     @param mol_1 the first molecule object.
     @param mol_2 the second molecule object.
     */
    CMolecule(const CMolecule& mol_1,
              const CMolecule& mol_2);

    /**
     Creates a sub-molecule object by slicing the molecule object.

     @param startIndex the starting index of the sub-molecule (0-based).
     @param numAtoms the number of atoms in the sub-molecule.
     */
    CMolecule getSubMolecule(int32_t startIndex,
                             int32_t numAtoms);
    
    
    /**
     Adds atom to molecule using given atom label and coordinates.

     @param atomLabel the label of atom.
     @param atomCoordinateX the coordinate X of atom.
     @param atomCoordinateY the coordinate Y of atom.
     @param atomCoordinateZ the coordinate Z of atom.
    */
    void addAtom(const std::string& atomLabel,
                 const double       atomCoordinateX,
                 const double       atomCoordinateY,
                 const double       atomCoordinateZ);

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
     Sets spin multiplicity of molecule object.

     @param multiplicity the multiplicity (2S+1) of molecule.
     */
    void setMultiplicity(const int32_t multiplicity);

    /**
     Gets charge of molecule.

     @return the charge of molecule.
     */
    double getCharge() const;

    /**
     Gets spin multiplicity of molecule.

     @return the multiplicity of molecules.
     */
    int32_t getMultiplicity() const;

    /**
     Gets total number of atoms in molecule.

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
     Gets masses of all atoms in molecule.

     @return the memory block object (all masses).
     */
    CMemBlock<double> getMasses() const;

    /**
     Computes vector of distances to closest neighbouring atom for each atom.

     @return the vector of distances.
     */
    CMemBlock<double> getMinDistances() const;
    
    /**
     Computes minimal distance from given external point (x,y,z) to closest atom
     in molecule.

     @param coordinateX the coordinate X of external point.
     @param coordinateY the coordinate Y of external point.
     @param coordinateZ the coordinate Z of external point.
     @return the minimal distance.
    */
    double getMinDistance(const double coordinateX,
                          const double coordinateY,
                          const double coordinateZ) const; 

    /**
     Gets nuclear repulsion energy for molecule assuming point charge model for
     nucleus.

     @return the nuclear repulsion energy.
     */
    double getNuclearRepulsionEnergy() const;

    /**
     Gets VDW radii of the atoms.

     @return the vector of VDW radii.
     */
    std::vector<double> getVdwRadii() const;

    /**
     Gets MK radii of the atoms.

     @return the vector of MK radii.
     */
    std::vector<double> getMkRadii() const;

    /**
     Gets covalent radii of the atoms.

     @return the vector of covalent radii.
     */
    std::vector<double> getCovalentRadii() const;

    /**
     Gets label of specific atom.

     @param iAtom the index of atom.
     @return the label of atom.
     */
    std::string getLabel(const int32_t iAtom) const;
    
    /**
     Gets indexes of atoms with given atomic label.

     @param atomLabel the label of requested atom type.
     @return the vector of atom indexes.
    */
    std::vector<int32_t> getAtomIndexes(const std::string& atomLabel) const;
    
    /**
     Gets vector of atom coordinates of specific atom.

     @param iAtom the index of atom.
     @return the vector of atom indexes.
    */
    std::vector<double> getAtomCoordinates(const int32_t iAtom) const;
    
    /**
     Gets index of nearest atom with given atom label to specific atom.

     @param iAtom the index of requested atom.
     @return the index of nearest atom to requested atom.
    */
    int32_t getIndexOfNearestAtom(const int32_t iAtom) const;
    
    
    /**
     Gets coordination number of specific atom in molecule.

     @param iAtom the index of requested atom.
     @param radius the effective coordination radius.
     @return the coordination number of atom.
    */
    int32_t getCoordinationNummber(const int32_t iAtom,
                                   const double  radius) const;

    /**
     Prints geometry of molecule as table to output stream.

     @return the output string.
     */
    std::string printGeometry() const;

    /**
     Checks if any pair of atoms in molecule is closer than given minimal
     distance. Prints error message to output stream for first pair of atoms,
     which are to close.

     @param minDistance the minimal distance.
     @return true if proximity condition is not violated, false otherwise.
     */
    bool checkProximity(const double minDistance) const;

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
