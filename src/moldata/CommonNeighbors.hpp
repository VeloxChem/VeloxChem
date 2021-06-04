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

#ifndef CommonNeighbors_hpp
#define CommonNeighbors_hpp

#include <cstdint>
#include <string>
#include <vector>
#include <set>

#include "MemBlock.hpp"
#include "DenseMatrix.hpp"
#include "TwoIndexes.hpp"
#include "FourIndexes.hpp"
#include "Molecule.hpp"

/**
 Class CCommonNeighbors stores data about molecular connectivity data.

 @author Z. Rinkevicius
 */
class CCommonNeighbors
{
    /**
     The bond cutoff radius.
     */
     double _cutRadius;
    
    /**
     The adjacency matrix.
     */
    CMemBlock<int32_t> _adjacencies;
    
    /**
     The bond distances matrix.
     */
    CDenseMatrix _bonds;
    
    /**
     The vector of unique  CNA signatures.
    */
    std::vector<CFourIndexes> _signatures;
    
    /**
     The vector of repetitions of CNA unique signatures.
    */
    std::vector<int32_t> _repetitions;
    
    /**
      The atomic identifiers of atoms.
     */
    std::vector<int32_t> _idAtomic;
    
    /**
      The elemental composition.
     */
    std::set<int32_t> _composition;
    
    /**
      The bond labels.
     */
    std::vector<std::string> _bondLabels;
    
    /**
     Computes the bonding matrix for the given molecule.
     
     @param molecule The molecule.
    */
    void _computeBonds(const CMolecule& molecule);
    
    /**
     Computes the adjacency matrix.
    */
    void _computeAdjacencies();
    
    /**
     Determines common surrounding atoms around the pair of atoms with the given radius.
        
     @param atomsPair the atoms pair.
     @param radius the radius around the atoms pair.
     @return the vector of common atoms.
    */
    std::vector<int32_t> _getCommonAtoms(const CTwoIndexes& atomsPair,
                                         const double       radius);
    
    /**
     Determines number of bonds between common atoms.
        
     @param atoms the vector of common atoms.
     @return the number of bonds between common atoms.
    */
    int32_t _getNumberOfCommonBonds(const std::vector<int32_t>& atoms);
    
    /**
     Determines longest consecutive bond between common atoms.
        
     @param atoms the vector of common atoms.
     @return the longest concecutive bond between common atoms.
    */
    int32_t _getLongestCommonBond(const std::vector<int32_t>& atoms);
    
    /**
     Determines longest consecutive bond for given path.
        
     @param path the atoms connection path.
     @return the longest concecutive bond for the given path.
    */
    int32_t _getLongestBondForPath(const std::vector<int32_t>& path);
    
    /**
     Determines bond identifier for the given pair of atoms.
        
     @param atomsPair the atoms pair.
     @return the unique identifier of bond.
    */
    int32_t _getBondIdentifier(const CTwoIndexes& atomsPair);
    
    /**
     Adds signature to vector of unique signatures.
        
     @param signature the signature to be added.
    */
    void _add(const CFourIndexes& signature);
    
    /**
     Sets bond labels.
    */
    void _setBondLabels();

   public:
    /**
     Creates an empty common neighbors object.
     */
    CCommonNeighbors();
    
    /**
     Creates a common neighbors object.
     
     @param cutRadius The bond cutoff radius.
     @param adjacencies The adjacencies matrix.
     @param bonds The bonds matrix.
     @param signatures The vector of unique CNA signatures.
     @param repetitions The vector of repetitions of unique CNA signatures.
     @param idAtomic The vector of elemental identifiers.
     @param composition The elemental composition of molecular system.
     @param bondLabels The vector of bond labels.
    */
    CCommonNeighbors(const double                     cutRadius,
                     const CMemBlock<int32_t>&        adjacencies,
                     const CDenseMatrix&              bonds,
                     const std::vector<CFourIndexes>& signatures,
                     const std::vector<int32_t>&      repetitions,
                     const std::vector<int32_t>&      idAtomic,
                     const std::set<int32_t>&         composition,
                     const std::vector<std::string>&  bondLabels);
    
    /**
     Creates a common neighbors object using molecule data.
     
     @param molecule The molecule.
     @param cutRadius The bond cutoff radius.
    */
    CCommonNeighbors(const CMolecule& molecule,
                     const double     cutRadius);

    /**
     Creates a common neighbors object by copying other common neighbors object.

     @param source the common neighbors object.
     */
    CCommonNeighbors(const CCommonNeighbors& source);

    /**
     Creates a common neighbors object by moving other common neighbors object.

     @param source the common neighbors object.
     */
    CCommonNeighbors(CCommonNeighbors&& source) noexcept;

    /**
     Destroys a common neighbors object.
     */
    ~CCommonNeighbors();

    /**
     Assigns a common neighbors object by copying other common neighbors object.

     @param source the common neighbors object.
     */
    CCommonNeighbors& operator=(const CCommonNeighbors& source);

    /**
     Assigns a common neighbors object by moving other common neighbors object.

     @param source the common neighbors object.
     */
    CCommonNeighbors& operator=(CCommonNeighbors&& source) noexcept;

    /**
     Compares common neighbors object with other common neighbors object.

     @param other the common neighbors object.
     @return true if common neighbors objects are equal, false otherwise.
     */
    bool operator==(const CCommonNeighbors& other) const;

    /**
     Compares common neighbors object with other common neighbors object.

     @param other the common neighbors object.
     @return true if common neighbors objects are not equal, false otherwise.
     */
    bool operator!=(const CCommonNeighbors& other) const;
    
    /**
     Generates signatures for the given radius.
    */
    void generate(const double radius);
    
    /**
     Gets bonds matrix.

     @return the dense matrix with pairwise bond distances.
    */
    CDenseMatrix getBonds() const;
    
    /**
     Gets  adjacencies matrix.

     @return the memory block with adjacencies matrix.
    */
    CMemBlock<int32_t> getAdjacencies() const;
    
    /**
     Gets  bond pairs vector.

     @return the vector of bond pairs.
    */
    std::vector<CTwoIndexes> getBondPairs() const;
    
    /**
     Gets unique signatures vector.

     @return the vector of signatures.
    */
    std::vector<CFourIndexes> getSignatures() const;

    /**
     Gets unique repetitions vector.

     @return the vector of repetitions.
    */
    std::vector<int32_t> getRepetitions() const;
    
    /**
     Computes Jaccard similarity index between this  object and other common neighbors object.

     @return the vector of repetitions.
    */
    double compJaccardIndex(const CCommonNeighbors& other);
    
    /**
     Finds requested signature in vector of unique signatures.
        
     @param signature the signature to be find.
     @return the index of requested signature in vector of unique signatures.
    */
    int32_t find(const CFourIndexes& signature) const;
    
    /**
     Gets signatures string representation.
        
     @return the string representation of signatures.
    */
    std::string getSignaturesRepr() const;
    
    /**
     Converts common neighbors object to text output and insert it into output text
     stream.

     @param output the output text stream.
     @param source the common neighbors object.
     */
    friend std::ostream& operator<<(      std::ostream&     output,
                                    const CCommonNeighbors& source);
};

#endif /* CommonNeighbors_hpp */
