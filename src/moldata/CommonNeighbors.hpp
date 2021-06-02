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
#include <vector>

#include "MemBlock.hpp"
#include "DenseMatrix.hpp"
#include "TwoIndexes.hpp"
#include "ThreeIndexes.hpp"
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
    std::vector<CThreeIndexes> _signatures;
    
    /**
     The vector of repetitions of CNA unique signatures.
    */
    std::vector<int32_t> _repetitions;
    
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
    */
    CCommonNeighbors(const double                      cutRadius,
                     const CMemBlock<int32_t>&         adjacencies,
                     const CDenseMatrix&               bonds,
                     const std::vector<CThreeIndexes>& signatures,
                     const std::vector<int32_t>&       repetitions);
    
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
     Converts common neighbors object to text output and insert it into output text
     stream.

     @param output the output text stream.
     @param source the common neighbors object.
     */
    friend std::ostream& operator<<(      std::ostream&     output,
                                    const CCommonNeighbors& source);
};

#endif /* CommonNeighbors_hpp */
