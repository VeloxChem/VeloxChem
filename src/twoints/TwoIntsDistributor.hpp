//
//                           VELOXCHEM 1.0-RC
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

#ifndef TwoIntsDistributor_hpp
#define TwoIntsDistributor_hpp

#include <cstdint>

#include "TwoIntsDistType.hpp"
#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"
#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "FockContainer.hpp"

/**
 Class CTwoIntsDistribution provides set of two electron integrals distribution
 methods.
 
 @author Z. Rinkevicius
 */
class CTwoIntsDistribution
{
    /**
     The two electron integrals distribution pattern.
     */
    dist2e _distPattern;
    
    /**
     The number of rows.
     */
    int32_t _nRows;
    
    /**
     The number of columns.
     */
    int32_t _nColumns;
    
    /**
     The identifier of GTOs pair (used for storing Q values).
     */
    int32_t _idGtoPair;
    
    /**
     The pointer to two electron integrals destination data buffer.
     */
    double* _intsData;
    
    /**
     The pointer to AO density matrix used to distribute integrals into Fock
     matrix.
     */
    const CAODensityMatrix* _aoDensity;
    
    /**
     The pointer to destination AO Fock matrix.
     */
    CAOFockMatrix* _aoFock;
    
    /**
     The pointer to Fock container object with partial Fock matrix data.
     */
    CFockContainer _fockContainer;
    
    /**
     Distributes two electron integrals into data batch.
     
     @param spherInts the spherical two electron integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
            bra and ket sides.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void _distSpherIntsIntoBatch(const CMemBlock2D<double>& spherInts,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const bool                 isBraEqualKet,
                                 const int32_t              iContrPair);
    
    /**
     Distributes two electron integrals into Q values vector. Only largest
     component from shell is stored.
     NOTE: GTOs pairs blocks on bra and ket sides must contain single
           contracted GTO.
     
     @param spherInts the spherical two electron integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
            bra and ket sides.
     @param iContrPair the index of contracted GTO pair being computed.
     */
    void _distSpherIntsIntoQValues(const CMemBlock2D<double>& spherInts,
                                   const CGtoPairsBlock&      braGtoPairsBlock,
                                   const CGtoPairsBlock&      ketGtoPairsBlock,
                                   const bool                 isBraEqualKet,
                                   const int32_t              iContrPair);
    
    /**
     Distributes two electron integrals into AO fock matrix.
     
     @param spherInts the spherical two electron integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
            bra and ket sides.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair being computed.
     */
    void _distSpherIntsIntoFock(const CMemBlock2D<double>& spherInts,
                                const CGtoPairsBlock&      braGtoPairsBlock,
                                const CGtoPairsBlock&      ketGtoPairsBlock,
                                const bool                 isBraEqualKet,
                                const int32_t              nKetContrPairs,
                                const int32_t              iContrPair);
    
    /**
     Gets starting index of spherical integrals vector in batch of integrals.

     @param nShellComponents the number of components in shell.
     @param iContrPair the index of contracted GTO pair on bra side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs
            blocks on bra and ket sides.
     @return the starting index.
     */
    int32_t _getStartIndexForBatch(const int32_t nShellComponents, 
                                   const int32_t iContrPair,
                                   const bool    isBraEqualKet) const;
    
public:
    
    /**
     Creates an empty two electron integrals distributor object.
     */
    CTwoIntsDistribution();
    
    /**
     Creates a two electron integrals distributor object.
     
     @param intsData the pointer to one electron integrals data buffer.
     @param nRows the number of rows in data buffer.
     @param nColumns the number of columns in data buffer.
     @param distPattern the two electron integrals distribution pattern.
     */
    CTwoIntsDistribution(      double* intsData,
                         const int32_t nRows,
                         const int32_t nColumns,
                         const dist2e  distPattern);
    
    /**
     Creates a two electron integrals distributor object.
     
     @param intsData the pointer to one electron integrals data buffer.
     @param nRows the number of rows in data buffer.
     @param nColumns the number of columns in data buffer.
     @param idGtoPair the index of GTOs pair.
     @param distPattern the two electron integrals distribution pattern.
     */
    CTwoIntsDistribution(      double* intsData,
                         const int32_t nRows,
                         const int32_t nColumns,
                         const int32_t idGtoPair,
                         const dist2e  distPattern);
    
    /**
     Creates a two electron integrals distributor object.
     
     @param aoFock the pointer AO Fock matrix.
     @param aoDensity the pointer AO Density matrix.
     */
    CTwoIntsDistribution(      CAOFockMatrix*    aoFock,
                         const CAODensityMatrix* aoDensity);
    
    /**
     Creates an two electron integrals distributor object by copying other
     two electron integrals distributor object.
     
     @param source the two electron integrals distributor object.
     */
    CTwoIntsDistribution(const CTwoIntsDistribution& source);
    
    /**
     Destroys an two electron integrals distributor object.
     */
    ~CTwoIntsDistribution();
    
    /**
     Assigns an two electron integrals distributor object by copying other
     two electron integrals distributor matrix object.
     
     @param source the two electron integrals distributor object.
     */
    CTwoIntsDistribution& operator=(const CTwoIntsDistribution& source);
    
    /**
     Compares two electron integrals distributor object with other two electron
     integrals distributor object.
     
     @param other the two electron integrals distributor object.
     @return true if two electron integrals distributor objects are equal,
     false otherwise.
     */
    bool operator==(const CTwoIntsDistribution& other) const;
    
    /**
     Compares two electron integrals distributor object with other two electron
     integrals distributor object.
     
     @param other the two electron integrals distributor object.
     @return true if two electron integrals distributor objects are not equal,
     false otherwise.
     */
    bool operator!=(const CTwoIntsDistribution& other) const;
    
    /**
     Allocates and initializes Fock container data for Fock matrices
     distribution method.

     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void setFockContainer(const CGtoPairsBlock&      braGtoPairsBlock,
                          const CGtoPairsBlock&      ketGtoPairsBlock);
    
    /**
     Distributes two electron integrals into data buffer.

     @param spherInts the spherical two electron integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
            bra and ket sides.
     @param nKetContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distribute(const CMemBlock2D<double>& spherInts,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const bool                 isBraEqualKet,
                    const int32_t              nKetContrPairs,
                    const int32_t              iContrPair);
    
    /**
     Accumulates AO Fock matrix from partial Fock matrix data computed inside
     single task. NOTE: Not threadsafe routine, must be allways called with
     appropiate guards to prevent raise conditions.
     */
    void accumulate();
    
    /**
     Determines maximum density elements over shell pair used for construction
     of Fock matrices for given set GTO pairs on ket side and GTO pair on bra
     side.

     @param maxDensityElements the vector of maximum density elements.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
             bra and ket sides.
     @param nKetContrPairs nKetContrPairs the number of contracted GTOs pairs
            on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void getMaxDensityElements(      CMemBlock<double>& maxDensityElements,
                               const CGtoPairsBlock&    braGtoPairsBlock,
                               const CGtoPairsBlock&    ketGtoPairsBlock,
                               const bool               isBraEqualKet,
                               const int32_t            nKetContrPairs,
                               const int32_t            iContrPair) const;
    
    /**
     Converts two electron integrals distributor object to text output and
     insert it into output text stream.
     
     @param output the output text stream.
     @param source the two electron integrals distributor object.
     */
    friend std::ostream& operator<<(      std::ostream&         output,
                                    const CTwoIntsDistribution& source);
};

#endif /* TwoIntsDistributor_hpp */
