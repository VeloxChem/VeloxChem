//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TwoIntsDistributor_hpp
#define TwoIntsDistributor_hpp

#include <cstdint>

#include "TwoIntsDistType.hpp"
#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"

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
     The flag indicating need of synchronization lock for updating integrals to
     integrals buffer. 
     */
    bool _needSyncLock;
    
    /**
     The number of rows.
     */
    int32_t _nRows;
    
    /**
     The number of columns.
     */
    int32_t _nColumns;
    
    /**
     The pointer to two electron integrals destination data buffer.
     */
    double* _intsData;
    
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
     Gets flag for synchronization lock

     @return true if synchronization lock is needed for distribution mode,
             false - otherwise.
     */
    bool needSyncLock() const;
    
    /**
     Distributes two electron integrals into data buffer.

     @param spherInts the spherical two electron integrals buffer.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs pairs blocks on
            bra and ket sides.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void distribute(const CMemBlock2D<double>& spherInts,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const bool                 isBraEqualKet,
                    const int32_t              iContrPair);
    
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
