//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef TwoIntsDistributor_hpp
#define TwoIntsDistributor_hpp

#include <cstdint>

#include "TwoIntsDistType.hpp"
#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

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
     Converts two electron integrals distributor object to text output and
     insert it into output text stream.
     
     @param output the output text stream.
     @param source the two electron integrals distributor object.
     */
    friend std::ostream& operator<<(      std::ostream&         output,
                                    const CTwoIntsDistribution& source);
};

#endif /* TwoIntsDistributor_hpp */
