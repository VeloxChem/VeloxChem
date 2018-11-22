//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef CauchySchwarzScreener_hpp
#define CauchySchwarzScreener_hpp

#include <cstdint>

#include "MemBlock.hpp"
#include "EriScreenerType.hpp"
#include "GtoPairsBlock.hpp"

/**
 Class CCauchySchwarzScreener class stores information required for screening of
 electron repulsion integrals using Cauchy-Schwarz (QQ) and distance dependent
 Cauchy-Schwarz (QQR) methods and provides set of methods for computing and
 handling of screening data.
 
 @author Z. Rinkevicius
 */
class CCauchySchwarzScreener
{
    /**
     The screening scheme used in Cauchy-Schwarz screener object. 
     */
    ericut _screeningScheme;
    
    /**
     The vector of Q values, sqrt([ab|ab]), for GTOs pairs batch on bra side.
     */
    CMemBlock<double> _braQValues;
    
    /**
     The vector of Q values, sqrt([cd|cd]), for GTOs pairs batch on ket side.
     */
    CMemBlock<double> _ketQValues;
    
    /**
     The vector of GTO pair extends for GTOs pairs batch on bra side.
     */
    CMemBlock<double> _braPairExtends;
    
    /**
     The vector of GTO pair extends for GTOs pairs batch on ket side.
     */
    CMemBlock<double> _ketPairExtends;
    
    /**
     The screening threshold of electron repulsion integrals.
     */
    double _threshold;
    
    /**
     Sets maximum GTOs pairs extents for given GTOs pairs block. Extents are
     computed using Eq. B4 from Maurer et al., J. Chem. Phys. 136,
     144107 (2012).

     @param gtoPairExtents the vector of GTOs pairs extents. 
     @param gtoPairsBlock the GTOs pairs block.
     */
    void _setPairsExtents(      CMemBlock<double>& gtoPairExtents,
                          const CGtoPairsBlock&    gtoPairsBlock);
    
public:
    
    /**
     Creates an empty Cauchy-Schwarz screener object.
     */
    CCauchySchwarzScreener();
    
    /**
     Creates a Cauchy-Schwarz screener object.
     
     @param braQValues the vector of maximum Q values on bra side.
     @param ketQValues the vector of maximum Q values on ket side.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param screeningScheme the screening scheme.
     @param threshold the cut-off threshold of electron repulsion integrals.
     */
    CCauchySchwarzScreener(const CMemBlock<double>& braQValues,
                           const CMemBlock<double>& ketQValues,
                           const CGtoPairsBlock&    braGtoPairsBlock,
                           const CGtoPairsBlock&    ketGtoPairsBlock,
                           const ericut             screeningScheme,
                           const double             threshold);
    
    /**
     Creates a Cauchy-Schwarz screener object by copying other Cauchy-Schwarz
     screener object.
     
     @param source the Cauchy-Schwarz screener object.
     */
    CCauchySchwarzScreener(const CCauchySchwarzScreener& source);
    
    /**
     Creates a Cauchy-Schwarz screener object by moving other Cauchy-Schwarz
     screener object.
     
     @param source the Cauchy-Schwarz screener object.
     */
    CCauchySchwarzScreener(CCauchySchwarzScreener&& source) noexcept;
    
    /**
     Destroys a Cauchy-Schwarz screener object.
     */
    ~CCauchySchwarzScreener();
    
    /**
     Assigns a Cauchy-Schwarz screener object by copying other Cauchy-Schwarz
     screener object.
     
     @param source the Cauchy-Schwarz screener object.
     */
    CCauchySchwarzScreener& operator=(const CCauchySchwarzScreener& source);
    
    /**
     Assigns a Cauchy-Schwarz screener object by moving other Cauchy-Schwarz
     screener object.
     
     @param source the Cauchy-Schwarz screener object.
     */
    CCauchySchwarzScreener& operator=(CCauchySchwarzScreener&& source) noexcept;
    
    /**
     Compares Cauchy-Schwarz screener object with other Cauchy-Schwarz screener
     object.
     
     @param other the Cauchy-Schwarz screener object.
     @return true if Cauchy-Schwarz screener objects are equal, false otherwise.
     */
    bool operator==(const CCauchySchwarzScreener& other) const;
    
    /**
     Compares Cauchy-Schwarz screener object with other Cauchy-Schwarz screener
     object.
     
     @param other the Cauchy-Schwarz screener object.
     @return true if Cauchy-Schwarz screener objects are not equal, false otherwise.
     */
    bool operator!=(const CCauchySchwarzScreener& other) const;
    
    /**
     Sets screening scheme for Cauchy-Schwarz screener object.

     @param screeningScheme the screening scheme.
     */
    void setScreeningScheme(const ericut screeningScheme);
    
    /**
     Sets threshold for screening of electron repulsion integrals.

     @param threshold the screening threshold.
     */
    void setThreshold(const double threshold);
    
    /**
     Gets screening scheme used by Cauchy-Schwarz screener object.

     @return the screening scheme.
     */
    ericut getScreeningScheme() const;
    
    /**
     Gets screening threshold of electron repulion integrals.

     @return the screening threshold.
     */
    double getThreshold() const;
    
    /**
     Gets constant pointer to first element of Q values vector on bra side.
     
     @return the constant pointer to Q values vector.
     */
    const double* getBraQValues() const;
    
    /**
     Gets constant pointer to first element of Q values vector on ket side.
     
     @return the constant pointer to Q values vector.
     */
    const double* getKetQValues() const;
    
    /**
     Checks if Cauchy-Schwarz screener object is empty.

     @return true if Cauchy-Schwarz screener object is empty, false - otherwise.
     */
    bool isEmpty() const;
    
    /**
     Sets screening vector with elements 0 (exclude GTOs pair) or 1 (include
     GTOs pair).

     @param qqVector the screening vector.
     @param pqDistances the vector of effective distances between contracted
            GTOs pairs on bra and ket sides.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void setScreeningVector(      CMemBlock<int32_t>& qqVector,
                            const CMemBlock<double>&  pqDistances,
                            const bool                isBraEqualKet,
                            const int32_t             iContrPair) const;
    
    /**
     Sets screening vector with elements 0 (exclude GTOs pair) or 1 (include
     GTOs pair) for AO based screening.

     @param qqVector the screening vector.
     @param qqIndexes the vector of indexes to QQ values.
     @param maxDensityElements the vector of maximum density elements.
     @param pqDistances the vector of effective distances between contracted
            GTOs pairs on bra and ket sides.
     @param nContrPairs the number of contracted GTOs pairs on ket side.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void setScreeningVector(      CMemBlock<int32_t>& qqVector,
                            const CMemBlock<int32_t>& qqIndexes,
                            const CMemBlock<double>&  maxDensityElements,
                            const CMemBlock<double>&  pqDistances,
                            const int32_t             nContrPairs,
                            const int32_t             iContrPair) const;
    
    /**
     Converts Cauchy-Schwarz screener object to text output and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the Cauchy-Schwarz screener object.
     */
    friend std::ostream& operator<<(      std::ostream&           output,
                                    const CCauchySchwarzScreener& source);
}; 

#endif /* CauchySchwarzScreener_hpp */
