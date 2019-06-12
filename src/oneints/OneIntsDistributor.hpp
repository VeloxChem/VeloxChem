//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef OneIntsDistributor_hpp
#define OneIntsDistributor_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "OneIntsDistType.hpp"

/**
 Class COneIntsDistribution provides set of one electron integrals distribution
 methods.

 @author Z. Rinkevicius
 */
class COneIntsDistribution
{
    /**
     The one electron integrals distribution pattern.
     */
    dist1e _distPattern;

    /**
     The number of rows.
     */
    int32_t _nRows;

    /**
     The number of columns.
     */
    int32_t _nColumns;

    /**
     The pointer to one electron integrals destination data buffer.
     */
    double* _intsData;

    /**
     Distributes one electron integrals into data batch.

     @param spherInts the spherical one electron integrals buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _distSpherIntsIntoBatch(const CMemBlock2D<double>& spherInts,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto);

    /**
     Distributes one electron integrals into symmetric matrix.

     @param spherInts the spherical one electron integrals buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs blocks on bra and
            ket sides.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _distSpherIntsIntoSymMatrix(const CMemBlock2D<double>& spherInts,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const bool                 isBraEqualKet,
                                     const int32_t              iContrGto);

    /**
     Distributes one electron integrals into anti-symmetric matrix.

     @param spherInts the spherical one electron integrals buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs blocks on bra and
            ket sides.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _distSpherIntsIntoAntiSymMatrix(const CMemBlock2D<double>& spherInts,
                                         const CGtoBlock&           braGtoBlock,
                                         const CGtoBlock&           ketGtoBlock,
                                         const bool                 isBraEqualKet,
                                         const int32_t              iContrGto);

    /**
     Distributes one electron integrals into rectangular matrix.

     @param spherInts the spherical one electron integrals buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _distSpherIntsIntoRectMatrix(const CMemBlock2D<double>& spherInts,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto);

   public:
    /**
     Creates an empty one electron integrals distributor object.
     */
    COneIntsDistribution();

    /**
     Creates an one electron integrals distributor object.

     @param intsData the pointer to one electron integrals data buffer.
     @param nRows the number of rows in data buffer.
     @param nColumns the number of columns in data buffer.
     @param distPattern the one electron integrals distribution pattern.
     */
    COneIntsDistribution(double* intsData, const int32_t nRows, const int32_t nColumns, const dist1e distPattern);

    /**
     Creates an one electron integrals distributor object by copying other
     one electron integrals distributor object.

     @param source the one electron integrals distributor object.
     */
    COneIntsDistribution(const COneIntsDistribution& source);

    /**
     Destroys an one electron integrals distributor object.
     */
    ~COneIntsDistribution();

    /**
     Assigns an one electron integrals distributor object by copying other
     one electron integrals distributor matrix object.

     @param source the one electron integrals distributor object.
     */
    COneIntsDistribution& operator=(const COneIntsDistribution& source);

    /**
     Compares one electron integrals distributor object with other one electron
     integrals distributor object.

     @param other the one electron integrals distributor object.
     @return true if one electron integrals distributor objects are equal,
             false otherwise.
     */
    bool operator==(const COneIntsDistribution& other) const;

    /**
     Compares one electron integrals distributor object with other one electron
     integrals distributor object.

     @param other the one electron integrals distributor object.
     @return true if one electron integrals distributor objects are not equal,
             false otherwise.
     */
    bool operator!=(const COneIntsDistribution& other) const;

    /**
     Distributes one electron integrals into data buffer.

     @param spherInts the spherical one electron integrals buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param isBraEqualKet the flag indicating equality of GTOs blocks on bra and
            ket sides.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void distribute(const CMemBlock2D<double>& spherInts,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const bool                 isBraEqualKet,
                    const int32_t              iContrGto);

    /**
     Converts one electron integrals distributor object to text output and
     insert it into output text stream.

     @param output the output text stream.
     @param source the one electron integrals distributor object.
     */
    friend std::ostream& operator<<(std::ostream& output, const COneIntsDistribution& source);
};

#endif /* OneIntsDistributor_hpp */
