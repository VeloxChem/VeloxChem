//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
#ifndef GtoPairsBlock_hpp
#define GtoPairsBlock_hpp

#include <cstdint>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"

/**
 Class CGtoPairsBlock stores data about basis function pairs and provides set of
 methods for manipulating with basis function pairs.
 
 @author Z. Rinkevicius
 */
class CGtoPairsBlock
{
    /**
     The angular momentum of bra side.
     */
    int32_t _braAngularMomentum;
    
    /**
     The angular momentum of ket side.
     */
    int32_t _ketAngularMomentum;
    
    /**
     The contraction scheme (contracted basis functions start/end position in
     primitive Gaussian functions vector, contracted basis function pair
     indexes).
     */
    CMemBlock2D<int32_t> _contrPattern;
    
    /**
     The primitives Gaussian function pair factors data (overlap, P coordinates,
     PA and PB distances, etc).
     */
    CMemBlock2D<double> _pairFactors;
    
    /**
     The GTO pairs screening threshold.
     */
    double _threshold;
    
    /**
     The number of unscreened primitive GTO pairs.
     */
    int32_t _nOriginalPrimPairs;
    
    /**
     The number of screened primitive GTO pairs.
     */
    int32_t _nScreenedPrimPairs;
    
    /**
     The number of unscreened contracted GTO pairs.
     */
    int32_t _nOriginalContrPairs;
    
    /**
     The number of screened contracted GTO pairs.
     */
    int32_t _nScreenedContrPairs;
    
public:
    
    /**
     Creates an empty GTOs pairs block object.
     */
    CGtoPairsBlock();
    
    /**
     Creates a GTOs pairs block object.
     
     @param contrPattern the contraction pattern (contracted basis functions
            start/end position in primitive Gaussian functions
            vector, basis functions indexes).
     @param pairFactors the various pair factors (overlap, P coordinates,
            PA and PB distances, etc)
     @param braAngularMomentum the angular momentum on bra side.
     @param ketAngularMomentum the angular momentum on ket side.
     @param threshold the primitive pairs screening threshold.
     */
    CGtoPairsBlock(const CMemBlock2D<int32_t>& contrPattern,
                   const CMemBlock2D<double>&  pairFactors,
                   const int32_t               braAngularMomentum,
                   const int32_t               ketAngularMomentum,
                   const double                threshold);
    
    /**
     Creates a GTOs pairs block object.
     
     @param braGtoBlock the GTO block object for bra side.
     @param ketGtoBlock the GTO block object for ket side.
     @param threshold the primitive pairs screening threshold.
     */
    CGtoPairsBlock(const CGtoBlock& braGtoBlock,
                   const CGtoBlock& ketGtoBlock,
                   const double     threshold);
    
    /**
     Creates a GTOs pairs block object.
     
     @param gtoBlock the GTO block object for bra and ket sides.
     @param threshold the primitive pairs screening threshold.
     */
    CGtoPairsBlock(const CGtoBlock& gtoBlock,
                   const double     threshold);
    
    /**
     Creates a GTOs pairs block object by copying other GTOs pairs block object.
     
     @param source the GTOs pairs block object.
     */
    CGtoPairsBlock(const CGtoPairsBlock& source);
    
    /**
     Creates a GTOs pairs block object by moving other GTOs pairs block object.
     
     @param source the GTOs pairs block object.
     */
    CGtoPairsBlock(CGtoPairsBlock&& source) noexcept;
    
    /**
     Destroys a GTOs pairs block object.
     */
    ~CGtoPairsBlock();
    
    /**
     Assigns a GTOs pairs block object by copying other GTOs pairs block object.
     
     @param source the GTOs pairs block object.
     */
    CGtoPairsBlock& operator=(const CGtoPairsBlock& source);
    
    /**
     Assigns a GTOs pairs block object by moving other GTOs pairs block object.
     
     @param source the GTOs pairs block object.
     */
    CGtoPairsBlock& operator=(CGtoPairsBlock&& source) noexcept;
    
    /**
     Compares GTOs pairs block object with other GTOs pairs block object.
     
     @param other the GTOs pairs block object.
     @return true if GTOs pairs block objects are equal, false otherwise.
     */
    bool operator==(const CGtoPairsBlock& other) const;
    
    /**
     Compares GTOs pairs block object with other GTOs pairs block object.
     
     @param other the GTOs pairs block object.
     @return true if GTOs pairs block objects are not equal, false otherwise.
     */
    bool operator!=(const CGtoPairsBlock& other) const;
    
    /**
     Gets angular momentum of bra side in GTOs pairs block object.
     
     @return the angular momentum of bra side.
     */
    int32_t getBraAngularMomentum() const;
    
    /**
     Gets angular momentum of ket side in GTOs pairs block object.
     
     @return the angular momentum of ket side.
     */
    int32_t getKetAngularMomentum() const;
    
    /**
     Checks if GTOs pairs block is empty.
     
     @return true if GTOs pairs block is empty, false otherwise.
     */
    bool empty() const;
    
    /**
     Gets constant vector of Obara-Saika factors: alpha_a + alpha_b.

     @return the vector of Obara-Saika factors.
     */
    const double* getFactorsXi() const;
    
    /**
     Gets constant vector of Obara-Saika factors: 1 / (alpha_a + alpha_b).
     
     @return the vector of Obara-Saika factors.
     */
    const double* getFactorsOneOverXi() const;
    
    /**
     Gets constant vector of Obara-Saika factors: alpha_a * alpha_b / (alpha_a + alpha_b).
     
     @return the vector of Obara-Saika factors.
     */
    const double* getFactorsZeta() const;
    
    /**
     Gets constant vector of primitive overlaps c_mu * c_nu  (s_mu | s_nu).
     
     @return the vector of primitive overlaps.
     */
    const double* getOverlaps() const;
    
    /**
     Gets constant vector of Cartesian X coordinates of Gaussian product center P.
     
     @return the vector of Cartesian X coordinates.
     */
    const double* getCoordinatesPX() const;
    
    /**
     Gets constant vector of Cartesian Y coordinates of Gaussian product center P.
     
     @return the vector of Cartesian Y coordinates.
     */
    const double* getCoordinatesPY() const;
    
    /**
     Gets constant vector of Cartesian Z coordinates of Gaussian product center P.
     
     @return the vector of Cartesian Z coordinates.
     */
    const double* getCoordinatesPZ() const;
    
    /**
     Gets constant vector of Cartesian X component of R(PA) = P - A distances.
     
     @return the vector of Cartesian R(PA)_X distances.
     */
    const double* getDistancesPAX() const;
    
    /**
     Gets constant vector of Cartesian Y component of R(PA) = P - A distances.
     
     @return the vector of Cartesian R(PA)_Y distances.
     */
    const double* getDistancesPAY() const;
    
    /**
     Gets constant vector of Cartesian Z component of R(PA) = P - A distances.
     
     @return the vector of Cartesian R(PA)_Y distances.
     */
    const double* getDistancesPAZ() const;
    
    /**
     Gets constant vector of Cartesian X component of R(PB) = P - B distances.
     
     @return the vector of Cartesian R(PB)_Z distances.
     */
    const double* getDistancesPBX() const;
    
    /**
     Gets vector of Cartesian Y component of R(PB) = P - B distances.
     
     @return the vector of Cartesian R(PB)_Y distances.
     */
    const double* getDistancesPBY() const;
    
    /**
     Gets constant vector of Cartesian Z component of R(PB) = P - B distances.
     
     @return the vector of Cartesian R(PB)_Z distances.
     */
    const double* getDistancesPBZ() const;
    
    /**
     Gets constant vector of Cartesian X coordinates of center A.
     
     @return the vector of Cartesian X coordinates.
     */
    const double* getCoordinatesAX() const;
    
    /**
     Gets constant vector of Cartesian Y coordinates of center A.
     
     @return the vector of Cartesian Y coordinates.
     */
    const double* getCoordinatesAY() const;
    
    /**
     Gets constant vector of Cartesian Z coordinates of center A.
     
     @return the vector of Cartesian Z coordinates.
     */
    const double* getCoordinatesAZ() const;
    
    /**
     Gets constant vector of Cartesian X coordinates of center B.
     
     @return the vector of Cartesian X coordinates.
     */
    const double* getCoordinatesBX() const;
    
    /**
     Gets constant vector of Cartesian Y coordinates of center B.
     
     @return the vector of Cartesian Y coordinates.
     */
    const double* getCoordinatesBY() const;
    
    /**
     Gets constant vector of Cartesian Z coordinates of center B.
     
     @return the vector of Cartesian Z coordinates.
     */
    const double* getCoordinatesBZ() const;
    
    /**
     Gets constant vector of Cartesian X component of R(AB) = A - B distances.
     
     @return the vector of Cartesian R(AB)_Z distances.
     */
    const double* getDistancesABX() const;
    
    /**
     Gets vector of Cartesian Y component of R(AB) = A - B distances.
     
     @return the vector of Cartesian R(AB)_Y distances.
     */
    const double* getDistancesABY() const;
    
    /**
     Gets constant vector of Cartesian Z component of R(AB) = A - B distances.
     
     @return the vector of Cartesian R(AB)_Z distances.
     */
    const double* getDistancesABZ() const;
    
    /**
     Gets constant pointer to contracted pair start positions in primitive
     pairs vector.
     
     @return the start positions of contracted pairs.
     */
    const int32_t* getStartPositions() const;
    
    /**
     Gets constant pointer to contracted pair end positions in primitive
     pairs vector.
     
     @return the end positions of contracted pairs.
     */
    const int32_t* getEndPositions() const;
    
    /**
     Gets constant pointer to contracted pairs bra indexes in full AO basis for
     specific angular momentum component.
     
     @param iComponent the component of angular momentum.
     @return the bra indexes in full AO basis.
     */
    const int32_t* getBraIdentifiers(const int32_t iComponent) const;
    
    /**
     Gets constant pointer to contracted pairs ket indexes in full AO basis for
     specific angular momentum component.
     
     @param iComponent the component of angular momentum.
     @return the ket indexes in full AO basis.
     */
    const int32_t* getKetIdentifiers(const int32_t iComponent) const;
    
    /**
     Gets number of initial primitive pairs generated from input data.

     @return the number of primitive pairs.
     */
    int32_t getNumberOfOriginalPrimPairs() const;
    
    /**
     Gets number of screened primitive pairs generated from input data.
     
     @return the number of primitive pairs.
     */
    int32_t getNumberOfScreenedPrimPairs() const;
    
    /**
     Gets number of initial contracted pairs generated from input data.
     
     @return the number of contractes pairs.
     */
    int32_t getNumberOfOriginalContrPairs() const;
    
    /**
     Gets number of screened contracted pairs generated from input data.
     
     @return the number of contracted pairs.
     */
    int32_t getNumberOfScreenedContrPairs() const;
    
    /**
     Gets pair type string for GTOs pairs object.

     @return the pair type string.
     */
    std::string getPairType() const;
    
    /**
     Gets raw size string for GTOs pairs object.
     
     @return the raw size string.
     */
    std::string getRawSizeString() const;
    
    /**
     Gets screened size string for GTOs pairs object.
     
     @return the raw size string.
     */
    std::string getScreenedSizeString() const;
    
    /**
     Converts GTOs pairs block object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the GTOs pairs block object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CGtoPairsBlock& source);
};

#endif /* GtoPairsBlock_hpp */
