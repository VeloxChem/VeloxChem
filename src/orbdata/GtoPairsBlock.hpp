//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
#ifndef GtoPairsBlock_hpp
#define GtoPairsBlock_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "CudaDevices.hpp"

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
     The contracted Gaussian function pair factors data (effective P coordinates).
     */
    CMemBlock2D<double> _contrFactors;

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

    /**
     Determines recommended dimensions of GTOs pairs block.

     @return the number of GTOs pairs.
     */
    int32_t _getBlockDimensions() const;
    
    /**
     Determines recommended dimensions of GTOs pairs block for GPU compute code.
     
     @return the number of GTOs pairs.
     */
    int32_t _getBlockDimensionsForGPU() const;

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
     @param contrFactors the various contracted pair factors (effective P
            coordintes)
     @param pairFactors the various primitive pair factors (overlap,
            P coordinates, PA and PB distances, etc)
     @param braAngularMomentum the angular momentum on bra side.
     @param ketAngularMomentum the angular momentum on ket side.
     @param threshold the primitive pairs screening threshold.
     */
    CGtoPairsBlock(const CMemBlock2D<int32_t>& contrPattern,
                   const CMemBlock2D<double>&  contrFactors,
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
    CGtoPairsBlock(const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock, const double threshold);

    /**
     Creates a GTOs pairs block object.

     @param gtoBlock the GTO block object for bra and ket sides.
     @param threshold the primitive pairs screening threshold.
     */
    CGtoPairsBlock(const CGtoBlock& gtoBlock, const double threshold);

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
     Creates vector of GTOs pairs objects by splitting GTOs pairs object.

     @param nodes the number of MPI ranks.
     @param cudaDevices the CUDA compute capable devices.
     @return the vector of GTOs pairs objects.
     */
    std::vector<CGtoPairsBlock> split(const int32_t       nodes,
                                      const CCudaDevices& cudaDevices = CCudaDevices()) const;

    /**
     Creates GTOs pairs object consisting from specific GTOs pairs object.

     @param iGtoPair the index of requested GTOs pair.
     @return the GTOs pairs object.
     */
    CGtoPairsBlock pick(const int32_t iGtoPair) const;

    /**
     Compresses other GTOs pairs block data into GTOs apirs block object without
     changing dimensions of GTOs pairs block object. Compression is performed
     using specified screening pattern.

     @param source the GTOs pairs block object.
     @param screeningPattern the screening pattern.
     @param nElements the number of elements in screening pattern.
     @return the number of contracted GTOs pairs in compressed GTOs pairs block
           object.
     */
    int32_t compress(const CGtoPairsBlock& source, const CMemBlock<int32_t>& screeningPattern, const int32_t nElements);

    /**
     Compresses other GTOs pairs block data into GTOs apirs block object without
     changing dimensions of GTOs pairs block object. Compression is performed
     using specified screening pattern.

     @param source the GTOs pairs block object.
     @param screeningPattern the screening pattern.
     @return the number of contracted GTOs pairs in compressed GTOs pairs block
             object.
     */
    int32_t compress(const CGtoPairsBlock& source, const CMemBlock<int32_t>& screeningPattern);

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
     Gets vector with Cartesian R(AB) = A - B distances for contracted GTOs
     pairs.

     @return vector of R(AB) distances.
     */
    CMemBlock2D<double> getDistancesAB() const;

    /**
     Gets vector with Cartesian R(AB) = A - B distances for specific set of
     contracted GTOs pairs.

     @param abDistances the vector of Cartesian distances.
     @param nContrPairs the number of contracted pairs.
     */
    void getDistancesAB(CMemBlock2D<double>& abDistances, const int32_t nContrPairs) const;

    /**
     Gets contants vector to Cartesian X coordinates of effective P center of
     contracted GTOs pairs.

     @return the vector of effective PX coordinates.
     */
    const double* getEffectiveCoordinatesPX() const;

    /**
     Gets contants vector to Cartesian Y coordinates of effective P center of
     contracted GTOs pairs.

     @return the vector of effective PY coordinates.
     */
    const double* getEffectiveCoordinatesPY() const;

    /**
     Gets contants vector to Cartesian Z coordinates of effective P center of
     contracted GTOs pairs.

     @return the vector of effective PZ coordinates.
     */
    const double* getEffectiveCoordinatesPZ() const;

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
     Gets maximum contraction depth of GTOs pair in GTOs pairs object.

     @return the maximum contraction depth.
     */
    int32_t getMaxContractionDepth() const;

    /**
     Gets number of primitive pairs in set of contracted pairs ([0..iContrPair]).

     @param iContrPair the index of lasr contracted pair.
     @return the number of primitive pairs.
     */
    int32_t getNumberOfPrimPairs(const int32_t iContrPair) const;

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
     Gets starting position of submatrix along bra side of full matrix for
     specific angular component of bra angular momentum.

     @param iComponent the angular momentum component.
     @return the starting position of submatrix.
     */
    int32_t getBraMatrixPosition(const int32_t iComponent) const;

    /**
     Gets starting position of submatrix along ket side of full matrix for
     specific angular component of ket angular momentum.

     @param iComponent the angular momentum component.
     @return the starting position of submatrix.
     */
    int32_t getKetMatrixPosition(const int32_t iComponent) const;

    /**
     Gets number of rows in submatrix along bra side of full matrix.

     @return the number of rows in submatrix.
     */
    int32_t getNumberOfRowsInBraMatrix() const;

    /**
     Gets number of rows in submatrix along ket side of full matrix.

     @return the number of rows in submatrix.
     */
    int32_t getNumberOfRowsInKetMatrix() const;

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
     Gets pair factors data from GTOs pair block object.

     @return <#return value description#>
     */
    CMemBlock2D<double> getPairFactors() const;

    /**
     Converts GTOs pairs block object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the GTOs pairs block object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CGtoPairsBlock& source);
};

#endif /* GtoPairsBlock_hpp */
