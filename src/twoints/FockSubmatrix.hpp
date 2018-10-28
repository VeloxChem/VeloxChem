//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef FockSubmatrix_hpp
#define FockSubmatrix_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "MemBlock.hpp"
#include "GtoPairsBlock.hpp"
#include "FockMatrixType.hpp"

/**
 Class CFockSubMatrix stores part of AO Fock matrix and provides set of methods
 for handling of partial AO Fock matrix data.
 
 @author Z. Rinkevicius
 */
class CFockSubMatrix
{
    /**
     The vector of partial AO Fock submatrices.
     */
    std::vector<CMemBlock2D<double>> _subFockMatrices;
    
    /**
     The starting positions of submatrix components on centre A on bra side.
     */
    CMemBlock<int32_t> _startPositionsA;
    
    /**
     The starting positions of submatrix components on centre B on bra side.
     */
    CMemBlock<int32_t> _startPositionsB;
    
    /**
     The starting positions of submatrix components on centre C on ket side.
     */
    CMemBlock<int32_t> _startPositionsC;
    
    /**
     The starting positions of submatrix components on centre D on ket side.
     */
    CMemBlock<int32_t> _startPositionsD;
    
    /**
     The dimensions of submatrix along centre A on bra side.
     */
    int32_t _dimSubMatrixA;
    
    /**
     The dimensions of submatrix along centre B on bra side.
     */
    int32_t _dimSubMatrixB;
    
    /**
     The dimensions of submatrix along centre C on ket side.
     */
    int32_t _dimSubMatrixC;
    
    /**
     The dimensions of submatrix along centre D on ket side.
     */
    int32_t _dimSubMatrixD;
    
    /**
     Allocates and initializes to zero Fock submatrices data according to type
     of associated AO Fock matrix.

     @param fockType the type of AO Fock matrix.
     @param nComponentsA the number of angular components on center A.
     @param nComponentsB the number of angular components on center B.
     @param nComponentsC the number of angular components on center C.
     @param nComponentsD the number of angular components on center D.
     */
    void _allocSubMatrices(const fockmat fockType,
                           const int32_t nComponentsA,
                           const int32_t nComponentsB,
                           const int32_t nComponentsC,
                           const int32_t nComponentsD);
    
public:
    
    /**
     Creates an empty AO Fock submatrix object.
     */
    CFockSubMatrix();
    
    /**
     Creates a AO Fock submatrix object.
     
     @param subFockMatrices the vector of partial AO Fock submatrices.
     @param startPositionsA the starting positions of centre A on bra side.
     @param startPositionsB the starting positions of centre B on bra side.
     @param startPositionsC the starting positions of centre C on ket side.
     @param startPositionsD the starting positions of centre D on ket side.
     @param dimSubMatrixA the dimensions of aling centre A along bra side.
     @param dimSubMatrixB the dimensions of aling centre B along bra side.
     @param dimSubMatrixC the dimensions of aling centre C along ket side.
     @param dimSubMatrixD the dimensions of aling centre D along ket side.
     */
    CFockSubMatrix(const std::vector<CMemBlock2D<double>>& subFockMatrices,
                   const CMemBlock<int32_t>&               startPositionsA,
                   const CMemBlock<int32_t>&               startPositionsB,
                   const CMemBlock<int32_t>&               startPositionsC,
                   const CMemBlock<int32_t>&               startPositionsD,
                   const int32_t                           dimSubMatrixA,
                   const int32_t                           dimSubMatrixB,
                   const int32_t                           dimSubMatrixC,
                   const int32_t                           dimSubMatrixD);
    
    /**
     Creates a AO Fock submatrix object.
     
     @param braGtoPairsBlock the GTOs pairsblock on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param fockType the type of Fock matrix associated with Fock submatrix.
     */
    CFockSubMatrix(const CGtoPairsBlock& braGtoPairsBlock,
                   const CGtoPairsBlock& ketGtoPairsBlock,
                   const fockmat         fockType);
    
    /**
     Creates a AO Fock submatrix object by copying other AO Fock submatrix
     object.
     
     @param source the AO Fock submatrix object.
     */
    CFockSubMatrix(const CFockSubMatrix& source);
    
    /**
     Creates a AO Fock submatrix object by moving other AO Fock submatrix
     object.
     
     @param source the AO Fock submatrix object.
     */
    CFockSubMatrix(CFockSubMatrix&& source) noexcept;
    
    /**
     Destroys a AO Fock submatrix object.
     */
    ~CFockSubMatrix();
    
    /**
     Assigns a AO Fock submatrix object by copying other AO Fock submatrix
     object.
     
     @param source the AO Fock submatrix object.
     */
    CFockSubMatrix& operator=(const CFockSubMatrix& source);
    
    /**
     Assigns a AO Fock submatrix object by moving other AO Fock submatrix
     object.
     
     @param source the AO Fock submatrix object.
     */
    CFockSubMatrix& operator=(CFockSubMatrix&& source) noexcept;
    
    /**
     Compares AO Fock submatrix object with other AO Fock submatrix object.
     
     @param other the AO Fock submatrix object.
     @return true if AO Fock submatrix objects are equal, false otherwise.
     */
    bool operator==(const CFockSubMatrix& other) const;
    
    /**
     Compares AO Fock submatrix object with other AO Fock submatrix object.
     
     @param other the AO Fock submatrix object.
     @return true if AO Fock submatrix objects are not equal, false otherwise.
     */
    bool operator!=(const CFockSubMatrix& other) const;
    
    /**
     Gets pointer to partial data of specific submatrix with requested angular
     component.

     @param iMatrix the index of submatrix.
     @param iComponent the angular component of submatrix.
     @return the pointer to partial submatrix data.
     */
    double* getSubMatrixData(const int32_t iMatrix,
                             const int32_t iComponent);
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center A on bra side.

     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsA() const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center B on bra side.
     
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsB() const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center C on ket side.
     
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsC() const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center D on ket side.
     
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsD() const;
    
    /**
     Gets submatrix dimensions along center A on bra side.

     @return the dimensions of submatrix.
     */
    int32_t getDimensionsA() const;
    
    /**
     Gets submatrix dimensions along center B on bra side.
     
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsB() const;
    
    /**
     Gets submatrix dimensions along center C on ket side.
     
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsC() const;
    
    /**
     Gets submatrix dimensions along center D on ket side.
     
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsD() const;
    
    
    /**
     Converts AO Fock submatrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the AO Fock submatrix object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CFockSubMatrix& source);
};

#endif /* FockSubmatrix_hpp */
