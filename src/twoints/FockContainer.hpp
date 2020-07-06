//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef FockContainer_hpp
#define FockContainer_hpp

#include <cstdint>
#include <vector>

#include "FockSubMatrix.hpp"
#include "AOFockMatrix.hpp"

/**
 Class CFockContainer stores vector of Fock submatrices and provides set of
 methods for handling of Fock submatrices data.
 
 @author Z. Rinkevicius
 */
class CFockContainer
{
    /**
     The vector of Fock submatrices.
     */
    std::vector<CFockSubMatrix> _subFockMatrices;
    
public:
    
    /**
     Creates an empty Fock container object.
     */
    CFockContainer();
    
    /**
     Creates a Fock container object.
     
     @param subFockMatrices the vector of Fock submatrices.
     */
    CFockContainer(const std::vector<CFockSubMatrix>& subFockMatrices);
    
    /**
     Creates a Fock container object.
     
     @param aoFockMatrix the pointer to AO Fock matrix.
     @param braGtoPairsBlock the GTOs pairsblock on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    CFockContainer(const CAOFockMatrix*  aoFockMatrix,
                   const CGtoPairsBlock& braGtoPairsBlock,
                   const CGtoPairsBlock& ketGtoPairsBlock);
    
    /**
     Creates a Fock container object by copying other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer(const CFockContainer& source);
    
    /**
     Creates a Fock container object by moving other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer(CFockContainer&& source) noexcept;
    
    /**
     Destroys a Fock container object.
     */
    ~CFockContainer();
    
    /**
     Assigns a Fock container object by copying other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer& operator=(const CFockContainer& source);
    
    /**
     Assigns a Fock container object by moving other Fock container object.
     
     @param source the Fock container object.
     */
    CFockContainer& operator=(CFockContainer&& source) noexcept;
    
    /**
     Compares Fock container object with other Fock container object.
     
     @param other the Fock container object.
     @return true if Fock container objects are equal, false otherwise.
     */
    bool operator==(const CFockContainer& other) const;
    
    /**
     Compares Fock container object with other Fock container object.
     
     @param other the Fock container object.
     @return true if Fock container objects are not equal, false otherwise.
     */
    bool operator!=(const CFockContainer& other) const;
    
    /**
     Accumulates partial Fock matrices data stored in Fock container into AO
     Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     */
    void accumulate(CAOFockMatrix* aoFockMatrix);
    
    /**
     Gets pointer to partial data of specific submatrix with requested angular
     component for requested AO Fock matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @param iSubMatrix the index of submatrix.
     @param iComponent the angular component of submatrix.
     @return the pointer to partial submatrix data.
     */
    double* getSubMatrixData(const int32_t iFockMatrix,
                             const int32_t iSubMatrix,
                             const int32_t iComponent);
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center A on bra side for requested AO Fock matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsA(const int32_t iFockMatrix) const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center B on bra side for requested AO Fock matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsB(const int32_t iFockMatrix) const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center C on ket side for requested AO Fock matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsC(const int32_t iFockMatrix) const;
    
    /**
     Get constant pointer to starting positions of submatrices data along
     center D on ket side for requested AO Fock matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the pointer to starting positions.
     */
    const int32_t* getStartPositionsD(const int32_t iFockMatrix) const;
    
    /**
     Gets submatrix dimensions along center A on bra side for requested AO Fock
     matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsA(const int32_t iFockMatrix) const;
    
    /**
     Gets submatrix dimensions along center B on bra side for requested AO Fock
     matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsB(const int32_t iFockMatrix) const;
    
    /**
     Gets submatrix dimensions along center C on ket side for requested AO Fock
     matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsC(const int32_t iFockMatrix) const;
    
    /**
     Gets submatrix dimensions along center D on ket side for requested AO Fock
     matrix.
     
     @param iFockMatrix the index of AO Fock matrix.
     @return the dimensions of submatrix.
     */
    int32_t getDimensionsD(const int32_t iFockMatrix) const;
    
    /**
     Converts Fock container object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the Fock container object.
     */
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CFockContainer& source);
};

#endif /* FockContainer_hpp */
