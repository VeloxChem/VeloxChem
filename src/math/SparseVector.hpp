//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SparseVector_hpp
#define SparseVector_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock.hpp"

/**
 Class CSparseVector stores sparse vector (zero-based indexing scheme) and
 provides set of methods for manipulating sparse vector data.
 
 @author Z. Rinkevicius
 */
class CSparseVector
{
    /**
     The values of vector elements.
     */
    CMemBlock<double> _values;
    
    /**
     The vector of element indexes.
     */
    CMemBlock<int32_t> _indexes;
    
    /**
     The maximum number of elements in values vector.
     */
    int32_t _nMaxElements;
    
    /**
     The current number of elements in values vector.
     */
    int32_t _nElements;
    
    /**
     The cut-off threshold for small elements in vector.
     */
    double _threshold;
    
public:
    
    /**
     Creates an empty sparse vector object.
     */
    CSparseVector();
    
    /**
     Creates a sparse vector object from vector of values and vector of indexes.
     
     @param values the vector of matrix elements.
     @param indexes the vector of indexes.
     @param threshold the cut-off threshold for small elements in vector.
     */
    CSparseVector(const std::vector<double>&  values,
                  const std::vector<int32_t>& indexes,
                  const double                threshold);
    
    /**
     Creates an empty sparse vector object with specific number of elements.
     
     @param nElements the number of elements in vector.
     @param threshold the cut-off thrshold for small matrix elements.
     */
    CSparseVector(const int32_t nElements,
                  const double  threshold);
    
    /**
     Creates a sparse vector object by copying other sparse vector object.
     
     @param source the sparse vector object.
     */
    CSparseVector(const CSparseVector& source);
    
    /**
     Creates a sparse vector object by moving other sparse vector object.
     
     @param source the sparse vector object.
     */
    CSparseVector(CSparseVector&& source) noexcept;
    
    /**
     Destroys a sparse vector object.
     */
    ~CSparseVector();
    
    /**
     Assigns a sparse vector object by copying other sparse vector object.
     
     @param source the sparse vector object.
     */
    CSparseVector& operator=(const CSparseVector& source);
    
    /**
     Assigns a sparse vector object by moving other sparse vector object.
     
     @param source the sparse vector object.
     */
    CSparseVector& operator=(CSparseVector&& source) noexcept;
    
    /**
     Compares sparse vector object with other sparse vector object.
     
     @param other the sparse vector object.
     @return true if sparse vector objects are equal, false otherwise.
     */
    bool operator==(const CSparseVector& other) const;
    
    /**
     Compares sparse vector object with other sparse vector object.
     
     @param other the sparse vector object.
     @return true if sparse vector objects are not equal, false otherwise.
     */
    bool operator!=(const CSparseVector& other) const;
    
    /**
     Converts sparse vector object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the sparse vector object.
     */
    friend std::ostream& operator<<(      std::ostream&  output,
                                    const CSparseVector& source);
};

#endif /* SparseVector_hpp */
