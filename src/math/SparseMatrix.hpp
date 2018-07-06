//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SparseMatrix_hpp
#define SparseMatrix_hpp

#include <cstdint>
#include <ostream>
#include <vector>

#include "MemBlock.hpp"

/**
 Class CSparseMatrix stores sparse matrix in coordinate format (zero-based
 indexing scheme) and provides set of methods for manipulating sparse matrix
 data.
 
 @author Z. Rinkevicius
 */
class CSparseMatrix
{
    /**
     The number of rows.
     */
    int32_t _nRows;
    
    /**
     The number of columns.
     */
    int32_t _nColumns;
    
    /**
     The vector of matrix element values.
     */
    CMemBlock<double> _values;
    
    /**
     The vector of matrix element row indexes.
     */
    CMemBlock<int32_t> _rows;
    
    /**
     The vector of matrix element column indexes.
     */
    CMemBlock<int32_t> _columns;
    
    /**
     The maximum number of elements in values vector.
     */
    int32_t _nMaxElements;
    
    /**
     The current number of elements in values vector.
     */
    int32_t _nElements;
    
    /**
     The vector of row starting positions in values vector.
     */
    CMemBlock<int32_t> _rowPositions;
    
    /**
     The vector of row sizes in values vector.
     */
    CMemBlock<int32_t> _rowSizes;
    
    /**
     The cut-off threshold for small matrix elements.
     */
    double _threshold;
    
    /**
     Sets rows access pattern (position and number of elements) for each row in
     matrix.
     */
    void _setAccessPattern();
    
public:
    
    /**
     Creates an empty sparse matrix object.
     */
    CSparseMatrix();
    
    /**
     Creates a sparse matrix object from vector of values and vectors of row and
     column indexes.
     
     @param values the vector of matrix elements (Coordinate format).
     @param rows the vector of row indexes (Coordinate format).
     @param columns the vector of column indexes (Coordinate format).
     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     @param threshold the cut-off thrshold for small matrix elements.
     */
    CSparseMatrix(const std::vector<double>&  values,
                  const std::vector<int32_t>& rows,
                  const std::vector<int32_t>& columns,
                  const int32_t               nRows,
                  const int32_t               nColumns,
                  const double                threshold);
    
    /**
     Creates a sparse matrix object by copying other sparse matrix object.
     
     @param source the sparse matrix object.
     */
    CSparseMatrix(const CSparseMatrix& source);
    
    /**
     Creates a sparse matrix object by moving other sparse matrix object.
     
     @param source the sparse matrix object.
     */
    CSparseMatrix(CSparseMatrix&& source) noexcept;
    
    /**
     Destroys a sparse matrix object.
     */
    ~CSparseMatrix();
    
    /**
     Assigns a sparse matrix object by copying other sparse matrix object.
     
     @param source the sparse matrix object.
     */
    CSparseMatrix& operator=(const CSparseMatrix& source);
    
    /**
     Assigns a sparse matrix object by moving other sparse matrix object.
     
     @param source the basis function object.
     */
    CSparseMatrix& operator=(CSparseMatrix&& source) noexcept;
    
    /**
     Compares sparse matrix object with other sparse matrix object.
     
     @param other the sparse matrix object.
     @return true if sparse matrix objects are equal, false otherwise.
     */
    bool operator==(const CSparseMatrix& other) const;
    
    /**
     Compares sparse matrix object with other sparse matrix object.
     
     @param other the sparse matrix object.
     @return true if sparse matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CSparseMatrix& other) const;
    
    /**
     Converts sparse matrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the sparse matrix object.
     */
    friend std::ostream& operator<<(      std::ostream&  output,
                                    const CSparseMatrix& source);
    
};

#endif /* SparseMatrix_hpp */
