//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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

    /**
     Determines if given row is last initialized row in sparse matrix.

     @param iRow the index of row.
     @return true if specific row is last initialized row in sparse matrix,
             false - otherwise.
     */
    bool _isLastRow(const int32_t iRow) const;

    /**
     Gets number of additional row added to sparse matrix.

     @param iRow the index of row.
     @return the number of additional rows.
     */
    int32_t _getAdditionalRows(const int32_t iRow) const;

    /**
     Sets initial maximum number of elements in data buffers of sparse matrix.

     @return the maximum number of elements.
     */
    int32_t _setMaxNumberOfElements() const;

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
     @param threshold the cut-off threshold for small matrix elements.
     */
    CSparseMatrix(const std::vector<double>&  values,
                  const std::vector<int32_t>& rows,
                  const std::vector<int32_t>& columns,
                  const int32_t               nRows,
                  const int32_t               nColumns,
                  const double                threshold);

    /**
     Creates an empty sparse matrix object with specific number of rows and
     columns.

     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     @param threshold the cut-off threshold for small matrix elements.
     */
    CSparseMatrix(const int32_t nRows, const int32_t nColumns, const double threshold);

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

     @param source the sparse matrix object.
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
     Appends new row (values and indexes) of specified size to sparse matrix.
     NOTE: The index of new row must be greater than index of last initialized
     row in sparse matrix.

     @param rowValues the vector of values.
     @param rowColumns the vector of column indexes.
     @param nElementsInRow the number of elements in added row.
     @param iRow the index of row.
     */
    void append(const CMemBlock<double>& rowValues, const CMemBlock<int32_t>& rowColumns, const int32_t nElementsInRow, const int32_t iRow);

    /**
     Appends new row (values and indexes) of specified size to sparse matrix.
     NOTE: The index of new row must be greater than index of last initialized
     row in sparse matrix.

     @param rowValues the vector of values.
     @param rowColumns the vector of column indexes.
     @param iRow the index of row.
     */
    void append(const CMemBlock<double>& rowValues, const CMemBlock<int32_t>& rowColumns, const int32_t iRow);

    /**
     Optimizes memory footprint of buffers used to store sparse matrix elements
     values and indexing information.
     */
    void optimize_storage();

    /**
     Check if sparse matrix data is stored im optimal way in memory.

     @return true if sparse matrix is stored in optimal way, false otherwise.
     */
    bool isOptimizedStorage() const;

    /**
     Gets number of rows in sparse matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;

    /**
     Gets number of columns in sparse matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;

    /**
     Gets number of elements in sparse matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;

    /**
     Gets number of elements in specific row of sparse matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements(const int32_t iRow) const;

    /**
     Gets cut-off threshold of matrix elements.

     @return the cut-off threshold.
     */
    double getThreshold() const;

    /**
     Computes sparsity number for sparse matrix.

     @return the sparsity number.
     */
    double getSparsity() const;

    /**
     Get constant pointer to selected row in sparse matrix.

     @param iRow the index of row.
     @return the constant pointer to first element in row.
     */
    const double* row(const int32_t iRow) const;

    /**
     Get pointer to selected row in sparse matrix.

     @param iRow the index of row.
     @return the pointer to first element in row.
     */
    double* row(const int32_t iRow);

    /**
     Get constant pointer to column indexes of selected row in sparse matrix.

     @param iRow the index of row.
     @return the constant pointer to first column index of selected row.
     */
    const int32_t* indexes(const int32_t iRow) const;

    /**
     Get pointer to column indexes of selected row in sparse matrix.

     @param iRow the index of row.
     @return the pointer to first column index of selected row.
     */
    int32_t* indexes(const int32_t iRow);

    /**
     Gets constant pointer to sparse matrix data.

     @return the constant pointer to first data element in sparse matrix.
     */
    const double* values() const;

    /**
     Gets pointer to sparse matrix data.

     @return the pointer to first data element in sparse matrix.
     */
    double* values();

    /**
     Gets constant pointer to indexes of rows in sparse matrix.

     @return the constant pointer to first index of rows.
     */
    const int32_t* rows() const;

    /**
     Gets constant pointer to indexes of columns in sparse matrix.

     @return the constant pointer to first index of columns.
     */
    const int32_t* columns() const;

    /**
     Converts sparse matrix object to text output and insert it into output
     text stream.

     @param output the output text stream.
     @param source the sparse matrix object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CSparseMatrix& source);
};

#endif /* SparseMatrix_hpp */
