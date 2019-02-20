//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DenseMatrix_hpp
#define DenseMatrix_hpp

#include <cstdint>
#include <ostream>
#include <vector>
#include <string>

#include "MpiFunc.hpp"
#include "MemBlock.hpp"

/**
 Class CDenseMatrix stores dense matrix in coordinate format (zero-based
 indexing scheme) and provides set of methods for manipulating dense matrix
 data.
 
 @author Z. Rinkevicius
 */
class CDenseMatrix
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
    
public:
    
    /**
     Creates an empty dense matrix object.
     */
    CDenseMatrix();
    
    /**
     Creates a dense matrix object from vector of values.
     
     @param values the vector of matrix elements.
     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const std::vector<double>& values,
                 const int32_t              nRows,
                 const int32_t              nColumns);
    
    /**
     Creates an empty dense matrix object with specific number of rows and
     columns.
     
     @param nRows the number of rows in matrix.
     @param nColumns the number of columns in matrix.
     */
    CDenseMatrix(const int32_t nRows,
                 const int32_t nColumns);
    
    /**
     Creates an empty square dense matrix object with specific number of rows.
     
     @param nRows the number of rows in matrix.
     */
    CDenseMatrix(const int32_t nRows);
    
    /**
     Creates a dense matrix object by copying other dense matrix object.
     
     @param source the dense matrix object.
     */
    CDenseMatrix(const CDenseMatrix& source);
    
    /**
     Creates a dense matrix object by moving other dense matrix object.
     
     @param source the dense matrix object.
     */
    CDenseMatrix(CDenseMatrix&& source) noexcept;
    
    /**
     Destroys a dense matrix object.
     */
    ~CDenseMatrix();
    
    /**
     Assigns a dense matrix object by copying other dense matrix object.
     
     @param source the dense matrix object.
     */
    CDenseMatrix& operator=(const CDenseMatrix& source);
    
    /**
     Assigns a dense matrix object by moving other dense matrix object.
     
     @param source the dense matrix object.
     */
    CDenseMatrix& operator=(CDenseMatrix&& source) noexcept;
    
    /**
     Compares dense matrix object with other dense matrix object.
     
     @param other the dense matrix object.
     @return true if dense matrix objects are equal, false otherwise.
     */
    bool operator==(const CDenseMatrix& other) const;
    
    /**
     Compares dense matrix object with other dense matrix object.
     
     @param other the dense matrix object.
     @return true if dense matrix objects are not equal, false otherwise.
     */
    bool operator!=(const CDenseMatrix& other) const;
    
    /**
     Sets all values in dense matrix to zero.
     */
    void zero();
    
    /**
     Creates transpose dense matrix.

     @return the transpose dense matrix.
     */
    CDenseMatrix transpose() const;
    
    /**
     Symmetrizes elements of square matrix: a_ij = a_ji = (a_ij + a_ji).
     */
    void symmetrize();
    
    /**
     Reduces dense matrix objects from all MPI process within domain of MPI
     communicator into dense matrix object on master node by summing them.
     
     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void reduce_sum(int32_t  rank,
                    int32_t  nodes,
                    MPI_Comm comm);
    
    /**
     Creates dense matrix object by slicing specified size submatrix at
     selected position from this dense matrix object.

     @param iRow the starting row of submatrix.
     @param iColumn the starting column of submatrix.
     @param nRows the number of rows in submatrix.
     @param nColumns the number of columns in submatrix.
     @return the dense matrix object.
     */
    CDenseMatrix slice(const int32_t iRow,
                       const int32_t iColumn,
                       const int32_t nRows,
                       const int32_t nColumns) const;
    
    /**
     Creates dense matrix from selected colums in this dense matrix object.

     @param iColumns the indexes of selected columns.
     @return the dense matrix.
     */
    CDenseMatrix selectByColumn(const std::vector<int32_t>& iColumns) const;
    
    /**
     Creates dense matrix from selected rows im this dense matrix object.

     @param iRows the indexes of selected rows.
     @return the dense matrix.
     */
    CDenseMatrix selectByRow(const std::vector<int32_t>& iRows) const;
    
    /**
     Gets number of rows in dense matrix.

     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in dense matrix.

     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets number of elements in dense matrix.

     @return the number of elements.
     */
    int32_t getNumberOfElements() const;
    
    /**
     Gets constant pointer to first element of dense matrix.

     @return the constant pointer to first element of dense matrix.
     */
    const double* values() const;
    
    /**
     Gets pointer to first element of dense matrix.
     
     @return the pointer to first element of dense matrix.
     */
    double* values();
    
    /**
     Gets constant pointer to first element of specific row in dense matrix.
     
     @return the constant pointer to first element of specific row.
     */
    const double* row(const int32_t iRow) const;
    
    /**
     Gets pointer to first element of specific row in dense matrix.
     
     @return the pointer to first element of specific row.
     */
    double* row(const int32_t iRow);

    /**
     Gets string representation of dense matrix object.

     @return the string representation.
     */
    std::string getString() const;
    
    /**
     Broadcasts density matrix object within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t  rank,
                   MPI_Comm comm);
    
    /**
     Converts dense matrix object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the dense matrix object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const CDenseMatrix& source);
};

#endif /* DenseMatrix_hpp */
