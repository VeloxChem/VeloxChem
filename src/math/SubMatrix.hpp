#ifndef SubMatrix_hpp
#define SubMatrix_hpp

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "T4Index.hpp"

/**
 Class CSubMatrix stores submatrix data and provides set of methods
 for handling of submatrix data.

 @author Z. Rinkevicius
 */
class CSubMatrix
{
    /**
     The pointer to memory containing submatrix values.
     */
    double* _values;

    /**
     The dimensions of submatrix: row and column offsets and dimensions.
     */
    T4Index _dimensions;

   public:
    /**
     Creates an empty submatrix.
     */
    CSubMatrix();

    /**
     Creates a submatrix.

     @param dimensions the dimensions of submatrix.
     */
    CSubMatrix(const T4Index& dimensions);

    /**
     Creates a submatrix.

     @param values the vector with submatrix values.
     @param dimensions the dimensions of submatrix.
     */
    CSubMatrix(const std::vector<double>& values, const T4Index& dimensions);

    /**
     Creates a submatrix.

     @param other the submatrix to be copied.
     */
    CSubMatrix(const CSubMatrix& other);

    /**
     Destroys a matrix.
     */
    ~CSubMatrix();

    /**
     Gets reference to specific submatirx element.

     @param irow the row index in supermatrix or submatrix indexing scheme.
     @param icol the column index in supermatrix or submatrix indexing scheme.
     @param supmat the flag to enable supermatrix indexing scheme.
     @return the reference to specific submatrix element.
     */
    inline auto
    at(const int64_t irow, const int64_t icol, const bool supmat) -> double&
    {
        if (supmat)
        {
            return _values[(irow - _dimensions[0]) * _dimensions[3] + icol - _dimensions[1]];
        }
        else
        {
            return _values[irow * _dimensions[3] + icol];
        }
    }

    /**
     Gets reference to specific submatirx element.

     @param irow the row index in supermatrix or submatrix indexing scheme.
     @param icol the column index in supermatrix or submatrix indexing scheme.
     @param supmat the flag to enable supermatrix indexing scheme.
     @return the reference to specific submatrix element.
     */
    inline auto
    at(const int64_t irow, const int64_t icol, const bool supmat) const -> const double&
    {
        if (supmat)
        {
            return _values[(irow - _dimensions[0]) * _dimensions[3] + icol - _dimensions[1]];
        }
        else
        {
            return _values[irow * _dimensions[3] + icol];
        }
    }

    /**
     Set offset of submatrix.

     @param row_offset the row offset.
     @param col_offset the column offset.
     */
    auto setOffsets(const int64_t row_offset, const int64_t col_offset) -> void;

    /**
     Set values of submatrix.

     @param values the vector with submatrix values.
     */
    auto setValues(const std::vector<double>& values) -> void;

    /**
     Set values of submatrix to zero.
     */
    auto zero() -> void;
    
    /**
     Gets dimensions of submatrix.

     @return The dimensions of submatrix.
     */
    auto getDimensions() const -> T4Index;

    /**
     Gets values of all elements in submatrix.

     @return the vector with submatrix values.
     */
    auto getValues() const -> std::vector<double>;

    /**
     Gets raw pointer to submatrix data.

     @return the raw pointer to submatrix data.
     */
    auto getData() -> double*;

    /**
     Gets raw constant pointer to submatrix data.

     @return the raw constant pointer to submatrix data.
     */
    auto getData() const -> const double*;

    /**
     Gets offset of rows in submatrix.

     @return The  offset of rows in submatrix.
     */
    inline auto
    getOffsetOfRows() const -> int64_t
    {
        return _dimensions[0];
    }

    /**
     Gets offset of columns in submatrix.

     @return The offset of columns in submatrix.
     */
    inline auto
    getOffsetOfColumns() const -> int64_t
    {
        return _dimensions[1];
    }

    /**
     Gets number of rows in submatrix.

     @return The number of rows in submatrix.
     */
    inline auto
    getNumberOfRows() const -> int64_t
    {
        return _dimensions[2];
    }

    /**
     Gets number of columns in submatrix.

     @return The number of columns in submatrix.
     */
    inline auto
    getNumberOfColumns() const -> int64_t
    {
        return _dimensions[3];
    }
};

#endif /* SubMatrix_hpp */
