#ifndef FockMatrix_hpp
#define FockMatrix_hpp

#include <cstdint>

#include "Matrix.hpp"
#include "FockType.hpp"

/**
 Class CFockMatrix stores Fock matrix and provides set of methods
 for handling of Fock matrix data.

 @author Z. Rinkevicius
 */
class CFockMatrix
{
    /**
     The matrix to store Fock matrix values.
     */
    CMatrix _matrix;
    
    /**
     The scaling factor of exchange contribution to Fock matrix.
     */
    double _exc_scale;
    
    /**
     The type of Fock matrix.
     */
    fock_t _ftype;
    
   public:
    /**
     Creates an empty Fock matrix.
     */
    CFockMatrix();

    /**
     Creates a Fock matrix.

     @param matrix  the matrix to store Fock matrix values.
     @param exc_scale the scaling factor of exchange contribution to Fock matrix.
     @param ftype the Fock matrix type.
     */
    CFockMatrix(const CMatrix& matrix,
                const double   exc_scale,
                const fock_t   ftype);

    /**
     Creates a Fock matrix.

     @param other the Fock matrix to copy.
     */
    CFockMatrix(const CFockMatrix& other);

    /**
     Sets the scaling factor of exchange contribution to Fock matrix type.

     @param exc_scale the scalling factor of exchange contribution to Fock matrix type.
     */
    auto setExchangeScale(const double exc_scale) -> void;
    
    /**
     Sets the type of Fock matrix.

     @param ftype the type of Fock matrix.
     */
    auto setFockType(const fock_t ftype) -> void;

    /**
     Sets Fock matrix to zero.
     */
    auto zero() -> void;
    
    /**
     Gets the scaling factor of exchange contribution to Fock matrix type.

     @return the scalling factor of exchange contribution to Fock matrix type.
     */
    auto getExchangeScale() const -> double;
    
    /**
     Gets the type of Fock matrix.

     @return the type of Fock matrix
     */
    auto getFockType() const -> fock_t;
    
    /**
     Get vector of angular pairs from map of Fock submatrices.

     @return the vector of angular pairs.
     */
    auto getAngularPairs() const -> std::vector<T2Pair>;

    /**
     Get storage type of Fock matrix.

     @return the storage type of Fock matrix.
     */
    auto getStorageType() const -> mat_t;

    /**
     Get pointer to specific Fock submatrix.

     @param angpair the angular pair of Fock submatrix.
     @return the pointer to requested Fock submatrix.
     */
    auto getSubMatrix(const T2Pair& angpair) -> CSubMatrix*;

    /**
     Get constant pointer to specific Fock submatrix.

     @param angpair the angular pair of Fock submatrix.
     @return the constant pointer to requested Fock submatrix.
     */
    auto getSubMatrix(const T2Pair& angpair) const -> const CSubMatrix*;

    /**
     Checks if Fock submatrix with requested angular pair is stored in Fock matrix.

     @param angpair the angular pair of Fock submatrix.
     @return True if the Fock submatrix with this angular pair is stored, False otherwise.
     */
    auto isAngularOrder(const T2Pair& angpair) const -> bool;

    /**
     Gets number of rows in Fock matrix.

     @return the number of rows in Fock matrix.
     */
    auto getNumberOfRows() const -> int64_t;

    /**
     Gets number of columns in Fock matrix.

     @return the number of columns in Fock matrix.
     */
    auto getNumberOfColumns() const -> int64_t;

    /**
     Reconstructs full Fock matrix as submatrix with dimensions (0, 0, naos, naos).

     @return the full Fock matrix.
     */
    auto getFullMatrix() const -> CSubMatrix;
};

#endif /* FockMatrix_hpp */
