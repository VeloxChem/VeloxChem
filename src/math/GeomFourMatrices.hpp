#ifndef GeomFourMatrices_hpp
#define GeomFourMatrices_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "GeomKeys.hpp"
#include "Matrix.hpp"

/**
 Class CGeomFourMatrices stores dictionary of matrices associated with fourth order geometrical derivatives and provides set of methods
 for handling of dictionary data.

 @author Z. Rinkevicius
 */
class CGeomFourMatrices
{
    /**
     The vector of matrices.
     */
    std::map<T4GeomKey, CMatrix*> _matrices;

   public:
    /**
     Creates an empty matrices.
     */
    CGeomFourMatrices();

    /**
     Creates a matrices.

     @param matrices the map of matrices.
     */
    CGeomFourMatrices(const std::map<T4GeomKey, CMatrix>& matrices);

    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param atoms the vector of atoms.
     */
    CGeomFourMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms);

    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param keys the vector of geometrical key pairs.
     */
    CGeomFourMatrices(const CMatrix& matrix, const std::vector<T4GeomKey>& keys);

    /**
     Creates a matrices.

     @param other the matrices to copy.
     */
    CGeomFourMatrices(const CGeomFourMatrices& other);

    /**
     Destroys a matrices.
     */
    ~CGeomFourMatrices();

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param key the key of matrix.
     */
    auto add(const CMatrix& matrix, const T4GeomKey& key) -> void;

    /**
     Get vector of keys from map  of matrices.

     @return the vector  of keys.
     */
    auto getKeys() const -> std::vector<T4GeomKey>;

    /**
     Get pointer to specific matrix.

     @param key the key of matrix.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const T4GeomKey& key) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param key the key of matrix.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const T4GeomKey& key) const -> const CMatrix*;
};

#endif /* GeomFourMatrices_hpp */
