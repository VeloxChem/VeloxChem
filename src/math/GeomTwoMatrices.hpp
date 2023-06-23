#ifndef GeomTwoMatrices_hpp
#define GeomTwoMatrices_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "Matrix.hpp"
#include "GeomKeys.hpp"

/**
 Class CGeomTwoMatrices stores dictionary of matrices associated with second order geometrical derivatives and provides set of methods
 for handling of dictionary data.

 @author Z. Rinkevicius
 */
class CGeomTwoMatrices
{
    /**
     The vector of matrices.
     */
    std::map<T2GeomKey, CMatrix*> _matrices;

   public:
    /**
     Creates an empty matrices.
     */
    CGeomTwoMatrices();

    /**
     Creates a matrices.

     @param matrices the map of matrices.
     */
    CGeomTwoMatrices(const std::map<T2GeomKey, CMatrix>& matrices);
    
    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param atoms the vector of atoms.
     */
    CGeomTwoMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms);
    
    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param keys the vector of geometrical key pairs.
     */
    CGeomTwoMatrices(const CMatrix& matrix, const std::vector<T2GeomKey>& keys);

    /**
     Creates a matrices.

     @param other the matrices to copy.
     */
    CGeomTwoMatrices(const CGeomTwoMatrices& other);

    /**
     Destroys a matrices.
     */
    ~CGeomTwoMatrices();

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param key the key of matrix.
     */
    auto add(const CMatrix& matrix, const T2GeomKey& key) -> void;

    /**
     Get vector of keys from map  of matrices.

     @return the vector  of keys.
     */
    auto getKeys() const -> std::vector<T2GeomKey>;

    /**
     Get pointer to specific matrix.

     @param key the key of matrix.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const T2GeomKey& key) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param key the key of matrix.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const T2GeomKey& key) const -> const CMatrix*;
};

#endif /* GeomTwoMatrices_hpp */
