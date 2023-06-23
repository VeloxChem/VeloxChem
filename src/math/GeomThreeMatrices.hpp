#ifndef GeomThreeMatrices_hpp
#define GeomThreeMatrices_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "Matrix.hpp"
#include "GeomKeys.hpp"

/**
 Class CGeomThreeMatrices stores dictionary of matrices associated with third order geometrical derivatives and provides set of methods
 for handling of dictionary data.

 @author Z. Rinkevicius
 */
class CGeomThreeMatrices
{
    /**
     The vector of matrices.
     */
    std::map<T3GeomKey, CMatrix*> _matrices;

   public:
    /**
     Creates an empty matrices.
     */
    CGeomThreeMatrices();

    /**
     Creates a matrices.

     @param matrices the map of matrices.
     */
    CGeomThreeMatrices(const std::map<T3GeomKey, CMatrix>& matrices);
    
    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param atoms the vector of atoms.
     */
    CGeomThreeMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms);
    
    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param keys the vector of geometrical key pairs.
     */
    CGeomThreeMatrices(const CMatrix& matrix, const std::vector<T3GeomKey>& keys);

    /**
     Creates a matrices.

     @param other the matrices to copy.
     */
    CGeomThreeMatrices(const CGeomThreeMatrices& other);

    /**
     Destroys a matrices.
     */
    ~CGeomThreeMatrices();

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param key the key of matrix.
     */
    auto add(const CMatrix& matrix, const T3GeomKey& key) -> void;

    /**
     Get vector of keys from map  of matrices.

     @return the vector  of keys.
     */
    auto getKeys() const -> std::vector<T3GeomKey>;

    /**
     Get pointer to specific matrix.

     @param key the key of matrix.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const T3GeomKey& key) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param key the key of matrix.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const T3GeomKey& key) const -> const CMatrix*;
};

#endif /* GeomThreeMatrices_hpp */
