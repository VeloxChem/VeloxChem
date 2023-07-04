#ifndef GeomOneMatrices_hpp
#define GeomOneMatrices_hpp

#include <cstdint>
#include <map>
#include <vector>

#include "GeomKeys.hpp"
#include "Matrix.hpp"

/**
 Class CGeomOneMatrices stores dictionary of matrices associated with first order geometrical derivatives and provides set of methods
 for handling of dictionary data.

 @author Z. Rinkevicius
 */
class CGeomOneMatrices
{
    /**
     The vector of matrices.
     */
    std::map<TGeomPair, CMatrix*> _matrices;

   public:
    /**
     Creates an empty matrices.
     */
    CGeomOneMatrices();

    /**
     Creates a matrices.

     @param matrices the map of matrices.
     */
    CGeomOneMatrices(const std::map<TGeomPair, CMatrix>& matrices);

    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param atoms the vector of atoms.
     */
    CGeomOneMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms);

    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param keys the vector of geometrical key pairs.
     */
    CGeomOneMatrices(const CMatrix& matrix, const std::vector<TGeomPair>& keys);

    /**
     Creates a matrices.

     @param other the matrices to copy.
     */
    CGeomOneMatrices(const CGeomOneMatrices& other);

    /**
     Destroys a matrices.
     */
    ~CGeomOneMatrices();

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param key the key of matrix.
     */
    auto add(const CMatrix& matrix, const TGeomPair& key) -> void;

    /**
     Get vector of keys from map  of matrices.

     @return the vector  of keys.
     */
    auto getKeys() const -> std::vector<TGeomPair>;

    /**
     Get pointer to specific matrix.

     @param key the key of matrix.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const TGeomPair& key) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param key the key of matrix.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const TGeomPair& key) const -> const CMatrix*;
};

#endif /* GeomOneMatrices_hpp */
