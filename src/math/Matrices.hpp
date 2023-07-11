#ifndef Matrices_hpp
#define Matrices_hpp

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "Matrix.hpp"

/**
 Class CMatrices stores dictionary of matrices and provides set of methods
 for handling of dictionary data.

 @author Z. Rinkevicius
 */
class CMatrices
{
    /**
     The vector of matrices.
     */
    std::map<int64_t, CMatrix*> _matrices;

    /**
     Converts string label to integer key of matrix.

     @param label the label of key.
     @return the integer key of matrix.
     */
    auto _to_key(const std::string& label) const -> int64_t;

   public:
    /**
     Creates an empty matrices.
     */
    CMatrices();

    /**
     Creates a matrices.

     @param matrices the map of matrices.
     */
    CMatrices(const std::map<int64_t, CMatrix>& matrices);

    /**
     Creates a matrices.

     @param matrix the matrix to be create map of matrices.
     @param keys the vector of matrix keys.
     */
    CMatrices(const CMatrix& matrix, const std::vector<int64_t>& keys);

    /**
     Creates a matrices.

     @param other the matrices to copy.
     */
    CMatrices(const CMatrices& other);

    /**
     Destroys a matrices.
     */
    ~CMatrices();

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param key the key of matrix.
     */
    auto add(const CMatrix& matrix, const int64_t key) -> void;

    /**
     Adds matrix to matrices.

     @param matrix the matrix to be added.
     @param label the label of key.
     */
    auto add(const CMatrix& matrix, const std::string& label) -> void;
    
    /**
    Adds matrix to matrices.

    @param matrix the matrix to be added.
    @param atom the atomic index.
    @param label the label of key.
    */
   auto add(const CMatrix& matrix, const int64_t atom, const std::string& label) -> void;

    /**
     Sets all matrices to zero.

     */
    auto zero() -> void;

    /**
     Get vector of keys from map  of matrices.

     @return the vector  of keys.
     */
    auto getKeys() const -> std::vector<int64_t>;

    /**
     Get pointer to specific matrix.

     @param key the key of matrix.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t key) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param key the key of matrix.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t key) const -> const CMatrix*;

    /**
     Get pointer to specific matrix.

     @param label the label of key.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const std::string& label) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param label the label of key.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const std::string& label) const -> const CMatrix*;
    
    /**
     Get pointer to specific matrix.

     @param atom the atomic index.
     @param label the label of key.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom, const std::string& label) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param atom the atomic index.
     @param label the label of key.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom, const std::string& label) const -> const CMatrix*;
};

#endif /* Matrices_hpp */
