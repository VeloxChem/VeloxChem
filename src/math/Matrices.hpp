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
    @param atom the atomic index of specific atom.
    @param label the label of specific atom.
    */
    auto add(const CMatrix& matrix, const int64_t atom, const std::string& label) -> void;

    /**
    Adds matrix to matrices.

    @param matrix the matrix to be added.
    @param atom_a the atomic index of atom A.
    @param label_a the label of atom A.
    @param atom_b the atomic index of atom B.
    @param label_b the label of atom B.
    */
    auto add(const CMatrix& matrix, const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) -> void;

    /**
    Adds matrix to matrices.

    @param matrix the matrix to be added.
    @param atom_a the atomic index of atom A.
    @param label_a the label of atom A.
    @param atom_b the atomic index of atom B.
    @param label_b the label of atom B.
    @param atom_c the atomic index of atom C.
    @param label_c the label of atom C.
    */
    auto add(const CMatrix&     matrix,
             const int64_t      atom_a,
             const std::string& label_a,
             const int64_t      atom_b,
             const std::string& label_b,
             const int64_t      atom_c,
             const std::string& label_c) -> void;

    /**
    Adds matrix to matrices.

    @param matrix the matrix to be added.
    @param atom_a the atomic index of atom A.
    @param label_a the label of atom A.
    @param atom_b the atomic index of atom B.
    @param label_b the label of atom B.
    @param atom_c the atomic index of atom C.
    @param label_c the label of atom C.
    @param atom_d the atomic index of atom D.
    @param label_d the label of atom D.
    */
    auto add(const CMatrix&     matrix,
             const int64_t      atom_a,
             const std::string& label_a,
             const int64_t      atom_b,
             const std::string& label_b,
             const int64_t      atom_c,
             const std::string& label_c,
             const int64_t      atom_d,
             const std::string& label_d) -> void;

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

     @param atom the atomic index of specific atom.
     @param label the label of specific atom.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom, const std::string& label) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param atom the atomic index of specific atom.
     @param label the label of specific atom.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom, const std::string& label) const -> const CMatrix*;

    /**
     Get pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t atom_a, const std::string& label_a, const int64_t atom_b, const std::string& label_b) const -> const CMatrix*;

    /**
     Get pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @param atom_c the atomic index of atom C.
     @param label_c the label of atom C.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t      atom_a,
                   const std::string& label_a,
                   const int64_t      atom_b,
                   const std::string& label_b,
                   const int64_t      atom_c,
                   const std::string& label_c) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @param atom_c the atomic index of atom C.
     @param label_c the label of atom C.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t      atom_a,
                   const std::string& label_a,
                   const int64_t      atom_b,
                   const std::string& label_b,
                   const int64_t      atom_c,
                   const std::string& label_c) const -> const CMatrix*;

    /**
     Get pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @param atom_c the atomic index of atom C.
     @param label_c the label of atom C.
     @param atom_d the atomic index of atom D.
     @param label_d the label of atom D.
     @return the pointer to requested matrix.
     */
    auto getMatrix(const int64_t      atom_a,
                   const std::string& label_a,
                   const int64_t      atom_b,
                   const std::string& label_b,
                   const int64_t      atom_c,
                   const std::string& label_c,
                   const int64_t      atom_d,
                   const std::string& label_d) -> CMatrix*;

    /**
     Get constant pointer to specific matrix.

     @param atom_a the atomic index of atom A.
     @param label_a the label of atom A.
     @param atom_b the atomic index of atom B.
     @param label_b the label of atom B.
     @param atom_c the atomic index of atom C.
     @param label_c the label of atom C.
     @param atom_d the atomic index of atom D.
     @param label_d the label of atom D.
     @return the constant pointer to requested matrix.
     */
    auto getMatrix(const int64_t      atom_a,
                   const std::string& label_a,
                   const int64_t      atom_b,
                   const std::string& label_b,
                   const int64_t      atom_c,
                   const std::string& label_c,
                   const int64_t      atom_d,
                   const std::string& label_d) const -> const CMatrix*;
};

#endif /* Matrices_hpp */
