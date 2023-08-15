#ifndef FockMatrices_hpp
#define FockMatrices_hpp

#include <cstdint>
#include <vector>

#include "FockMatrix.hpp"

/**
 Class CFockMatrices stores vector of Fock matrices and provides set of methods
 for handling of Fock matrices data.

 @author Z. Rinkevicius
 */
class CFockMatrices
{
    /**
     The vector of Fock matrices.
     */
    std::vector<CFockMatrix*> _matrices;

   public:
    /**
     Creates an empty Fock matrices.
     */
    CFockMatrices();

    /**
     Creates a Fock matrices.

     @param matrices the vector of Fock matrices.
     */
    CFockMatrices(const std::vector<CFockMatrix>& matrices);

    /**
     Creates a Fock matrices.

     @param other the Fock matrices to copy.
     */
    CFockMatrices(const CFockMatrices& other);

    /**
     Destroys a Fock matrices.
     */
    ~CFockMatrices();

    /**
     Adds Fock matrix to matrices.

     @param matrix the Fock matrix to be added.
     */
    auto add(const CFockMatrix& matrix) -> void;

    /**
     Sets all Fock matrices to zero.

     */
    auto zero() -> void;

    /**
     Get pointer to specific Fock matrix.

     @param index the index of Fock matrix.
     @return the pointer to requested Fock matrix.
     */
    auto getMatrix(const int64_t index) -> CFockMatrix*;

    /**
     Get constant pointer to specific Fock matrix.

     @param index the index of Fock matrix.
     @return the constant pointer to requested Fock matrix.
     */
    auto getMatrix(const int64_t index) const -> const CFockMatrix*;

    /**
     Gets number of Fock matrices.

     @return the number of Fock matrices.
    */
    auto getNumberOfMatrices() const -> int64_t;
};

#endif /* FockMatrices_hpp */
