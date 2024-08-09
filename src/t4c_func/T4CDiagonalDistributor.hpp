#ifndef T4CDiagonalDistributor_hpp
#define T4CDiagonalDistributor_hpp

#include <vector>
#include <array>

#include "SimdArray.hpp"

/// Class CT4CDiagonalDistributor provides methods for distributing vector of maximum values of integral shells.
class CT4CDiagonalDistributor
{
    /// The vector of maximum values of integral shells.
    double* _max_values;
    
   public:
    
    /// Creates a maximum values of integral shells distributor.
    /// - Parameter max_values : the vector of maximum values.
    CT4CDiagonalDistributor(double* max_values);
    
    /// Distributes maximum values of integral shells  into vector.
    /// - Parameter max_integrals: the vector of integral values.
    /// - Parameter gto_range: the index of the range [ket_first, ket_last) of integrals.
    auto distribute(      std::vector<double>& max_integrals,
                    const std::array<int, 2>&  gto_range) -> void;
};


#endif /* T4CDiagonalDistributor_hpp */
