#ifndef T4CDiagonalDistributor_hpp
#define T4CDiagonalDistributor_hpp

#include <cstddef>
#include <vector>
#include <utility>

#include "SimdArray.hpp"

/// @brief Class CT4CDiagonalDistributor provides methods for distributing vector of maximum values of integral shells.
class CT4CDiagonalDistributor
{
   public:
    
    /// @brief Creates a maximum values of integral shells distributor.
    /// @param max_values  The pointer to vector of maximum values.
    CT4CDiagonalDistributor(double* max_values);
    
    /// @brief The default copy constructor.
    /// @param other The distributor to be copied.
    CT4CDiagonalDistributor(const CT4CDiagonalDistributor& other) = delete;

    /// @brief The default move constructor.
    /// @param other The distributor to be moved.
    CT4CDiagonalDistributor(CT4CDiagonalDistributor&& other) noexcept = delete;

    /// @brief The default destructor.
    ~CT4CDiagonalDistributor() = default;

    /// @brief The default copy assignment operator.
    /// @param other The distributor to be copy assigned.
    /// @return The assigned distributor.
    auto operator=(const CT4CDiagonalDistributor& other) -> CT4CDiagonalDistributor& = delete;

    /// @brief The default move assignment operator.
    /// @param other The distributor to be move assigned.
    /// @return The assigned distributor.
    auto operator=(CT4CDiagonalDistributor&& other) noexcept -> CT4CDiagonalDistributor& = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are equal, False otherwise.
    auto operator==(const CT4CDiagonalDistributor& other) const -> bool = delete;

    /// @brief The equality operator.
    /// @param other The distributor to be compared.
    /// @return True if distributors are not equal, False otherwise.
    auto operator!=(const CT4CDiagonalDistributor& other) const -> bool = delete;
    
    /// @brief Distributes maximum values of integral shells  into vector.
    /// @param max_integrals The vector of integral values.
    /// @param gto_range The index of the range [ket_first, ket_last) of integrals.
    auto distribute(const std::vector<double>&        max_integrals,
                    const std::pair<size_t, size_t>&  gto_range) -> void;
    
private:
    
    /// @brief The vector of maximum values of integral shells.
    double* _max_values;
};


#endif /* T4CDiagonalDistributor_hpp */
