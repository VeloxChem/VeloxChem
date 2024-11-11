#ifndef TensorComponents_hpp
#define TensorComponents_hpp

#include <array>
#include <cstddef>
#include <numeric>

namespace tensor {  // tensor

/// @brief Gets compound number of Cartesian components of canonical tensors array.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @return The number of Cartesian components.
template <size_t N>
auto
number_of_cartesian_components(const std::array<int, N> &orders) -> int
{
    return std::accumulate(orders.begin(), orders.end(), 1, [](const int prod, const int order) { return prod * (order + 1) * (order + 2) / 2; });
}

/// @brief Gets compound number of spherical components of canonical tensors array.
/// @tparam N The size of canonical tensors array.
/// @param orders The array of cononical tensor orders.
/// @return The number of spherical components.
template <size_t N>
auto
number_of_spherical_components(const std::array<int, N> &orders) -> int
{
    return std::accumulate(orders.begin(), orders.end(), 1, [](const int prod, const int order) { return prod * (2 * order + 1); });
}

}  // namespace tensor

#endif /* TensorComponents_hpp */
