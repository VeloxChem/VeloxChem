#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <array>
#include <cmath>
#include <vector>
#include <ranges>

#include "CustomConstrains.hpp"

namespace mathfunc {  // mathfunc

/// @brief Compares two floating point values.
/// @tparam T The floating point type.
/// @param flhs The left-hand-side floating point value to compare.
/// @param frhs The rigth-hand-side floating point value to compare.
/// @param rtol The relative tollerance.
/// @param atol The absolute tollerance.
/// @return True if floating point values are equal, False otherwise.
template <FloatingPoint T>
inline auto
equal(const T flhs, const T frhs, const T rtol, const T atol) -> bool
{
    return std::fabs(flhs - frhs) < std::max(atol, rtol * std::max(std::fabs(flhs), std::fabs(frhs)));
}

/// @brief Gets upper triangular matrix linearized index (column wise scheme).
/// @param i The index of row in matrix.
/// @param j The index of collumn in matrix.
/// @return The linearized index.
template <Integral T>
inline auto
uplo_index(const T i, const T j) -> T
{
    return i + j * (j + 1) / 2;
}

/// @brief Counts number of elements in vector matching given selector.
/// @param values  The vector of values.
/// @param selector  The selector to march vector values.
/// @return The number of elements in vector matching given selector.
template <Integral T>
inline auto
count_elements_by_values(const std::vector<T>& values,
                         const T               selector) -> T
{
    return static_cast<T>(std::ranges::count(values, selector));
}

inline auto
batch_sizes(const int nElements, const int nodes) -> std::vector<int>
{
    int ave = nElements / nodes;

    int rem = nElements % nodes;

    std::vector<int> counts;

    for (int p = 0; p < nodes; p++)
    {
        counts.push_back((p < rem) ? (ave + 1) : ave);
    }

    return counts;
}

inline auto
batch_offsets(const int nElements, const int nodes) -> std::vector<int>
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    std::vector<int> displs;

    int index = 0;

    for (int p = 0; p < nodes; p++)
    {
        displs.push_back(index);

        index += counts[p];
    }

    return displs;
}

inline auto
batch_size(const int nElements, const int rank, const int nodes) -> int
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    return counts[rank];
}

inline auto
batch_offset(const int nElements, const int rank, const int nodes) -> int
{
    auto displs = mathfunc::batch_offsets(nElements, nodes);

    return displs[rank];
}

}  // namespace mathfunc

#endif /* MathFunc_hpp */
