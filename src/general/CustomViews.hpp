#ifndef CustomViews_hpp
#define CustomViews_hpp

#include <ranges>
#include <utility>

#include "CustomConstrains.hpp"

namespace views {  // views

/// @brief Creates indices (i,j) view of rectangular matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @param ncolumns The number of columns in matrix.
/// @return The indices (i,j) view of rectangular matrix.
template <Integral T>
inline auto
rectangular(const T nrows, const T ncolumns) -> decltype(auto)
{
    auto indices = [=](const T index) { return std::make_pair(index / ncolumns, index % ncolumns); };

    return std::views::iota(T{0}, nrows * ncolumns) | std::views::transform(indices);
}

/// @brief Creates triangular indices (i,j) view of square matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @return The triangular indices (i,j) view of square matrix.
template <Integral T>
inline auto
triangular(const T nrows) -> decltype(auto)
{
    auto rec_indices = [](const auto &index) { return index.first <= index.second; };

    return views::rectangular(nrows, nrows) | std::views::filter(rec_indices);
}

/// @brief Creates upper triangular indices (i,j) view of square matrix.
/// @tparam T The integral indexing type.
/// @param nrows The number of rows in matrix.
/// @return The upper triangular indices (i,j) view of square matrix.
template <Integral T>
inline auto
upper_triangular(const T nrows) -> decltype(auto)
{
    auto rec_up_indices = [](const auto &index) { return index.first < index.second; };

    return views::rectangular(nrows, nrows) | std::views::filter(rec_up_indices);
}

}  // namespace views

#endif /* CustomViews_hpp */
