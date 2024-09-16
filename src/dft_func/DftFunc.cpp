#include "DftFunc.hpp"

#include <algorithm>
#include <ranges>

namespace gtoval {  // gtoval namespace

auto
distribute(CSubMatrix* matrix, const std::vector<double>& values, const size_t irow) -> void
{
    if (const auto ncols = values.size(); ncols > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, irow), [&](const auto i) {
            matrix->at({irow, i}) += values[i]; 
        });
    }
}

auto
distribute(CSubMatrix* matrix, const std::vector<double>& values, const double factor, const size_t irow) -> void
{
    if (const auto ncols = values.size(); ncols > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, irow), [&](const auto i) {
            matrix->at({irow, i}) += factor * values[i];
        });
    }
}

}  // namespace gtoval
