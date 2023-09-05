#include "DftFunc.hpp"

namespace dft {  // dft namespace

auto
distribute(      CSubMatrix*          matrix,
           const std::vector<double>& values,
           const int64_t              irow) -> void
{
    const auto ncols = static_cast<int64_t>(values.size());
    
    for (int64_t i = 0; i < ncols; i++)
    {
        matrix->at(irow, i, false) = values[i]; 
    }
}

auto
distribute(      CSubMatrix*          matrix,
           const std::vector<double>& values,
           const double               factor,
           const int64_t              irow) -> void
{
    const auto ncols = static_cast<int64_t>(values.size());
    
    for (int64_t i = 0; i < ncols; i++)
    {
        matrix->at(irow, i, false) = factor * values[i];
    }
}

}  // namespace dft
