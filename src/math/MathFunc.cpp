#include "MathFunc.hpp"

namespace mathfunc {  // mathfunc namespace

auto
countSignificantElements(const std::vector<int64_t>& mask) -> int64_t
{
    int64_t nelems = 0;
    
    for (auto mvalue : mask) if (mvalue == 1) nelems++;
    
    return nelems;
}

auto
zero(std::vector<double>& values) -> void
{
    const auto ndim = values.size();
    
    auto ptr_values = values.data();
    
    #pragma omp simd
    for (size_t i = 0; i < ndim; i++)
    {
        ptr_values[i] = 0.0;
    }
}

}  // namespace mathfunc
