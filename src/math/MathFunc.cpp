#include "MathFunc.hpp"

#include <cmath>
#include <cstdlib>

#include "MathConst.hpp"

namespace mathfunc {  // mathfunc namespace

auto
countSignificantElements(const std::vector<int64_t>& mask) -> int64_t
{
    int64_t nelems = 0;

    for (auto mvalue : mask)
    {
        if (mvalue == 1) nelems++;
    }

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

auto
quadChebyshevOfKindTwo(double* coordinates, double* weights, const int32_t nPoints) -> void
{
    // prefactor

    auto fstep = mathconst::getPiValue() / (static_cast<double>(nPoints) + 1.0);

    // loop over grid points

    for (int32_t i = 1; i < nPoints + 1; i++)
    {
        auto farg = static_cast<double>(i) * fstep;

        coordinates[i - 1] = std::cos(farg);

        auto warg = std::sin(farg);

        weights[i - 1] = fstep * warg * warg;
    }
}

}  // namespace mathfunc
