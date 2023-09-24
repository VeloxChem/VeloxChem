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
quadChebyshevOfKindTwo(double* coordinates, double* weights, const int64_t nPoints) -> void
{
    // prefactor

    auto fstep = mathconst::getPiValue() / (static_cast<double>(nPoints) + 1.0);

    // loop over grid points

    for (int64_t i = 1; i < nPoints + 1; i++)
    {
        auto farg = static_cast<double>(i) * fstep;

        coordinates[i - 1] = std::cos(farg);

        auto warg = std::sin(farg);

        weights[i - 1] = fstep * warg * warg;
    }
}

auto
batch_sizes(const int64_t nElements, const int64_t nodes) -> std::vector<int64_t>
{
    int64_t ave = nElements / nodes;

    int64_t rem = nElements % nodes;

    std::vector<int64_t> counts;

    for (int64_t p = 0; p < nodes; p++)
    {
        counts.push_back((p < rem) ? (ave + 1) : ave);
    }

    return counts;
}

auto
batch_offsets(const int64_t nElements, const int64_t nodes) -> std::vector<int64_t>
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    std::vector<int64_t> displs;

    int64_t index = 0;

    for (int64_t p = 0; p < nodes; p++)
    {
        displs.push_back(index);

        index += counts[p];
    }

    return displs;
}

auto
batch_size(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t
{
    auto counts = mathfunc::batch_sizes(nElements, nodes);

    return counts[rank];
}

auto
batch_offset(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t
{
    auto displs = mathfunc::batch_offsets(nElements, nodes);

    return displs[rank];
}

}  // namespace mathfunc
