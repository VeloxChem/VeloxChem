#include "T4CDiagonalDistributor.hpp"

CT4CDiagonalDistributor::CT4CDiagonalDistributor(double* max_values)
{
    _max_values = max_values;
}

auto
CT4CDiagonalDistributor::distribute(const std::vector<double>& max_integrals, const std::pair<size_t, size_t>& gto_range) -> void
{
    for (auto i = gto_range.first; i < gto_range.second; i++)
    {
        _max_values[i] = max_integrals[i - gto_range.first];
    }
}
