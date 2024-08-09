#include "T4CDiagonalDistributor.hpp"

CT4CDiagonalDistributor::CT4CDiagonalDistributor(double* max_values)
{
    _max_values = max_values;
}

auto
CT4CDiagonalDistributor::distribute(      std::vector<double>& max_integrals,
                                    const std::array<int, 2>&  gto_range) -> void
{
    for (int i = gto_range[0]; i < gto_range[1]; i++)
    {
        _max_values[i] = max_integrals[i - gto_range[0]]; 
    }
}
