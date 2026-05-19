//
//  Tabula — custom-recursion molecular-integral machinery.
//  Per-CGTO screening bounds of a basis-function block.
//

#include "TabulaGtoBlockScreeningData.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace tabula {  // tabula namespace

auto
make_screening_data(const CGtoBlock& block) -> GtoBlockScreeningData
{
    const auto ncgtos = block.number_of_basis_functions();
    const auto npgtos = block.number_of_primitives();

    // primitive-major SoA: element (primitive p, CGTO c) is at p * ncgtos + c
    const auto exponents = block.exponents();
    const auto norms     = block.normalization_factors();

    GtoBlockScreeningData data;
    data.maxCoefficient.assign(static_cast<std::size_t>(ncgtos), 0.0);
    data.minExponent.assign(static_cast<std::size_t>(ncgtos), 0.0);
    data.maxExponent.assign(static_cast<std::size_t>(ncgtos), 0.0);

    for (int c = 0; c < ncgtos; c++)
    {
        double maxCoef = 0.0;
        double minExp  = 0.0;
        double maxExp  = 0.0;

        for (int p = 0; p < npgtos; p++)
        {
            const auto idx = static_cast<std::size_t>(p * ncgtos + c);

            maxCoef = std::max(maxCoef, std::abs(norms[idx]));

            const auto exponent = exponents[idx];

            if (p == 0)
            {
                minExp = exponent;
                maxExp = exponent;
            }
            else
            {
                minExp = std::min(minExp, exponent);
                maxExp = std::max(maxExp, exponent);
            }
        }

        data.maxCoefficient[static_cast<std::size_t>(c)] = maxCoef;
        data.minExponent[static_cast<std::size_t>(c)]    = minExp;
        data.maxExponent[static_cast<std::size_t>(c)]    = maxExp;
    }

    return data;
}

}  // namespace tabula
