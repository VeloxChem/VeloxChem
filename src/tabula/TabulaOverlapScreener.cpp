//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap screening predicate.
//

#include "TabulaOverlapScreener.hpp"

#include <cmath>
#include <limits>

#include "MathConst.hpp"

namespace tabula {  // tabula namespace

auto
make_overlap_screening_estimator(const int l_a, const int l_c)
    -> std::function<double(const CGtoBlockScreeningData &, const CGtoBlockScreeningData &, const double)>
{
    // m — the max seed order of the overlap (l_a, l_c) recursion: the total
    // angular momentum, the seed-ladder depth.
    const int m = l_a + l_c;

    return [l_a, l_c, m](const CGtoBlockScreeningData &bra, const CGtoBlockScreeningData &ket, const double r) -> double {
        // same-center pairs are always kept — the r^(l_a+l_c) factor below is
        // an asymptotic (large-R) bound and is not valid at R = 0
        if (r == 0.0)
        {
            return std::numeric_limits<double>::infinity();
        }

        const auto a_min = bra.min_exponent;
        const auto a_max = bra.max_exponent;
        const auto b_min = ket.min_exponent;
        const auto b_max = ket.max_exponent;

        // reduced exponents — rho_min from the slowest-decaying primitives,
        // rho_max from the fastest
        const auto rho_min = a_min * b_min / (a_min + b_min);
        const auto rho_max = a_max * b_max / (a_max + b_max);

        // estimate = cMaxA·cMaxB · r^(l_a+l_c) · (2·rho_max)^m · exp(-rho_min·r²)
        //          · (pi/(a_min+b_min))^(3/2) · (1/2·a_min)^l_a · (1/2·b_min)^l_c
        auto estimate = bra.max_coefficient * ket.max_coefficient;

        estimate *= std::pow(r, l_a + l_c);
        estimate *= std::pow(2.0 * rho_max, m);
        estimate *= std::exp(-rho_min * r * r);

        const auto fact = mathconst::pi_value() / (a_min + b_min);
        estimate *= fact * std::sqrt(fact);

        if (l_a > 0) estimate *= std::pow(0.5 / a_min, l_a);
        if (l_c > 0) estimate *= std::pow(0.5 / b_min, l_c);

        return estimate;
    };
}

}  // namespace tabula
