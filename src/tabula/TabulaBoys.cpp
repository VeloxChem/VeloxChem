//
//  Tabula — custom-recursion molecular-integral machinery.
//  Boys function — three-region rational-minimax evaluation.
//

#include "TabulaBoys.hpp"

#include <cmath>

#include "TabulaBoysTable.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

constexpr double boys_half_sqrt_pi = 0.886226925452758014;  // √π / 2
constexpr double boys_region_a_end = 11.899848152108484;    // x_0
constexpr double boys_region_b_end = 28.989337738820740;    // x_1

/// @brief Horner evaluation of the degree-increasing coefficient run
/// `coeff[begin, end)` at `x`.
inline auto
horner(const double *coeff, const int begin, const int end, const double x) -> double
{
    double result = coeff[end - 1];
    for (int i = end - 2; i >= begin; i--)
    {
        result = result * x + coeff[i];
    }
    return result;
}

}  // namespace

auto
boys(const int order, const double x, double *results) -> void
{
    if (x < boys_region_a_end)
    {
        // region A — a per-order rational minimax for F_order, then the
        // downward recursion F_{l-1} = (2x·F_l + e^{-x}) / (2l-1)
        const int num_begin = detail::boys_region_a_numerator_offset[order];
        const int num_end   = detail::boys_region_a_numerator_offset[order + 1];
        const int den_begin = detail::boys_region_a_denominator_offset[order];
        const int den_end   = detail::boys_region_a_denominator_offset[order + 1];

        const double numerator   = horner(detail::boys_region_a_numerator, num_begin, num_end, x);
        const double denominator = horner(detail::boys_region_a_denominator, den_begin, den_end, x);

        results[order] = numerator / denominator;

        if (order > 0)
        {
            const double ex    = std::exp(-x);
            const double two_x = 2.0 * x;
            for (int l = order; l >= 1; l--)
            {
                results[l - 1] = (two_x * results[l] + ex) / static_cast<double>(2 * l - 1);
            }
        }
    }
    else if (x < boys_region_b_end)
    {
        // region B — a rational minimax for F_0, then the upward recursion
        // F_l = ((2l-1)·F_{l-1} - e^{-x}) / (2x)
        constexpr int b_num_end = static_cast<int>(sizeof(detail::boys_region_b_numerator) / sizeof(double));
        constexpr int b_den_end = static_cast<int>(sizeof(detail::boys_region_b_denominator) / sizeof(double));

        const double numerator   = horner(detail::boys_region_b_numerator, 0, b_num_end, x);
        const double denominator = horner(detail::boys_region_b_denominator, 0, b_den_end, x);

        results[0] = numerator / denominator;

        if (order > 0)
        {
            const double ex     = std::exp(-x);
            const double inv_2x = 1.0 / (2.0 * x);
            for (int l = 1; l <= order; l++)
            {
                results[l] = (static_cast<double>(2 * l - 1) * results[l - 1] - ex) * inv_2x;
            }
        }
    }
    else
    {
        // region C — the asymptotic F_0 = (√π/2)/√x, then the upward
        // recursion F_l = (2l-1)·F_{l-1} / (2x); the e^{-x} tail is negligible
        results[0] = boys_half_sqrt_pi / std::sqrt(x);

        const double inv_2x = 1.0 / (2.0 * x);
        for (int l = 1; l <= order; l++)
        {
            results[l] = static_cast<double>(2 * l - 1) * results[l - 1] * inv_2x;
        }
    }
}

}  // namespace tabula
