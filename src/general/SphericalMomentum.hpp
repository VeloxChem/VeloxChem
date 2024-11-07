#ifndef SphericalMomentum_hpp
#define SphericalMomentum_hpp

#include <cmath>
#include <utility>
#include <vector>

namespace spher_mom {

/// @brief Creates vector of Cartesian to spherical transformation factors.
/// @tparam N The order of angular momentum tensor.
/// @param component The component of spherical momentum to generate transformation factors for.
template <int N>
auto
transformation_factors(const int component) -> std::vector<std::pair<int, double>>
{
    // S type real solid harmonics
    if constexpr (N == 0)
    {
        if (component == 0)
            return {
                {0, 1.0},
            };

        return std::vector<std::pair<int, double>>();
    }

    // P type real solid harmonics
    if constexpr (N == 1)
    {
        if (component == 0)
            return {
                {1, 1.0},
            };
        if (component == 1)
            return {
                {2, 1.0},
            };
        if (component == 2)
            return {
                {0, 1.0},
            };

        return std::vector<std::pair<int, double>>();
    }

    // D type real solid harmonics
    if constexpr (N == 2)
    {
        const double f3 = 2.0 * std::sqrt(3.0);

        if (component == 0)
            return {
                {1, f3},
            };
        if (component == 1)
            return {
                {4, f3},
            };
        if (component == 2)
            return {
                {0, -1.0},
                {3, -1.0},
                {5, 2.0},
            };
        if (component == 3)
            return {
                {2, f3},
            };
        if (component == 4)
            return {
                {0, 0.5 * f3},
                {3, -0.5 * f3},
            };

        return std::vector<std::pair<int, double>>();
    }

    // F type real solid harmonics
    if constexpr (N == 3)
    {
        const double f5 = std::sqrt(2.5);

        const double f15 = 2.0 * std::sqrt(15.0);

        const double f3 = std::sqrt(1.5);

        if (component == 0)
            return {
                {1, 3.0 * f5},
                {6, -f5},
            };
        if (component == 1)
            return {
                {4, f15},
            };
        if (component == 2)
            return {
                {8, 4.0 * f3},
                {1, -f3},
                {6, -f3},
            };
        if (component == 3)
            return {
                {9, 2.0},
                {2, -3.0},
                {7, -3.0},
            };
        if (component == 4)
            return {
                {5, 4.0 * f3},
                {0, -f3},
                {3, -f3},
            };
        if (component == 5)
            return {
                {2, 0.5 * f15},
                {7, -0.5 * f15},
            };
        if (component == 6)
            return {
                {0, f5},
                {3, -3.0 * f5},
            };

        return std::vector<std::pair<int, double>>();
    }

    // G type real solid harmonics
    if constexpr (N == 4)
    {
        const double f35 = 4.0 * std::sqrt(35);

        const double f17 = 4.0 * std::sqrt(17.5);

        const double f5 = 4.0 * std::sqrt(5.0);

        const double f2 = 4.0 * std::sqrt(2.5);

        if (component == 0)
            return {
                {1, f35},
                {6, -f35},
            };
        if (component == 1)
            return {
                {4, 3.0 * f17},
                {11, -f17},
            };
        if (component == 2)
            return {
                {8, 6.0 * f5},
                {1, -f5},
                {6, -f5},
            };
        if (component == 3)
            return {
                {13, 4.0 * f2},
                {4, -3.0 * f2},
                {11, -3.0 * f2},
            };
        if (component == 4)
            return {
                {14, 8.0},
                {0, 3.0},
                {10, 3.0},
                {3, 6.0},
                {5, -24.0},
                {12, -24.0},
            };
        if (component == 5)
            return {
                {9, 4.0 * f2},
                {2, -3.0 * f2},
                {7, -3.0 * f2},
            };
        if (component == 6)
            return {
                {5, 3.0 * f5},
                {12, -3.0 * f5},
                {0, -0.5 * f5},
                {10, 0.5 * f5},
            };
        if (component == 7)
            return {
                {2, f17},
                {7, -3.0 * f17},
            };
        if (component == 8)
            return {
                {0, 0.25 * f35},
                {10, 0.25 * f35},
                {3, -1.5 * f35},
            };

        return std::vector<std::pair<int, double>>();
    }

    // TODO: Add higher order transformation factors l > 4

    return std::vector<std::pair<int, double>>();
}

}  // namespace spher_mom

#endif /* SphericalMomentum_hpp */
