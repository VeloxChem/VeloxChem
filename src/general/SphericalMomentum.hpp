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
    // NOTE: Cartesian solid harmonics expansion factors are generated
    // recursively using Eqs. A2-A6 (see Supporting Info).
    // J. Chem. Theory Comput. 2020, https://doi.org/10.1021/acs.jctc.9b01296
    
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
        const double f3 = std::sqrt(3.0);

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
                {0, -0.5},
                {3, -0.5},
                {5, 1.0},
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
        const double f10 = 0.25 * std::sqrt(10.0);

        const double f15 = std::sqrt(15.0);

        const double f6 = std::sqrt(6.0);

        if (component == 0)
            return {
                {1, 3.0 * f10},
                {6, -f10},
            };
        if (component == 1)
            return {
                {4, f15},
            };
        if (component == 2)
            return {
                {8, f6},
                {1, -0.25 * f6},
                {6, -0.25 * f6},
            };
        if (component == 3)
            return {
                {9, 1.0},
                {2, -1.5},
                {7, -1.5},
            };
        if (component == 4)
            return {
                {5,  f6},
                {0, -0.25 * f6},
                {3, -0.25 * f6},
            };
        if (component == 5)
            return {
                {2, 0.5 * f15},
                {7, -0.5 * f15},
            };
        if (component == 6)
            return {
                {0, f10},
                {3, -3.0 * f10},
            };

        return std::vector<std::pair<int, double>>();
    }

    // G type real solid harmonics
    if constexpr (N == 4)
    {
        const double f35 = 0.5 * std::sqrt(35.0);

        const double f70 = 0.25 * std::sqrt(70.0);

        const double f10 = std::sqrt(10.0);

        const double f5 = 0.5 * std::sqrt(5.0);

        if (component == 0)
            return {
                {1, f35},
                {6, -f35},
            };
        if (component == 1)
            return {
                {4, 3.0 * f70},
                {11, -f70},
            };
        if (component == 2)
            return {
                {8, 6.0 * f5},
                {1, -f5},
                {6, -f5},
            };
        if (component == 3)
            return {
                {13, f10},
                {4, -0.75 * f10},
                {11, -0.75 * f10},
            };
        if (component == 4)
            return {
                {14, 1.0},
                {0, 0.375},
                {10, 0.375},
                {3, 0.75},
                {5, -3.0},
                {12, -3.0},
            };
        if (component == 5)
            return {
                {9, f10},
                {2, -0.75 * f10},
                {7, -0.75 * f10},
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
                {2, f70},
                {7, -3.0 * f70},
            };
        if (component == 8)
            return {
                {0, 0.25 * f35},
                {10, 0.25 * f35},
                {3, -1.5 * f35},
            };

        return std::vector<std::pair<int, double>>();
    }
    
    // H type real solid harmonics
    if constexpr (N == 5)
    {
        const double f14 = 0.3750 * std::sqrt(14.0);
        
        const double f35 = 0.3750 * std::sqrt(35.0);
        
        const double f70 = 0.0625 * std::sqrt(70.0);
        
        const double f105 = 0.25 * std::sqrt(105.0);
        
        const double f15 = 0.125 * std::sqrt(15.0);
        
        if (component == 0)
            return {
                {15, 0.5 * f14},
                {6, -5.0 * f14},
                {1,  2.5 * f14},
            };
        if (component == 1)
            return {
                {11, -4.0 * f35},
                {4,   4.0 * f35},
            };
        if (component == 2)
            return {
                {17, -8.0 * f70},
                {15, f70},
                {8,  24.0 * f70},
                {6,  -2.0 * f70},
                {1,  -3.0 * f70},
            };
        if (component == 3)
            return {
                {13,  4.0 * f105},
                {11, -2.0 * f105},
                {4,  -2.0 * f105},
            };
        if (component == 4)
            return {
                {19,  8.0 * f15},
                {17, -12.0 * f15},
                {15,  f15},
                {8,  -12.0 * f15},
                {6,   2.0 * f15},
                {1,   f15},
            };
        if (component == 5)
            return {
                {20,  1.0},
                {18, -5.0},
                {16,  1.875},
                {9,  -5.0},
                {7,   3.75},
                {2,   1.875},
            };
        if (component == 6)
            return {
                {14,  8.0 * f15},
                {12, -12.0 * f15},
                {10,  f15},
                {5,  -12.0 * f15},
                {3,   2.0 * f15},
                {0,   f15},
            };
        if (component == 7)
            return {
                {18, -2.0 * f105},
                {16,  f105},
                {9,   2.0 * f105},
                {2,  -f105},
            };
        if (component == 8)
            return {
                {12, -24.0 * f70},
                {10,  3.0 * f70},
                {5,  8.0 * f70},
                {3,   2.0 * f70},
                {0,  -f70},
            };
        if (component == 9)
            return {
                {16, f35},
                {7, -6.0 * f35},
                {2,  f35},
            };
        if (component == 10)
            return {
                {10, 2.5 * f14},
                {3, -5.0 * f14},
                {0,  0.5 * f14},
            };
        
        return std::vector<std::pair<int, double>>();
    }
    
    // I type real solid harmonics
    if constexpr (N == 6)
    {
        const double f462 = 0.03125 * std::sqrt(462.0);
        
        const double f154 = 0.18750 * std::sqrt(154.0);
        
        const double f7 = 0.18750 * std::sqrt(7.0);
       
        if (component == 0)
            return {
                {15, 6.0 * f462},
                {6, -20.0 * f462},
                {1,  6.0 * f462},
            };
        if (component == 1)
            return {
                {22, f154},
                {11, -10.0 * f154},
                {4,  5.0 * f154},
            };
        if (component == 2)
            return {
                {17, -40.0 * f7},
                {15,   4.0 * f7},
                {4,   40.0 * f7},
                {4,   -4.0 * f7},
            };
        
        S(l=6,m=-4) = -3*sqrt(7)*x**5*y/4 + 15*sqrt(7)*x**3*y*z**2/2 + 3*sqrt(7)*x*y**5/4 - 15*sqrt(7)*x*y**3*z**2/2
           |--> s^{6,-4}_{1,3,2} = -15*sqrt(7)/2
           |--> s^{6,-4}_{1,5,0} = 3*sqrt(7)/4
           |--> s^{6,-4}_{3,1,2} = 15*sqrt(7)/2
           |--> s^{6,-4}_{5,1,0} = -3*sqrt(7)/4
        
        S(l=6,m=+4) = -3*sqrt(7)*x**6/16 + 15*sqrt(7)*x**4*y**2/16 + 15*sqrt(7)*x**4*z**2/8 + 15*sqrt(7)*x**2*y**4/16 - 45*sqrt(7)*x**2*y**2*z**2/4 - 3*sqrt(7)*y**6/16 + 15*sqrt(7)*y**4*z**2/8
          |--> s^{6,4}_{0,4,2} = 15*sqrt(7)/8
          |--> s^{6,4}_{0,6,0} = -3*sqrt(7)/16
          |--> s^{6,4}_{2,2,2} = -45*sqrt(7)/4
          |--> s^{6,4}_{2,4,0} = 15*sqrt(7)/16
          |--> s^{6,4}_{4,0,2} = 15*sqrt(7)/8
          |--> s^{6,4}_{4,2,0} = 15*sqrt(7)/16
          |--> s^{6,4}_{6,0,0} = -3*sqrt(7)/16
        
        
        if (component == 11)
            return {
                {16, 5.0 *f154},
                {7, -10.0 * f154},
                {2,  f154},
            };
        if (component == 12)
            return {
                {21, -f462},
                {10,  15.0 * f462},
                {3,  -15.0 * f462},
                {0,  f462},
            };
    }
    
    // TODO: Add higher order transformation factors l > 6

    return std::vector<std::pair<int, double>>();
}

}  // namespace spher_mom

#endif /* SphericalMomentum_hpp */
