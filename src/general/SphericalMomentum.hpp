//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
        
        const double f210 = 0.0625 * std::sqrt(210.0);
        
        const double f21 = 0.1250 * std::sqrt(21.0);
       
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
                {8,   40.0 * f7},
                {1,   -4.0 * f7},
            };
        if (component == 3)
            return {
                {24, -8.0 * f210},
                {22,  3.0 * f210},
                {13,  24.0 * f210},
                {11, -6.0 * f210},
                {4,  -9.0 * f210},
            };
        if (component == 4)
            return {
                {19,  16.0 * f210},
                {17, -16.0 * f210},
                {15,  f210},
                {8,  -16.0 * f210},
                {6,   2.0 * f210},
                {1,  f210},
            };
        if (component == 5)
            return {
                {26,  8.0 * f21},
                {24, -20.0 * f21},
                {22,  5.0  * f21},
                {13, -20.0 * f21},
                {11,  10.0 * f21},
                {4,   5.0 * f21},
            };
        if (component == 6)
            return {
                {27,  1.0000},
                {25, -7.5000},
                {23,  5.6250},
                {21, -0.3125},
                {14, -7.5000},
                {12,  11.250},
                {10,  -0.9375},
                {5,   5.6250},
                {3,  -0.9375},
                {0, -0.3125},
            };
        if (component == 7)
            return {
                {20,  8.0 * f21},
                {18, -20.0 * f21},
                {16,  5.0  * f21},
                {9,  -20.0 * f21},
                {7,   10.0 * f21},
                {2,   5.0 * f21},
            };
        if (component == 8)
            return {
                {25, -8.0 * f210},
                {23,  8.0 * f210},
                {21, -0.5 * f210},
                {14,   8.0 * f210},
                {10,  -0.5 * f210},
                {5, -8.0 * f210},
                {3,   0.5 * f210},
                {0,   0.5 * f210},
            };
        if (component == 9)
            return {
                {18, -24.0 * f210},
                {16,  9.0 * f210},
                {9,   8.0 * f210},
                {7,  6.0 * f210},
                {2,  -3.0 * f210},
            };
        if (component == 10)
            return {
                {23,  10.0 * f7},
                {21, -f7},
                {12, -60.0 * f7},
                {10,  5.0 * f7},
                {5,  10.0 * f7},
                {3,   5.0 * f7},
                {0, -f7},
            };
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
