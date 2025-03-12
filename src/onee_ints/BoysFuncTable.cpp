//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "BoysFuncTable.hpp"

#include <cmath>

#include "BoysFunc.hpp"
#include "MathConst.hpp"

#define MATH_CONST_PI 3.14159265358979323846

#define MATH_CONST_HALF_SQRT_PI 0.88622692545275794096

namespace onee {  // onee namespace

auto
getBoysFuncFactors() -> std::vector<double>
{
    return std::vector<double>({1.0,        1.0 / 3.0,  1.0 / 5.0,  1.0 / 7.0,  1.0 / 9.0,  1.0 / 11.0, 1.0 / 13.0,
                                1.0 / 15.0, 1.0 / 17.0, 1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0, 1.0 / 27.0,
                                1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0, 1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                                1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0, 1.0 / 51.0, 1.0 / 53.0, 1.0 / 55.0});
}

auto
getFullBoysFuncTable() -> std::vector<double>
{
    // Boys function (tabulated for order 0-28)

    std::vector<double> boys_func_table((28 + 1) * 121 * 7);

    for (int bf_order = 0; bf_order <= 28; bf_order++)
    {
        const auto bf_table = onee::getBoysFuncTable(bf_order);

        auto bf_data = boys_func_table.data() + bf_order * 121 * 7;

        for (int r = 0; r < 121; r++)
        {
            for (int c = 0; c < 7; c++)
            {
                bf_data[r * 7 + c] = bf_table[r][c];
            }
        }
    }

    return boys_func_table;
}

auto
getBoysFuncTable(const int N) -> std::array<std::array<double, 7>, 121>
{
    switch (N)
    {
        case 0:
            {
                const CBoysFunc<0> bf_table;
                return bf_table.get_table();
            }
        case 1:
            {
                const CBoysFunc<1> bf_table;
                return bf_table.get_table();
            }
        case 2:
            {
                const CBoysFunc<2> bf_table;
                return bf_table.get_table();
            }
        case 3:
            {
                const CBoysFunc<3> bf_table;
                return bf_table.get_table();
            }
        case 4:
            {
                const CBoysFunc<4> bf_table;
                return bf_table.get_table();
            }
        case 5:
            {
                const CBoysFunc<5> bf_table;
                return bf_table.get_table();
            }
        case 6:
            {
                const CBoysFunc<6> bf_table;
                return bf_table.get_table();
            }
        case 7:
            {
                const CBoysFunc<7> bf_table;
                return bf_table.get_table();
            }
        case 8:
            {
                const CBoysFunc<8> bf_table;
                return bf_table.get_table();
            }
        case 9:
            {
                const CBoysFunc<9> bf_table;
                return bf_table.get_table();
            }
        case 10:
            {
                const CBoysFunc<10> bf_table;
                return bf_table.get_table();
            }
        case 11:
            {
                const CBoysFunc<11> bf_table;
                return bf_table.get_table();
            }
        case 12:
            {
                const CBoysFunc<12> bf_table;
                return bf_table.get_table();
            }
        case 13:
            {
                const CBoysFunc<13> bf_table;
                return bf_table.get_table();
            }
        case 14:
            {
                const CBoysFunc<14> bf_table;
                return bf_table.get_table();
            }
        case 15:
            {
                const CBoysFunc<15> bf_table;
                return bf_table.get_table();
            }
        case 16:
            {
                const CBoysFunc<16> bf_table;
                return bf_table.get_table();
            }
        case 17:
            {
                const CBoysFunc<17> bf_table;
                return bf_table.get_table();
            }
        case 18:
            {
                const CBoysFunc<18> bf_table;
                return bf_table.get_table();
            }
        case 19:
            {
                const CBoysFunc<19> bf_table;
                return bf_table.get_table();
            }
        case 20:
            {
                const CBoysFunc<20> bf_table;
                return bf_table.get_table();
            }
        case 21:
            {
                const CBoysFunc<21> bf_table;
                return bf_table.get_table();
            }
        case 22:
            {
                const CBoysFunc<22> bf_table;
                return bf_table.get_table();
            }
        case 23:
            {
                const CBoysFunc<23> bf_table;
                return bf_table.get_table();
            }
        case 24:
            {
                const CBoysFunc<24> bf_table;
                return bf_table.get_table();
            }
        case 25:
            {
                const CBoysFunc<25> bf_table;
                return bf_table.get_table();
            }
        case 26:
            {
                const CBoysFunc<26> bf_table;
                return bf_table.get_table();
            }
        case 27:
            {
                const CBoysFunc<27> bf_table;
                return bf_table.get_table();
            }
        case 28:
            {
                const CBoysFunc<28> bf_table;
                return bf_table.get_table();
            }
        default:
            {
                return std::array<std::array<double, 7>, 121>();
            }
    }
}

auto
computeBoysFunction(double* values, const double fa, const uint32_t N, const double* bf_table, const double* ft) -> void
{
    // Note: 847 = 121 * 7
    const double* bf_data = bf_table + N * 847;

    uint32_t pnt = (fa > 1.0e5) ? 1000000 : static_cast<uint32_t>(10.0 * fa + 0.5);

    if (pnt < 121)
    {
        const double w = fa - 0.1 * pnt;

        const double w2 = w * w;

        const double w4 = w2 * w2;

        values[N] = bf_data[pnt * 7 + 0] + bf_data[pnt * 7 + 1] * w + bf_data[pnt * 7 + 2] * w2 + bf_data[pnt * 7 + 3] * w2 * w

                    + bf_data[pnt * 7 + 4] * w4 + bf_data[pnt * 7 + 5] * w4 * w + bf_data[pnt * 7 + 6] * w4 * w2;

        const double f2a = fa + fa;

        const double fx = exp(-fa);

        for (uint32_t j = 0; j < N; j++)
        {
            values[N - j - 1] = ft[N - j - 1] * (f2a * values[N - j] + fx);
        }
    }
    else
    {
        const double fia = 1.0 / fa;

        double pf = 0.5 * fia;

        values[0] = MATH_CONST_HALF_SQRT_PI * sqrt(fia);

        if (pnt < 921)
        {
            const double fia2 = fia * fia;

            const double f = 0.4999489092 * fia - 0.2473631686 * fia2 + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

            const double fx = exp(-fa);

            values[0] -= f * fx;

            const double rterm = pf * fx;

            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1] - rterm;

                pf += fia;
            }
        }
        else
        {
            for (uint32_t j = 1; j <= N; j++)
            {
                values[j] = pf * values[j - 1];

                pf += fia;
            }
        }
    }
}

}  // namespace onee
