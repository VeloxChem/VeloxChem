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

#ifndef RealSolidHarmonicAB_hpp
#define RealSolidHarmonicAB_hpp

#include <cmath>

namespace harm {

inline double Y_ll_0_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 1.0;
}

inline double Y_ll_1_m_n1(double AB_x, double AB_y, double AB_z)
{
    return AB_y;
}

inline double Y_ll_1_m_p0(double AB_x, double AB_y, double AB_z)
{
    return AB_z;
}

inline double Y_ll_1_m_p1(double AB_x, double AB_y, double AB_z)
{
    return AB_x;
}

inline double Y_ll_2_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3 = std::sqrt(3.0);

    return sqrt3 * AB_x * AB_y;
}

inline double Y_ll_2_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3 = std::sqrt(3.0);

    return sqrt3 * AB_y * AB_z;
}

inline double Y_ll_2_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -0.5 * AB_x * AB_x - 0.5 * AB_y * AB_y + AB_z * AB_z;
}

inline double Y_ll_2_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3 = std::sqrt(3.0);

    return sqrt3 * AB_x * AB_z;
}

inline double Y_ll_2_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3 = std::sqrt(3.0);

    return 0.5 * sqrt3 * AB_x * AB_x - 0.5 * sqrt3 * AB_y * AB_y;
}

inline double Y_ll_3_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt10 = std::sqrt(10.0);

    return 0.75 * sqrt10 * AB_x * AB_x * AB_y - 0.25 * sqrt10 * AB_y * AB_y * AB_y;
}

inline double Y_ll_3_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt15 = std::sqrt(15.0);

    return sqrt15 * AB_x * AB_y * AB_z;
}

inline double Y_ll_3_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt6 = std::sqrt(6.0);

    return -0.25 * sqrt6 * AB_x * AB_x * AB_y - 0.25 * sqrt6 * AB_y * AB_y * AB_y + sqrt6 * AB_y * AB_z * AB_z;
}

inline double Y_ll_3_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -1.5 * AB_x * AB_x * AB_z - 1.5 * AB_y * AB_y * AB_z + AB_z * AB_z * AB_z;
}

inline double Y_ll_3_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt6 = std::sqrt(6.0);

    return -0.25 * sqrt6 * AB_x * AB_x * AB_x - 0.25 * sqrt6 * AB_x * AB_y * AB_y + sqrt6 * AB_x * AB_z * AB_z;
}

inline double Y_ll_3_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt15 = std::sqrt(15.0);

    return 0.5 * sqrt15 * AB_x * AB_x * AB_z - 0.5 * sqrt15 * AB_y * AB_y * AB_z;
}

inline double Y_ll_3_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt10 = std::sqrt(10.0);

    return 0.25 * sqrt10 * AB_x * AB_x * AB_x - 0.75 * sqrt10 * AB_x * AB_y * AB_y;
}

inline double Y_ll_4_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt35 = std::sqrt(35.0);

    return 0.5 * sqrt35 * AB_x * AB_x * AB_x * AB_y - 0.5 * sqrt35 * AB_x * AB_y * AB_y * AB_y;
}

inline double Y_ll_4_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return 0.75 * sqrt70 * AB_x * AB_x * AB_y * AB_z - 0.25 * sqrt70 * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_4_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5 = std::sqrt(5.0);

    return -0.5 * sqrt5 * AB_x * AB_x * AB_x * AB_y - 0.5 * sqrt5 * AB_x * AB_y * AB_y * AB_y + 3.0 * sqrt5 * AB_x * AB_y * AB_z * AB_z;
}

inline double Y_ll_4_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt10 = std::sqrt(10.0);

    return -0.75 * sqrt10 * AB_x * AB_x * AB_y * AB_z - 0.75 * sqrt10 * AB_y * AB_y * AB_y * AB_z + sqrt10 * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_4_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 0.375 * AB_x * AB_x * AB_x * AB_x + 0.75 * AB_x * AB_x * AB_y * AB_y - 3.0 * AB_x * AB_x * AB_z * AB_z + 0.375 * AB_y * AB_y * AB_y * AB_y - 3.0 * AB_y * AB_y * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_4_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt10 = std::sqrt(10.0);

    return -0.75 * sqrt10 * AB_x * AB_x * AB_x * AB_z - 0.75 * sqrt10 * AB_x * AB_y * AB_y * AB_z + sqrt10 * AB_x * AB_z * AB_z * AB_z;
}

inline double Y_ll_4_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5 = std::sqrt(5.0);

    return -0.25 * sqrt5 * AB_x * AB_x * AB_x * AB_x + 1.5 * sqrt5 * AB_x * AB_x * AB_z * AB_z + 0.25 * sqrt5 * AB_y * AB_y * AB_y * AB_y - 1.5 * sqrt5 * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_4_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return 0.25 * sqrt70 * AB_x * AB_x * AB_x * AB_z - 0.75 * sqrt70 * AB_x * AB_y * AB_y * AB_z;
}

inline double Y_ll_4_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt35 = std::sqrt(35.0);

    return 0.125 * sqrt35 * AB_x * AB_x * AB_x * AB_x - 0.75 * sqrt35 * AB_x * AB_x * AB_y * AB_y + 0.125 * sqrt35 * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_5_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt14 = std::sqrt(14.0);

    return 0.9375 * sqrt14 * AB_x * AB_x * AB_x * AB_x * AB_y - 1.875 * sqrt14 * AB_x * AB_x * AB_y * AB_y * AB_y + 0.1875 * sqrt14 * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_5_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt35 = std::sqrt(35.0);

    return 1.5 * sqrt35 * AB_x * AB_x * AB_x * AB_y * AB_z - 1.5 * sqrt35 * AB_x * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_5_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return -0.1875 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_y - 0.125 * sqrt70 * AB_x * AB_x * AB_y * AB_y * AB_y + 1.5 * sqrt70 * AB_x * AB_x * AB_y * AB_z * AB_z + 0.0625 * sqrt70 * AB_y * AB_y * AB_y * AB_y * AB_y - 0.5 * sqrt70 * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_5_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt105 = std::sqrt(105.0);

    return -0.5 * sqrt105 * AB_x * AB_x * AB_x * AB_y * AB_z - 0.5 * sqrt105 * AB_x * AB_y * AB_y * AB_y * AB_z + sqrt105 * AB_x * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_5_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt15 = std::sqrt(15.0);

    return 0.125 * sqrt15 * AB_x * AB_x * AB_x * AB_x * AB_y + 0.25 * sqrt15 * AB_x * AB_x * AB_y * AB_y * AB_y - 1.5 * sqrt15 * AB_x * AB_x * AB_y * AB_z * AB_z + 0.125 * sqrt15 * AB_y * AB_y * AB_y * AB_y * AB_y - 1.5 * sqrt15 * AB_y * AB_y * AB_y * AB_z * AB_z + sqrt15 * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_5_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 1.875 * AB_x * AB_x * AB_x * AB_x * AB_z + 3.75 * AB_x * AB_x * AB_y * AB_y * AB_z - 5.0 * AB_x * AB_x * AB_z * AB_z * AB_z + 1.875 * AB_y * AB_y * AB_y * AB_y * AB_z - 5.0 * AB_y * AB_y * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_5_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt15 = std::sqrt(15.0);

    return 0.125 * sqrt15 * AB_x * AB_x * AB_x * AB_x * AB_x + 0.25 * sqrt15 * AB_x * AB_x * AB_x * AB_y * AB_y - 1.5 * sqrt15 * AB_x * AB_x * AB_x * AB_z * AB_z + 0.125 * sqrt15 * AB_x * AB_y * AB_y * AB_y * AB_y - 1.5 * sqrt15 * AB_x * AB_y * AB_y * AB_z * AB_z + sqrt15 * AB_x * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_5_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt105 = std::sqrt(105.0);

    return -0.25 * sqrt105 * AB_x * AB_x * AB_x * AB_x * AB_z + 0.5 * sqrt105 * AB_x * AB_x * AB_z * AB_z * AB_z + 0.25 * sqrt105 * AB_y * AB_y * AB_y * AB_y * AB_z - 0.5 * sqrt105 * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_5_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return -0.0625 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x + 0.125 * sqrt70 * AB_x * AB_x * AB_x * AB_y * AB_y + 0.5 * sqrt70 * AB_x * AB_x * AB_x * AB_z * AB_z + 0.1875 * sqrt70 * AB_x * AB_y * AB_y * AB_y * AB_y - 1.5 * sqrt70 * AB_x * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_5_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt35 = std::sqrt(35.0);

    return 0.375 * sqrt35 * AB_x * AB_x * AB_x * AB_x * AB_z - 2.25 * sqrt35 * AB_x * AB_x * AB_y * AB_y * AB_z + 0.375 * sqrt35 * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_5_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt14 = std::sqrt(14.0);

    return 0.1875 * sqrt14 * AB_x * AB_x * AB_x * AB_x * AB_x - 1.875 * sqrt14 * AB_x * AB_x * AB_x * AB_y * AB_y + 0.9375 * sqrt14 * AB_x * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_6_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt462 = std::sqrt(462.0);

    return 0.1875 * sqrt462 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.625 * sqrt462 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.1875 * sqrt462 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_6_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt154 = std::sqrt(154.0);

    return 0.9375 * sqrt154 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 1.875 * sqrt154 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.1875 * sqrt154 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_6_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt7 = std::sqrt(7.0);

    return -0.75 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 7.5 * sqrt7 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.75 * sqrt7 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 7.5 * sqrt7 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_6_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt210 = std::sqrt(210.0);

    return -0.5625 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.375 * sqrt210 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.5 * sqrt210 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.1875 * sqrt210 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.5 * sqrt210 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt210 = std::sqrt(210.0);

    return 0.0625 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.125 * sqrt210 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - sqrt210 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.0625 * sqrt210 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - sqrt210 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + sqrt210 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt21 = std::sqrt(21.0);

    return 0.625 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 1.25 * sqrt21 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 2.5 * sqrt21 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.625 * sqrt21 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 2.5 * sqrt21 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + sqrt21 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -0.3125 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.9375 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 5.625 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.9375 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 11.25 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 7.5 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.3125 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 5.625 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt21 = std::sqrt(21.0);

    return 0.625 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 1.25 * sqrt21 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 2.5 * sqrt21 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.625 * sqrt21 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 2.5 * sqrt21 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + sqrt21 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt210 = std::sqrt(210.0);

    return 0.03125 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.03125 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.5 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.03125 * sqrt210 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 0.5 * sqrt210 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.03125 * sqrt210 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.5 * sqrt210 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.5 * sqrt210 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt210 = std::sqrt(210.0);

    return -0.1875 * sqrt210 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.375 * sqrt210 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.5 * sqrt210 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.5625 * sqrt210 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 1.5 * sqrt210 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_6_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt7 = std::sqrt(7.0);

    return -0.1875 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.9375 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 1.875 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.9375 * sqrt7 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 11.25 * sqrt7 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.1875 * sqrt7 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.875 * sqrt7 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_6_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt154 = std::sqrt(154.0);

    return 0.1875 * sqrt154 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 1.875 * sqrt154 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.9375 * sqrt154 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_6_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt462 = std::sqrt(462.0);

    return 0.03125 * sqrt462 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.46875 * sqrt462 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.46875 * sqrt462 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.03125 * sqrt462 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_7_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt429 = std::sqrt(429.0);

    return 0.21875 * sqrt429 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 1.09375 * sqrt429 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.65625 * sqrt429 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.03125 * sqrt429 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_7_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt6006 = std::sqrt(6006.0);

    return 0.1875 * sqrt6006 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.625 * sqrt6006 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.1875 * sqrt6006 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_7_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt231 = std::sqrt(231.0);

    return -0.15625 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.15625 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.875 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.28125 * sqrt231 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 3.75 * sqrt231 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.03125 * sqrt231 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.375 * sqrt231 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_7_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt231 = std::sqrt(231.0);

    return -0.75 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 2.5 * sqrt231 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.75 * sqrt231 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 2.5 * sqrt231 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt21 = std::sqrt(21.0);

    return 0.28125 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.46875 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 5.625 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.09375 * sqrt21 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 3.75 * sqrt21 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 7.5 * sqrt21 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.09375 * sqrt21 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.875 * sqrt21 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 2.5 * sqrt21 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt42 = std::sqrt(42.0);

    return 0.9375 * sqrt42 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 1.875 * sqrt42 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 5.0 * sqrt42 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.9375 * sqrt42 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 5.0 * sqrt42 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 3.0 * sqrt42 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt7 = std::sqrt(7.0);

    return -0.15625 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.46875 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 3.75 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.46875 * sqrt7 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 7.5 * sqrt7 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt7 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.15625 * sqrt7 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.75 * sqrt7 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt7 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 2.0 * sqrt7 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -2.1875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 6.5625 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 13.125 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 6.5625 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 26.25 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 10.5 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 2.1875 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 13.125 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 10.5 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt7 = std::sqrt(7.0);

    return -0.15625 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.46875 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 3.75 * sqrt7 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.46875 * sqrt7 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 7.5 * sqrt7 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt7 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.15625 * sqrt7 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.75 * sqrt7 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt7 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 2.0 * sqrt7 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt42 = std::sqrt(42.0);

    return 0.46875 * sqrt42 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.46875 * sqrt42 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 2.5 * sqrt42 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.46875 * sqrt42 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 1.5 * sqrt42 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 0.46875 * sqrt42 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 2.5 * sqrt42 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 1.5 * sqrt42 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt21 = std::sqrt(21.0);

    return 0.09375 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.09375 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 1.875 * sqrt21 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.46875 * sqrt21 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 3.75 * sqrt21 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 2.5 * sqrt21 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.28125 * sqrt21 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 5.625 * sqrt21 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt21 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt231 = std::sqrt(231.0);

    return -0.1875 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.9375 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.625 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.9375 * sqrt231 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 3.75 * sqrt231 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.1875 * sqrt231 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.625 * sqrt231 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_7_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt231 = std::sqrt(231.0);

    return -0.03125 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.28125 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.375 * sqrt231 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.15625 * sqrt231 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 3.75 * sqrt231 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.15625 * sqrt231 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.875 * sqrt231 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_7_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt6006 = std::sqrt(6006.0);

    return 0.03125 * sqrt6006 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.46875 * sqrt6006 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.46875 * sqrt6006 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.03125 * sqrt6006 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_7_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt429 = std::sqrt(429.0);

    return 0.03125 * sqrt429 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.65625 * sqrt429 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 1.09375 * sqrt429 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.21875 * sqrt429 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_8_m_n8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt715 = std::sqrt(715.0);

    return 0.1875 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 1.3125 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.3125 * sqrt715 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.1875 * sqrt715 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_8_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt715 = std::sqrt(715.0);

    return 0.65625 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 3.28125 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.96875 * sqrt715 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.09375 * sqrt715 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_8_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt858 = std::sqrt(858.0);

    return -0.09375 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.21875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.3125 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.21875 * sqrt858 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 4.375 * sqrt858 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.09375 * sqrt858 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.3125 * sqrt858 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_8_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1001 = std::sqrt(1001.0);

    return -0.46875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.46875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.84375 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.75 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.09375 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.375 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt77 = std::sqrt(77.0);

    return 0.1875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.1875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 4.5 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.1875 * sqrt77 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 7.5 * sqrt77 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.1875 * sqrt77 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 4.5 * sqrt77 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt77 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1155 = std::sqrt(1155.0);

    return 0.28125 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.46875 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 1.875 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.09375 * sqrt1155 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.25 * sqrt1155 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 1.5 * sqrt1155 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.09375 * sqrt1155 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.625 * sqrt1155 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.5 * sqrt1155 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return -0.09375 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.28125 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 2.8125 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.28125 * sqrt70 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 5.625 * sqrt70 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt70 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.09375 * sqrt70 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 2.8125 * sqrt70 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 7.5 * sqrt70 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt70 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_n1(double AB_x, double AB_y, double AB_z)
{
    return -3.28125 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 9.84375 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 26.25 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 9.84375 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 52.5 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 31.5 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 3.28125 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 26.25 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 31.5 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 6.0 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 0.2734375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 1.09375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 8.75 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 1.640625 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 26.25 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 26.25 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 1.09375 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 26.25 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 52.5 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 14.0 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.2734375 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 8.75 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 26.25 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 14.0 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p1(double AB_x, double AB_y, double AB_z)
{
    return -3.28125 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 9.84375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 26.25 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 9.84375 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 52.5 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 31.5 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 3.28125 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 26.25 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 31.5 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 6.0 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt70 = std::sqrt(70.0);

    return -0.046875 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.09375 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 1.40625 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 1.40625 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 3.75 * sqrt70 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.09375 * sqrt70 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt70 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 1.5 * sqrt70 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.046875 * sqrt70 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt70 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.75 * sqrt70 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 1.5 * sqrt70 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1155 = std::sqrt(1155.0);

    return 0.09375 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.09375 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 0.625 * sqrt1155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.46875 * sqrt1155 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 1.25 * sqrt1155 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 0.5 * sqrt1155 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 0.28125 * sqrt1155 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.875 * sqrt1155 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 1.5 * sqrt1155 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt77 = std::sqrt(77.0);

    return 0.046875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.1875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 1.125 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.46875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 5.625 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 1.875 * sqrt77 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.1875 * sqrt77 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 5.625 * sqrt77 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 11.25 * sqrt77 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.046875 * sqrt77 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.125 * sqrt77 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 1.875 * sqrt77 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1001 = std::sqrt(1001.0);

    return -0.09375 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.84375 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.375 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.46875 * sqrt1001 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 3.75 * sqrt1001 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.46875 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.875 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_8_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt858 = std::sqrt(858.0);

    return -0.015625 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.21875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.21875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 3.28125 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.21875 * sqrt858 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.28125 * sqrt858 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.015625 * sqrt858 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.21875 * sqrt858 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_8_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt715 = std::sqrt(715.0);

    return 0.09375 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 1.96875 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 3.28125 * sqrt715 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt715 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_8_m_p8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt715 = std::sqrt(715.0);

    return 0.0234375 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.65625 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 1.640625 * sqrt715 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.65625 * sqrt715 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.0234375 * sqrt715 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_9_m_n9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt24310 = std::sqrt(24310.0);

    return 0.03515625 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.328125 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.4921875 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.140625 * sqrt24310 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.00390625 * sqrt24310 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_9_m_n8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt12155 = std::sqrt(12155.0);

    return 0.1875 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 1.3125 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt12155 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.1875 * sqrt12155 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_9_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1430 = std::sqrt(1430.0);

    return -0.08203125 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.328125 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.3125 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.1640625 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 6.5625 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.234375 * sqrt1430 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.9375 * sqrt1430 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.01171875 * sqrt1430 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.1875 * sqrt1430 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_9_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return -0.28125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.65625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.65625 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 4.375 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.28125 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt286 = std::sqrt(286.0);

    return 0.1171875 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 3.28125 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.328125 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 3.28125 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.1875 * sqrt286 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 5.90625 * sqrt286 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt286 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.0234375 * sqrt286 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.65625 * sqrt286 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 1.3125 * sqrt286 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5005 = std::sqrt(5005.0);

    return 0.1875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.1875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 1.5 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 0.1875 * sqrt5005 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.5 * sqrt5005 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.1875 * sqrt5005 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.5 * sqrt5005 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 1.5 * sqrt5005 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2310 = std::sqrt(2310.0);

    return -0.0234375 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.0625 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.84375 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.046875 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 1.40625 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 2.8125 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.28125 * sqrt2310 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 1.875 * sqrt2310 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 1.5 * sqrt2310 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0078125 * sqrt2310 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.28125 * sqrt2310 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.9375 * sqrt2310 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.5 * sqrt2310 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt110 = std::sqrt(110.0);

    return -0.65625 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 1.96875 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 6.5625 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 1.96875 * sqrt110 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 13.125 * sqrt110 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 10.5 * sqrt110 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.65625 * sqrt110 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 6.5625 * sqrt110 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 10.5 * sqrt110 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt110 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5 = std::sqrt(5.0);

    return 0.1640625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.65625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 6.5625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.984375 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 26.25 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.65625 * sqrt5 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt5 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 52.5 * sqrt5 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 21.0 * sqrt5 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.1640625 * sqrt5 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 6.5625 * sqrt5 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 26.25 * sqrt5 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 21.0 * sqrt5 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt5 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 2.4609375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 9.84375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 26.25 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 14.765625 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 78.75 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 47.25 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 9.84375 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 78.75 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 94.5 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 18.0 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 2.4609375 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 26.25 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 47.25 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 18.0 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5 = std::sqrt(5.0);

    return 0.1640625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.65625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 6.5625 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.984375 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 26.25 * sqrt5 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.65625 * sqrt5 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt5 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 52.5 * sqrt5 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 21.0 * sqrt5 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.1640625 * sqrt5 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 6.5625 * sqrt5 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 26.25 * sqrt5 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 21.0 * sqrt5 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt5 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt110 = std::sqrt(110.0);

    return -0.328125 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.65625 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 3.28125 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 3.28125 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 5.25 * sqrt110 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 0.65625 * sqrt110 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt110 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 1.5 * sqrt110 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.328125 * sqrt110 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt110 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 5.25 * sqrt110 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 1.5 * sqrt110 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2310 = std::sqrt(2310.0);

    return -0.0078125 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.28125 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.046875 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.28125 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.9375 * sqrt2310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.0625 * sqrt2310 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt2310 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 1.875 * sqrt2310 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.5 * sqrt2310 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0234375 * sqrt2310 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.84375 * sqrt2310 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 2.8125 * sqrt2310 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 1.5 * sqrt2310 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt5005 = std::sqrt(5005.0);

    return 0.046875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.1875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 0.375 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.46875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 1.875 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 0.375 * sqrt5005 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 0.1875 * sqrt5005 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.875 * sqrt5005 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 2.25 * sqrt5005 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.046875 * sqrt5005 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.375 * sqrt5005 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.375 * sqrt5005 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt286 = std::sqrt(286.0);

    return 0.0234375 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.1875 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.65625 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.328125 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 5.90625 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 1.3125 * sqrt286 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 3.28125 * sqrt286 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt286 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.1171875 * sqrt286 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 3.28125 * sqrt286 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt286 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return -0.046875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.65625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.21875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 3.28125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.65625 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 3.28125 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.046875 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.21875 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_9_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1430 = std::sqrt(1430.0);

    return -0.01171875 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.234375 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.1875 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.1640625 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 3.9375 * sqrt1430 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.328125 * sqrt1430 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 6.5625 * sqrt1430 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.08203125 * sqrt1430 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.3125 * sqrt1430 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_9_m_p8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt12155 = std::sqrt(12155.0);

    return 0.0234375 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.65625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 1.640625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt12155 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.0234375 * sqrt12155 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_9_m_p9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt24310 = std::sqrt(24310.0);

    return 0.00390625 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.140625 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.4921875 * sqrt24310 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.328125 * sqrt24310 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.03515625 * sqrt24310 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_10_m_n10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt92378 = std::sqrt(92378.0);

    return 0.01953125 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.234375 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.4921875 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.234375 * sqrt92378 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.01953125 * sqrt92378 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_10_m_n9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt461890 = std::sqrt(461890.0);

    return 0.03515625 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.328125 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.4921875 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.140625 * sqrt461890 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.00390625 * sqrt461890 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_10_m_n8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt12155 = std::sqrt(12155.0);

    return -0.03125 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.1875 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.5625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 3.9375 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.1875 * sqrt12155 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.9375 * sqrt12155 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.03125 * sqrt12155 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.5625 * sqrt12155 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_10_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt72930 = std::sqrt(72930.0);

    return -0.08203125 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.328125 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.4375 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.1640625 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 2.1875 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.234375 * sqrt72930 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt72930 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.01171875 * sqrt72930 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.0625 * sqrt72930 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return 0.03515625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.046875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 1.125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.1640625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 2.625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 2.625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.046875 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 2.625 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 8.75 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.03515625 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.125 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 2.625 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt858 = std::sqrt(858.0);

    return 0.5859375 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 5.46875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 1.640625 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 5.46875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.9375 * sqrt858 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 9.84375 * sqrt858 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 13.125 * sqrt858 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.1171875 * sqrt858 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.09375 * sqrt858 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 1.3125 * sqrt858 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2145 = std::sqrt(2145.0);

    return -0.03125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.0625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 5.25 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.0625 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.5 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.03125 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.3125 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 5.25 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 3.5 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return -0.1640625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.4375 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 1.96875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 0.328125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 3.28125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.65625 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 2.625 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 1.5 * sqrt4290 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0546875 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 1.3125 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.5 * sqrt4290 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt165 = std::sqrt(165.0);

    return 0.0546875 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.21875 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 2.625 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.328125 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 7.875 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 13.125 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.21875 * sqrt165 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 7.875 * sqrt165 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 26.25 * sqrt165 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 14.0 * sqrt165 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0546875 * sqrt165 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 2.625 * sqrt165 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 13.125 * sqrt165 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 14.0 * sqrt165 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt165 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt55 = std::sqrt(55.0);

    return 0.4921875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 1.96875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 2.953125 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 19.6875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 1.96875 * sqrt55 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 19.6875 * sqrt55 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 31.5 * sqrt55 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 9.0 * sqrt55 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.4921875 * sqrt55 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt55 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt55 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 9.0 * sqrt55 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt55 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -0.24609375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 1.23046875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 12.3046875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 2.4609375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 49.21875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 65.625 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 2.4609375 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 73.828125 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 196.875 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 78.75 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 1.23046875 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 49.21875 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 196.875 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 157.5 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 22.5 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.24609375 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 12.3046875 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 65.625 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 78.75 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 22.5 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt55 = std::sqrt(55.0);

    return 0.4921875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 1.96875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 6.5625 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 2.953125 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 19.6875 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt55 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 1.96875 * sqrt55 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 19.6875 * sqrt55 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 31.5 * sqrt55 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 9.0 * sqrt55 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.4921875 * sqrt55 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt55 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt55 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 9.0 * sqrt55 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt55 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt165 = std::sqrt(165.0);

    return 0.02734375 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.08203125 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 1.3125 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.0546875 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 2.625 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.0546875 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 6.5625 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 7.0 * sqrt165 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.08203125 * sqrt165 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 2.625 * sqrt165 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt165 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 1.5 * sqrt165 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.02734375 * sqrt165 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.3125 * sqrt165 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt165 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 7.0 * sqrt165 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 1.5 * sqrt165 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return -0.0546875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.65625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.328125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 1.3125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 0.4375 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 2.625 * sqrt4290 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.5 * sqrt4290 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.1640625 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.96875 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 3.9375 * sqrt4290 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 1.5 * sqrt4290 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2145 = std::sqrt(2145.0);

    return -0.0078125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.0234375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.328125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.109375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.109375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 3.28125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.875 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0234375 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.3125 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 5.25 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.0078125 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.328125 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 1.3125 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.875 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt858 = std::sqrt(858.0);

    return 0.1171875 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.9375 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 1.09375 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 1.640625 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 9.84375 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 1.3125 * sqrt858 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 5.46875 * sqrt858 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 13.125 * sqrt858 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.5859375 * sqrt858 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 5.46875 * sqrt858 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt858 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4290 = std::sqrt(4290.0);

    return 0.005859375 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.076171875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.1875 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.08203125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 2.625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 0.4375 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.08203125 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 6.5625 * sqrt4290 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.076171875 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 2.625 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt4290 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.005859375 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.1875 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.4375 * sqrt4290 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt72930 = std::sqrt(72930.0);

    return -0.01171875 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.234375 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.0625 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.1640625 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 1.3125 * sqrt72930 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.328125 * sqrt72930 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 2.1875 * sqrt72930 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.08203125 * sqrt72930 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.4375 * sqrt72930 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_10_m_p8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt12155 = std::sqrt(12155.0);

    return -0.00390625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.10546875 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.0703125 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.1640625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 1.96875 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.1640625 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 4.921875 * sqrt12155 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.10546875 * sqrt12155 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.96875 * sqrt12155 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.00390625 * sqrt12155 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.0703125 * sqrt12155 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_10_m_p9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt461890 = std::sqrt(461890.0);

    return 0.00390625 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.140625 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.4921875 * sqrt461890 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.328125 * sqrt461890 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.03515625 * sqrt461890 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_10_m_p10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt92378 = std::sqrt(92378.0);

    return 0.001953125 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.087890625 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.41015625 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.41015625 * sqrt92378 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.087890625 * sqrt92378 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.001953125 * sqrt92378 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_11_m_n11(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt88179 = std::sqrt(88179.0);

    return 0.021484375 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.322265625 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.90234375 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.64453125 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.107421875 * sqrt88179 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.001953125 * sqrt88179 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_11_m_n10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1939938 = std::sqrt(1939938.0);

    return 0.01953125 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.234375 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.4921875 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.234375 * sqrt1939938 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.01953125 * sqrt1939938 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_11_m_n9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt46189 = std::sqrt(46189.0);

    return -0.017578125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.146484375 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.3515625 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.08203125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 3.28125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.17578125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 4.921875 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.068359375 * sqrt46189 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt46189 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.001953125 * sqrt46189 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.0390625 * sqrt46189 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_11_m_n8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt692835 = std::sqrt(692835.0);

    return -0.03125 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.1875 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.1875 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 1.3125 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.1875 * sqrt692835 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt692835 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.03125 * sqrt692835 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.1875 * sqrt692835 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt36465 = std::sqrt(36465.0);

    return 0.013671875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.041015625 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 0.4921875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.08203125 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 1.96875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 1.3125 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.01171875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.984375 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.037109375 * sqrt36465 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt36465 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.9375 * sqrt36465 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.001953125 * sqrt36465 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.0703125 * sqrt36465 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.1875 * sqrt36465 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt14586 = std::sqrt(14586.0);

    return 0.17578125 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.234375 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 1.875 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 0.8203125 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 4.375 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 2.625 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.234375 * sqrt14586 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 4.375 * sqrt14586 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 8.75 * sqrt14586 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.17578125 * sqrt14586 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.875 * sqrt14586 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 2.625 * sqrt14586 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt143 = std::sqrt(143.0);

    return -0.146484375 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.146484375 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 7.03125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.41015625 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 32.8125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.64453125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 32.8125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.205078125 * sqrt143 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 11.25 * sqrt143 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 59.0625 * sqrt143 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 52.5 * sqrt143 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.029296875 * sqrt143 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.40625 * sqrt143 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt143 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 5.25 * sqrt143 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1001 = std::sqrt(1001.0);

    return -0.46875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.9375 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 6.5625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 15.75 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.9375 * sqrt1001 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt1001 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 7.5 * sqrt1001 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.46875 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt1001 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt30030 = std::sqrt(30030.0);

    return 0.005859375 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.021484375 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 0.328125 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.02734375 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.875 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 1.96875 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.01171875 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.65625 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.28125 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 2.625 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.001953125 * sqrt30030 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.65625 * sqrt30030 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 1.75 * sqrt30030 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.75 * sqrt30030 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.001953125 * sqrt30030 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.109375 * sqrt30030 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.65625 * sqrt30030 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.875 * sqrt30030 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.25 * sqrt30030 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2145 = std::sqrt(2145.0);

    return 0.1640625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.65625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 2.625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.984375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 7.875 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 7.875 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.65625 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 7.875 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 15.75 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 6.0 * sqrt2145 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.1640625 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 2.625 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 7.875 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 6.0 * sqrt2145 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt2145 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt66 = std::sqrt(66.0);

    return -0.041015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.205078125 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 2.4609375 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.41015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 9.84375 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 16.40625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.41015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 14.765625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 49.21875 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.205078125 * sqrt66 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 9.84375 * sqrt66 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 49.21875 * sqrt66 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 52.5 * sqrt66 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt66 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.041015625 * sqrt66 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 2.4609375 * sqrt66 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 16.40625 * sqrt66 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt66 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt66 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt66 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p0(double AB_x, double AB_y, double AB_z)
{
    return -2.70703125 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 13.53515625 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 45.1171875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 27.0703125 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 180.46875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 144.375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 27.0703125 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 270.703125 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 433.125 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 123.75 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 13.53515625 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 180.46875 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 433.125 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 247.5 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 27.5 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 2.70703125 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 45.1171875 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 144.375 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 123.75 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 27.5 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt66 = std::sqrt(66.0);

    return -0.041015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.205078125 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 2.4609375 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.41015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 9.84375 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 16.40625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.41015625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 14.765625 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 49.21875 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt66 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.205078125 * sqrt66 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 9.84375 * sqrt66 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 49.21875 * sqrt66 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 52.5 * sqrt66 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt66 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.041015625 * sqrt66 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 2.4609375 * sqrt66 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 16.40625 * sqrt66 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt66 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt66 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt66 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2145 = std::sqrt(2145.0);

    return 0.08203125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.24609375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 1.3125 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.1640625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 2.625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 3.9375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 0.1640625 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 3.9375 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 3.0 * sqrt2145 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.24609375 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 2.625 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt2145 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.5 * sqrt2145 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.08203125 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.3125 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 3.0 * sqrt2145 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.5 * sqrt2145 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt30030 = std::sqrt(30030.0);

    return 0.001953125 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.001953125 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.109375 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.01171875 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 0.65625 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.02734375 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.65625 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.65625 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.875 * sqrt30030 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.021484375 * sqrt30030 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.875 * sqrt30030 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 3.28125 * sqrt30030 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 1.75 * sqrt30030 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.25 * sqrt30030 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.005859375 * sqrt30030 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.328125 * sqrt30030 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 1.96875 * sqrt30030 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 2.625 * sqrt30030 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.75 * sqrt30030 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1001 = std::sqrt(1001.0);

    return -0.1171875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.3515625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 1.640625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 1.640625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 1.640625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 16.40625 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 19.6875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 1.875 * sqrt1001 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.3515625 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 19.6875 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt1001 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.1171875 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.640625 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 1.875 * sqrt1001 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt143 = std::sqrt(143.0);

    return -0.029296875 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.205078125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 1.40625 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.64453125 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 11.25 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.41015625 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 19.6875 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 59.0625 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 5.25 * sqrt143 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.146484375 * sqrt143 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 32.8125 * sqrt143 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 52.5 * sqrt143 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.146484375 * sqrt143 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 7.03125 * sqrt143 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 32.8125 * sqrt143 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt143 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt14586 = std::sqrt(14586.0);

    return 0.029296875 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.380859375 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 0.3125 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.41015625 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 4.375 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 0.4375 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 0.41015625 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt14586 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.380859375 * sqrt14586 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 4.375 * sqrt14586 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt14586 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.029296875 * sqrt14586 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.3125 * sqrt14586 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.4375 * sqrt14586 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt36465 = std::sqrt(36465.0);

    return 0.001953125 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.037109375 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.0703125 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.01171875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 1.40625 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 0.1875 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.08203125 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.984375 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 3.9375 * sqrt36465 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.041015625 * sqrt36465 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.96875 * sqrt36465 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt36465 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.013671875 * sqrt36465 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.4921875 * sqrt36465 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 1.3125 * sqrt36465 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt692835 = std::sqrt(692835.0);

    return -0.00390625 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.10546875 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.0234375 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.1640625 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.1640625 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.640625 * sqrt692835 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.10546875 * sqrt692835 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.65625 * sqrt692835 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.00390625 * sqrt692835 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.0234375 * sqrt692835 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_11_m_p9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt46189 = std::sqrt(46189.0);

    return -0.001953125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.068359375 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.0390625 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.17578125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 1.40625 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 0.08203125 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 4.921875 * sqrt46189 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.146484375 * sqrt46189 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 3.28125 * sqrt46189 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.017578125 * sqrt46189 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.3515625 * sqrt46189 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_11_m_p10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1939938 = std::sqrt(1939938.0);

    return 0.001953125 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.087890625 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.41015625 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.41015625 * sqrt1939938 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.087890625 * sqrt1939938 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.001953125 * sqrt1939938 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_11_m_p11(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt88179 = std::sqrt(88179.0);

    return 0.001953125 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.107421875 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.64453125 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.90234375 * sqrt88179 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.322265625 * sqrt88179 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.021484375 * sqrt88179 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_12_m_n12(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1352078 = std::sqrt(1352078.0);

    return 0.005859375 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.107421875 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.38671875 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 0.38671875 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.107421875 * sqrt1352078 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.005859375 * sqrt1352078 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

inline double Y_ll_12_m_n11(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2028117 = std::sqrt(2028117.0);

    return 0.021484375 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.322265625 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.90234375 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.64453125 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.107421875 * sqrt2028117 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.001953125 * sqrt2028117 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_12_m_n10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt176358 = std::sqrt(176358.0);

    return -0.009765625 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.107421875 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.21484375 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.12890625 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 2.578125 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 0.12890625 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 5.4140625 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.107421875 * sqrt176358 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 2.578125 * sqrt176358 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.009765625 * sqrt176358 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.21484375 * sqrt176358 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_12_m_n9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt323323 = std::sqrt(323323.0);

    return -0.052734375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.439453125 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 0.3515625 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 0.24609375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.52734375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 4.921875 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.205078125 * sqrt323323 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.40625 * sqrt323323 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.005859375 * sqrt323323 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.0390625 * sqrt323323 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt138567 = std::sqrt(138567.0);

    return 0.0078125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.0390625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 0.3125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.046875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 1.875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 0.9375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.046875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 6.5625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.0390625 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.875 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.0078125 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.3125 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.9375 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt138567 = std::sqrt(138567.0);

    return 0.068359375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.205078125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 0.8203125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 0.41015625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 3.28125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 1.3125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.05859375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.640625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 6.5625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.185546875 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 2.34375 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 3.9375 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.009765625 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.1171875 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.1875 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4862 = std::sqrt(4862.0);

    return -0.029296875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.009765625 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 1.58203125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.17578125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 2.109375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 8.4375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z + 0.17578125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 7.3828125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 19.6875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 7.875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.009765625 * sqrt4862 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 2.109375 * sqrt4862 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 19.6875 * sqrt4862 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 26.25 * sqrt4862 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.029296875 * sqrt4862 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.58203125 * sqrt4862 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 8.4375 * sqrt4862 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 7.875 * sqrt4862 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt17017 = std::sqrt(17017.0);

    return -0.146484375 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 0.146484375 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 2.34375 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 0.41015625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.64453125 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 3.75 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.205078125 * sqrt17017 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.75 * sqrt17017 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 11.8125 * sqrt17017 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt17017 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.029296875 * sqrt17017 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.46875 * sqrt17017 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 1.3125 * sqrt17017 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.75 * sqrt17017 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2002 = std::sqrt(2002.0);

    return 0.029296875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y + 0.087890625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y - 1.875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z + 0.05859375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y - 3.75 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z + 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.05859375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 21.0 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.087890625 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.75 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 7.5 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.029296875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 21.0 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2002 = std::sqrt(2002.0);

    return 0.263671875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z + 0.966796875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z - 4.921875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z + 1.23046875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 17.71875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.52734375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 9.84375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 29.53125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 16.875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.087890625 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 5.90625 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 3.75 * sqrt2002 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.087890625 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 1.640625 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 5.90625 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 5.625 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 1.25 * sqrt2002 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3003 = std::sqrt(3003.0);

    return -0.01171875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y - 0.05859375 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y + 0.8203125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z - 0.1171875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y + 3.28125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z - 0.1171875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 4.921875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 19.6875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 13.125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.05859375 * sqrt3003 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.28125 * sqrt3003 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 19.6875 * sqrt3003 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt3003 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt3003 * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.01171875 * sqrt3003 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.8203125 * sqrt3003 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt3003 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 13.125 * sqrt3003 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt3003 * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt3003 * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_n1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt78 = std::sqrt(78.0);

    return -0.451171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z - 2.255859375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z + 9.0234375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z - 4.51171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 36.09375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 36.09375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 4.51171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 54.140625 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 108.28125 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 41.25 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 2.255859375 * sqrt78 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 36.09375 * sqrt78 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 108.28125 * sqrt78 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 82.5 * sqrt78 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 13.75 * sqrt78 * AB_x * AB_x * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.451171875 * sqrt78 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 9.0234375 * sqrt78 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 36.09375 * sqrt78 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 41.25 * sqrt78 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 13.75 * sqrt78 * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt78 * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p0(double AB_x, double AB_y, double AB_z)
{
    return 0.2255859375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 1.353515625 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 16.2421875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + (3465.0 / 1024.0) * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 81.2109375 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 135.3515625 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 4.51171875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 162.421875 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 541.40625 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 288.75 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + (3465.0 / 1024.0) * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 162.421875 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 812.109375 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 866.25 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 185.625 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 1.353515625 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 81.2109375 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 541.40625 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 866.25 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 371.25 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 33.0 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.2255859375 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 16.2421875 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 135.3515625 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 288.75 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 185.625 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 33.0 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p1(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt78 = std::sqrt(78.0);

    return -0.451171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 2.255859375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 9.0234375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 4.51171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 36.09375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 36.09375 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 4.51171875 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 54.140625 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 108.28125 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 41.25 * sqrt78 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 2.255859375 * sqrt78 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 36.09375 * sqrt78 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 108.28125 * sqrt78 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 82.5 * sqrt78 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 13.75 * sqrt78 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.451171875 * sqrt78 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 9.0234375 * sqrt78 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 36.09375 * sqrt78 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 41.25 * sqrt78 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 13.75 * sqrt78 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + sqrt78 * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p2(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt3003 = std::sqrt(3003.0);

    return -0.005859375 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.0234375 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.41015625 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.029296875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 1.23046875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 3.28125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.8203125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 6.5625 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 6.5625 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.029296875 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.8203125 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 3.75 * sqrt3003 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0234375 * sqrt3003 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.23046875 * sqrt3003 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 6.5625 * sqrt3003 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 6.5625 * sqrt3003 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.5 * sqrt3003 * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.005859375 * sqrt3003 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.41015625 * sqrt3003 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.28125 * sqrt3003 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 6.5625 * sqrt3003 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 3.75 * sqrt3003 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.5 * sqrt3003 * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p3(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2002 = std::sqrt(2002.0);

    return 0.087890625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.087890625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 1.640625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.52734375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 5.90625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z - 1.23046875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 9.84375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 5.90625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 5.625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.966796875 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 29.53125 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 11.25 * sqrt2002 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 1.25 * sqrt2002 * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.263671875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 4.921875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 17.71875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 16.875 * sqrt2002 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 3.75 * sqrt2002 * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p4(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2002 = std::sqrt(2002.0);

    return 0.00732421875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.0146484375 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.46875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + (-255.0 / 2048.0) * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 1.40625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 3.28125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 0.205078125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 6.5625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 5.25 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + (-255.0 / 2048.0) * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 6.5625 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 32.8125 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 1.875 * sqrt2002 * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.0146484375 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.40625 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 13.125 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 26.25 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 11.25 * sqrt2002 * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.00732421875 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.46875 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 3.28125 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 5.25 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 1.875 * sqrt2002 * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p5(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt17017 = std::sqrt(17017.0);

    return -0.029296875 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.205078125 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.46875 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z + 0.64453125 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 3.75 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 1.3125 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 0.41015625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 6.5625 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 11.8125 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.75 * sqrt17017 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.146484375 * sqrt17017 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 6.5625 * sqrt17017 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 7.5 * sqrt17017 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.146484375 * sqrt17017 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 2.34375 * sqrt17017 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 6.5625 * sqrt17017 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 3.75 * sqrt17017 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p6(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt4862 = std::sqrt(4862.0);

    return -0.0048828125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.05859375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.263671875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.1318359375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 3.427734375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z - 1.40625 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z - 3.69140625 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 19.6875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 1.3125 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.1318359375 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.69140625 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 19.6875 * sqrt4862 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z - 0.05859375 * sqrt4862 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 3.427734375 * sqrt4862 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 19.6875 * sqrt4862 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 19.6875 * sqrt4862 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z + 0.0048828125 * sqrt4862 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.263671875 * sqrt4862 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 1.40625 * sqrt4862 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 1.3125 * sqrt4862 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p7(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt138567 = std::sqrt(138567.0);

    return 0.009765625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.185546875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z - 0.1171875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.05859375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z + 2.34375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z + 0.1875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z * AB_z + 0.41015625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 1.640625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 3.9375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z + 0.205078125 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 6.5625 * sqrt138567 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z - 0.068359375 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.8203125 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 1.3125 * sqrt138567 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p8(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt138567 = std::sqrt(138567.0);

    return 0.0009765625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.025390625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y - 0.0390625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z + 0.0146484375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y + 1.0546875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 0.1171875 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z * AB_z + 0.08203125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.640625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 3.28125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.0146484375 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 1.640625 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 8.203125 * sqrt138567 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z - 0.025390625 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 1.0546875 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 3.28125 * sqrt138567 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z + 0.0009765625 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.0390625 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.1171875 * sqrt138567 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p9(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt323323 = std::sqrt(323323.0);

    return -0.005859375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z + 0.205078125 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.0390625 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z * AB_z - 0.52734375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 1.40625 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z * AB_z - 0.24609375 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 4.921875 * sqrt323323 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z + 0.439453125 * sqrt323323 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 3.28125 * sqrt323323 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z - 0.052734375 * sqrt323323 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.3515625 * sqrt323323 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z * AB_z;
}

inline double Y_ll_12_m_p10(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt176358 = std::sqrt(176358.0);

    return -0.0009765625 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x + 0.04296875 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + 0.021484375 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z * AB_z - 0.1611328125 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.966796875 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z * AB_z + 4.51171875 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.1611328125 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 4.51171875 * sqrt176358 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z - 0.04296875 * sqrt176358 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.966796875 * sqrt176358 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z + 0.0009765625 * sqrt176358 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.021484375 * sqrt176358 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z * AB_z;
}

inline double Y_ll_12_m_p11(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt2028117 = std::sqrt(2028117.0);

    return 0.001953125 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_z - 0.107421875 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_z + 0.64453125 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_z - 0.90234375 * sqrt2028117 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z + 0.322265625 * sqrt2028117 * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z - 0.021484375 * sqrt2028117 * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_z;
}

inline double Y_ll_12_m_p12(double AB_x, double AB_y, double AB_z)
{
    const auto sqrt1352078 = std::sqrt(1352078.0);

    return 0.00048828125 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x - 0.0322265625 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y + (495.0 / 2048.0) * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y - 0.451171875 * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + (495.0 / 2048.0) * sqrt1352078 * AB_x * AB_x * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y - 0.0322265625 * sqrt1352078 * AB_x * AB_x * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y + 0.00048828125 * sqrt1352078 * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y * AB_y;
}

}  // namespace harm

#endif /* RealSolidHarmonicAB_hpp */