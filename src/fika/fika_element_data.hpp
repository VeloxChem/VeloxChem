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

// Element data from AME2020 and NUBASE2020.
#ifndef fika_element_data_hpp
#define fika_element_data_hpp

#include <array>
#include <string_view>

namespace fika::detail {

struct ElementData {
  std::string_view label;
  double mass;  // Da
};

// Indexed by atomic number - 1.
inline constexpr std::array<ElementData, 118> element_data{{
    {"H", 1.007825031898},    //   1 1H, 99.9855% natural abundance
    {"He", 4.00260325413},    //   2 4He, 99.9998% natural abundance
    {"Li", 7.01600343426},    //   3 7Li, 95.15% natural abundance
    {"Be", 9.012183062},      //   4 9Be, 100% natural abundance
    {"B", 11.009305166},      //   5 11B, 80.35% natural abundance
    {"C", 12.0},              //   6 12C, 98.94% natural abundance
    {"N", 14.00307400425},    //   7 14N, 99.6205% natural abundance
    {"O", 15.99491461926},    //   8 16O, 99.757% natural abundance
    {"F", 18.99840316207},    //   9 19F, 100% natural abundance
    {"Ne", 19.99244017525},   //  10 20Ne, 90.48% natural abundance
    {"Na", 22.98976928195},   //  11 23Na, 100% natural abundance
    {"Mg", 23.985041689},     //  12 24Mg, 78.965% natural abundance
    {"Al", 26.981538408},     //  13 27Al, 100% natural abundance
    {"Si", 27.97692653442},   //  14 28Si, 92.2545% natural abundance
    {"P", 30.97376199768},    //  15 31P, 100% natural abundance
    {"S", 31.97207117354},    //  16 32S, 94.85% natural abundance
    {"Cl", 34.968852694},     //  17 35Cl, 75.8% natural abundance
    {"Ar", 39.96238312204},   //  18 40Ar, 99.6035% natural abundance
    {"K", 38.96370648482},    //  19 39K, 93.2581% natural abundance
    {"Ca", 39.96259085},      //  20 40Ca, 96.941% natural abundance
    {"Sc", 44.955907051},     //  21 45Sc, 100% natural abundance
    {"Ti", 47.947940677},     //  22 48Ti, 73.72% natural abundance
    {"V", 50.943957664},      //  23 51V, 99.75% natural abundance
    {"Cr", 51.940504714},     //  24 52Cr, 83.789% natural abundance
    {"Mn", 54.93804304},      //  25 55Mn, 100% natural abundance
    {"Fe", 55.934935537},     //  26 56Fe, 91.754% natural abundance
    {"Co", 58.933193524},     //  27 59Co, 100% natural abundance
    {"Ni", 57.93534165},      //  28 58Ni, 68.0769% natural abundance
    {"Cu", 62.929597119},     //  29 63Cu, 69.15% natural abundance
    {"Zn", 63.929141776},     //  30 64Zn, 49.17% natural abundance
    {"Ga", 68.925573528},     //  31 69Ga, 60.108% natural abundance
    {"Ge", 73.92117776},      //  32 74Ge, 36.52% natural abundance
    {"As", 74.921594562},     //  33 75As, 100% natural abundance
    {"Se", 79.916521761},     //  34 80Se, 49.8% natural abundance
    {"Br", 78.918337574},     //  35 79Br, 50.65% natural abundance
    {"Kr", 83.91149772708},   //  36 84Kr, 56.987% natural abundance
    {"Rb", 84.91178973604},   //  37 85Rb, 72.17% natural abundance
    {"Sr", 87.905612253},     //  38 88Sr, 82.58% natural abundance
    {"Y", 88.905838156},      //  39 89Y, 100% natural abundance
    {"Zr", 89.904698755},     //  40 90Zr, 51.45% natural abundance
    {"Nb", 92.90637317},      //  41 93Nb, 100% natural abundance
    {"Mo", 97.905403609},     //  42 98Mo, 24.292% natural abundance
    {"Tc", 96.90636072},      //  43 97Tc, longest lived, t1/2 = 4.21e+06 y
    {"Ru", 101.904340312},    //  44 102Ru, 31.55% natural abundance
    {"Rh", 102.905494081},    //  45 103Rh, 100% natural abundance
    {"Pd", 105.903480287},    //  46 106Pd, 27.33% natural abundance
    {"Ag", 106.905091509},    //  47 107Ag, 51.839% natural abundance
    {"Cd", 113.903364998},    //  48 114Cd, 28.754% natural abundance
    {"In", 114.903878772},    //  49 115In, 95.719% natural abundance
    {"Sn", 119.902202557},    //  50 120Sn, 32.58% natural abundance
    {"Sb", 120.903811353},    //  51 121Sb, 57.21% natural abundance
    {"Te", 129.906222745},    //  52 130Te, 34.08% natural abundance
    {"I", 126.904472592},     //  53 127I, 100% natural abundance
    {"Xe", 131.90415508346},  //  54 132Xe, 26.909% natural abundance
    {"Cs", 132.905451958},    //  55 133Cs, 100% natural abundance
    {"Ba", 137.905247059},    //  56 138Ba, 71.7% natural abundance
    {"La", 138.906362927},    //  57 139La, 99.9112% natural abundance
    {"Ce", 139.905448433},    //  58 140Ce, 88.449% natural abundance
    {"Pr", 140.907659604},    //  59 141Pr, 100% natural abundance
    {"Nd", 141.907728824},    //  60 142Nd, 27.153% natural abundance
    {"Pm", 144.912755748},    //  61 145Pm, longest lived, t1/2 = 17.7 y
    {"Sm", 151.919738646},    //  62 152Sm, 26.74% natural abundance
    {"Eu", 152.921236789},    //  63 153Eu, 52.19% natural abundance
    {"Gd", 157.9241112},      //  64 158Gd, 24.84% natural abundance
    {"Tb", 158.925353707},    //  65 159Tb, 100% natural abundance
    {"Dy", 163.929180819},    //  66 164Dy, 28.26% natural abundance
    {"Ho", 164.930329116},    //  67 165Ho, 100% natural abundance
    {"Er", 165.930301067},    //  68 166Er, 33.503% natural abundance
    {"Tm", 168.934218956},    //  69 169Tm, 100% natural abundance
    {"Yb", 173.938867545},    //  70 174Yb, 32.025% natural abundance
    {"Lu", 174.940777211},    //  71 175Lu, 97.401% natural abundance
    {"Hf", 179.946559537},    //  72 180Hf, 35.08% natural abundance
    {"Ta", 180.947998528},    //  73 181Ta, 99.988% natural abundance
    {"W", 183.95093318},      //  74 184W, 30.64% natural abundance
    {"Re", 186.955752217},    //  75 187Re, 62.6% natural abundance
    {"Os", 191.961478765},    //  76 192Os, 40.78% natural abundance
    {"Ir", 192.962923753},    //  77 193Ir, 62.7% natural abundance
    {"Pt", 194.964794325},    //  78 195Pt, 33.775% natural abundance
    {"Au", 196.966570103},    //  79 197Au, 100% natural abundance
    {"Hg", 201.970643604},    //  80 202Hg, 29.74% natural abundance
    {"Tl", 204.974427318},    //  81 205Tl, 70.485% natural abundance
    {"Pb", 207.976652005},    //  82 208Pb, 52.4% natural abundance
    {"Bi", 208.980398599},    //  83 209Bi, 100% natural abundance
    {"Po", 208.982430361},    //  84 209Po, longest lived, t1/2 = 124 y
    {"At", 209.987147423},    //  85 210At, longest lived, t1/2 = 0.000924 y
    {"Rn", 222.017576017},    //  86 222Rn, longest lived, t1/2 = 0.0105 y
    {"Fr", 223.019734241},    //  87 223Fr, longest lived, t1/2 = 4.18e-05 y
    {"Ra", 226.025408186},    //  88 226Ra, longest lived, t1/2 = 1.6e+03 y
    {"Ac", 227.027750594},    //  89 227Ac, longest lived, t1/2 = 21.8 y
    {"Th", 232.038053606},    //  90 232Th, 99.98% natural abundance
    {"Pa", 231.0358825},      //  91 231Pa, 100% natural abundance
    {"U", 238.050786936},     //  92 238U, 99.2742% natural abundance
    {"Np", 237.04817164},     //  93 237Np, longest lived, t1/2 = 2.14e+06 y
    {"Pu", 244.064204401},    //  94 244Pu, longest lived, t1/2 = 8.13e+07 y
    {"Am", 243.061379889},    //  95 243Am, longest lived, t1/2 = 7.35e+03 y
    {"Cm", 247.070352678},    //  96 247Cm, longest lived, t1/2 = 1.56e+07 y
    {"Bk", 247.070305889},    //  97 247Bk, longest lived, t1/2 = 1.38e+03 y
    {"Cf", 251.079587171},    //  98 251Cf, longest lived, t1/2 = 898 y
    {"Es", 252.082979173},    //  99 252Es, longest lived, t1/2 = 1.29 y
    {"Fm", 257.095105419},    // 100 257Fm, longest lived, t1/2 = 0.275 y
    {"Md", 258.098433634},    // 101 258Md, longest lived, t1/2 = 0.141 y
    {"No", 259.100998364},    // 102 259No, longest lived, t1/2 = 0.00011 y
    {"Lr", 266.119874},       // 103 266Lr, longest lived, t1/2 = 0.00251 y
    {"Rf", 267.121787},       // 104 267Rf, longest lived, t1/2 = 0.000285 y
    {"Db", 268.125669},       // 105 268Db, longest lived, t1/2 = 0.00331 y
    {"Sg", 269.128495},       // 106 269Sg, longest lived, t1/2 = 9.51e-06 y
    {"Bh", 270.133366},       // 107 270Bh, longest lived, t1/2 = 7.23e-06 y
    {"Hs", 269.133649},       // 108 269Hs, longest lived, t1/2 = 4.75e-07 y
    {"Mt", 277.153525},       // 109 277Mt, longest lived, t1/2 = 2.85e-07 y
    {"Ds", 282.166174},       // 110 282Ds, longest lived, t1/2 = 7.99e-06 y
    {"Rg", 282.169343},       // 111 282Rg, longest lived, t1/2 = 4.12e-06 y
    {"Cn", 285.177227},       // 112 285Cn, longest lived, t1/2 = 9.51e-07 y
    {"Nh", 286.182456},       // 113 286Nh, longest lived, t1/2 = 3.8e-07 y
    {"Fl", 290.191875},       // 114 290Fl, longest lived, t1/2 = 2.54e-06 y
    {"Mc", 290.196235},       // 115 290Mc, longest lived, t1/2 = 2.66e-08 y
    {"Lv", 293.204583},       // 116 293Lv, longest lived, t1/2 = 2.22e-09 y
    {"Ts", 294.21084},        // 117 294Ts, longest lived, t1/2 = 2.22e-09 y
    {"Og", 295.216178},       // 118 295Og, longest lived, t1/2 = 2.15e-08 y
}};

}  // namespace fika::detail

#endif  // fika_element_data_hpp
