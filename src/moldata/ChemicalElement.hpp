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

#ifndef ChemicalElement_hpp
#define ChemicalElement_hpp

#include <array>
#include <string>

namespace chem_elem {
namespace {
/// @brief Table of chemical element names.
static const std::array<std::string, 103> _names = {
    "BQ", "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE", "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR", "K",  "CA", "SC",
    "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y",  "ZR", "NB", "MO", "TC",
    "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",  "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR",
    "RA", "AC", "TH", "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO"};

/// @brief Table of chemical element masses (most abundant isotope).
static const std::array<double, 103> _masses = {
    0.000000,   1.007825,   4.002603,   7.016005,   9.012182,   11.009305,  12.000000,  14.003074,  15.994915,  18.998403,  19.992440,
    22.989769,  23.985042,  26.981539,  27.976927,  30.973762,  31.972071,  34.968853,  39.962383,  38.963707,  39.962591,  44.955912,
    47.947946,  50.943960,  51.940508,  54.938045,  55.934938,  58.933195,  57.935343,  62.929598,  63.929142,  68.925574,  73.921178,
    74.921597,  79.916521,  78.918337,  83.911507,  84.911790,  87.905612,  88.905848,  89.904704,  92.906378,  97.905408,  97.907216,
    101.904349, 102.905504, 105.903486, 106.905097, 113.903359, 114.903878, 119.902195, 120.903816, 129.906224, 126.904473, 131.904153,
    132.905452, 137.905247, 138.906353, 139.905439, 140.907653, 141.907723, 146.915139, 151.919732, 152.921230, 157.924104, 158.925347,
    163.929175, 164.930322, 165.930293, 168.934213, 173.938862, 174.940772, 179.946550, 180.947996, 183.950931, 186.955753, 191.961481,
    192.962926, 194.964791, 196.966569, 201.970643, 204.974428, 207.976652, 208.980399, 208.982430, 209.987148, 222.017578, 223.019731,
    226.025403, 227.027747, 232.038050, 231.035879, 238.050783, 237.048167, 244.064198, 243.061373, 247.070347, 247.070299, 251.079580,
    252.082972, 257.095099, 258.098425, 259.101024};
}  // namespace

/// @brief Checks if identifier is valid chemical element number.
/// @param id The identifier of chemical element to validate.
/// @return True if identifier is valid chemical element number, False
/// otherwise.
inline auto
valid_identifier(const int id) -> bool
{
    return (id >= 0) && (id <= 102);
}

/// @brief Gets name of chemical element.
/// @param id The identifier of chemicla element.
/// @return The name of chemical element.
auto name(const int id) -> std::string;

/// @brief Gets label of chemical element.
/// @param id The identifier of chemicla element.
/// @return The standart label of chemical element.
auto label(const int id) -> std::string;

/// @brief Gets identifier of chemical element.
/// @param name The upper cased name of chemical element.
/// @return The identifier of chemical element.
auto identifier(const std::string &name) -> int;

/// @brief Gets mass of chemical element's most abundant isotope.
/// @param id The identifier of chemicla element.
/// @return The mass of chemical element.
auto mass(const int id) -> double;

/// @brief Gets maximum angular momentum of atomic shell in chemical element.
/// @param id The identifier of chemicla element.
/// @return The maximum angular momentum of atomic shell.
auto max_angular_momentum(const int id) -> int;

/// @brief Gets maximum valid value of chemical element number.
/// @return The maximum value of chemicla element number.
auto max_identifier() -> int;

}  // namespace chem_elem

#endif /* ChemicalElement_hpp */
