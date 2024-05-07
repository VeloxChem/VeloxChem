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

#include "StringFormat.hpp"

#include <algorithm>
#include <cctype>
#include <iterator>
#include <sstream>

namespace fstr {  // fstr namespace

auto
upcase(const std::string& source) -> std::string
{
    std::string str;

    for (size_t i = 0; i < source.size(); i++)
    {
        str.push_back(std::toupper(source[i]));
    }

    return str;
}

auto
format(const std::string& source, const size_t width, const fmt_t aligment) -> std::string
{
    auto str = source;

    auto strwidth = source.size();

    if (strwidth > width)
    {
        str.erase(width, strwidth - width);
    }
    else
    {
        if (aligment == fmt_t::left) str.append(width - strwidth, ' ');

        if (aligment == fmt_t::center)
        {
            strwidth = (width - strwidth) / 2;

            str.insert(0, strwidth, ' ');

            strwidth = str.size();

            str.append(width - strwidth, ' ');
        }

        if (aligment == fmt_t::right) str.insert(0, width - strwidth, ' ');
    }

    return str;
}

auto
to_string(const double source, const size_t presicion, const size_t width, const fmt_t aligment) -> std::string
{
    std::stringstream ss;

    ss.setf(std::ios::fixed);

    ss.precision(presicion);

    ss << source;

    std::string str(ss.str());

    if (!(source < 0.0)) str.insert(0, " ");

    return fstr::format(str, width, aligment);
}

auto
to_string(const double source, const size_t presicion) -> std::string
{
    std::stringstream ss;

    ss.setf(std::ios::fixed);

    ss.precision(presicion);

    ss << source;

    return std::string(ss.str());
}

auto
to_string(const int64_t source, const size_t width, const fmt_t aligment) -> std::string
{
    auto str = std::to_string(source);

    if (source >= 0) str.insert(0, " ");

    return fstr::format(str, width, aligment);
}

auto
to_string(const bool source) -> std::string
{
    if (source) return std::string("True");

    return std::string("False");
}

auto
to_AngularMomentum(const std::string& label) -> int64_t
{
    auto key = fstr::upcase(label);

    if (key == std::string("S")) return 0;

    if (key == std::string("P")) return 1;

    if (key == std::string("D")) return 2;

    if (key == std::string("F")) return 3;

    if (key == std::string("G")) return 4;

    if (key == std::string("H")) return 5;

    if (key == std::string("I")) return 6;

    return -1;
}

auto
to_AngularMomentum(const int64_t angmom) -> std::string
{
    if (angmom == 0) return std::string("S");

    if (angmom == 1) return std::string("P");

    if (angmom == 2) return std::string("D");

    if (angmom == 3) return std::string("F");

    if (angmom == 4) return std::string("G");

    if (angmom == 5) return std::string("H");

    if (angmom == 6) return std::string("I");

    return std::string();
}

auto
to_TensorComponents(const int64_t torder) -> std::vector<std::string>
{
    if (torder == 1)
    {
        return {"X", "Y", "Z"};
    }

    if (torder == 2)
    {
        return {"XX", "XY", "XZ", "YY", "YZ", "ZZ"};
    }

    if (torder == 3)
    {
        return {"XXX", "XXY", "XXZ", "XYY", "XYZ", "XZZ", "YYY", "YYZ", "YZZ", "ZZZ"};
    }

    if (torder == 4)
    {
        return {"XXXX", "XXXY", "XXXZ", "XXYY", "XXYZ", "XXZZ", "XYYY", "XYYZ", "XYZZ", "XZZZ", "YYYY", "YYYZ", "YYZZ", "YZZZ", "ZZZZ"};
    }

    return std::vector<std::string>();
}

auto
to_TensorComponent(const std::string& tlabel) -> int64_t
{
    const auto index = fstr::upcase(tlabel);

    if (index.size() == 1)
    {
        if (index == "X") return 0;

        if (index == "Y") return 1;

        if (index == "Z") return 2;
    }

    if (index.size() == 2)
    {
        if (index == "XX") return 0;

        if (index == "XY") return 1;

        if (index == "XZ") return 2;

        if (index == "YY") return 3;

        if (index == "YZ") return 4;

        if (index == "ZZ") return 5;
    }

    if (index.size() == 3)
    {
        if (index == "XXX") return 0;

        if (index == "XXY") return 1;

        if (index == "XXZ") return 2;

        if (index == "XYY") return 3;

        if (index == "XYZ") return 4;

        if (index == "XZZ") return 5;

        if (index == "YYY") return 6;

        if (index == "YYZ") return 7;

        if (index == "YZZ") return 8;

        if (index == "ZZZ") return 9;
    }

    if (index.size() == 4)
    {
        if (index == "XXXX") return 0;

        if (index == "XXXY") return 1;

        if (index == "XXXZ") return 2;

        if (index == "XXYY") return 3;

        if (index == "XXYZ") return 4;

        if (index == "XXZZ") return 5;

        if (index == "XYYY") return 6;

        if (index == "XYYZ") return 7;

        if (index == "XYZZ") return 8;

        if (index == "XZZZ") return 9;

        if (index == "YYYY") return 10;

        if (index == "YYYZ") return 11;

        if (index == "YYZZ") return 12;

        if (index == "YZZZ") return 13;

        if (index == "ZZZZ") return 14;
    }

    return -1;
}

}  // namespace fstr
