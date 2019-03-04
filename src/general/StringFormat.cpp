//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "StringFormat.hpp"

#include <sstream>

namespace fstr { // fstr namespace

std::string
upcase(const std::string& source)
{
    auto str = source;

    for (size_t i = 0; i < str.size(); i++)
    {
        str[i] = static_cast<char>(std::toupper(str[i]));
    }

    return str;
}

std::string
format(const std::string& source,
       const size_t       width,
       const fmt          aligment)
{
    auto str = source;

    auto strwidth = source.size();

    if (strwidth > width)
    {
        str.erase(width, strwidth - width);
    }
    else
    {
        if (aligment == fmt::left) str.append(width - strwidth, ' ');

        if (aligment == fmt::center)
        {
            strwidth = (width - strwidth) / 2;

            str.insert(0, strwidth, ' ');

            strwidth = str.size();

            str.append(width - strwidth, ' ');
        }

        if (aligment == fmt::right) str.insert(0, width - strwidth, ' ');
    }

    return str;
}

std::string
to_string(const char*  source,
          const size_t width,
          const fmt    aligment)
{
    return fstr::format(std::string(source), width, aligment);
}

std::string
to_string(const char*  source,
          const size_t width)
{
    return fstr::to_string(source, width, fmt::center);
}

std::string
to_string(const double source,
          const size_t presicion,
          const size_t width,
          const fmt    aligment)
{
    std::stringstream sStream;

    sStream.setf(std::ios::fixed);

    sStream.precision(presicion);

    sStream << source;

    std::string str(sStream.str());

    if (!(source < 0.0)) str.insert(0, " ");

    return fstr::format(str, width, aligment);
}
    
std::string
to_string(const double source,
          const size_t presicion)
{
    std::stringstream sStream;
    
    sStream.setf(std::ios::fixed);
        
    sStream.precision(presicion);
        
    sStream << source;
        
    return std::string(sStream.str());
}

std::string
to_string(const int32_t source,
          const size_t  width,
          const fmt     aligment)
{
    auto str = std::to_string(source);

    if (source >= 0) str.insert(0, " ");

    return fstr::format(str, width, aligment);
}

std::string
to_string(const size_t source,
          const size_t width,
          const fmt    aligment)
{
    auto str = std::to_string(source);

    return fstr::format(str, width, aligment);
}

std::string
to_string(const bool source)
{
    if (source) return std::string("True");
        
    return std::string("False");
}
    
int32_t
to_AngularMomentum(const std::string& label)
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

std::string
to_AngularMomentum(const int32_t angularmomentum)
{
    if (angularmomentum == 0) return std::string("S");

    if (angularmomentum == 1) return std::string("P");

    if (angularmomentum == 2) return std::string("D");

    if (angularmomentum == 3) return std::string("F");

    if (angularmomentum == 4) return std::string("G");

    if (angularmomentum == 5) return std::string("H");

    if (angularmomentum == 6) return std::string("I");

    return std::string();
}

} // fstr namespace
