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

}  // namespace fstr
