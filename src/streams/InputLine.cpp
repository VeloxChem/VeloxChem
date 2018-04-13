//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "InputLine.hpp"

#include <utility>
#include <regex>

#include "StringFormat.hpp"

CInputLine::CInputLine()
{

}

CInputLine::CInputLine(const std::string& string)

    : _originalString(string)

    , _parsedString(string)
{
    _trimParsedString();
}

CInputLine::CInputLine(const CInputLine& source)

    : _originalString(source._originalString)

    , _parsedString(source._parsedString)
{

}

CInputLine::CInputLine(CInputLine&& source) noexcept

    : _originalString(std::move(source._originalString))

    , _parsedString(std::move(source._parsedString))
{

}

CInputLine::~CInputLine()
{

}

CInputLine& CInputLine::operator=(const CInputLine& source)
{
    if (this == &source) return *this;

    _originalString = source._originalString;

    _parsedString = source._parsedString;

    return *this;
}

CInputLine& CInputLine::operator=(CInputLine&& source) noexcept
{
    if (this == &source) return *this;

    _originalString = std::move(source._originalString);

    _parsedString = std::move(source._parsedString);

    return *this;
}

bool CInputLine::operator==(const CInputLine& other) const
{
    if (_originalString != other._originalString) return false;

    if (_parsedString != other._parsedString) return false;

    return true;
}

bool CInputLine::operator!=(const CInputLine& other) const
{
    return !(*this == other);
}

std::string CInputLine::getKeyword(const size_t iKeyword) const
{
    auto str = _parsedString;

    for (size_t i = 0; i < iKeyword; i++)
    {
        auto pos = str.find_first_of(" ");

        if (pos == std::string::npos) return std::string();

        str.erase(0, pos);

        pos = str.find_first_not_of(" ");

        str.erase(0, pos);
    }

    auto pos = str.find_first_of(" ");

    if (pos != std::string::npos) str.erase(pos, std::string::npos);

    return str;
}

std::string CInputLine::getUpcasedKeyword(const size_t iKeyword) const
{
    return fstr::upcase(getKeyword(iKeyword));
}

size_t CInputLine::getNumberOfKeywords() const
{
    size_t cnt = 0;

    for (size_t i = 0; ; i++)
    {
        auto str = getKeyword(i);

        if (str.empty()) break;

        cnt++;
    }

    return cnt;
}

bool CInputLine::isRealNumber(const size_t iKeyword) const
{
    auto str = getKeyword(iKeyword);

    if (!str.empty())
    {
        auto pos = str.find_first_not_of("1234567890.+-Ee");

        if (pos == std::string::npos) return true;
    }

    return false;
}

double CInputLine::getRealNumber(const size_t iKeyword) const
{
    if (isRealNumber(iKeyword)) return std::stod(getKeyword(iKeyword));

    return std::numeric_limits<double>::quiet_NaN();
}

bool CInputLine::isIntegerNumber(const size_t iKeyword) const
{
    auto str = getKeyword(iKeyword);

    if (!str.empty())
    {
        auto pos = str.find_first_not_of("1234567890+-");

        if (pos == std::string::npos) return true;
    }

    return false;
}

int32_t CInputLine::getIntegerNumber(const size_t iKeyword) const
{
    if (isIntegerNumber(iKeyword))
    {
        return static_cast<int32_t>(std::stol(getKeyword(iKeyword)));
    }

    return std::numeric_limits<int32_t>::quiet_NaN();
}

bool CInputLine::isKeyword(const size_t iKeyword,
                           const std::string& keyLabel) const
{
    auto str = fstr::upcase(getKeyword(iKeyword));

    if (!str.empty())
    {
        if (str == fstr::upcase(keyLabel)) return true;
    }

    return false;
}

bool CInputLine::isKeyword(const size_t iKeyword, const char* keyLabel) const
{
    return isKeyword(iKeyword, std::string(keyLabel));
}

bool CInputLine::isControlKeyword(const std::string& keyLabel) const
{
    auto str = fstr::upcase(getKeyword(0));

    if (!str.empty())
    {
        auto rkey = fstr::upcase(keyLabel);

        rkey.insert(0, "@");

        if (str == rkey) return true;
    }

    return false;
}

bool CInputLine::isControlKeyword(const char* keyLabel) const
{
    return isControlKeyword(std::string(keyLabel));
}

bool CInputLine::isControlLine()const
{
    auto  str = getKeyword(0);

    if (!str.empty())
    {
        if (str[0] == '@') return true;
    }

    return false;
}

std::string CInputLine::getParsedString() const
{
    return _parsedString;
}

std::string CInputLine::getOriginalString() const
{
    return _originalString;
}

bool CInputLine::isEmpty() const
{
    return _parsedString.empty();
}

void CInputLine::clear()
{
    _parsedString.clear();

    _originalString.clear();
}

void CInputLine::_trimParsedString()
{
    std::smatch match;

    if (std::regex_search(_parsedString, match, std::regex("([^!#]*)[!#].*")))
    {
        _parsedString = match[1].str();
    }

    if (std::regex_search(_parsedString, match,
                          std::regex("\\s*(\\S.*\\S|\\S)\\s*")))
    {
        _parsedString = match[1].str();
    }
    else
    {
        _parsedString.clear();
    }
}

std::ostream& operator<<(std::ostream& output, const CInputLine& source)
{
    output << std::endl;

    output << "[CInputLine (Object):" << &source << "]" << std::endl;

    output << "_originalString: \"" << source._originalString << "\"";

    output <<  std::endl;

    output << "_parsedString: \"" << source._parsedString << "\"";

    output <<  std::endl;

    return output;
}
