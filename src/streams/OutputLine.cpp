//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OutputLine.hpp"

#include <utility>

COutputLine::COutputLine()

    : _offset(0)

    , _leftSymbol(' ')

    , _rightSymbol(' ')

    , _fillSymbol(' ')

    , _width(0)
{

}

COutputLine::COutputLine(const std::string& line,
                         const size_t       offset,
                         const char         leftSymol,
                         const char         rightSymbol,
                         const char         fillSymbol,
                         const size_t       width)

    : _line(line)

    , _offset(offset)

    , _leftSymbol(leftSymol)

    , _rightSymbol(rightSymbol)

    , _fillSymbol(fillSymbol)

    , _width(width)
{

}

COutputLine::COutputLine(const COutputLine& source)

    : _line(source._line)

    , _offset(source._offset)

    , _leftSymbol(source._leftSymbol)

    , _rightSymbol(source._rightSymbol)

    , _fillSymbol(source._fillSymbol)

    , _width(source._width)
{

}

COutputLine::COutputLine(COutputLine&& source) noexcept

    : _line(std::move(source._line))

    , _offset(std::move(source._offset))

    , _leftSymbol(std::move(source._leftSymbol))

    , _rightSymbol(std::move(source._rightSymbol))

    , _fillSymbol(std::move(source._fillSymbol))

    , _width(std::move(source._width))
{

}

COutputLine::~COutputLine()
{

}

COutputLine&
COutputLine::operator=(const COutputLine& source)
{
    if (this == &source) return *this;

    _line = source._line;

    _offset = source._offset;

    _leftSymbol = source._leftSymbol;

    _rightSymbol = source._rightSymbol;

    _fillSymbol = source._fillSymbol;

    _width = source._width;

    return *this;
}

COutputLine&
COutputLine::operator=(COutputLine&& source) noexcept
{
    if (this == &source) return *this;

    _line = std::move(source._line);

    _offset = std::move(source._offset);

    _leftSymbol = std::move(source._leftSymbol);

    _rightSymbol = std::move(source._rightSymbol);

    _fillSymbol = std::move(source._fillSymbol);

    _width = std::move(source._width);

    return *this;
}

bool
COutputLine::operator==(const COutputLine& other) const
{
    if (_line != other._line) return false;

    if (_offset != other._offset) return false;

    if (_leftSymbol != other._leftSymbol) return false;

    if (_rightSymbol != other._rightSymbol) return false;

    if (_fillSymbol != other._fillSymbol) return false;

    if (_width != other._width) return false;

    return true;
}

bool
COutputLine::operator!=(const COutputLine& other) const
{
    return !(*this == other);
}

std::ofstream&
operator<<(      std::ofstream& output,
           const COutputLine&   source)
{
    std::string str(source._offset, source._fillSymbol);

    str.append(source._line);

    auto strwidth = str.size();

    if (strwidth <= source._width)
    {
        str.append(source._width - strwidth, source._fillSymbol);
    }
    else
    {
        str.erase(str.begin() + source._width, str.end());
    }

    str.insert(str.begin(), 1, source._leftSymbol);

    str.insert(str.end(), 1, source._rightSymbol);

    output << str << std::endl;

    return output;
}

std::ostream&
operator<<(      std::ostream& output,
           const COutputLine&  source)
{
    output << std::endl;

    output << "[COutputLine (Object):" << &source << "]" << std::endl;

    output << "_line: \"" << source._line << "\"" <<  std::endl;

    output << "_offset: " << source._offset << std::endl;

    output << "_leftSymbol: '" << source._leftSymbol << "'" << std::endl;

    output << "_rightSymbol: '" << source._rightSymbol <<"'"<< std::endl;

    output << "_fillSymbol: '" << source._fillSymbol << "'"<< std::endl;

    output << "_width: " << source._width << std::endl;

    return output;
}
