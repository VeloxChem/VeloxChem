//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OutputStream.hpp"

#include <fstream>
#include <iostream>
#include <utility>

COutputStream::COutputStream(const std::string& filename)

    : _state(true)

    , _filename(filename)

    , _maxWidth(120)

    , _leftSymbol(' ')

    , _rightSymbol(' ')

    , _fillSymbol(' ')

    , _alignment(fmt::left)
{
    if (!_filename.empty())
    {
        std::ofstream ostream;

        ostream.open(_filename.c_str(), std::ios_base::trunc);

        _state = ostream.good();

        ostream.close();

        if (!_state) _errorFileOpen();
    }
}

COutputStream::~COutputStream()
{

}

bool COutputStream::getState() const
{
    return _state;
}

void COutputStream::flush()
{
    if (!_filename.empty())
    {
        std::ofstream ostream;

        ostream.open(_filename.c_str(), std::ios_base::app);

        _state = ostream.good();

        if (_state)
        {
            for (size_t i = 0;  i < _buffer.size(); i++) ostream << _buffer[i];
        }
        else
        {
            _errorFileOpen();
        }

        _buffer.clear();

        ostream.close();
    }
}

void COutputStream::_errorFileOpen() const
{
    std::cerr << "*** ERROR @ COutputStream @: ";

    std::cerr << "failed to open file " << _filename  << std::endl;
}

void COutputStream::_addLineToBuffer()
{
    if (!_filename.empty())
    {
        if (!_prefix.empty()) _line.insert(0, _prefix);

        _buffer.push_back(COutputLine(_line, _getLinePosition(),
                                      _leftSymbol, _rightSymbol,
                                      _fillSymbol, _maxWidth));

        _line.clear();
    }
}

size_t COutputStream::_getLinePosition() const
{
    auto linwidth = _line.size();

    if (_maxWidth <= linwidth) return 0;

    if (_alignment == fmt::center) return (_maxWidth - linwidth) / 2;

    if (_alignment == fmt::right) return (_maxWidth - linwidth);

    return 0;
}

void COutputStream::_appendToLine(const std::string& line)
{
    _line.append(line);
}


COutputStream& operator<<(COutputStream& output, const fmt& source)
{
    switch (source)
    {
        case fmt::end:

            output._addLineToBuffer();

            break;

        case fmt::left:

            output._alignment = fmt::left;

            break;

        case fmt::center:

            output._alignment = fmt::center;

            break;

        case fmt::right:

            output._alignment = fmt::right;

            break;

        case fmt::blank:

            if (!output._line.empty()) output._addLineToBuffer();

            output._leftSymbol = ' ';

            output._rightSymbol = ' ';

            output._fillSymbol = ' ';

            output._prefix.clear();

            output._addLineToBuffer();

            break;

        case fmt::title:

            output._leftSymbol = '!';

            output._rightSymbol = '!';

            output._fillSymbol = ' ';

            output._alignment = fmt::center;

            output._prefix.clear();

            break;

        case fmt::tsep:

            if (!output._line.empty()) output._addLineToBuffer();

            output._fillSymbol = '=';

            output._addLineToBuffer();

            output._fillSymbol = ' ';

            break;

        case fmt::cerror:

            output._leftSymbol = '*';

            output._rightSymbol = ' ';

            output._fillSymbol = ' ';

            output._alignment = fmt::left;

            output._prefix.assign("** CRITICAL ERROR *** ");

            break;

        case fmt::info:

            output._leftSymbol = '*';

            output._rightSymbol = ' ';

            output._fillSymbol = ' ';

            output._alignment = fmt::left;

            output._prefix.assign(" Info * ");

            break;

        case fmt::error:

            output._leftSymbol = '*';

            output._rightSymbol = ' ';

            output._fillSymbol = ' ';

            output._alignment = fmt::left;

            output._prefix.assign("** SYNTAX ERROR *** ");

            break;

        case fmt::header:

            output._leftSymbol = ' ';

            output._rightSymbol = ' ';

            output._fillSymbol = ' ';

            output._alignment = fmt::center;

            output._prefix.clear();

            break;

        default:

            break;
    }

    return output;
}

COutputStream& operator<<(COutputStream& output, const std::string& source)
{
    output._appendToLine(source);

    return output;
}

COutputStream& operator<<(COutputStream& output, const char* source)
{
    output._appendToLine(std::string(source));

    return output;
}
