//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ControlGroup.hpp"

#include <utility>

CControlGroup::CControlGroup()
{

}

CControlGroup::CControlGroup(const CControlGroup& source)

    : _header(source._header)

    , _commands(source._commands)
{

}

CControlGroup::CControlGroup(CControlGroup&& source) noexcept

    : _header(std::move(source._header))

    , _commands(std::move(source._commands))
{

}

CControlGroup::~CControlGroup()
{

}

CControlGroup& CControlGroup::operator=(const CControlGroup& source)
{
    if (this == &source) return *this;

    _header = source._header;

    _commands = source._commands;

    return *this;
}

CControlGroup& CControlGroup::operator=(CControlGroup&& source) noexcept
{
    if (this == &source) return *this;

    _header = std::move(source._header);

    _commands = std::move(source._commands);

    return *this;
}

bool CControlGroup::operator==(const CControlGroup& other) const
{
    if (_header != other._header) return false;

    if (_commands.size() != other._commands.size()) return false;

    for (size_t i = 0; i < _commands.size(); i++)
    {
        if (_commands[i] != other._commands[i]) return false;
    }

    return true;
}

bool CControlGroup::operator!=(const CControlGroup& other) const
{
    return !(*this == other);
}

void CControlGroup::setHeader(const CInputLine& header)
{
    _header = header;
}

void CControlGroup::addCommand(const CInputLine& command)
{
    _commands.push_back(command);
}

void CControlGroup::clear()
{
    _header.clear();

    _commands.clear();
}

bool CControlGroup::isEmpty() const
{
    return _commands.empty();
}

bool CControlGroup::isNameOfControlGroup(const std::string& name) const
{
    return _header.isControlKeyword(name);
}

bool CControlGroup::isNameOfControlGroup(const char* name) const
{
    return isNameOfControlGroup(std::string(name));
}

size_t CControlGroup::getNumberOfCommands() const
{
    return _commands.size();
}

size_t CControlGroup::getNumberOfCommands(const std::string& keyword) const
{
    size_t nkeys = 0;

    for (size_t i = 0; i < _commands.size(); i++)
    {
        if (_commands[i].isKeyword(0, keyword)) nkeys++;
    }

    return nkeys;
}

size_t CControlGroup::getNumberOfCommands(const char* keyword) const
{
    return getNumberOfCommands(std::string(keyword));
}

CInputLine CControlGroup::getCommand(const size_t index) const
{
    return _commands[index];
}

CInputLine CControlGroup::getCommand(const std::string& keyword) const
{
    for (size_t i = 0; i < _commands.size(); i++)
    {
        if (_commands[i].isKeyword(0, keyword)) return _commands[i];
    }

    return CInputLine();
}

CInputLine CControlGroup::getCommand(const char* keyword) const
{
    return getCommand(std::string(keyword));
}

std::ostream& operator<<(std::ostream&  output, const CControlGroup& source)
{
    output << std::endl;

    output << "[CControlGroup (Object):" << &source << "]" << std::endl;

    output << "_header: ";

    output << source._header <<  std::endl;

    output << "_commands: " << std::endl;

    for (size_t i = 0; i < source._commands.size(); i++)
    {
        output << "_commands[" << i << "]:" << std::endl;

        output << source._commands[i] << std::endl;
    }

    return output;
}
