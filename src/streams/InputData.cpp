//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "InputData.hpp"

CInputData::CInputData()
{

}

CInputData::~CInputData()
{

}

bool
CInputData::operator==(const CInputData& other) const
{
    if (_controlGroups.size() != other._controlGroups.size()) return false;

    for (size_t i = 0; i < _controlGroups.size(); i++)
    {
        if (_controlGroups[i] != other._controlGroups[i]) return false;
    }

    return true;
}

bool
CInputData::operator!=(const CInputData& other) const
{
    return !(*this == other);
}

void
CInputData::addControlGroup(const CControlGroup& controlGroup)
{
    _controlGroups.push_back(controlGroup);
}

int32_t
CInputData::getNumberOfControlGroups() const
{
    return static_cast<int32_t>(_controlGroups.size());
}

int32_t
CInputData::getNumberOfControlGroups(const std::string& nameOfControlGroup) const
{
    int32_t ngroups = 0;

    for (size_t i = 0; i < _controlGroups.size(); i++)
    {
        if (_controlGroups[i].isNameOfControlGroup(nameOfControlGroup))
        {
            ngroups++;
        }
    }

    return ngroups;
}

int32_t
CInputData::getNumberOfControlGroups(const char* nameOfControlGroup) const
{
    return getNumberOfControlGroups(std::string(nameOfControlGroup));
}

CControlGroup
CInputData::getControlGroup(const size_t indexOfControlGroup) const
{
    return _controlGroups[indexOfControlGroup];
}

CControlGroup
CInputData::getControlGroup(const size_t       indexOfControlGroup,
                            const std::string& nameOfControlGroup) const
{
    auto ngrps = static_cast<size_t>(getNumberOfControlGroups(nameOfControlGroup));

    if (indexOfControlGroup < ngrps)
    {
        size_t curgroup = 0;

        for (size_t i = 0; i < _controlGroups.size(); i++)
        {
            if (_controlGroups[i].isNameOfControlGroup(nameOfControlGroup))
            {
                if (curgroup == indexOfControlGroup) return _controlGroups[i];

                curgroup++;
            }
        }
    }

    return CControlGroup();
}

CControlGroup
CInputData::getControlGroup(const size_t indexOfControlGroup,
                            const char*  nameOfControlGroup) const
{
    return getControlGroup(indexOfControlGroup,
                           std::string(nameOfControlGroup));
}

std::ostream&
operator<<(      std::ostream& output,
           const CInputData&   source)
{
    output << std::endl;

    output << "[CInputData (Object):" << &source << "]" << std::endl;

    output << " _controlGroups: " << std::endl;

    for (size_t i = 0; i < source._controlGroups.size(); i++)
    {
        output << "_controlGroups[" << i << "]: " << std::endl;

        output << source._controlGroups[i] << std::endl;
    }

    return output;
}
