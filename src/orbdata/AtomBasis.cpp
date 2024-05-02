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

#include "AtomBasis.hpp"

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"

CAtomBasis::CAtomBasis(const std::vector<CBasisFunction>& functions, const int64_t identifier, const std::string& name, const std::string& ecp_label)

    : _functions(functions)

    , _identifier(identifier)

    , _name(name)

    , _ecp_label(ecp_label)
{
}

auto
CAtomBasis::setIdentifier(const int64_t identifier) -> void
{
    _identifier = identifier;
}

auto
CAtomBasis::setName(const std::string& name) -> void
{
    _name = name;
}

auto
CAtomBasis::setEffectiveCorePotentialLabel(const std::string& label) -> void
{
    _ecp_label = label;
}

auto
CAtomBasis::add(const CBasisFunction& function) -> void
{
    _functions.push_back(function);
}

auto
CAtomBasis::reduceToValenceBasis() const -> CAtomBasis
{
    CAtomBasis vbasis;

    vbasis.setIdentifier(_identifier);

    vbasis.setName(_name + "(Valence)");

    vbasis.setEffectiveCorePotentialLabel(_ecp_label);

    if (CChemicalElement elem; elem.setAtomType(_identifier))
    {
        if (const auto mang = elem.getMaxAngularMomentum(); mang >= 0)
        {
            if (const auto ncgtos = _functions.size(); ncgtos > 0)
            {
                for (size_t i = 0; i < ncgtos; i++)
                {
                    if (_functions[i].getAngularMomentum() <= mang)
                    {
                        vbasis.add(_functions[i]);
                    }
                }
            }
        }
    }

    return vbasis;
}

auto
CAtomBasis::getBasisFunctions() const -> std::vector<CBasisFunction>
{
    return _functions;
}

auto
CAtomBasis::getIdentifier() const -> int64_t
{
    return _identifier;
}

auto
CAtomBasis::getName() const -> std::string
{
    return _name;
}

auto
CAtomBasis::getEffectiveCorePotentialLabel() const -> std::string
{
    return _ecp_label;
}

auto
CAtomBasis::needEffectiveCorePotential() const -> bool
{
    return _ecp_label.size() > 0;
}

auto
CAtomBasis::getMaxAngularMomentum() const -> int64_t
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        auto mang = _functions[0].getAngularMomentum();

        for (size_t i = 1; i < ncgtos; i++)
        {
            if (const auto angmom = _functions[i].getAngularMomentum(); angmom > mang)
            {
                mang = angmom;
            }
        }

        return mang;
    }
    else
    {
        return -1;
    }
}

auto
CAtomBasis::getNumberOfBasisFunctions(const int64_t angmom) const -> int64_t
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        auto ngtos = 0;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if (_functions[i].getAngularMomentum() == angmom) ngtos++;
        }

        return ngtos;
    }
    else
    {
        return 0;
    }
}

auto
CAtomBasis::getNumberOfBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> int64_t
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        auto ngtos = 0;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if ((_functions[i].getAngularMomentum() == angmom) && (_functions[i].getNumberOfPrimitiveFunctions() == npgtos)) ngtos++;
        }

        return ngtos;
    }
    else
    {
        return 0;
    }
}

auto
CAtomBasis::getNumberOfPrimitiveFunctions(const int64_t angmom) const -> int64_t
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        auto npgtos = 0;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if (_functions[i].getAngularMomentum() == angmom)
            {
                npgtos += _functions[i].getNumberOfPrimitiveFunctions();
            }
        }

        return npgtos;
    }
    else
    {
        return 0;
    }
}

auto
CAtomBasis::getContractionDepths(const int64_t angmom) const -> std::set<int64_t>
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        std::set<int64_t> depths;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if (_functions[i].getAngularMomentum() == angmom)
            {
                depths.insert(_functions[i].getNumberOfPrimitiveFunctions());
            }
        }

        return depths;
    }
    else
    {
        return std::set<int64_t>();
    }
}

auto
CAtomBasis::getContractionString() const -> std::string
{
    if (const auto mang = getMaxAngularMomentum(); mang >= 0)
    {
        std::string str("(");

        for (int64_t i = 0; i <= mang; i++)
        {
            if (const auto ncgtos = getNumberOfBasisFunctions(i); ncgtos > 0)
            {
                str.append(std::to_string(ncgtos));

                str.append(fstr::to_AngularMomentum(i));

                if (i != mang) str.append(",");
            }
        }

        str.append(")");

        return str;
    }
    else
    {
        return std::string();
    }
}

auto
CAtomBasis::getPrimitivesString() const -> std::string
{
    if (const auto mang = getMaxAngularMomentum(); mang >= 0)
    {
        std::string str("(");

        for (int64_t i = 0; i <= mang; i++)
        {
            if (const auto npgtos = getNumberOfPrimitiveFunctions(i); npgtos > 0)
            {
                str.append(std::to_string(npgtos));

                str.append(fstr::to_AngularMomentum(i));

                if (i != mang) str.append(",");
            }
        }

        str.append(")");

        return str;
    }
    else
    {
        return std::string();
    }
}

auto
CAtomBasis::getBasisFunctions(const int64_t angmom) const -> std::vector<CBasisFunction>
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if (_functions[i].getAngularMomentum() == angmom)
            {
                gtos.push_back(_functions[i]);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}

auto
CAtomBasis::getBasisFunctions(const int64_t angmom, const int64_t npgtos) const -> std::vector<CBasisFunction>
{
    if (const auto ncgtos = _functions.size(); ncgtos > 0)
    {
        std::vector<CBasisFunction> gtos;

        for (size_t i = 0; i < ncgtos; i++)
        {
            if ((_functions[i].getAngularMomentum() == angmom) && (_functions[i].getNumberOfPrimitiveFunctions() == npgtos))
            {
                gtos.push_back(_functions[i]);
            }
        }

        return gtos;
    }
    else
    {
        return std::vector<CBasisFunction>();
    }
}
