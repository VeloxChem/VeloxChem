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
