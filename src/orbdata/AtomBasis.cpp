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

#include <algorithm>
#include <numeric>
#include <ranges>

#include "ChemicalElement.hpp"
#include "TensorLabels.hpp"

CAtomBasis::CAtomBasis()

    : _functions{}

    , _name{}

    , _ecp_label{}

    , _identifier{-1}
{
}

CAtomBasis::CAtomBasis(const std::vector<CBasisFunction> &functions, const std::string &name, const std::string &ecp_label, const int identifier)

    : _functions(functions)

    , _name(name)

    , _ecp_label(ecp_label)

    , _identifier(identifier)
{
}

CAtomBasis::CAtomBasis(const CAtomBasis &other)

    : _functions(other._functions)

    , _name(other._name)

    , _ecp_label(other._ecp_label)

    , _identifier(other._identifier)
{
}

CAtomBasis::CAtomBasis(CAtomBasis &&other) noexcept

    : _functions(std::move(other._functions))

    , _name(std::move(other._name))

    , _ecp_label(std::move(other._ecp_label))

    , _identifier(std::move(other._identifier))
{
}

auto
CAtomBasis::operator=(const CAtomBasis &other) -> CAtomBasis &
{
    _functions = other._functions;

    _name = other._name;

    _ecp_label = other._ecp_label;

    _identifier = other._identifier;

    return *this;
}

auto
CAtomBasis::operator=(CAtomBasis &&other) noexcept -> CAtomBasis &
{
    if (this != &other)
    {
        _functions = std::move(other._functions);

        _name = std::move(other._name);

        _ecp_label = std::move(other._ecp_label);

        _identifier = std::move(other._identifier);
    }

    return *this;
}

auto
CAtomBasis::operator==(const CAtomBasis &other) const -> bool
{
    if (_identifier != other._identifier)
    {
        return false;
    }
    else if (_name != other._name)
    {
        return false;
    }
    else if (_ecp_label != other._ecp_label)
    {
        return false;
    }
    else
    {
        return _functions == other._functions;
    }
}

auto
CAtomBasis::set_identifier(const int identifier) -> void
{
    _identifier = identifier;
}

auto
CAtomBasis::set_name(const std::string &name) -> void
{
    _name = name;
}

auto
CAtomBasis::set_ecp_label(const std::string &label) -> void
{
    _ecp_label = label;
}

auto
CAtomBasis::add(const CBasisFunction &function) -> void
{
    _functions.push_back(function);
}

auto
CAtomBasis::reduce_to_valence_basis() const -> CAtomBasis
{
    CAtomBasis vbasis;

    vbasis.set_identifier(_identifier);

    vbasis.set_name(_name + "(Valence)");

    vbasis.set_ecp_label(_ecp_label);

    const auto mang = chem_elem::max_angular_momentum(_identifier);

    std::ranges::for_each(_functions, [&](const auto &bf) {
        if (bf.get_angular_momentum() <= mang) vbasis.add(bf);
    });

    return vbasis;
}

auto
CAtomBasis::basis_functions() const -> std::vector<CBasisFunction>
{
    return _functions;
}

auto
CAtomBasis::basis_functions(const int angular_momentum) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    bfs.reserve(_functions.size());

    std::ranges::copy_if(_functions, std::back_inserter(bfs), [=](const auto &bf) { return bf.get_angular_momentum() == angular_momentum; });

    return bfs;
}

auto
CAtomBasis::basis_functions(const int angular_momentum, const size_t npgtos) const -> std::vector<CBasisFunction>
{
    std::vector<CBasisFunction> bfs;

    bfs.reserve(_functions.size());

    std::ranges::copy_if(_functions, std::back_inserter(bfs), [=](const auto &bf) {
        return (bf.get_angular_momentum() == angular_momentum) && (bf.number_of_primitive_functions() == npgtos);
    });

    return bfs;
}

auto
CAtomBasis::get_identifier() const -> int
{
    return _identifier;
}

auto
CAtomBasis::get_name() const -> std::string
{
    return _name;
}

auto
CAtomBasis::get_ecp_label() const -> std::string
{
    return _ecp_label;
}

auto
CAtomBasis::need_ecp() const -> bool
{
    return !(_ecp_label.empty());
}

auto
CAtomBasis::max_angular_momentum() const -> int
{
    auto pos = std::ranges::max_element(_functions,
                                        [&](const auto &lbf, const auto &rbf) { return lbf.get_angular_momentum() < rbf.get_angular_momentum(); });

    return (pos == _functions.end()) ? -1 : pos->get_angular_momentum();
}

auto
CAtomBasis::number_of_basis_functions(const int angular_momentum) const -> size_t
{
    return static_cast<size_t>(std::ranges::count_if(_functions, [=](const auto &bf) { return bf.get_angular_momentum() == angular_momentum; }));
}

auto
CAtomBasis::number_of_basis_functions(const int angular_momentum, const size_t npgtos) const -> size_t
{
    return static_cast<size_t>(std::ranges::count_if(_functions, [=](const auto &bf) {
        return (bf.get_angular_momentum() == angular_momentum) && (bf.number_of_primitive_functions() == npgtos);
    }));
}

auto
CAtomBasis::number_of_primitive_functions(const int angular_momentum) const -> size_t
{
    return std::accumulate(_functions.begin(), _functions.end(), size_t{0}, [=](const size_t &sum, const auto &bf) {
        return (bf.get_angular_momentum() == angular_momentum) ? sum + bf.number_of_primitive_functions() : sum;
    });
}

auto
CAtomBasis::contraction_depths(const int angular_momentum) const -> std::set<size_t>
{
    std::set<size_t> depths;

    std::ranges::for_each(_functions, [&](const auto &bf) {
        if (bf.get_angular_momentum() == angular_momentum) depths.insert(bf.number_of_primitive_functions());
    });

    return depths;
}

auto
CAtomBasis::contraction_string() const -> std::string
{
    auto str = std::string("(");

    auto mang = max_angular_momentum();

    std::ranges::for_each(std::views::iota(0, mang + 1), [&](const int i) {
        if (const auto ncgtos = number_of_basis_functions(i); ncgtos > 0)
        {
            str.append(std::to_string(ncgtos));
            str.append(1, tensor::label(i));
            if (i != mang) str.append(",");
        }
    });

    str.append(")");

    return str;
}

auto
CAtomBasis::primitives_string() const -> std::string
{
    auto str = std::string("(");

    const auto mang = max_angular_momentum();

    std::ranges::for_each(std::views::iota(0, mang + 1), [&](const int i) {
        if (const auto npgtos = number_of_primitive_functions(i); npgtos > 0)
        {
            str.append(std::to_string(npgtos));
            str.append(1, tensor::label(i));
            if (i != mang) str.append(",");
        }
    });

    str.append(")");

    return str;
}
