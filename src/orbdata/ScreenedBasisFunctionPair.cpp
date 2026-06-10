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

#include "ScreenedBasisFunctionPair.hpp"

CScreenedBasisFunctionPair::CScreenedBasisFunctionPair()

    : _bra_function()

    , _ket_function()

    , _bra_x{}

    , _bra_y{}

    , _bra_z{}

    , _ket_x{}

    , _ket_y{}

    , _ket_z{}

    , _orb_indices{}
{
}

CScreenedBasisFunctionPair::CScreenedBasisFunctionPair(const CBasisFunction                         &bra_function,
                                                       const CBasisFunction                         &ket_function,
                                                       const std::vector<double>                    &bra_x,
                                                       const std::vector<double>                    &bra_y,
                                                       const std::vector<double>                    &bra_z,
                                                       const std::vector<double>                    &ket_x,
                                                       const std::vector<double>                    &ket_y,
                                                       const std::vector<double>                    &ket_z,
                                                       const std::vector<std::pair<size_t, size_t>> &orb_indices)

    : _bra_function(bra_function)

    , _ket_function(ket_function)

    , _bra_x(bra_x)

    , _bra_y(bra_y)

    , _bra_z(bra_z)

    , _ket_x(ket_x)

    , _ket_y(ket_y)

    , _ket_z(ket_z)

    , _orb_indices(orb_indices)
{
}

CScreenedBasisFunctionPair::CScreenedBasisFunctionPair(const CScreenedBasisFunctionPair &other)

    : _bra_function(other._bra_function)

    , _ket_function(other._ket_function)

    , _bra_x(other._bra_x)

    , _bra_y(other._bra_y)

    , _bra_z(other._bra_z)

    , _ket_x(other._ket_x)

    , _ket_y(other._ket_y)

    , _ket_z(other._ket_z)

    , _orb_indices(other._orb_indices)
{
}

CScreenedBasisFunctionPair::CScreenedBasisFunctionPair(CScreenedBasisFunctionPair &&other) noexcept

    : _bra_function(std::move(other._bra_function))

    , _ket_function(std::move(other._ket_function))

    , _bra_x(std::move(other._bra_x))

    , _bra_y(std::move(other._bra_y))

    , _bra_z(std::move(other._bra_z))

    , _ket_x(std::move(other._ket_x))

    , _ket_y(std::move(other._ket_y))

    , _ket_z(std::move(other._ket_z))

    , _orb_indices(std::move(other._orb_indices))
{
}

auto
CScreenedBasisFunctionPair::operator=(const CScreenedBasisFunctionPair &other) -> CScreenedBasisFunctionPair &
{
    _bra_function = other._bra_function;

    _ket_function = other._ket_function;

    _bra_x = other._bra_x;

    _bra_y = other._bra_y;

    _bra_z = other._bra_z;

    _ket_x = other._ket_x;

    _ket_y = other._ket_y;

    _ket_z = other._ket_z;

    _orb_indices = other._orb_indices;

    return *this;
}

auto
CScreenedBasisFunctionPair::operator=(CScreenedBasisFunctionPair &&other) noexcept -> CScreenedBasisFunctionPair &
{
    if (this != &other)
    {
        _bra_function = std::move(other._bra_function);

        _ket_function = std::move(other._ket_function);

        _bra_x = std::move(other._bra_x);

        _bra_y = std::move(other._bra_y);

        _bra_z = std::move(other._bra_z);

        _ket_x = std::move(other._ket_x);

        _ket_y = std::move(other._ket_y);

        _ket_z = std::move(other._ket_z);

        _orb_indices = std::move(other._orb_indices);
    }

    return *this;
}

auto
CScreenedBasisFunctionPair::operator==(const CScreenedBasisFunctionPair &other) const -> bool
{
    return (_bra_function == other._bra_function) && (_ket_function == other._ket_function) && (_bra_x == other._bra_x) && (_bra_y == other._bra_y) &&
           (_bra_z == other._bra_z) && (_ket_x == other._ket_x) && (_ket_y == other._ket_y) && (_ket_z == other._ket_z) &&
           (_orb_indices == other._orb_indices);
}

auto
CScreenedBasisFunctionPair::operator!=(const CScreenedBasisFunctionPair &other) const -> bool
{
    return !(*this == other);
}

auto
CScreenedBasisFunctionPair::_append(const TPoint<double> &bra_center,
                                    const TPoint<double> &ket_center,
                                    const size_t          bra_index,
                                    const size_t          ket_index) -> void
{
    const auto r_bra = bra_center.coordinates();

    const auto r_ket = ket_center.coordinates();

    _bra_x.push_back(r_bra[0]);

    _bra_y.push_back(r_bra[1]);

    _bra_z.push_back(r_bra[2]);

    _ket_x.push_back(r_ket[0]);

    _ket_y.push_back(r_ket[1]);

    _ket_z.push_back(r_ket[2]);

    _orb_indices.push_back({bra_index, ket_index});
}

auto
CScreenedBasisFunctionPair::bra_function() const -> CBasisFunction
{
    return _bra_function;
}

auto
CScreenedBasisFunctionPair::ket_function() const -> CBasisFunction
{
    return _ket_function;
}

auto
CScreenedBasisFunctionPair::bra_x() const -> std::vector<double>
{
    return _bra_x;
}

auto
CScreenedBasisFunctionPair::bra_y() const -> std::vector<double>
{
    return _bra_y;
}

auto
CScreenedBasisFunctionPair::bra_z() const -> std::vector<double>
{
    return _bra_z;
}

auto
CScreenedBasisFunctionPair::ket_x() const -> std::vector<double>
{
    return _ket_x;
}

auto
CScreenedBasisFunctionPair::ket_y() const -> std::vector<double>
{
    return _ket_y;
}

auto
CScreenedBasisFunctionPair::ket_z() const -> std::vector<double>
{
    return _ket_z;
}

auto
CScreenedBasisFunctionPair::orbital_indices() const -> std::vector<std::pair<size_t, size_t>>
{
    return _orb_indices;
}

auto
CScreenedBasisFunctionPair::number_of_pairs() const -> size_t
{
    return _orb_indices.size();
}
