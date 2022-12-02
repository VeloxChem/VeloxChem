//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "PrimitiveFunctional.hpp"

#include <xc.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>

#include "ErrorHandler.hpp"
#include "MemAlloc.hpp"

CPrimitiveFunctional::CPrimitiveFunctional()

    : _label(std::string())

    , _xcFuncType(xcfun::undefined)

    , _abFirstOrderFunction(nullptr)

    , _aFirstOrderFunction(nullptr)

    , _bFirstOrderFunction(nullptr)

    , _abSecondOrderFunction(nullptr)

    , _aSecondOrderFunction(nullptr)

    , _bSecondOrderFunction(nullptr)

    , _abThirdOrderFunction(nullptr)

    , _aThirdOrderFunction(nullptr)

    , _bThirdOrderFunction(nullptr)
{
}

CPrimitiveFunctional::CPrimitiveFunctional(const std::string&                      label,
                                           const xcfun                             xcFuncType,
                                           const std::function<def_vxc_func_typ>&  abFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>&  aFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>&  bFirstOrderFunction,
                                           const std::function<def_vxc2_func_typ>& abSecondOrderFunction,
                                           const std::function<def_vxc2_func_typ>& aSecondOrderFunction,
                                           const std::function<def_vxc2_func_typ>& bSecondOrderFunction,
                                           const std::function<def_vxc3_func_typ>& abThirdOrderFunction,
                                           const std::function<def_vxc3_func_typ>& aThirdOrderFunction,
                                           const std::function<def_vxc3_func_typ>& bThirdOrderFunction)

    : _label(label)

    , _xcFuncType(xcFuncType)

    , _abFirstOrderFunction(abFirstOrderFunction)

    , _aFirstOrderFunction(aFirstOrderFunction)

    , _bFirstOrderFunction(bFirstOrderFunction)

    , _abSecondOrderFunction(abSecondOrderFunction)

    , _aSecondOrderFunction(aSecondOrderFunction)

    , _bSecondOrderFunction(bSecondOrderFunction)

    , _abThirdOrderFunction(abThirdOrderFunction)

    , _aThirdOrderFunction(aThirdOrderFunction)

    , _bThirdOrderFunction(bThirdOrderFunction)
{
}

CPrimitiveFunctional::CPrimitiveFunctional(const CPrimitiveFunctional& source)

    : _label(source._label)

    , _xcFuncType(source._xcFuncType)

    , _abFirstOrderFunction(source._abFirstOrderFunction)

    , _aFirstOrderFunction(source._aFirstOrderFunction)

    , _bFirstOrderFunction(source._bFirstOrderFunction)

    , _abSecondOrderFunction(source._abSecondOrderFunction)

    , _aSecondOrderFunction(source._aSecondOrderFunction)

    , _bSecondOrderFunction(source._bSecondOrderFunction)

    , _abThirdOrderFunction(source._abThirdOrderFunction)

    , _aThirdOrderFunction(source._aThirdOrderFunction)

    , _bThirdOrderFunction(source._bThirdOrderFunction)
{
}

CPrimitiveFunctional::CPrimitiveFunctional(CPrimitiveFunctional&& source) noexcept

    : _label(std::move(source._label))

    , _xcFuncType(std::move(source._xcFuncType))

    , _abFirstOrderFunction(std::move(source._abFirstOrderFunction))

    , _aFirstOrderFunction(std::move(source._aFirstOrderFunction))

    , _bFirstOrderFunction(std::move(source._bFirstOrderFunction))

    , _abSecondOrderFunction(std::move(source._abSecondOrderFunction))

    , _aSecondOrderFunction(std::move(source._aSecondOrderFunction))

    , _bSecondOrderFunction(std::move(source._bSecondOrderFunction))

    , _abThirdOrderFunction(std::move(source._abThirdOrderFunction))

    , _aThirdOrderFunction(std::move(source._aThirdOrderFunction))

    , _bThirdOrderFunction(std::move(source._bThirdOrderFunction))
{
}

CPrimitiveFunctional::~CPrimitiveFunctional()
{
}

CPrimitiveFunctional&
CPrimitiveFunctional::operator=(const CPrimitiveFunctional& source)
{
    if (this == &source) return *this;

    _label = source._label;

    _xcFuncType = source._xcFuncType;

    _abFirstOrderFunction = source._abFirstOrderFunction;

    _aFirstOrderFunction = source._aFirstOrderFunction;

    _bFirstOrderFunction = source._bFirstOrderFunction;

    _abSecondOrderFunction = source._abSecondOrderFunction;

    _aSecondOrderFunction = source._aSecondOrderFunction;

    _bSecondOrderFunction = source._bSecondOrderFunction;

    _abThirdOrderFunction = source._abThirdOrderFunction;

    _aThirdOrderFunction = source._aThirdOrderFunction;

    _bThirdOrderFunction = source._bThirdOrderFunction;

    return *this;
}

CPrimitiveFunctional&
CPrimitiveFunctional::operator=(CPrimitiveFunctional&& source) noexcept
{
    if (this == &source) return *this;

    _label = std::move(source._label);

    _xcFuncType = std::move(source._xcFuncType);

    _abFirstOrderFunction = std::move(source._abFirstOrderFunction);

    _aFirstOrderFunction = std::move(source._aFirstOrderFunction);

    _bFirstOrderFunction = std::move(source._bFirstOrderFunction);

    _abSecondOrderFunction = std::move(source._abSecondOrderFunction);

    _aSecondOrderFunction = std::move(source._aSecondOrderFunction);

    _bSecondOrderFunction = std::move(source._bSecondOrderFunction);

    _abThirdOrderFunction = std::move(source._abThirdOrderFunction);

    _aThirdOrderFunction = std::move(source._aThirdOrderFunction);

    _bThirdOrderFunction = std::move(source._bThirdOrderFunction);

    return *this;
}

bool
CPrimitiveFunctional::operator==(const CPrimitiveFunctional& other) const
{
    if (this == &other) return true;

    if (_label != other._label) return false;

    if (_xcFuncType != other._xcFuncType) return false;

    return true;
}

bool
CPrimitiveFunctional::operator!=(const CPrimitiveFunctional& other) const
{
    return !((*this) == other);
}

void
CPrimitiveFunctional::compute(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid) const
{
    if (densityGrid.getDensityGridType() == dengrid::ab) _abFirstOrderFunction(xcGradientGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::lima) _aFirstOrderFunction(xcGradientGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::limb) _bFirstOrderFunction(xcGradientGrid, factor, densityGrid);
}

void
CPrimitiveFunctional::compute(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid) const
{
    if (densityGrid.getDensityGridType() == dengrid::ab) _abSecondOrderFunction(xcHessianGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::lima) _aSecondOrderFunction(xcHessianGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::limb) _bSecondOrderFunction(xcHessianGrid, factor, densityGrid);
}

void
CPrimitiveFunctional::compute(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid) const
{
    if (densityGrid.getDensityGridType() == dengrid::ab) _abThirdOrderFunction(xcCubicHessianGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::lima) _aThirdOrderFunction(xcCubicHessianGrid, factor, densityGrid);

    if (densityGrid.getDensityGridType() == dengrid::limb) _bThirdOrderFunction(xcCubicHessianGrid, factor, densityGrid);
}

std::string
CPrimitiveFunctional::getLabel() const
{
    return _label;
}

xcfun
CPrimitiveFunctional::getFunctionalType() const
{
    return _xcFuncType;
}

std::ostream&
operator<<(std::ostream& output, const CPrimitiveFunctional& source)
{
    output << std::endl;

    output << "[CPrimitiveFunctional (Object):" << &source << "]" << std::endl;

    output << "_label: " << source._label << std::endl;

    output << "_xcFuncType:" << to_string(source._xcFuncType);

    output << "_abFirstOrderFunction: " << &(source._abFirstOrderFunction) << std::endl;

    output << "_aFirstOrderFunction: " << &(source._aFirstOrderFunction) << std::endl;

    output << "_bFirstOrderFunction: " << &(source._bFirstOrderFunction) << std::endl;

    return output;
}
