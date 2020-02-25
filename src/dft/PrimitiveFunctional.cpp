//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "PrimitiveFunctional.hpp"

CPrimitiveFunctional::CPrimitiveFunctional()

    : _label(std::string())

    , _xcFuncType(xcfun::undefined)

    , _abFirstOrderFunction(nullptr)

    , _aFirstOrderFunction(nullptr)

    , _bFirstOrderFunction(nullptr)

    , _abSecondOrderFunction(nullptr)

    , _aSecondOrderFunction(nullptr)

    , _bSecondOrderFunction(nullptr)
{
}

CPrimitiveFunctional::CPrimitiveFunctional(const std::string&                      label,
                                           const xcfun                             xcFuncType,
                                           const std::function<def_vxc_func_typ>&  abFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>&  aFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>&  bFirstOrderFunction,
                                           const std::function<def_vxc2_func_typ>& abSecondOrderFunction,
                                           const std::function<def_vxc2_func_typ>& aSecondOrderFunction,
                                           const std::function<def_vxc2_func_typ>& bSecondOrderFunction)

    : _label(label)

    , _xcFuncType(xcFuncType)

    , _abFirstOrderFunction(abFirstOrderFunction)

    , _aFirstOrderFunction(aFirstOrderFunction)

    , _bFirstOrderFunction(bFirstOrderFunction)

    , _abSecondOrderFunction(abSecondOrderFunction)

    , _aSecondOrderFunction(aSecondOrderFunction)

    , _bSecondOrderFunction(bSecondOrderFunction)
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
CPrimitiveFunctional::compute(      CXCGradientGrid& xcGradientGrid,
                              const double           factor,
                              const CDensityGrid&    densityGrid) const
{
    if (densityGrid.getDensityGridType() == dengrid::ab) _abFirstOrderFunction(xcGradientGrid, factor, densityGrid);
    
    if (densityGrid.getDensityGridType() == dengrid::lima) _aFirstOrderFunction(xcGradientGrid, factor, densityGrid);
    
    if (densityGrid.getDensityGridType() == dengrid::limb) _bFirstOrderFunction(xcGradientGrid, factor, densityGrid);
}

void
CPrimitiveFunctional::compute(      CXCHessianGrid& xcHessianGrid,
                              const double          factor,
                              const CDensityGrid&   densityGrid) const
{
    if (densityGrid.getDensityGridType() == dengrid::ab) _abSecondOrderFunction(xcHessianGrid, factor, densityGrid);
    
    if (densityGrid.getDensityGridType() == dengrid::lima) _aSecondOrderFunction(xcHessianGrid, factor, densityGrid);
    
    if (densityGrid.getDensityGridType() == dengrid::limb) _bSecondOrderFunction(xcHessianGrid, factor, densityGrid);
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
