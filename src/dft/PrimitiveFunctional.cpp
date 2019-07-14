//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "PrimitiveFunctional.hpp"

CPrimitiveFunctional::CPrimitiveFunctional()

    : _label(std::string())

    , _xcFuncType(xcfun::undefined)

    , _abFirstOrderFunction(nullptr)

    , _aFirstOrderFunction(nullptr)

    , _bFirstOrderFunction(nullptr)
{
}

CPrimitiveFunctional::CPrimitiveFunctional(const std::string&                     label,
                                           const xcfun                            xcFuncType,
                                           const std::function<def_vxc_func_typ>& abFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>& aFirstOrderFunction,
                                           const std::function<def_vxc_func_typ>& bFirstOrderFunction)

    : _label(label)

    , _xcFuncType(xcFuncType)

    , _abFirstOrderFunction(abFirstOrderFunction)

    , _aFirstOrderFunction(aFirstOrderFunction)

    , _bFirstOrderFunction(bFirstOrderFunction)
{
    
}

CPrimitiveFunctional::CPrimitiveFunctional(const CPrimitiveFunctional& source)

    : _label(source._label)

    , _xcFuncType(source._xcFuncType)

    , _abFirstOrderFunction(source._abFirstOrderFunction)

    , _aFirstOrderFunction(source._aFirstOrderFunction)

    , _bFirstOrderFunction(source._bFirstOrderFunction)
{
    
}

CPrimitiveFunctional::CPrimitiveFunctional(CPrimitiveFunctional&& source) noexcept

    : _label(std::move(source._label))

    , _xcFuncType(std::move(source._xcFuncType))

    , _abFirstOrderFunction(std::move(source._abFirstOrderFunction))

    , _aFirstOrderFunction(std::move(source._aFirstOrderFunction))

    , _bFirstOrderFunction(std::move(source._bFirstOrderFunction))
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
