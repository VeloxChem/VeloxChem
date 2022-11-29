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
#include <iostream>
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

Functional::Functional(const std::vector<std::string>& labels, const std::vector<double>& coeffs)
{
    errors::assertMsgCritical(labels.size() == coeffs.size(),
                              std::string(__func__) +
                                  ": exchange-correlation functional labels and mixing coefficients must have the same size! labels.size() = " +
                                  std::to_string(labels.size()) + ", coeffs.size() = " + std::to_string(coeffs.size()));

    auto is_in_family = [](const auto& options, auto family) {
        return std::any_of(options.begin(), options.end(), [family](auto a) { return a == family; });
    };

    // sum of coefficients of exchange components
    auto total_x = 0.0;

    // sum of coefficients of correlation components
    auto total_c = 0.0;

    std::vector<Component> x_funcs, c_funcs;

    // TODO write function to compute this number based on functional family and order of derivatives available
    auto n_xc_outputs = 0;

    for (auto i = 0; i < labels.size(); ++i)
    {
        auto label = labels[i];
        auto coeff = coeffs[i];

        auto funcID = xc_functional_get_number(label.c_str());

        xc_func_type func;
        if (xc_func_init(&func, funcID, XC_POLARIZED))
        {
            errors::msgCritical(std::string("Could not find required LibXC functional ") + label);
        }

        auto kind = func.info->kind;

        if (kind == XC_EXCHANGE)
        {
            total_x += coeff;
            x_funcs.push_back({coeff, func});
        }
        else if (kind == XC_CORRELATION)
        {
            total_c += coeff;
            c_funcs.push_back({coeff, func});
        }
        else
        {
            // clean up
            xc_func_end(&func);

            // error out if the kind is XC_EXCHANGE_CORRELATION ("canned" x-c mix in LibXC) or XC_KINETIC
            errors::msgCritical("We cannot handle functionals of kind XC_EXCHANGE_CORRELATION or XC_KINETIC");
        }

        auto flags = func.info->flags;

        // which derivative orders do we have for this x-c mixture?
        _hasExc = _hasExc && (flags & XC_FLAGS_HAVE_EXC);
        _hasVxc = _hasVxc && (flags & XC_FLAGS_HAVE_VXC);
        _hasFxc = _hasFxc && (flags & XC_FLAGS_HAVE_FXC);
        _hasKxc = _hasKxc && (flags & XC_FLAGS_HAVE_KXC);
        _hasLxc = _hasLxc && (flags & XC_FLAGS_HAVE_LXC);

        // which family does this x-c mixture belong to?
        auto family = func.info->family;

        // LDA?
        if (is_in_family(std::array{XC_FAMILY_LDA, XC_FAMILY_HYB_LDA}, family))
        {
            _isLDA       = true;
            n_xc_outputs = 15;
        }

        // GGA?
        if (is_in_family(std::array{XC_FAMILY_GGA, XC_FAMILY_HYB_GGA}, family))
        {
            _isGGA       = true;
            n_xc_outputs = 126;
        }

        // metaGGA?
        if (is_in_family(std::array{XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA}, family))
        {
            _isMetaGGA   = true;
            n_xc_outputs = 767;
        }

        // hybrid?
        // FIXME possibly can be removed: I believe these checks are only
        // relevant if we allow usage of any of the "prebaked" mixes in LibXC
        // hybrid-ness should be established from passed coefficients for global/range-separation
        if (is_in_family(std::array{XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA}, family))
        {
            // range separation using the error function or the Yukawa function
            if ((flags & XC_FLAGS_HYB_CAM) || (flags & XC_FLAGS_HYB_CAMY))
            {
                _isRangeSeparatedHybrid = true;
            }

            // global hybrid
            if (!_isRangeSeparatedHybrid)
            {
                _isGlobalHybrid = true;
            }
        }
    }

    _numberOfExchangeFunctionals    = x_funcs.size();
    _numberOfCorrelationFunctionals = c_funcs.size();
    // join the two lists together: first all exchange, then all correlation functionals.
    _components.reserve(_numberOfExchangeFunctionals + _numberOfCorrelationFunctionals);
    _components.insert(_components.end(), x_funcs.cbegin(), x_funcs.cend());
    _components.insert(_components.end(), c_funcs.cbegin(), c_funcs.cend());

    // warn if the sum of X coefficients and C coefficients is greater than 1
    errors::assertMsgWarning(total_x <= 1.0, "Sum of coefficients for exchange functionals is greater than 1. We hope you know what your are doing!");
    errors::assertMsgWarning(total_c <= 1.0,
                             "Sum of coefficients for correlation functionals is greater than 1. We hope you know what your are doing!");

    _xcName   = detail::getXCName(x_funcs, c_funcs);
    _citation = detail::getCitation(x_funcs, c_funcs);

    // allocate _stagingBuffer
    _stagingBuffer = mem::malloc<double>(n_xc_outputs * _ldStaging);
}

Functional::~Functional()
{
    // clean up allocated LibXC objects
    std::for_each(std::begin(_components), std::end(_components), [](Component& c) { xc_func_end(&std::get<1>(c)); });
    _components.clear();

    // clean up staging buffer
    mem::free(_stagingBuffer);
}

auto
Functional::repr() const -> std::string
{
    std::ostringstream os;

    os << std::endl;

    os << "[Functional (Object):" << this << "]" << std::endl;

    os << "Functional name: " << getXCName() << std::endl;

    os << "Functional family: ";
    if (_isMetaGGA)
    {
        os << "metaGGA." << std::endl;
    }
    else if (_isGGA)
    {
        os << "GGA." << std::endl;
    }
    else
    {
        os << "LDA." << std::endl;
    }
    os << std::endl;

    os << "Available derivatives up to and including ";
    if (_hasLxc)
    {
        os << "4th order." << std::endl;
    }
    else if (_hasKxc)
    {
        os << "3rd order." << std::endl;
    }
    else if (_hasFxc)
    {
        os << "2nd order." << std::endl;
    }
    else if (_hasVxc)
    {
        os << "1st order." << std::endl;
    }
    else
    {
        os << "0th order." << std::endl;
    }

    os << std::endl;

    os << getCitation();

    return os.str();
}

auto
Functional::compute_exc(int32_t np, const double* rho, double* exc) const -> void
{
    errors::assertMsgCritical(_hasExc, std::string(__func__) + ": exchange-correlation functional does not provide an evaluator for Exc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];

    for (const auto& [coeff, func] : _components)
    {
        xc_lda_exc(&func, np, rho, stage_exc);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
    }
}

auto
Functional::compute_vxc(int32_t np, const double* rho, double* vrho) const -> void
{
    errors::assertMsgCritical(_hasVxc, std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_vrho = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_lda_vxc(&func, np, rho, stage_vrho);

        const auto c = coeff;

#pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
        }
    }

    if (alloc)
    {
        mem::free(stage_vrho);
    }
}

auto
Functional::compute_exc_vxc(int32_t np, const double* rho, double* exc, double* vrho) const -> void
{
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc  = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];
    auto stage_vrho = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_lda_exc_vxc(&func, np, rho, stage_exc, stage_vrho);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];

            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
        mem::free(stage_vrho);
    }
}

auto
Functional::compute_exc(int32_t np, const double* rho, const double* sigma, double* exc) const -> void
{
    errors::assertMsgCritical(_hasExc, std::string(__func__) + ": exchange-correlation functional does not provide an evaluator for Exc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];

    for (const auto& [coeff, func] : _components)
    {
        xc_gga_exc(&func, np, rho, sigma, stage_exc);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
    }
}

auto
Functional::compute_vxc(int32_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_hasVxc, std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_gga_vxc(&func, np, rho, sigma, stage_vrho, stage_vsigma);

        const auto c = coeff;

#pragma omp simd aligned(vrho, stage_vrho, vsigma, stage_vsigma : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];

            vsigma[3 * g + 0] += c * stage_vsigma[3 * g + 0];
            vsigma[3 * g + 1] += c * stage_vsigma[3 * g + 1];
            vsigma[3 * g + 2] += c * stage_vsigma[3 * g + 2];
        }
    }

    if (alloc)
    {
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
    }
}

auto
Functional::compute_exc_vxc(int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_gga_exc_vxc(&func, np, rho, sigma, stage_exc, stage_vrho, stage_vsigma);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho, vsigma, stage_vsigma : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];

            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];

            vsigma[3 * g + 0] += c * stage_vsigma[3 * g + 0];
            vsigma[3 * g + 1] += c * stage_vsigma[3 * g + 1];
            vsigma[3 * g + 2] += c * stage_vsigma[3 * g + 2];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
    }
}

auto
Functional::compute_exc(int32_t np, const double* rho, const double* sigma, const double* lapl, const double* tau, double* exc) const -> void
{
    errors::assertMsgCritical(_hasExc, std::string(__func__) + ": exchange-correlation functional does not provide an evaluator for Exc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];

    for (const auto& [coeff, func] : _components)
    {
        xc_mgga_exc(&func, np, rho, sigma, lapl, tau, stage_exc);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
    }
}

auto
Functional::compute_vxc(int32_t       np,
                        const double* rho,
                        const double* sigma,
                        const double* lapl,
                        const double* tau,
                        double*       vrho,
                        double*       vsigma,
                        double*       vlapl,
                        double*       vtau) const -> void
{
    errors::assertMsgCritical(_hasVxc, std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_vlapl  = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_vtau   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_mgga_vxc(&func, np, rho, sigma, lapl, tau, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

        const auto c = coeff;

#pragma omp simd aligned(vrho, stage_vrho, vsigma, stage_vsigma, vlapl, stage_vlapl, vtau, stage_vtau : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];

            vsigma[3 * g + 0] += c * stage_vsigma[3 * g + 0];
            vsigma[3 * g + 1] += c * stage_vsigma[3 * g + 1];
            vsigma[3 * g + 2] += c * stage_vsigma[3 * g + 2];

            vlapl[2 * g + 0] += c * stage_vlapl[2 * g + 0];
            vlapl[2 * g + 1] += c * stage_vlapl[2 * g + 1];

            vtau[2 * g + 0] += c * stage_vtau[2 * g + 0];
            vtau[2 * g + 1] += c * stage_vtau[2 * g + 1];
        }
    }

    if (alloc)
    {
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
        mem::free(stage_vlapl);
        mem::free(stage_vtau);
    }
}

auto
Functional::compute_exc_vxc(int32_t       np,
                            const double* rho,
                            const double* sigma,
                            const double* lapl,
                            const double* tau,
                            double*       exc,
                            double*       vrho,
                            double*       vsigma,
                            double*       vlapl,
                            double*       vtau) const -> void
{
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(np) : &_stagingBuffer[0];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_vlapl  = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_vtau   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];

    for (const auto& [coeff, func] : _components)
    {
        xc_mgga_exc_vxc(&func, np, rho, sigma, lapl, tau, stage_exc, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

        const auto c = coeff;

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho, vsigma, stage_vsigma, vlapl, stage_vlapl, vtau, stage_vtau : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            exc[g] += c * stage_exc[g];

            vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
            vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];

            vsigma[3 * g + 0] += c * stage_vsigma[3 * g + 0];
            vsigma[3 * g + 1] += c * stage_vsigma[3 * g + 1];
            vsigma[3 * g + 2] += c * stage_vsigma[3 * g + 2];

            vlapl[2 * g + 0] += c * stage_vlapl[2 * g + 0];
            vlapl[2 * g + 1] += c * stage_vlapl[2 * g + 1];

            vtau[2 * g + 0] += c * stage_vtau[2 * g + 0];
            vtau[2 * g + 1] += c * stage_vtau[2 * g + 1];
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
        mem::free(stage_vlapl);
        mem::free(stage_vtau);
    }
}

std::ostream&
operator<<(std::ostream& output, const Functional& source)
{
    return (output << source.repr());
}

auto
getLibXCDescription() -> std::string
{
    std::ostringstream os;

    os << "LibXC version " << std::string(xc_version_string()) << std::endl;

    os << xc_reference() << " (" << xc_reference_doi() << ")";

    return os.str();
}

namespace detail {
auto
getXCName(const std::vector<Component>& x_funcs, const std::vector<Component>& c_funcs) -> std::string
{
    // create label for this mix of functionals
    std::ostringstream os;

    auto i = 1;
    for (const auto& [coeff, func] : x_funcs)
    {
        auto name = std::string(xc_functional_get_name(func.info->number));

        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

        os << std::setprecision(4) << coeff << "*" << name << (i != x_funcs.size() ? " + " : "");
        ++i;
    }

    os << ", ";

    i = 1;
    for (const auto& [coeff, func] : c_funcs)
    {
        auto name = std::string(xc_functional_get_name(func.info->number));

        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::toupper(c); });

        os << std::setprecision(4) << coeff << "*" << name << (i != c_funcs.size() ? " + " : "");
        ++i;
    }

    return os.str();
}

auto
getFunctionalCitation(const xc_func_type& func) -> std::string
{
    std::ostringstream os;

    for (auto i = 0; i < XC_MAX_REFERENCES; ++i)
    {
        if (func.info->refs[i])
        {
            if (i != 0)
            {
                os << std::endl;
            }
            os << "    ";
            os << func.info->refs[i]->ref;
            if (std::strlen(func.info->refs[i]->doi) > 0)
            {
                os << " (";
                os << func.info->refs[i]->doi;
                os << ")";
            }
        }
    }

    return os.str();
}

auto
getCitation(const std::vector<Component>& x_funcs, const std::vector<Component>& c_funcs) -> std::string
{
    std::ostringstream os;

    os << "Exchange components: " << std::endl;

    for (const auto& [_, func] : x_funcs)
    {
        os << " * " << std::string(func.info->name) << std::endl;
        os << detail::getFunctionalCitation(func) << std::endl;
    }

    os << std::endl;

    os << "Correlation components: " << std::endl;

    for (const auto& [_, func] : c_funcs)
    {
        os << " * " << std::string(func.info->name) << std::endl;
        os << detail::getFunctionalCitation(func) << std::endl;
    }

    return os.str();
}
}  // namespace detail
