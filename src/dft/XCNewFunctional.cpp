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

#include "XCNewFunctional.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ErrorHandler.hpp"
#include "MemAlloc.hpp"
#include "StringFormat.hpp"

CXCNewFunctional::CXCNewFunctional(const std::string&              nameOfFunctional,
                                   const std::vector<std::string>& labels,
                                   const std::vector<double>&      coeffs,
                                   const double                    fractionOfExactExchange,
                                   const double                    rangeSeparationParameter)

    : _nameOfFunctional(fstr::upcase(nameOfFunctional))

    , _fractionOfExactExchange(fractionOfExactExchange)

    , _rangeSeparationParameter(rangeSeparationParameter)
{
    std::string errmsg("XCNewFunctional: Inconsistent sizes of functional labels and coefficients");

    errors::assertMsgCritical(labels.size() == coeffs.size(), errmsg);

    _components.clear();

    bool hasExc = true, hasVxc = true, hasFxc = true, hasKxc = true, hasLxc = true;

    bool isLDA = false, isGGA = false, isMGGA = false;

    for (int32_t i = 0; i < static_cast<int32_t>(labels.size()); i++)
    {
        auto label = labels[i];

        auto coeff = coeffs[i];

        auto xccomp = CXCComponent(label, coeff);

        auto funcptr = xccomp.getFunctionalPointer();

        auto kind = funcptr->info->kind;

        if ((kind == XC_EXCHANGE) || (kind == XC_CORRELATION))
        {
            _components.push_back(xccomp);
        }
        else if (kind == XC_EXCHANGE_CORRELATION)
        {
            errors::assertMsgCritical(labels.size() == 1,
                                      std::string("XCNewFunctional: Cannot mix ") + label + std::string(" with other functionals"));

            _components.push_back(xccomp);
        }
        else
        {
            errors::assertMsgCritical(false, std::string("XCNewFunctional: Unsupported functional ") + label);
        }

        auto flags = funcptr->info->flags;

        // which derivative orders do we have for this x-c mixture?
        hasExc = hasExc && (flags & XC_FLAGS_HAVE_EXC);
        hasVxc = hasVxc && (flags & XC_FLAGS_HAVE_VXC);
        hasFxc = hasFxc && (flags & XC_FLAGS_HAVE_FXC);
        hasKxc = hasKxc && (flags & XC_FLAGS_HAVE_KXC);
        hasLxc = hasLxc && (flags & XC_FLAGS_HAVE_LXC);

        // which family does this x-c mixture belong to?
        auto family = funcptr->info->family;

        // LDA, GGA, metaGGA
        isMGGA = (isMGGA || (family == XC_FAMILY_MGGA));
        isGGA  = ((!isMGGA) && (isGGA || (family == XC_FAMILY_GGA)));
        isLDA  = ((!isMGGA) && (!isGGA) && (isLDA || (family == XC_FAMILY_LDA)));
    }

    if (hasExc) _maxDerivOrder = 0;
    if (hasVxc) _maxDerivOrder = 1;
    if (hasFxc) _maxDerivOrder = 2;
    if (hasKxc) _maxDerivOrder = 3;
    if (hasLxc) _maxDerivOrder = 4;

    if (isLDA) _familyOfFunctional = std::string("LDA");
    if (isGGA) _familyOfFunctional = std::string("GGA");
    if (isMGGA) _familyOfFunctional = std::string("MGGA");

    // figure out fraction of exact exchange from "prebaked" functional (e.g. HYB_GGA_XC_B3LYP)
    if (_components.size() == 1)
    {
        auto funcptr = _components[0].getFunctionalPointer();
        if (funcptr->info->kind == XC_EXCHANGE_CORRELATION)
        {
            _fractionOfExactExchange = xc_hyb_exx_coef(funcptr);
        }
    }

    // TODO figure out whether a functional is range-separated (e.g. HYB_GGA_XC_LRC_WPBEH)

    _allocateStagingBuffer();
}

CXCNewFunctional::CXCNewFunctional(const CXCNewFunctional& source)

    : _nameOfFunctional(source._nameOfFunctional)

    , _fractionOfExactExchange(source._fractionOfExactExchange)

    , _rangeSeparationParameter(source._rangeSeparationParameter)

    , _maxDerivOrder(source._maxDerivOrder)

    , _familyOfFunctional(source._familyOfFunctional)

    , _ldStaging(source._ldStaging)

    , _components(source._components)
{
    _allocateStagingBuffer();
}

CXCNewFunctional::CXCNewFunctional(CXCNewFunctional&& source) noexcept

    : _nameOfFunctional(std::move(source._nameOfFunctional))

    , _fractionOfExactExchange(std::move(source._fractionOfExactExchange))

    , _rangeSeparationParameter(std::move(source._rangeSeparationParameter))

    , _maxDerivOrder(std::move(source._maxDerivOrder))

    , _familyOfFunctional(std::move(source._familyOfFunctional))

    , _ldStaging(std::move(source._ldStaging))

    , _components(std::move(source._components))
{
    _allocateStagingBuffer();

    source._freeStagingBuffer();
}

CXCNewFunctional::~CXCNewFunctional()
{
    _components.clear();

    _freeStagingBuffer();
}

void
CXCNewFunctional::_allocateStagingBuffer()
{
    if (_stagingBuffer == nullptr)
    {
        // TODO write function to compute this number based on functional
        // family and order of derivatives available
        int32_t n_xc_outputs = 0;

        if (_familyOfFunctional == std::string("LDA")) n_xc_outputs = 15;

        if (_familyOfFunctional == std::string("GGA")) n_xc_outputs = 126;

        if (_familyOfFunctional == std::string("MGGA")) n_xc_outputs = 767;

        _stagingBuffer = mem::malloc<double>(n_xc_outputs * _ldStaging);
    }
}

void
CXCNewFunctional::_freeStagingBuffer()
{
    if (_stagingBuffer != nullptr)
    {
        mem::free(_stagingBuffer);

        _stagingBuffer = nullptr;
    }
}

CXCNewFunctional&
CXCNewFunctional::operator=(const CXCNewFunctional& source)
{
    if (this == &source) return *this;

    _nameOfFunctional = source._nameOfFunctional;

    _fractionOfExactExchange = source._fractionOfExactExchange;

    _rangeSeparationParameter = source._rangeSeparationParameter;

    _maxDerivOrder = source._maxDerivOrder;

    _familyOfFunctional = source._familyOfFunctional;

    _ldStaging = source._ldStaging;

    _components = source._components;

    _freeStagingBuffer();

    _allocateStagingBuffer();

    return *this;
}

CXCNewFunctional&
CXCNewFunctional::operator=(CXCNewFunctional&& source) noexcept
{
    if (this == &source) return *this;

    _nameOfFunctional = std::move(source._nameOfFunctional);

    _fractionOfExactExchange = std::move(source._fractionOfExactExchange);

    _rangeSeparationParameter = std::move(source._rangeSeparationParameter);

    _maxDerivOrder = std::move(source._maxDerivOrder);

    _familyOfFunctional = std::move(source._familyOfFunctional);

    _ldStaging = std::move(source._ldStaging);

    _components = std::move(source._components);

    _freeStagingBuffer();

    _allocateStagingBuffer();

    source._freeStagingBuffer();

    return *this;
}

bool
CXCNewFunctional::operator==(const CXCNewFunctional& other) const
{
    if (_nameOfFunctional != other._nameOfFunctional) return false;

    if (_fractionOfExactExchange != other._fractionOfExactExchange) return false;

    if (_rangeSeparationParameter != other._rangeSeparationParameter) return false;

    if (_maxDerivOrder != other._maxDerivOrder) return false;

    if (_familyOfFunctional != other._familyOfFunctional) return false;

    if (_ldStaging != other._ldStaging) return false;

    if (_components != other._components) return false;

    return true;
}

bool
CXCNewFunctional::operator!=(const CXCNewFunctional& other) const
{
    return !(*this == other);
}

std::string
CXCNewFunctional::getFunctionalLabel() const
{
    return _nameOfFunctional;
}

xcfun
CXCNewFunctional::getFunctionalType() const
{
    return to_xcfun(_familyOfFunctional);
}

bool
CXCNewFunctional::isUndefined() const
{
    return (fstr::upcase(_nameOfFunctional) == "UNDEFINED");
}

bool
CXCNewFunctional::isHybrid() const
{
    return (std::fabs(_fractionOfExactExchange) > 1.0e-13);
}

double
CXCNewFunctional::getFractionOfExactExchange() const
{
    return _fractionOfExactExchange;
}

auto
CXCNewFunctional::compute_exc_vxc_for_lda(int32_t np, const double* rho, double* exc, double* vrho) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc  = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];

#pragma omp simd aligned(exc, vrho : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        exc[g] = 0.0;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        xc_lda_exc_vxc(funcptr, np, rho, stage_exc, stage_vrho);

        const auto c = xccomp.getScalingFactor();

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
CXCNewFunctional::compute_fxc_for_lda(int32_t np, const double* rho, double* v2rho2) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc      (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho     (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_v2rho2 = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

#pragma omp simd aligned(v2rho2 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v2rho2[3 * g + 0] = 0.0;
        v2rho2[3 * g + 1] = 0.0;
        v2rho2[3 * g + 2] = 0.0;
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        xc_lda_fxc(funcptr, np, rho, stage_v2rho2);

        const auto c = xccomp.getScalingFactor();

#pragma omp simd aligned(v2rho2, stage_v2rho2 : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            v2rho2[3 * g + 0] += c * stage_v2rho2[3 * g + 0];
            v2rho2[3 * g + 1] += c * stage_v2rho2[3 * g + 1];
            v2rho2[3 * g + 2] += c * stage_v2rho2[3 * g + 2];
        }
    }

    if (alloc)
    {
        mem::free(stage_v2rho2);
    }
}

auto
CXCNewFunctional::compute_kxc_for_lda(int32_t np, const double* rho, double* v3rho3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    // stage_exc        (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    // stage_vrho       (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    // stage_v2rho2     (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_v3rho3 = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[6 * _ldStaging];

#pragma omp simd aligned(v3rho3 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v3rho3[4 * g + 0] = 0.0;
        v3rho3[4 * g + 1] = 0.0;
        v3rho3[4 * g + 2] = 0.0;
        v3rho3[4 * g + 3] = 0.0;
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        xc_lda_kxc(funcptr, np, rho, stage_v3rho3);

        const auto c = xccomp.getScalingFactor();

#pragma omp simd aligned(v3rho3, stage_v3rho3 : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            v3rho3[4 * g + 0] += c * stage_v3rho3[4 * g + 0];
            v3rho3[4 * g + 1] += c * stage_v3rho3[4 * g + 1];
            v3rho3[4 * g + 2] += c * stage_v3rho3[4 * g + 2];
            v3rho3[4 * g + 3] += c * stage_v3rho3[4 * g + 3];
        }
    }

    if (alloc)
    {
        mem::free(stage_v3rho3);
    }
}

auto
CXCNewFunctional::compute_exc_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

#pragma omp simd aligned(exc, vrho, vsigma : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        exc[g] = 0.0;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_exc_vxc(funcptr, np, rho, stage_exc, stage_vrho);

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                exc[g] += c * stage_exc[g];

                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_exc_vxc(funcptr, np, rho, sigma, stage_exc, stage_vrho, stage_vsigma);

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
    }

    if (alloc)
    {
        mem::free(stage_exc);
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
    }
}

auto
CXCNewFunctional::compute_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

#pragma omp simd aligned(vrho, vsigma : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc      (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

#pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_vxc(funcptr, np, rho, sigma, stage_vrho, stage_vsigma);

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
    }

    if (alloc)
    {
        mem::free(stage_vrho);
        mem::free(stage_vsigma);
    }
}

auto
CXCNewFunctional::compute_fxc_for_gga(int32_t np, const double* rho, const double* sigma, double* v2rho2, double* v2rhosigma, double* v2sigma2) const
    -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

#pragma omp simd aligned(v2rho2, v2rhosigma, v2sigma2 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v2rho2[3 * g + 0] = 0.0;
        v2rho2[3 * g + 1] = 0.0;
        v2rho2[3 * g + 2] = 0.0;

        v2rhosigma[6 * g + 0] = 0.0;
        v2rhosigma[6 * g + 1] = 0.0;
        v2rhosigma[6 * g + 2] = 0.0;
        v2rhosigma[6 * g + 3] = 0.0;
        v2rhosigma[6 * g + 4] = 0.0;
        v2rhosigma[6 * g + 5] = 0.0;

        v2sigma2[6 * g + 0] = 0.0;
        v2sigma2[6 * g + 1] = 0.0;
        v2sigma2[6 * g + 2] = 0.0;
        v2sigma2[6 * g + 3] = 0.0;
        v2sigma2[6 * g + 4] = 0.0;
        v2sigma2[6 * g + 5] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc          (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho         (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    //   stage_vsigma       (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_v2rho2     = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_v2rhosigma = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[9 * _ldStaging];
    auto stage_v2sigma2   = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[15 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_fxc(funcptr, np, rho, stage_v2rho2);

#pragma omp simd aligned(v2rho2, stage_v2rho2 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v2rho2[3 * g + 0] += c * stage_v2rho2[3 * g + 0];
                v2rho2[3 * g + 1] += c * stage_v2rho2[3 * g + 1];
                v2rho2[3 * g + 2] += c * stage_v2rho2[3 * g + 2];
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_fxc(funcptr, np, rho, sigma, stage_v2rho2, stage_v2rhosigma, stage_v2sigma2);

#pragma omp simd aligned(v2rho2, stage_v2rho2, v2rhosigma, stage_v2rhosigma, v2sigma2, stage_v2sigma2 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v2rho2[3 * g + 0] += c * stage_v2rho2[3 * g + 0];
                v2rho2[3 * g + 1] += c * stage_v2rho2[3 * g + 1];
                v2rho2[3 * g + 2] += c * stage_v2rho2[3 * g + 2];

                v2rhosigma[6 * g + 0] += c * stage_v2rhosigma[6 * g + 0];
                v2rhosigma[6 * g + 1] += c * stage_v2rhosigma[6 * g + 1];
                v2rhosigma[6 * g + 2] += c * stage_v2rhosigma[6 * g + 2];
                v2rhosigma[6 * g + 3] += c * stage_v2rhosigma[6 * g + 3];
                v2rhosigma[6 * g + 4] += c * stage_v2rhosigma[6 * g + 4];
                v2rhosigma[6 * g + 5] += c * stage_v2rhosigma[6 * g + 5];

                v2sigma2[6 * g + 0] += c * stage_v2sigma2[6 * g + 0];
                v2sigma2[6 * g + 1] += c * stage_v2sigma2[6 * g + 1];
                v2sigma2[6 * g + 2] += c * stage_v2sigma2[6 * g + 2];
                v2sigma2[6 * g + 3] += c * stage_v2sigma2[6 * g + 3];
                v2sigma2[6 * g + 4] += c * stage_v2sigma2[6 * g + 4];
                v2sigma2[6 * g + 5] += c * stage_v2sigma2[6 * g + 5];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v2rho2);
        mem::free(stage_v2rhosigma);
        mem::free(stage_v2sigma2);
    }
}

auto
CXCNewFunctional::compute_kxc_for_gga(int32_t       np,
                                      const double* rho,
                                      const double* sigma,
                                      double*       v3rho3,
                                      double*       v3rho2sigma,
                                      double*       v3rhosigma2,
                                      double*       v3sigma3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");

#pragma omp simd aligned(v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v3rho3[4 * g + 0] = 0.0;
        v3rho3[4 * g + 1] = 0.0;
        v3rho3[4 * g + 2] = 0.0;
        v3rho3[4 * g + 3] = 0.0;

        v3rho2sigma[9 * g + 0] = 0.0;
        v3rho2sigma[9 * g + 1] = 0.0;
        v3rho2sigma[9 * g + 2] = 0.0;
        v3rho2sigma[9 * g + 3] = 0.0;
        v3rho2sigma[9 * g + 4] = 0.0;
        v3rho2sigma[9 * g + 5] = 0.0;
        v3rho2sigma[9 * g + 6] = 0.0;
        v3rho2sigma[9 * g + 7] = 0.0;
        v3rho2sigma[9 * g + 8] = 0.0;

        v3rhosigma2[12 * g + 0]  = 0.0;
        v3rhosigma2[12 * g + 1]  = 0.0;
        v3rhosigma2[12 * g + 2]  = 0.0;
        v3rhosigma2[12 * g + 3]  = 0.0;
        v3rhosigma2[12 * g + 4]  = 0.0;
        v3rhosigma2[12 * g + 5]  = 0.0;
        v3rhosigma2[12 * g + 6]  = 0.0;
        v3rhosigma2[12 * g + 7]  = 0.0;
        v3rhosigma2[12 * g + 8]  = 0.0;
        v3rhosigma2[12 * g + 9]  = 0.0;
        v3rhosigma2[12 * g + 10] = 0.0;
        v3rhosigma2[12 * g + 11] = 0.0;

        v3sigma3[10 * g + 0] = 0.0;
        v3sigma3[10 * g + 1] = 0.0;
        v3sigma3[10 * g + 2] = 0.0;
        v3sigma3[10 * g + 3] = 0.0;
        v3sigma3[10 * g + 4] = 0.0;
        v3sigma3[10 * g + 5] = 0.0;
        v3sigma3[10 * g + 6] = 0.0;
        v3sigma3[10 * g + 7] = 0.0;
        v3sigma3[10 * g + 8] = 0.0;
        v3sigma3[10 * g + 9] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc           (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho          (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    //   stage_vsigma        (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    //   stage_v2rho2        (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[6 * _ldStaging];
    //   stage_v2rhosigma    (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[9 * _ldStaging];
    //   stage_v2sigma2      (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[15 * _ldStaging];
    auto stage_v3rho3      = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[21 * _ldStaging];
    auto stage_v3rho2sigma = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[25 * _ldStaging];
    auto stage_v3rhosigma2 = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[34 * _ldStaging];
    auto stage_v3sigma3    = (alloc) ? mem::malloc<double>(10 * np) : &_stagingBuffer[46 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_kxc(funcptr, np, rho, stage_v3rho3);

#pragma omp simd aligned(v3rho3, stage_v3rho3 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v3rho3[4 * g + 0] += c * stage_v3rho3[4 * g + 0];
                v3rho3[4 * g + 1] += c * stage_v3rho3[4 * g + 1];
                v3rho3[4 * g + 2] += c * stage_v3rho3[4 * g + 2];
                v3rho3[4 * g + 3] += c * stage_v3rho3[4 * g + 3];
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_kxc(funcptr, np, rho, sigma, stage_v3rho3, stage_v3rho2sigma, stage_v3rhosigma2, stage_v3sigma3);

#pragma omp simd aligned(v3rho3, stage_v3rho3, v3rho2sigma, stage_v3rho2sigma, v3rhosigma2, stage_v3rhosigma2, v3sigma3, stage_v3sigma3 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v3rho3[4 * g + 0] += c * stage_v3rho3[4 * g + 0];
                v3rho3[4 * g + 1] += c * stage_v3rho3[4 * g + 1];
                v3rho3[4 * g + 2] += c * stage_v3rho3[4 * g + 2];
                v3rho3[4 * g + 3] += c * stage_v3rho3[4 * g + 3];

                v3rho2sigma[9 * g + 0] += c * stage_v3rho2sigma[9 * g + 0];
                v3rho2sigma[9 * g + 1] += c * stage_v3rho2sigma[9 * g + 1];
                v3rho2sigma[9 * g + 2] += c * stage_v3rho2sigma[9 * g + 2];
                v3rho2sigma[9 * g + 3] += c * stage_v3rho2sigma[9 * g + 3];
                v3rho2sigma[9 * g + 4] += c * stage_v3rho2sigma[9 * g + 4];
                v3rho2sigma[9 * g + 5] += c * stage_v3rho2sigma[9 * g + 5];
                v3rho2sigma[9 * g + 6] += c * stage_v3rho2sigma[9 * g + 6];
                v3rho2sigma[9 * g + 7] += c * stage_v3rho2sigma[9 * g + 7];
                v3rho2sigma[9 * g + 8] += c * stage_v3rho2sigma[9 * g + 8];

                v3rhosigma2[12 * g + 0] += c * stage_v3rhosigma2[12 * g + 0];
                v3rhosigma2[12 * g + 1] += c * stage_v3rhosigma2[12 * g + 1];
                v3rhosigma2[12 * g + 2] += c * stage_v3rhosigma2[12 * g + 2];
                v3rhosigma2[12 * g + 3] += c * stage_v3rhosigma2[12 * g + 3];
                v3rhosigma2[12 * g + 4] += c * stage_v3rhosigma2[12 * g + 4];
                v3rhosigma2[12 * g + 5] += c * stage_v3rhosigma2[12 * g + 5];
                v3rhosigma2[12 * g + 6] += c * stage_v3rhosigma2[12 * g + 6];
                v3rhosigma2[12 * g + 7] += c * stage_v3rhosigma2[12 * g + 7];
                v3rhosigma2[12 * g + 8] += c * stage_v3rhosigma2[12 * g + 8];
                v3rhosigma2[12 * g + 9] += c * stage_v3rhosigma2[12 * g + 9];
                v3rhosigma2[12 * g + 10] += c * stage_v3rhosigma2[12 * g + 10];
                v3rhosigma2[12 * g + 11] += c * stage_v3rhosigma2[12 * g + 11];

                v3sigma3[10 * g + 0] += c * stage_v3sigma3[10 * g + 0];
                v3sigma3[10 * g + 1] += c * stage_v3sigma3[10 * g + 1];
                v3sigma3[10 * g + 2] += c * stage_v3sigma3[10 * g + 2];
                v3sigma3[10 * g + 3] += c * stage_v3sigma3[10 * g + 3];
                v3sigma3[10 * g + 4] += c * stage_v3sigma3[10 * g + 4];
                v3sigma3[10 * g + 5] += c * stage_v3sigma3[10 * g + 5];
                v3sigma3[10 * g + 6] += c * stage_v3sigma3[10 * g + 6];
                v3sigma3[10 * g + 7] += c * stage_v3sigma3[10 * g + 7];
                v3sigma3[10 * g + 8] += c * stage_v3sigma3[10 * g + 8];
                v3sigma3[10 * g + 9] += c * stage_v3sigma3[10 * g + 9];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v3rho3);
        mem::free(stage_v3rho2sigma);
        mem::free(stage_v3rhosigma2);
        mem::free(stage_v3sigma3);
    }
}

auto
CXCNewFunctional::compute_exc_vxc_for_mgga(int32_t       np,
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
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

#pragma omp simd aligned(exc, vrho, vsigma, vlapl, vtau : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        exc[g] = 0.0;

        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;

        vsigma[3 * g + 0] = 0.0;
        vsigma[3 * g + 1] = 0.0;
        vsigma[3 * g + 2] = 0.0;

        vlapl[2 * g + 0] = 0.0;
        vlapl[2 * g + 1] = 0.0;

        vtau[2 * g + 0] = 0.0;
        vtau[2 * g + 1] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_vlapl  = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_vtau   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_exc_vxc(funcptr, np, rho, stage_exc, stage_vrho);

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                exc[g] += c * stage_exc[g];

                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_exc_vxc(funcptr, np, rho, sigma, stage_exc, stage_vrho, stage_vsigma);

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
        else if (family == XC_FAMILY_MGGA)
        {
            xc_mgga_exc_vxc(funcptr, np, rho, sigma, lapl, tau, stage_exc, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

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
