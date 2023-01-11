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
                                   const double                    fractionOfExactExchange)

    : _nameOfFunctional(fstr::upcase(nameOfFunctional))

    , _fractionOfExactExchange(fractionOfExactExchange)
{
    std::string errmsg("XCNewFunctional: Inconsistent sizes of functional labels and coefficients");

    errors::assertMsgCritical(labels.size() == coeffs.size(), errmsg);

    _components.clear();

    bool hasExc = true, hasVxc = true, hasFxc = true, hasKxc = true, hasLxc = true;

    bool isLDA = false, isGGA = false, isMGGA = false;

    bool needLaplacian = false;

    for (int32_t i = 0; i < static_cast<int32_t>(labels.size()); i++)
    {
        auto label = labels[i];

        auto coeff = coeffs[i];

        auto xccomp = CXCComponent(label, coeff);

        auto funcptr = xccomp.getFunctionalPointer();

        // check functional kind

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

        // check functional flags

        auto flags = funcptr->info->flags;

        // which derivative orders do we have for this x-c mixture?
        hasExc = (hasExc && (flags & XC_FLAGS_HAVE_EXC));
        hasVxc = (hasVxc && (flags & XC_FLAGS_HAVE_VXC));
        hasFxc = (hasFxc && (flags & XC_FLAGS_HAVE_FXC));
        hasKxc = (hasKxc && (flags & XC_FLAGS_HAVE_KXC));
        hasLxc = (hasLxc && (flags & XC_FLAGS_HAVE_LXC));

        // whether Laplacian is needed
        needLaplacian = (needLaplacian || (flags & XC_FLAGS_NEEDS_LAPLACIAN));

        // check functional family

        // which family does this x-c mixture belong to?
        auto family = funcptr->info->family;

        // LDA, GGA, metaGGA
        isMGGA = (isMGGA || (family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA));
        isGGA  = ((!isMGGA) && (isGGA || (family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA)));
        isLDA  = ((!isMGGA) && (!isGGA) && (isLDA || (family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA)));

        // check functional hybrid type
        // TODO use xc_hyb_type when it is available

        auto hyb_exx_coeff = xc_hyb_exx_coef(funcptr);

        if (hyb_exx_coeff != 0.0) _fractionOfExactExchange += coeff * hyb_exx_coeff;

        // TODO figure out range-separation parameters
    }

    if (hasExc) _maxDerivOrder = 0;
    if (hasVxc) _maxDerivOrder = 1;
    if (hasFxc) _maxDerivOrder = 2;
    if (hasKxc) _maxDerivOrder = 3;
    if (hasLxc) _maxDerivOrder = 4;

    if (isLDA) _familyOfFunctional = std::string("LDA");
    if (isGGA) _familyOfFunctional = std::string("GGA");
    if (isMGGA) _familyOfFunctional = std::string("MGGA");

    errors::assertMsgCritical(!needLaplacian, std::string("XCNewFunctional: Density Laplacian is not supported"));

    _allocateStagingBuffer();
}

CXCNewFunctional::CXCNewFunctional(const CXCNewFunctional& source)

    : _nameOfFunctional(source._nameOfFunctional)

    , _fractionOfExactExchange(source._fractionOfExactExchange)

    , _rangeSeparationParameters(source._rangeSeparationParameters)

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

    , _rangeSeparationParameters(std::move(source._rangeSeparationParameters))

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

    _rangeSeparationParameters = source._rangeSeparationParameters;

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

    _rangeSeparationParameters = std::move(source._rangeSeparationParameters);

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

    if (_rangeSeparationParameters != other._rangeSeparationParameters) return false;

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
CXCNewFunctional::compute_lxc_for_lda(int32_t np, const double* rho, double* v4rho4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    // stage_exc        (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    // stage_vrho       (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    // stage_v2rho2     (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    // stage_v3rho3     (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_v4rho4 = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[10 * _ldStaging];

#pragma omp simd aligned(v4rho4 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v4rho4[5 * g + 0] = 0.0;
        v4rho4[5 * g + 1] = 0.0;
        v4rho4[5 * g + 2] = 0.0;
        v4rho4[5 * g + 3] = 0.0;
        v4rho4[5 * g + 4] = 0.0;
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

        const auto c = xccomp.getScalingFactor();

#pragma omp simd aligned(v4rho4, stage_v4rho4 : VLX_ALIGN)
        for (auto g = 0; g < np; ++g)
        {
            v4rho4[5 * g + 0] += c * stage_v4rho4[5 * g + 0];
            v4rho4[5 * g + 1] += c * stage_v4rho4[5 * g + 1];
            v4rho4[5 * g + 2] += c * stage_v4rho4[5 * g + 2];
            v4rho4[5 * g + 3] += c * stage_v4rho4[5 * g + 3];
            v4rho4[5 * g + 4] += c * stage_v4rho4[5 * g + 4];
        }
    }

    if (alloc)
    {
        mem::free(stage_v4rho4);
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

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

#pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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
CXCNewFunctional::compute_lxc_for_gga(int32_t       np,
                                      const double* rho,
                                      const double* sigma,
                                      double*       v4rho4,
                                      double*       v4rho3sigma,
                                      double*       v4rho2sigma2,
                                      double*       v4rhosigma3,
                                      double*       v4sigma4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

#pragma omp simd aligned(v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v4rho4[5 * g + 0] = 0.0;
        v4rho4[5 * g + 1] = 0.0;
        v4rho4[5 * g + 2] = 0.0;
        v4rho4[5 * g + 3] = 0.0;
        v4rho4[5 * g + 4] = 0.0;

        v4rho3sigma[12 * g + 0]  = 0.0;
        v4rho3sigma[12 * g + 1]  = 0.0;
        v4rho3sigma[12 * g + 2]  = 0.0;
        v4rho3sigma[12 * g + 3]  = 0.0;
        v4rho3sigma[12 * g + 4]  = 0.0;
        v4rho3sigma[12 * g + 5]  = 0.0;
        v4rho3sigma[12 * g + 6]  = 0.0;
        v4rho3sigma[12 * g + 7]  = 0.0;
        v4rho3sigma[12 * g + 8]  = 0.0;
        v4rho3sigma[12 * g + 9]  = 0.0;
        v4rho3sigma[12 * g + 10] = 0.0;
        v4rho3sigma[12 * g + 11] = 0.0;

        v4rho2sigma2[18 * g + 0]  = 0.0;
        v4rho2sigma2[18 * g + 1]  = 0.0;
        v4rho2sigma2[18 * g + 2]  = 0.0;
        v4rho2sigma2[18 * g + 3]  = 0.0;
        v4rho2sigma2[18 * g + 4]  = 0.0;
        v4rho2sigma2[18 * g + 5]  = 0.0;
        v4rho2sigma2[18 * g + 6]  = 0.0;
        v4rho2sigma2[18 * g + 7]  = 0.0;
        v4rho2sigma2[18 * g + 8]  = 0.0;
        v4rho2sigma2[18 * g + 9]  = 0.0;
        v4rho2sigma2[18 * g + 10] = 0.0;
        v4rho2sigma2[18 * g + 11] = 0.0;
        v4rho2sigma2[18 * g + 12] = 0.0;
        v4rho2sigma2[18 * g + 13] = 0.0;
        v4rho2sigma2[18 * g + 14] = 0.0;
        v4rho2sigma2[18 * g + 15] = 0.0;
        v4rho2sigma2[18 * g + 16] = 0.0;
        v4rho2sigma2[18 * g + 17] = 0.0;

        v4rhosigma3[20 * g + 0]  = 0.0;
        v4rhosigma3[20 * g + 1]  = 0.0;
        v4rhosigma3[20 * g + 2]  = 0.0;
        v4rhosigma3[20 * g + 3]  = 0.0;
        v4rhosigma3[20 * g + 4]  = 0.0;
        v4rhosigma3[20 * g + 5]  = 0.0;
        v4rhosigma3[20 * g + 6]  = 0.0;
        v4rhosigma3[20 * g + 7]  = 0.0;
        v4rhosigma3[20 * g + 8]  = 0.0;
        v4rhosigma3[20 * g + 9]  = 0.0;
        v4rhosigma3[20 * g + 10] = 0.0;
        v4rhosigma3[20 * g + 11] = 0.0;
        v4rhosigma3[20 * g + 12] = 0.0;
        v4rhosigma3[20 * g + 13] = 0.0;
        v4rhosigma3[20 * g + 14] = 0.0;
        v4rhosigma3[20 * g + 15] = 0.0;
        v4rhosigma3[20 * g + 16] = 0.0;
        v4rhosigma3[20 * g + 17] = 0.0;
        v4rhosigma3[20 * g + 18] = 0.0;
        v4rhosigma3[20 * g + 19] = 0.0;

        v4sigma4[15 * g + 0]  = 0.0;
        v4sigma4[15 * g + 1]  = 0.0;
        v4sigma4[15 * g + 2]  = 0.0;
        v4sigma4[15 * g + 3]  = 0.0;
        v4sigma4[15 * g + 4]  = 0.0;
        v4sigma4[15 * g + 5]  = 0.0;
        v4sigma4[15 * g + 6]  = 0.0;
        v4sigma4[15 * g + 7]  = 0.0;
        v4sigma4[15 * g + 8]  = 0.0;
        v4sigma4[15 * g + 9]  = 0.0;
        v4sigma4[15 * g + 10] = 0.0;
        v4sigma4[15 * g + 11] = 0.0;
        v4sigma4[15 * g + 12] = 0.0;
        v4sigma4[15 * g + 13] = 0.0;
        v4sigma4[15 * g + 14] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc            (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho           (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    //   stage_vsigma         (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    //   stage_v2rho2         (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[6 * _ldStaging];
    //   stage_v2rhosigma     (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[9 * _ldStaging];
    //   stage_v2sigma2       (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[15 * _ldStaging];
    //   stage_v3rho3         (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[21 * _ldStaging];
    //   stage_v3rho2sigma    (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[25 * _ldStaging];
    //   stage_v3rhosigma2    (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[34 * _ldStaging];
    //   stage_v3sigma3       (alloc) ? mem::malloc<double>(10 * np) : &_stagingBuffer[46 * _ldStaging];
    auto stage_v4rho4       = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[56 * _ldStaging];
    auto stage_v4rho3sigma  = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[61 * _ldStaging];
    auto stage_v4rho2sigma2 = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[73 * _ldStaging];
    auto stage_v4rhosigma3  = (alloc) ? mem::malloc<double>(20 * np) : &_stagingBuffer[91 * _ldStaging];
    auto stage_v4sigma4     = (alloc) ? mem::malloc<double>(15 * np) : &_stagingBuffer[111 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
        {
            xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

#pragma omp simd aligned(v4rho4, stage_v4rho4 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v4rho4[5 * g + 0] += c * stage_v4rho4[5 * g + 0];
                v4rho4[5 * g + 1] += c * stage_v4rho4[5 * g + 1];
                v4rho4[5 * g + 2] += c * stage_v4rho4[5 * g + 2];
                v4rho4[5 * g + 3] += c * stage_v4rho4[5 * g + 3];
                v4rho4[5 * g + 4] += c * stage_v4rho4[5 * g + 4];
            }
        }
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
        {
            xc_gga_lxc(funcptr, np, rho, sigma, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4);

#pragma omp simd aligned(                                                                                                                            \
    v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4 \
    : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v4rho4[5 * g + 0] += c * stage_v4rho4[5 * g + 0];
                v4rho4[5 * g + 1] += c * stage_v4rho4[5 * g + 1];
                v4rho4[5 * g + 2] += c * stage_v4rho4[5 * g + 2];
                v4rho4[5 * g + 3] += c * stage_v4rho4[5 * g + 3];
                v4rho4[5 * g + 4] += c * stage_v4rho4[5 * g + 4];

                v4rho3sigma[12 * g + 0] += c * stage_v4rho3sigma[12 * g + 0];
                v4rho3sigma[12 * g + 1] += c * stage_v4rho3sigma[12 * g + 1];
                v4rho3sigma[12 * g + 2] += c * stage_v4rho3sigma[12 * g + 2];
                v4rho3sigma[12 * g + 3] += c * stage_v4rho3sigma[12 * g + 3];
                v4rho3sigma[12 * g + 4] += c * stage_v4rho3sigma[12 * g + 4];
                v4rho3sigma[12 * g + 5] += c * stage_v4rho3sigma[12 * g + 5];
                v4rho3sigma[12 * g + 6] += c * stage_v4rho3sigma[12 * g + 6];
                v4rho3sigma[12 * g + 7] += c * stage_v4rho3sigma[12 * g + 7];
                v4rho3sigma[12 * g + 8] += c * stage_v4rho3sigma[12 * g + 8];
                v4rho3sigma[12 * g + 9] += c * stage_v4rho3sigma[12 * g + 9];
                v4rho3sigma[12 * g + 10] += c * stage_v4rho3sigma[12 * g + 10];
                v4rho3sigma[12 * g + 11] += c * stage_v4rho3sigma[12 * g + 11];

                v4rho2sigma2[18 * g + 0] += c * stage_v4rho2sigma2[18 * g + 0];
                v4rho2sigma2[18 * g + 1] += c * stage_v4rho2sigma2[18 * g + 1];
                v4rho2sigma2[18 * g + 2] += c * stage_v4rho2sigma2[18 * g + 2];
                v4rho2sigma2[18 * g + 3] += c * stage_v4rho2sigma2[18 * g + 3];
                v4rho2sigma2[18 * g + 4] += c * stage_v4rho2sigma2[18 * g + 4];
                v4rho2sigma2[18 * g + 5] += c * stage_v4rho2sigma2[18 * g + 5];
                v4rho2sigma2[18 * g + 6] += c * stage_v4rho2sigma2[18 * g + 6];
                v4rho2sigma2[18 * g + 7] += c * stage_v4rho2sigma2[18 * g + 7];
                v4rho2sigma2[18 * g + 8] += c * stage_v4rho2sigma2[18 * g + 8];
                v4rho2sigma2[18 * g + 9] += c * stage_v4rho2sigma2[18 * g + 9];
                v4rho2sigma2[18 * g + 10] += c * stage_v4rho2sigma2[18 * g + 10];
                v4rho2sigma2[18 * g + 11] += c * stage_v4rho2sigma2[18 * g + 11];
                v4rho2sigma2[18 * g + 12] += c * stage_v4rho2sigma2[18 * g + 12];
                v4rho2sigma2[18 * g + 13] += c * stage_v4rho2sigma2[18 * g + 13];
                v4rho2sigma2[18 * g + 14] += c * stage_v4rho2sigma2[18 * g + 14];
                v4rho2sigma2[18 * g + 15] += c * stage_v4rho2sigma2[18 * g + 15];
                v4rho2sigma2[18 * g + 16] += c * stage_v4rho2sigma2[18 * g + 16];
                v4rho2sigma2[18 * g + 17] += c * stage_v4rho2sigma2[18 * g + 17];

                v4rhosigma3[20 * g + 0] += c * stage_v4rhosigma3[20 * g + 0];
                v4rhosigma3[20 * g + 1] += c * stage_v4rhosigma3[20 * g + 1];
                v4rhosigma3[20 * g + 2] += c * stage_v4rhosigma3[20 * g + 2];
                v4rhosigma3[20 * g + 3] += c * stage_v4rhosigma3[20 * g + 3];
                v4rhosigma3[20 * g + 4] += c * stage_v4rhosigma3[20 * g + 4];
                v4rhosigma3[20 * g + 5] += c * stage_v4rhosigma3[20 * g + 5];
                v4rhosigma3[20 * g + 6] += c * stage_v4rhosigma3[20 * g + 6];
                v4rhosigma3[20 * g + 7] += c * stage_v4rhosigma3[20 * g + 7];
                v4rhosigma3[20 * g + 8] += c * stage_v4rhosigma3[20 * g + 8];
                v4rhosigma3[20 * g + 9] += c * stage_v4rhosigma3[20 * g + 9];
                v4rhosigma3[20 * g + 10] += c * stage_v4rhosigma3[20 * g + 10];
                v4rhosigma3[20 * g + 11] += c * stage_v4rhosigma3[20 * g + 11];
                v4rhosigma3[20 * g + 12] += c * stage_v4rhosigma3[20 * g + 12];
                v4rhosigma3[20 * g + 13] += c * stage_v4rhosigma3[20 * g + 13];
                v4rhosigma3[20 * g + 14] += c * stage_v4rhosigma3[20 * g + 14];
                v4rhosigma3[20 * g + 15] += c * stage_v4rhosigma3[20 * g + 15];
                v4rhosigma3[20 * g + 16] += c * stage_v4rhosigma3[20 * g + 16];
                v4rhosigma3[20 * g + 17] += c * stage_v4rhosigma3[20 * g + 17];
                v4rhosigma3[20 * g + 18] += c * stage_v4rhosigma3[20 * g + 18];
                v4rhosigma3[20 * g + 19] += c * stage_v4rhosigma3[20 * g + 19];

                v4sigma4[15 * g + 0] += c * stage_v4sigma4[15 * g + 0];
                v4sigma4[15 * g + 1] += c * stage_v4sigma4[15 * g + 1];
                v4sigma4[15 * g + 2] += c * stage_v4sigma4[15 * g + 2];
                v4sigma4[15 * g + 3] += c * stage_v4sigma4[15 * g + 3];
                v4sigma4[15 * g + 4] += c * stage_v4sigma4[15 * g + 4];
                v4sigma4[15 * g + 5] += c * stage_v4sigma4[15 * g + 5];
                v4sigma4[15 * g + 6] += c * stage_v4sigma4[15 * g + 6];
                v4sigma4[15 * g + 7] += c * stage_v4sigma4[15 * g + 7];
                v4sigma4[15 * g + 8] += c * stage_v4sigma4[15 * g + 8];
                v4sigma4[15 * g + 9] += c * stage_v4sigma4[15 * g + 9];
                v4sigma4[15 * g + 10] += c * stage_v4sigma4[15 * g + 10];
                v4sigma4[15 * g + 11] += c * stage_v4sigma4[15 * g + 11];
                v4sigma4[15 * g + 12] += c * stage_v4sigma4[15 * g + 12];
                v4sigma4[15 * g + 13] += c * stage_v4sigma4[15 * g + 13];
                v4sigma4[15 * g + 14] += c * stage_v4sigma4[15 * g + 14];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v4rho4);
        mem::free(stage_v4rho3sigma);
        mem::free(stage_v4rho2sigma2);
        mem::free(stage_v4rhosigma3);
        mem::free(stage_v4sigma4);
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

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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
        else if ((family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA))
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

auto
CXCNewFunctional::compute_vxc_for_mgga(int32_t       np,
                                       const double* rho,
                                       const double* sigma,
                                       const double* lapl,
                                       const double* tau,
                                       double*       vrho,
                                       double*       vsigma,
                                       double*       vlapl,
                                       double*       vtau) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

#pragma omp simd aligned(vrho, vsigma, vlapl, vtau : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
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

    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    auto stage_vlapl  = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_vtau   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

#pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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
        else if ((family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA))
        {
            xc_mgga_vxc(funcptr, np, rho, sigma, lapl, tau, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

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
CXCNewFunctional::compute_fxc_for_mgga(int32_t       np,
                                       const double* rho,
                                       const double* sigma,
                                       const double* lapl,
                                       const double* tau,
                                       double*       v2rho2,
                                       double*       v2rhosigma,
                                       double*       v2rholapl,
                                       double*       v2rhotau,
                                       double*       v2sigma2,
                                       double*       v2sigmalapl,
                                       double*       v2sigmatau,
                                       double*       v2lapl2,
                                       double*       v2lapltau,
                                       double*       v2tau2) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

#pragma omp simd aligned(v2rho2, v2rhosigma, v2sigma2,v2rholapl,\
                         v2rhotau, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2: VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v2rho2[3 * g + 0] = 0.0;
        v2rho2[3 * g + 1] = 0.0;
        v2rho2[3 * g + 2] = 0.0;

        v2rholapl[4 * g + 0] = 0.0;
        v2rholapl[4 * g + 1] = 0.0;
        v2rholapl[4 * g + 2] = 0.0;
        v2rholapl[4 * g + 3] = 0.0;

        v2rhotau[4 * g + 0] = 0.0;
        v2rhotau[4 * g + 1] = 0.0;
        v2rhotau[4 * g + 2] = 0.0;
        v2rhotau[4 * g + 3] = 0.0;

        v2lapl2[3 * g + 0] = 0.0;
        v2lapl2[3 * g + 1] = 0.0;
        v2lapl2[3 * g + 2] = 0.0;

        v2lapltau[4 * g + 0] = 0.0;
        v2lapltau[4 * g + 1] = 0.0;
        v2lapltau[4 * g + 2] = 0.0;
        v2lapltau[4 * g + 3] = 0.0;

        v2tau2[3 * g + 0] = 0.0;
        v2tau2[3 * g + 1] = 0.0;
        v2tau2[3 * g + 2] = 0.0;

        v2sigmalapl[6 * g + 0] = 0.0;
        v2sigmalapl[6 * g + 1] = 0.0;
        v2sigmalapl[6 * g + 2] = 0.0;
        v2sigmalapl[6 * g + 3] = 0.0;
        v2sigmalapl[6 * g + 4] = 0.0;
        v2sigmalapl[6 * g + 5] = 0.0;

        v2sigmatau[6 * g + 0] = 0.0;
        v2sigmatau[6 * g + 1] = 0.0;
        v2sigmatau[6 * g + 2] = 0.0;
        v2sigmatau[6 * g + 3] = 0.0;
        v2sigmatau[6 * g + 4] = 0.0;
        v2sigmatau[6 * g + 5] = 0.0;

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

    //   stage_exc           (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho          (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    //   stage_vsigma        (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    //   stage_vlapl         (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    //   stage_vtau          (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];
    auto stage_v2rho2      = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[10 * _ldStaging];
    auto stage_v2rhosigma  = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[13 * _ldStaging];
    auto stage_v2rholapl   = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[19 * _ldStaging];
    auto stage_v2rhotau    = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[23 * _ldStaging];
    auto stage_v2sigma2    = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[27 * _ldStaging];
    auto stage_v2sigmalapl = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[33 * _ldStaging];
    auto stage_v2sigmatau  = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[39 * _ldStaging];
    auto stage_v2lapl2     = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[45 * _ldStaging];
    auto stage_v2lapltau   = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[48 * _ldStaging];
    auto stage_v2tau2      = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[52 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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
        else if ((family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA))
        {
            xc_mgga_fxc(funcptr, np, rho, sigma, lapl, tau, stage_v2rho2, stage_v2rhosigma, stage_v2rholapl, stage_v2rhotau, stage_v2sigma2, stage_v2sigmalapl, stage_v2sigmatau, stage_v2lapl2, stage_v2lapltau,  stage_v2tau2);

#pragma omp simd aligned(v2rho2, v2rhosigma, v2sigma2,v2rholapl,\
                         v2rhotau, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,\
                         stage_v2rho2, stage_v2rhosigma, stage_v2sigma2,stage_v2rholapl,\
                         stage_v2rhotau, stage_v2sigmalapl, stage_v2sigmatau, stage_v2lapl2,\
                         stage_v2lapltau, stage_v2tau2 : VLX_ALIGN)
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

                v2rholapl[4 * g + 0] += c * stage_v2rholapl[4 * g + 0];
                v2rholapl[4 * g + 1] += c * stage_v2rholapl[4 * g + 1];
                v2rholapl[4 * g + 2] += c * stage_v2rholapl[4 * g + 2];
                v2rholapl[4 * g + 3] += c * stage_v2rholapl[4 * g + 3];

                v2rhotau[4 * g + 0] += c * stage_v2rhotau[4 * g + 0];
                v2rhotau[4 * g + 1] += c * stage_v2rhotau[4 * g + 1];
                v2rhotau[4 * g + 2] += c * stage_v2rhotau[4 * g + 2];
                v2rhotau[4 * g + 3] += c * stage_v2rhotau[4 * g + 3];

                v2lapl2[3 * g + 0] += c * stage_v2lapl2[3 * g + 0];
                v2lapl2[3 * g + 1] += c * stage_v2lapl2[3 * g + 1];
                v2lapl2[3 * g + 2] += c * stage_v2lapl2[3 * g + 2];

                v2lapltau[4 * g + 0] += c * stage_v2lapltau[4 * g + 0];
                v2lapltau[4 * g + 1] += c * stage_v2lapltau[4 * g + 1];
                v2lapltau[4 * g + 2] += c * stage_v2lapltau[4 * g + 2];
                v2lapltau[4 * g + 3] += c * stage_v2lapltau[4 * g + 3];

                v2tau2[3 * g + 0] += c * stage_v2tau2[3 * g + 0];
                v2tau2[3 * g + 1] += c * stage_v2tau2[3 * g + 1];
                v2tau2[3 * g + 2] += c * stage_v2tau2[3 * g + 2];

                v2sigmalapl[6 * g + 0] += c * stage_v2sigmalapl[6 * g + 0];
                v2sigmalapl[6 * g + 1] += c * stage_v2sigmalapl[6 * g + 1];
                v2sigmalapl[6 * g + 2] += c * stage_v2sigmalapl[6 * g + 2];
                v2sigmalapl[6 * g + 3] += c * stage_v2sigmalapl[6 * g + 3];
                v2sigmalapl[6 * g + 4] += c * stage_v2sigmalapl[6 * g + 4];
                v2sigmalapl[6 * g + 5] += c * stage_v2sigmalapl[6 * g + 5];

                v2sigmatau[6 * g + 0] += c * stage_v2sigmatau[6 * g + 0];
                v2sigmatau[6 * g + 1] += c * stage_v2sigmatau[6 * g + 1];
                v2sigmatau[6 * g + 2] += c * stage_v2sigmatau[6 * g + 2];
                v2sigmatau[6 * g + 3] += c * stage_v2sigmatau[6 * g + 3];
                v2sigmatau[6 * g + 4] += c * stage_v2sigmatau[6 * g + 4];
                v2sigmatau[6 * g + 5] += c * stage_v2sigmatau[6 * g + 5];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v2rho2);
        mem::free(stage_v2rhosigma);
        mem::free(stage_v2sigma2);
        mem::free(stage_v2lapltau);
        mem::free(stage_v2rhotau);
        mem::free(stage_v2rholapl);
        mem::free(stage_v2lapl2);
        mem::free(stage_v2sigmatau);
        mem::free(stage_v2sigmalapl);
        mem::free(stage_v2tau2);
    }
}

auto
CXCNewFunctional::compute_kxc_for_mgga(int32_t       np,
                                       const double* rho,
                                       const double* sigma,
                                       const double* lapl,
                                       const double* tau,
                                       double*       v3rho3,
                                       double*       v3rho2sigma,
                                       double*       v3rho2lapl,
                                       double*       v3rho2tau,
                                       double*       v3rhosigma2,
                                       double*       v3rhosigmalapl,
                                       double*       v3rhosigmatau,
                                       double*       v3rholapl2,
                                       double*       v3rholapltau,
                                       double*       v3rhotau2,
                                       double*       v3sigma3,
                                       double*       v3sigma2lapl,
                                       double*       v3sigma2tau,
                                       double*       v3sigmalapl2,
                                       double*       v3sigmalapltau,
                                       double*       v3sigmatau2,
                                       double*       v3lapl3,
                                       double*       v3lapl2tau,
                                       double*       v3lapltau2,
                                       double*       v3tau3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");
    
#pragma omp simd aligned(v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2,v3rhosigmalapl,v3rhosigmatau, \
                         v3rholapl2, v3rholapltau,v3rhotau2, v3sigma3, v3sigma2lapl, v3sigma2tau,v3sigmalapl2,\
                         v3sigmalapltau, v3sigmatau2, v3lapl3,v3lapl2tau,v3lapltau2,v3tau3 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v3rhosigma2[12 * g + 0] = 0.0;
        v3rhosigma2[12 * g + 1] = 0.0;
        v3rhosigma2[12 * g + 2] = 0.0;
        v3rhosigma2[12 * g + 3] = 0.0;
        v3rhosigma2[12 * g + 4] = 0.0;
        v3rhosigma2[12 * g + 5] = 0.0;
        v3rhosigma2[12 * g + 6] = 0.0;
        v3rhosigma2[12 * g + 7] = 0.0;
        v3rhosigma2[12 * g + 8] = 0.0;
        v3rhosigma2[12 * g + 9] = 0.0;
        v3rhosigma2[12 * g + 10] = 0.0;
        v3rhosigma2[12 * g + 11] = 0.0;
        
        v3rhosigmalapl[12 * g + 0] = 0.0;
        v3rhosigmalapl[12 * g + 1] = 0.0;
        v3rhosigmalapl[12 * g + 2] = 0.0;
        v3rhosigmalapl[12 * g + 3] = 0.0;
        v3rhosigmalapl[12 * g + 4] = 0.0;
        v3rhosigmalapl[12 * g + 5] = 0.0;
        v3rhosigmalapl[12 * g + 6] = 0.0;
        v3rhosigmalapl[12 * g + 7] = 0.0;
        v3rhosigmalapl[12 * g + 8] = 0.0;
        v3rhosigmalapl[12 * g + 9] = 0.0;
        v3rhosigmalapl[12 * g + 10] = 0.0;
        v3rhosigmalapl[12 * g + 11] = 0.0;
        
        v3rhosigmatau[12 * g + 0] = 0.0;
        v3rhosigmatau[12 * g + 1] = 0.0;
        v3rhosigmatau[12 * g + 2] = 0.0;
        v3rhosigmatau[12 * g + 3] = 0.0;
        v3rhosigmatau[12 * g + 4] = 0.0;
        v3rhosigmatau[12 * g + 5] = 0.0;
        v3rhosigmatau[12 * g + 6] = 0.0;
        v3rhosigmatau[12 * g + 7] = 0.0;
        v3rhosigmatau[12 * g + 8] = 0.0;
        v3rhosigmatau[12 * g + 9] = 0.0;
        v3rhosigmatau[12 * g + 10] = 0.0;
        v3rhosigmatau[12 * g + 11] = 0.0;

        v3sigma2lapl[12 * g + 0] = 0.0;
        v3sigma2lapl[12 * g + 1] = 0.0;
        v3sigma2lapl[12 * g + 2] = 0.0;
        v3sigma2lapl[12 * g + 3] = 0.0;
        v3sigma2lapl[12 * g + 4] = 0.0;
        v3sigma2lapl[12 * g + 5] = 0.0;
        v3sigma2lapl[12 * g + 6] = 0.0;
        v3sigma2lapl[12 * g + 7] = 0.0;
        v3sigma2lapl[12 * g + 8] = 0.0;
        v3sigma2lapl[12 * g + 9] = 0.0;
        v3sigma2lapl[12 * g + 10] = 0.0;
        v3sigma2lapl[12 * g + 11] = 0.0;
        
        v3sigma2tau[12 * g + 0] = 0.0;
        v3sigma2tau[12 * g + 1] = 0.0;
        v3sigma2tau[12 * g + 2] = 0.0;
        v3sigma2tau[12 * g + 3] = 0.0;
        v3sigma2tau[12 * g + 4] = 0.0;
        v3sigma2tau[12 * g + 5] = 0.0;
        v3sigma2tau[12 * g + 6] = 0.0;
        v3sigma2tau[12 * g + 7] = 0.0;
        v3sigma2tau[12 * g + 8] = 0.0;
        v3sigma2tau[12 * g + 9] = 0.0;
        v3sigma2tau[12 * g + 10] = 0.0;
        v3sigma2tau[12 * g + 11] = 0.0;

        v3sigmalapltau[12 * g + 0] = 0.0;
        v3sigmalapltau[12 * g + 1] = 0.0;
        v3sigmalapltau[12 * g + 2] = 0.0;
        v3sigmalapltau[12 * g + 3] = 0.0;
        v3sigmalapltau[12 * g + 4] = 0.0;
        v3sigmalapltau[12 * g + 5] = 0.0;
        v3sigmalapltau[12 * g + 6] = 0.0;
        v3sigmalapltau[12 * g + 7] = 0.0;
        v3sigmalapltau[12 * g + 8] = 0.0;
        v3sigmalapltau[12 * g + 9] = 0.0;
        v3sigmalapltau[12 * g + 10] = 0.0;
        v3sigmalapltau[12 * g + 11] = 0.0;

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

        v3rho2sigma[9 * g + 0] = 0.0;
        v3rho2sigma[9 * g + 1] = 0.0;
        v3rho2sigma[9 * g + 2] = 0.0;
        v3rho2sigma[9 * g + 3] = 0.0;
        v3rho2sigma[9 * g + 4] = 0.0;
        v3rho2sigma[9 * g + 5] = 0.0;
        v3rho2sigma[9 * g + 6] = 0.0;
        v3rho2sigma[9 * g + 7] = 0.0;
        v3rho2sigma[9 * g + 8] = 0.0;

        v3sigmalapl2[9 * g + 0] = 0.0;
        v3sigmalapl2[9 * g + 1] = 0.0;
        v3sigmalapl2[9 * g + 2] = 0.0;
        v3sigmalapl2[9 * g + 3] = 0.0;
        v3sigmalapl2[9 * g + 4] = 0.0;
        v3sigmalapl2[9 * g + 5] = 0.0;
        v3sigmalapl2[9 * g + 6] = 0.0;
        v3sigmalapl2[9 * g + 7] = 0.0;
        v3sigmalapl2[9 * g + 8] = 0.0;
        
        v3sigmatau2[9 * g + 0] = 0.0;
        v3sigmatau2[9 * g + 1] = 0.0;
        v3sigmatau2[9 * g + 2] = 0.0;
        v3sigmatau2[9 * g + 3] = 0.0;
        v3sigmatau2[9 * g + 4] = 0.0;
        v3sigmatau2[9 * g + 5] = 0.0;
        v3sigmatau2[9 * g + 6] = 0.0;
        v3sigmatau2[9 * g + 7] = 0.0;
        v3sigmatau2[9 * g + 8] = 0.0;

        v3lapl3[4 * g + 0] = 0.0;
        v3lapl3[4 * g + 1] = 0.0;
        v3lapl3[4 * g + 2] = 0.0;
        v3lapl3[4 * g + 3] = 0.0;

        v3tau3[4 * g + 0] = 0.0;
        v3tau3[4 * g + 1] = 0.0;
        v3tau3[4 * g + 2] = 0.0;
        v3tau3[4 * g + 3] = 0.0;

        v3rho3[4 * g + 0] = 0.0;
        v3rho3[4 * g + 1] = 0.0;
        v3rho3[4 * g + 2] = 0.0;
        v3rho3[4 * g + 3] = 0.0;
        
        v3rho2lapl[6 * g + 0] = 0.0;
        v3rho2lapl[6 * g + 1] = 0.0;
        v3rho2lapl[6 * g + 2] = 0.0;
        v3rho2lapl[6 * g + 3] = 0.0;
        v3rho2lapl[6 * g + 4] = 0.0;
        v3rho2lapl[6 * g + 5] = 0.0;
        
        v3rho2tau[6 * g + 0] = 0.0;
        v3rho2tau[6 * g + 1] = 0.0;
        v3rho2tau[6 * g + 2] = 0.0;
        v3rho2tau[6 * g + 3] = 0.0;
        v3rho2tau[6 * g + 4] = 0.0;
        v3rho2tau[6 * g + 5] = 0.0;
        
        v3rholapl2[6 * g + 0] = 0.0;
        v3rholapl2[6 * g + 1] = 0.0;
        v3rholapl2[6 * g + 2] = 0.0;
        v3rholapl2[6 * g + 3] = 0.0;
        v3rholapl2[6 * g + 4] = 0.0;
        v3rholapl2[6 * g + 5] = 0.0;
        
        v3rholapltau[8 * g + 0] = 0.0;
        v3rholapltau[8 * g + 1] = 0.0;
        v3rholapltau[8 * g + 2] = 0.0;
        v3rholapltau[8 * g + 3] = 0.0;
        v3rholapltau[8 * g + 4] = 0.0;
        v3rholapltau[8 * g + 5] = 0.0;
        v3rholapltau[8 * g + 6] = 0.0;
        v3rholapltau[8 * g + 7] = 0.0;
        
        v3rhotau2[6 * g + 0] = 0.0;
        v3rhotau2[6 * g + 1] = 0.0;
        v3rhotau2[6 * g + 2] = 0.0;
        v3rhotau2[6 * g + 3] = 0.0;
        v3rhotau2[6 * g + 4] = 0.0;
        v3rhotau2[6 * g + 5] = 0.0;

        v3lapl2tau[6 * g + 0] = 0.0;
        v3lapl2tau[6 * g + 1] = 0.0;
        v3lapl2tau[6 * g + 2] = 0.0;
        v3lapl2tau[6 * g + 3] = 0.0;
        v3lapl2tau[6 * g + 4] = 0.0;
        v3lapl2tau[6 * g + 5] = 0.0;
        
        v3lapltau2[6 * g + 0] = 0.0;
        v3lapltau2[6 * g + 1] = 0.0;
        v3lapltau2[6 * g + 2] = 0.0;
        v3lapltau2[6 * g + 3] = 0.0;
        v3lapltau2[6 * g + 4] = 0.0;
        v3lapltau2[6 * g + 5] = 0.0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    //   stage_exc              (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    //   stage_vrho             (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    //   stage_vsigma           (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];
    //   stage_vlapl            (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[6 * _ldStaging];
    //   stage_vtau             (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[8 * _ldStaging];
    //   stage_v2rho2           (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[10 * _ldStaging];
    //   stage_v2rhosigma       (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[13 * _ldStaging];
    //   stage_v2rholapl        (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[19 * _ldStaging];
    //   stage_v2rhotau         (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[23 * _ldStaging];
    //   stage_v2sigma2         (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[27 * _ldStaging];
    //   stage_v2sigmalapl      (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[33 * _ldStaging];
    //   stage_v2sigmatau       (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[39 * _ldStaging];
    //   stage_v2lapl2          (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[45 * _ldStaging];
    //   stage_v2lapltau        (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[48 * _ldStaging];
    //   stage_v2tau2           (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[52 * _ldStaging];
    auto stage_v3rho3         = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[55 * _ldStaging];
    auto stage_v3rho2sigma    = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[59 * _ldStaging];
    auto stage_v3rho2lapl     = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[68 * _ldStaging];
    auto stage_v3rho2tau      = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[74 * _ldStaging];
    auto stage_v3rhosigma2    = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[80 * _ldStaging];
    auto stage_v3rhosigmalapl = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[92 * _ldStaging];
    auto stage_v3rhosigmatau  = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[104 * _ldStaging];
    auto stage_v3rholapl2     = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[116 * _ldStaging];
    auto stage_v3rholapltau   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[122 * _ldStaging];
    auto stage_v3rhotau2      = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[130 * _ldStaging];
    auto stage_v3sigma3       = (alloc) ? mem::malloc<double>(10 * np) : &_stagingBuffer[136 * _ldStaging];
    auto stage_v3sigma2lapl   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[146 * _ldStaging];
    auto stage_v3sigma2tau    = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[158 * _ldStaging];
    auto stage_v3sigmalapl2   = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[170 * _ldStaging];
    auto stage_v3sigmalapltau = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[179 * _ldStaging];
    auto stage_v3sigmatau2    = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[191 * _ldStaging];
    auto stage_v3lapl3        = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[200 * _ldStaging];
    auto stage_v3lapl2tau     = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[204 * _ldStaging];
    auto stage_v3lapltau2     = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[210 * _ldStaging];
    auto stage_v3tau3         = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[216 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if ((family == XC_FAMILY_LDA) || (family == XC_FAMILY_HYB_LDA))
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
        else if ((family == XC_FAMILY_GGA) || (family == XC_FAMILY_HYB_GGA))
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
        else if ((family == XC_FAMILY_MGGA) || (family == XC_FAMILY_HYB_MGGA))
        {
            xc_mgga_kxc(funcptr, np, rho, sigma, lapl, tau,
                        stage_v3rho3, 
                        stage_v3rho2sigma, 
                        stage_v3rho2lapl, 
                        stage_v3rho2tau, 
                        stage_v3rhosigma2,
                        stage_v3rhosigmalapl, 
                        stage_v3rhosigmatau, 
                        stage_v3rholapl2, 
                        stage_v3rholapltau,
                        stage_v3rhotau2, 
                        stage_v3sigma3, 
                        stage_v3sigma2lapl, 
                        stage_v3sigma2tau,
                        stage_v3sigmalapl2, 
                        stage_v3sigmalapltau, 
                        stage_v3sigmatau2, 
                        stage_v3lapl3,
                        stage_v3lapl2tau, 
                        stage_v3lapltau2, 
                        stage_v3tau3);

#pragma omp simd aligned(v3rho3,stage_v3rho3, \
                         v3rho2sigma,stage_v3rho2sigma, \
                         v3rho2lapl,stage_v3rho2lapl, \
                         v3rho2tau,stage_v3rho2tau, \
                         v3rhosigma2,stage_v3rhosigma2, \
                         v3rhosigmalapl,stage_v3rhosigmalapl, \
                         v3rhosigmatau,stage_v3rhosigmatau, \
                         v3rholapl2,stage_v3rholapl2, \
                         v3rholapltau,stage_v3rholapltau, \
                         v3rhotau2,stage_v3rhotau2, \
                         v3sigma3,stage_v3sigma3, \
                         v3sigma2lapl,stage_v3sigma2lapl, \
                         v3sigma2tau,stage_v3sigma2tau, \
                         v3sigmalapl2,stage_v3sigmalapl2, \
                         v3sigmalapltau,stage_v3sigmalapltau, \
                         v3sigmatau2,stage_v3sigmatau2, \
                         v3lapl3,stage_v3lapl3, \
                         v3lapl2tau,stage_v3lapl2tau, \
                         v3lapltau2,stage_v3lapltau2, \
                         v3tau3,stage_v3tau3: VLX_ALIGN) 
            for (auto g = 0; g < np; ++g)
            {
                v3sigmalapltau[12 * g + 0] += c * stage_v3sigmalapltau[12 * g + 0];
                v3sigmalapltau[12 * g + 1] += c * stage_v3sigmalapltau[12 * g + 1];
                v3sigmalapltau[12 * g + 2] += c * stage_v3sigmalapltau[12 * g + 2];
                v3sigmalapltau[12 * g + 3] += c * stage_v3sigmalapltau[12 * g + 3];
                v3sigmalapltau[12 * g + 4] += c * stage_v3sigmalapltau[12 * g + 4];
                v3sigmalapltau[12 * g + 5] += c * stage_v3sigmalapltau[12 * g + 5];
                v3sigmalapltau[12 * g + 6] += c * stage_v3sigmalapltau[12 * g + 6];
                v3sigmalapltau[12 * g + 7] += c * stage_v3sigmalapltau[12 * g + 7];
                v3sigmalapltau[12 * g + 8] += c * stage_v3sigmalapltau[12 * g + 8];
                v3sigmalapltau[12 * g + 9] += c * stage_v3sigmalapltau[12 * g + 9];
                v3sigmalapltau[12 * g + 10] += c * stage_v3sigmalapltau[12 * g + 10];
                v3sigmalapltau[12 * g + 11] += c * stage_v3sigmalapltau[12 * g + 11];

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
                
                v3rhosigmalapl[12 * g + 0] += c * stage_v3rhosigmalapl[12 * g + 0];
                v3rhosigmalapl[12 * g + 1] += c * stage_v3rhosigmalapl[12 * g + 1];
                v3rhosigmalapl[12 * g + 2] += c * stage_v3rhosigmalapl[12 * g + 2];
                v3rhosigmalapl[12 * g + 3] += c * stage_v3rhosigmalapl[12 * g + 3];
                v3rhosigmalapl[12 * g + 4] += c * stage_v3rhosigmalapl[12 * g + 4];
                v3rhosigmalapl[12 * g + 5] += c * stage_v3rhosigmalapl[12 * g + 5];
                v3rhosigmalapl[12 * g + 6] += c * stage_v3rhosigmalapl[12 * g + 6];
                v3rhosigmalapl[12 * g + 7] += c * stage_v3rhosigmalapl[12 * g + 7];
                v3rhosigmalapl[12 * g + 8] += c * stage_v3rhosigmalapl[12 * g + 8];
                v3rhosigmalapl[12 * g + 9] += c * stage_v3rhosigmalapl[12 * g + 9];
                v3rhosigmalapl[12 * g + 10] += c * stage_v3rhosigmalapl[12 * g + 10];
                v3rhosigmalapl[12 * g + 11] += c * stage_v3rhosigmalapl[12 * g + 11];
                
                v3rhosigmatau[12 * g + 0] += c * stage_v3rhosigmatau[12 * g + 0];
                v3rhosigmatau[12 * g + 1] += c * stage_v3rhosigmatau[12 * g + 1];
                v3rhosigmatau[12 * g + 2] += c * stage_v3rhosigmatau[12 * g + 2];
                v3rhosigmatau[12 * g + 3] += c * stage_v3rhosigmatau[12 * g + 3];
                v3rhosigmatau[12 * g + 4] += c * stage_v3rhosigmatau[12 * g + 4];
                v3rhosigmatau[12 * g + 5] += c * stage_v3rhosigmatau[12 * g + 5];
                v3rhosigmatau[12 * g + 6] += c * stage_v3rhosigmatau[12 * g + 6];
                v3rhosigmatau[12 * g + 7] += c * stage_v3rhosigmatau[12 * g + 7];
                v3rhosigmatau[12 * g + 8] += c * stage_v3rhosigmatau[12 * g + 8];
                v3rhosigmatau[12 * g + 9] += c * stage_v3rhosigmatau[12 * g + 9];
                v3rhosigmatau[12 * g + 10] += c * stage_v3rhosigmatau[12 * g + 10];
                v3rhosigmatau[12 * g + 11] += c * stage_v3rhosigmatau[12 * g + 11];

                v3sigma2lapl[12 * g + 0] += c * stage_v3sigma2lapl[12 * g + 0];
                v3sigma2lapl[12 * g + 1] += c * stage_v3sigma2lapl[12 * g + 1];
                v3sigma2lapl[12 * g + 2] += c * stage_v3sigma2lapl[12 * g + 2];
                v3sigma2lapl[12 * g + 3] += c * stage_v3sigma2lapl[12 * g + 3];
                v3sigma2lapl[12 * g + 4] += c * stage_v3sigma2lapl[12 * g + 4];
                v3sigma2lapl[12 * g + 5] += c * stage_v3sigma2lapl[12 * g + 5];
                v3sigma2lapl[12 * g + 6] += c * stage_v3sigma2lapl[12 * g + 6];
                v3sigma2lapl[12 * g + 7] += c * stage_v3sigma2lapl[12 * g + 7];
                v3sigma2lapl[12 * g + 8] += c * stage_v3sigma2lapl[12 * g + 8];
                v3sigma2lapl[12 * g + 9] += c * stage_v3sigma2lapl[12 * g + 9];
                v3sigma2lapl[12 * g + 10] += c * stage_v3sigma2lapl[12 * g + 10];
                v3sigma2lapl[12 * g + 11] += c * stage_v3sigma2lapl[12 * g + 11];
                
                v3sigma2tau[12 * g + 0] += c * stage_v3sigma2tau[12 * g + 0];
                v3sigma2tau[12 * g + 1] += c * stage_v3sigma2tau[12 * g + 1];
                v3sigma2tau[12 * g + 2] += c * stage_v3sigma2tau[12 * g + 2];
                v3sigma2tau[12 * g + 3] += c * stage_v3sigma2tau[12 * g + 3];
                v3sigma2tau[12 * g + 4] += c * stage_v3sigma2tau[12 * g + 4];
                v3sigma2tau[12 * g + 5] += c * stage_v3sigma2tau[12 * g + 5];
                v3sigma2tau[12 * g + 6] += c * stage_v3sigma2tau[12 * g + 6];
                v3sigma2tau[12 * g + 7] += c * stage_v3sigma2tau[12 * g + 7];
                v3sigma2tau[12 * g + 8] += c * stage_v3sigma2tau[12 * g + 8];
                v3sigma2tau[12 * g + 9] += c * stage_v3sigma2tau[12 * g + 9];
                v3sigma2tau[12 * g + 10] += c * stage_v3sigma2tau[12 * g + 10];
                v3sigma2tau[12 * g + 11] += c * stage_v3sigma2tau[12 * g + 11];
  
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
                
                v3rho2lapl[6 * g + 0] += c * stage_v3rho2lapl[6 * g + 0];
                v3rho2lapl[6 * g + 1] += c * stage_v3rho2lapl[6 * g + 1];
                v3rho2lapl[6 * g + 2] += c * stage_v3rho2lapl[6 * g + 2];
                v3rho2lapl[6 * g + 3] += c * stage_v3rho2lapl[6 * g + 3];
                v3rho2lapl[6 * g + 4] += c * stage_v3rho2lapl[6 * g + 4];
                v3rho2lapl[6 * g + 5] += c * stage_v3rho2lapl[6 * g + 5];
                
                v3rho2tau[6 * g + 0] += c * stage_v3rho2tau[6 * g + 0];
                v3rho2tau[6 * g + 1] += c * stage_v3rho2tau[6 * g + 1];
                v3rho2tau[6 * g + 2] += c * stage_v3rho2tau[6 * g + 2];
                v3rho2tau[6 * g + 3] += c * stage_v3rho2tau[6 * g + 3];
                v3rho2tau[6 * g + 4] += c * stage_v3rho2tau[6 * g + 4];
                v3rho2tau[6 * g + 5] += c * stage_v3rho2tau[6 * g + 5];
                
                v3rholapl2[6 * g + 0] += c * stage_v3rholapl2[6 * g + 0];
                v3rholapl2[6 * g + 1] += c * stage_v3rholapl2[6 * g + 1];
                v3rholapl2[6 * g + 2] += c * stage_v3rholapl2[6 * g + 2];
                v3rholapl2[6 * g + 3] += c * stage_v3rholapl2[6 * g + 3];
                v3rholapl2[6 * g + 4] += c * stage_v3rholapl2[6 * g + 4];
                v3rholapl2[6 * g + 5] += c * stage_v3rholapl2[6 * g + 5];
                
                v3rholapltau[8 * g + 0] += c * stage_v3rholapltau[8 * g + 0];
                v3rholapltau[8 * g + 1] += c * stage_v3rholapltau[8 * g + 1];
                v3rholapltau[8 * g + 2] += c * stage_v3rholapltau[8 * g + 2];
                v3rholapltau[8 * g + 3] += c * stage_v3rholapltau[8 * g + 3];
                v3rholapltau[8 * g + 4] += c * stage_v3rholapltau[8 * g + 4];
                v3rholapltau[8 * g + 5] += c * stage_v3rholapltau[8 * g + 5];
                v3rholapltau[8 * g + 6] += c * stage_v3rholapltau[8 * g + 6];
                v3rholapltau[8 * g + 7] += c * stage_v3rholapltau[8 * g + 7];
                
                v3rhotau2[6 * g + 0] += c * stage_v3rhotau2[6 * g + 0];
                v3rhotau2[6 * g + 1] += c * stage_v3rhotau2[6 * g + 1];
                v3rhotau2[6 * g + 2] += c * stage_v3rhotau2[6 * g + 2];
                v3rhotau2[6 * g + 3] += c * stage_v3rhotau2[6 * g + 3];
                v3rhotau2[6 * g + 4] += c * stage_v3rhotau2[6 * g + 4];
                v3rhotau2[6 * g + 5] += c * stage_v3rhotau2[6 * g + 5];
                
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

                v3sigmalapl2[9 * g + 0] += c * stage_v3sigmalapl2[9 * g + 0];
                v3sigmalapl2[9 * g + 1] += c * stage_v3sigmalapl2[9 * g + 1];
                v3sigmalapl2[9 * g + 2] += c * stage_v3sigmalapl2[9 * g + 2];
                v3sigmalapl2[9 * g + 3] += c * stage_v3sigmalapl2[9 * g + 3];
                v3sigmalapl2[9 * g + 4] += c * stage_v3sigmalapl2[9 * g + 4];
                v3sigmalapl2[9 * g + 5] += c * stage_v3sigmalapl2[9 * g + 5];
                v3sigmalapl2[9 * g + 6] += c * stage_v3sigmalapl2[9 * g + 6];
                v3sigmalapl2[9 * g + 7] += c * stage_v3sigmalapl2[9 * g + 7];
                v3sigmalapl2[9 * g + 8] += c * stage_v3sigmalapl2[9 * g + 8];
                
                v3sigmatau2[9 * g + 0] += c * stage_v3sigmatau2[9 * g + 0];
                v3sigmatau2[9 * g + 1] += c * stage_v3sigmatau2[9 * g + 1];
                v3sigmatau2[9 * g + 2] += c * stage_v3sigmatau2[9 * g + 2];
                v3sigmatau2[9 * g + 3] += c * stage_v3sigmatau2[9 * g + 3];
                v3sigmatau2[9 * g + 4] += c * stage_v3sigmatau2[9 * g + 4];
                v3sigmatau2[9 * g + 5] += c * stage_v3sigmatau2[9 * g + 5];
                v3sigmatau2[9 * g + 6] += c * stage_v3sigmatau2[9 * g + 6];
                v3sigmatau2[9 * g + 7] += c * stage_v3sigmatau2[9 * g + 7];
                v3sigmatau2[9 * g + 8] += c * stage_v3sigmatau2[9 * g + 8];
                 
                v3lapl3[4 * g + 0] += c * stage_v3lapl3[4 * g + 0];
                v3lapl3[4 * g + 1] += c * stage_v3lapl3[4 * g + 1];
                v3lapl3[4 * g + 2] += c * stage_v3lapl3[4 * g + 2];
                v3lapl3[4 * g + 3] += c * stage_v3lapl3[4 * g + 3];
                
                v3lapl2tau[6 * g + 0] += c * stage_v3lapl2tau[6 * g + 0];
                v3lapl2tau[6 * g + 1] += c * stage_v3lapl2tau[6 * g + 1];
                v3lapl2tau[6 * g + 2] += c * stage_v3lapl2tau[6 * g + 2];
                v3lapl2tau[6 * g + 3] += c * stage_v3lapl2tau[6 * g + 3];
                v3lapl2tau[6 * g + 4] += c * stage_v3lapl2tau[6 * g + 4];
                v3lapl2tau[6 * g + 5] += c * stage_v3lapl2tau[6 * g + 5];
                
                v3lapltau2[6 * g + 0] += c * stage_v3lapltau2[6 * g + 0];
                v3lapltau2[6 * g + 1] += c * stage_v3lapltau2[6 * g + 1];
                v3lapltau2[6 * g + 2] += c * stage_v3lapltau2[6 * g + 2];
                v3lapltau2[6 * g + 3] += c * stage_v3lapltau2[6 * g + 3];
                v3lapltau2[6 * g + 4] += c * stage_v3lapltau2[6 * g + 4];
                v3lapltau2[6 * g + 5] += c * stage_v3lapltau2[6 * g + 5];
                
                v3tau3[4 * g + 0] += c * stage_v3tau3[4 * g + 0];
                v3tau3[4 * g + 1] += c * stage_v3tau3[4 * g + 1];
                v3tau3[4 * g + 2] += c * stage_v3tau3[4 * g + 2];
                v3tau3[4 * g + 3] += c * stage_v3tau3[4 * g + 3];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v3rho3);
        mem::free(stage_v3rho2sigma);
        mem::free(stage_v3rho2lapl);
        mem::free(stage_v3rho2tau);
        mem::free(stage_v3rhosigma2);
        mem::free(stage_v3rhosigmalapl);
        mem::free(stage_v3rhosigmatau);
        mem::free(stage_v3rholapl2);
        mem::free(stage_v3rholapltau);
        mem::free(stage_v3rhotau2);
        mem::free(stage_v3sigma3);
        mem::free(stage_v3sigma2lapl);
        mem::free(stage_v3sigma2tau);
        mem::free(stage_v3sigmalapl2);
        mem::free(stage_v3sigmalapltau);
        mem::free(stage_v3sigmatau2);
        mem::free(stage_v3lapl3);
        mem::free(stage_v3lapl2tau);
        mem::free(stage_v3lapltau2);
        mem::free(stage_v3tau3);
    }
}

auto
CXCNewFunctional::compute_lxc_for_mgga(int32_t       np,
                                       const double* rho,
                                       const double* sigma,
                                       const double* lapl,
                                       const double* tau,
                                       double*      v4rho4,
                                       double*      v4rho3sigma,
                                       double*      v4rho3lapl,
                                       double*      v4rho3tau,
                                       double*      v4rho2sigma2,
                                       double*      v4rho2sigmalapl,
                                       double*      v4rho2sigmatau,
                                       double*      v4rho2lapl2,
                                       double*      v4rho2lapltau,
                                       double*      v4rho2tau2,
                                       double*      v4rhosigma3,
                                       double*      v4rhosigma2lapl,
                                       double*      v4rhosigma2tau,
                                       double*      v4rhosigmalapl2,
                                       double*      v4rhosigmalapltau,
                                       double*      v4rhosigmatau2,
                                       double*      v4rholapl3,
                                       double*      v4rholapl2tau,
                                       double*      v4rholapltau2,
                                       double*      v4rhotau3,
                                       double*      v4sigma4,
                                       double*      v4sigma3lapl,
                                       double*      v4sigma3tau,
                                       double*      v4sigma2lapl2,
                                       double*      v4sigma2lapltau,
                                       double*      v4sigma2tau2,
                                       double*      v4sigmalapl3,
                                       double*      v4sigmalapl2tau,
                                       double*      v4sigmalapltau2,
                                       double*      v4sigmatau3,
                                       double*      v4lapl4,
                                       double*      v4lapl3tau,
                                       double*      v4lapl2tau2,
                                       double*      v4lapltau3,
                                       double*       v4tau4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

#pragma omp simd aligned(v4rho4,v4rho3sigma,v4rho3lapl,v4rho3tau,v4rho2sigma2,v4rho2sigmalapl,\
                         v4rho2sigmatau,v4rho2lapl2,v4rho2lapltau,v4rho2tau2,v4rhosigma3,v4rhosigma2lapl,\
                         v4rhosigma2tau,v4rhosigmalapl2,v4rhosigmalapltau,v4rhosigmatau2,v4rholapl3,v4rholapl2tau,\
                         v4rholapltau2,v4rhotau3,v4sigma4,v4sigma3lapl,v4sigma3tau,v4sigma2lapl2,v4sigma2lapltau,v4sigma2tau2,\
                         v4sigmalapl3,v4sigmalapl2tau,v4sigmalapltau2,v4sigmatau3,v4lapl4,v4lapl3tau,v4lapl2tau2,v4lapltau3, v4tau4 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v4rho4[5 * g + 0] = 0;
        v4rho4[5 * g + 1] = 0;
        v4rho4[5 * g + 2] = 0;
        v4rho4[5 * g + 3] = 0;
        v4rho4[5 * g + 4] = 0;

        v4rho3sigma[12 * g + 0] = 0;
        v4rho3sigma[12 * g + 1] = 0;
        v4rho3sigma[12 * g + 2] = 0;
        v4rho3sigma[12 * g + 3] = 0;
        v4rho3sigma[12 * g + 4] = 0;
        v4rho3sigma[12 * g + 5] = 0;
        v4rho3sigma[12 * g + 6] = 0;
        v4rho3sigma[12 * g + 7] = 0;
        v4rho3sigma[12 * g + 8] = 0;
        v4rho3sigma[12 * g + 9] = 0;
        v4rho3sigma[12 * g + 10] = 0;
        v4rho3sigma[12 * g + 11] = 0;

        v4rho3lapl[8 * g + 0] = 0;
        v4rho3lapl[8 * g + 1] = 0;
        v4rho3lapl[8 * g + 2] = 0;
        v4rho3lapl[8 * g + 3] = 0;
        v4rho3lapl[8 * g + 4] = 0;
        v4rho3lapl[8 * g + 5] = 0;
        v4rho3lapl[8 * g + 6] = 0;
        v4rho3lapl[8 * g + 7] = 0;

        v4rho3tau[8 * g + 0] = 0;
        v4rho3tau[8 * g + 1] = 0;
        v4rho3tau[8 * g + 2] = 0;
        v4rho3tau[8 * g + 3] = 0;
        v4rho3tau[8 * g + 4] = 0;
        v4rho3tau[8 * g + 5] = 0;
        v4rho3tau[8 * g + 6] = 0;
        v4rho3tau[8 * g + 7] = 0;

        v4rho2sigma2[18 * g + 0] = 0;
        v4rho2sigma2[18 * g + 1] = 0;
        v4rho2sigma2[18 * g + 2] = 0;
        v4rho2sigma2[18 * g + 3] = 0;
        v4rho2sigma2[18 * g + 4] = 0;
        v4rho2sigma2[18 * g + 5] = 0;
        v4rho2sigma2[18 * g + 6] = 0;
        v4rho2sigma2[18 * g + 7] = 0;
        v4rho2sigma2[18 * g + 8] = 0;
        v4rho2sigma2[18 * g + 9] = 0;
        v4rho2sigma2[18 * g + 10] = 0;
        v4rho2sigma2[18 * g + 11] = 0;
        v4rho2sigma2[18 * g + 12] = 0;
        v4rho2sigma2[18 * g + 13] = 0;
        v4rho2sigma2[18 * g + 14] = 0;
        v4rho2sigma2[18 * g + 15] = 0;
        v4rho2sigma2[18 * g + 16] = 0;
        v4rho2sigma2[18 * g + 17] = 0;

        v4rho2sigmalapl[18 * g + 0] = 0;
        v4rho2sigmalapl[18 * g + 1] = 0;
        v4rho2sigmalapl[18 * g + 2] = 0;
        v4rho2sigmalapl[18 * g + 3] = 0;
        v4rho2sigmalapl[18 * g + 4] = 0;
        v4rho2sigmalapl[18 * g + 5] = 0;
        v4rho2sigmalapl[18 * g + 6] = 0;
        v4rho2sigmalapl[18 * g + 7] = 0;
        v4rho2sigmalapl[18 * g + 8] = 0;
        v4rho2sigmalapl[18 * g + 9] = 0;
        v4rho2sigmalapl[18 * g + 10] = 0;
        v4rho2sigmalapl[18 * g + 11] = 0;
        v4rho2sigmalapl[18 * g + 12] = 0;
        v4rho2sigmalapl[18 * g + 13] = 0;
        v4rho2sigmalapl[18 * g + 14] = 0;
        v4rho2sigmalapl[18 * g + 15] = 0;
        v4rho2sigmalapl[18 * g + 16] = 0;
        v4rho2sigmalapl[18 * g + 17] = 0;

        v4rho2sigmatau[18 * g + 0] = 0;
        v4rho2sigmatau[18 * g + 1] = 0;
        v4rho2sigmatau[18 * g + 2] = 0;
        v4rho2sigmatau[18 * g + 3] = 0;
        v4rho2sigmatau[18 * g + 4] = 0;
        v4rho2sigmatau[18 * g + 5] = 0;
        v4rho2sigmatau[18 * g + 6] = 0;
        v4rho2sigmatau[18 * g + 7] = 0;
        v4rho2sigmatau[18 * g + 8] = 0;
        v4rho2sigmatau[18 * g + 9] = 0;
        v4rho2sigmatau[18 * g + 10] = 0;
        v4rho2sigmatau[18 * g + 11] = 0;
        v4rho2sigmatau[18 * g + 12] = 0;
        v4rho2sigmatau[18 * g + 13] = 0;
        v4rho2sigmatau[18 * g + 14] = 0;
        v4rho2sigmatau[18 * g + 15] = 0;
        v4rho2sigmatau[18 * g + 16] = 0;
        v4rho2sigmatau[18 * g + 17] = 0;

        v4rho2lapl2[9 * g + 0] = 0;
        v4rho2lapl2[9 * g + 1] = 0;
        v4rho2lapl2[9 * g + 2] = 0;
        v4rho2lapl2[9 * g + 3] = 0;
        v4rho2lapl2[9 * g + 4] = 0;
        v4rho2lapl2[9 * g + 5] = 0;
        v4rho2lapl2[9 * g + 6] = 0;
        v4rho2lapl2[9 * g + 7] = 0;
        v4rho2lapl2[9 * g + 8] = 0;

        v4rho2lapltau[12 * g + 0] = 0;
        v4rho2lapltau[12 * g + 1] = 0;
        v4rho2lapltau[12 * g + 2] = 0;
        v4rho2lapltau[12 * g + 3] = 0;
        v4rho2lapltau[12 * g + 4] = 0;
        v4rho2lapltau[12 * g + 5] = 0;
        v4rho2lapltau[12 * g + 6] = 0;
        v4rho2lapltau[12 * g + 7] = 0;
        v4rho2lapltau[12 * g + 8] = 0;
        v4rho2lapltau[12 * g + 9] = 0;
        v4rho2lapltau[12 * g + 10] = 0;
        v4rho2lapltau[12 * g + 11] = 0;

        v4rho2tau2[9 * g + 0] = 0;
        v4rho2tau2[9 * g + 1] = 0;
        v4rho2tau2[9 * g + 2] = 0;
        v4rho2tau2[9 * g + 3] = 0;
        v4rho2tau2[9 * g + 4] = 0;
        v4rho2tau2[9 * g + 5] = 0;
        v4rho2tau2[9 * g + 6] = 0;
        v4rho2tau2[9 * g + 7] = 0;
        v4rho2tau2[9 * g + 8] = 0;

        v4rhosigma3[20 * g + 0] = 0;
        v4rhosigma3[20 * g + 1] = 0;
        v4rhosigma3[20 * g + 2] = 0;
        v4rhosigma3[20 * g + 3] = 0;
        v4rhosigma3[20 * g + 4] = 0;
        v4rhosigma3[20 * g + 5] = 0;
        v4rhosigma3[20 * g + 6] = 0;
        v4rhosigma3[20 * g + 7] = 0;
        v4rhosigma3[20 * g + 8] = 0;
        v4rhosigma3[20 * g + 9] = 0;
        v4rhosigma3[20 * g + 10] = 0;
        v4rhosigma3[20 * g + 11] = 0;
        v4rhosigma3[20 * g + 12] = 0;
        v4rhosigma3[20 * g + 13] = 0;
        v4rhosigma3[20 * g + 14] = 0;
        v4rhosigma3[20 * g + 15] = 0;
        v4rhosigma3[20 * g + 16] = 0;
        v4rhosigma3[20 * g + 17] = 0;
        v4rhosigma3[20 * g + 18] = 0;
        v4rhosigma3[20 * g + 19] = 0;

        v4rhosigma2lapl[36 * g + 0] = 0;
        v4rhosigma2lapl[36 * g + 1] = 0;
        v4rhosigma2lapl[36 * g + 2] = 0;
        v4rhosigma2lapl[36 * g + 3] = 0;
        v4rhosigma2lapl[36 * g + 4] = 0;
        v4rhosigma2lapl[36 * g + 5] = 0;
        v4rhosigma2lapl[36 * g + 6] = 0;
        v4rhosigma2lapl[36 * g + 7] = 0;
        v4rhosigma2lapl[36 * g + 8] = 0;
        v4rhosigma2lapl[36 * g + 9] = 0;
        v4rhosigma2lapl[36 * g + 10] = 0;
        v4rhosigma2lapl[36 * g + 11] = 0;
        v4rhosigma2lapl[36 * g + 12] = 0;
        v4rhosigma2lapl[36 * g + 13] = 0;
        v4rhosigma2lapl[36 * g + 14] = 0;
        v4rhosigma2lapl[36 * g + 15] = 0;
        v4rhosigma2lapl[36 * g + 16] = 0;
        v4rhosigma2lapl[36 * g + 17] = 0;
        v4rhosigma2lapl[36 * g + 18] = 0;
        v4rhosigma2lapl[36 * g + 19] = 0;
        v4rhosigma2lapl[36 * g + 20] = 0;
        v4rhosigma2lapl[36 * g + 21] = 0;
        v4rhosigma2lapl[36 * g + 22] = 0;
        v4rhosigma2lapl[36 * g + 23] = 0;

        v4rhosigma2tau[36 * g + 0] = 0;
        v4rhosigma2tau[36 * g + 1] = 0;
        v4rhosigma2tau[36 * g + 2] = 0;
        v4rhosigma2tau[36 * g + 3] = 0;
        v4rhosigma2tau[36 * g + 4] = 0;
        v4rhosigma2tau[36 * g + 5] = 0;
        v4rhosigma2tau[36 * g + 6] = 0;
        v4rhosigma2tau[36 * g + 7] = 0;
        v4rhosigma2tau[36 * g + 8] = 0;
        v4rhosigma2tau[36 * g + 9] = 0;
        v4rhosigma2tau[36 * g + 10] = 0;
        v4rhosigma2tau[36 * g + 11] = 0;
        v4rhosigma2tau[36 * g + 12] = 0;
        v4rhosigma2tau[36 * g + 13] = 0;
        v4rhosigma2tau[36 * g + 14] = 0;
        v4rhosigma2tau[36 * g + 15] = 0;
        v4rhosigma2tau[36 * g + 16] = 0;
        v4rhosigma2tau[36 * g + 17] = 0;
        v4rhosigma2tau[36 * g + 18] = 0;
        v4rhosigma2tau[36 * g + 19] = 0;
        v4rhosigma2tau[36 * g + 20] = 0;
        v4rhosigma2tau[36 * g + 21] = 0;
        v4rhosigma2tau[36 * g + 22] = 0;
        v4rhosigma2tau[36 * g + 23] = 0;

        v4rhosigmalapl2[18 * g + 0] = 0;
        v4rhosigmalapl2[18 * g + 1] = 0;
        v4rhosigmalapl2[18 * g + 2] = 0;
        v4rhosigmalapl2[18 * g + 3] = 0;
        v4rhosigmalapl2[18 * g + 4] = 0;
        v4rhosigmalapl2[18 * g + 5] = 0;
        v4rhosigmalapl2[18 * g + 6] = 0;
        v4rhosigmalapl2[18 * g + 7] = 0;
        v4rhosigmalapl2[18 * g + 8] = 0;
        v4rhosigmalapl2[18 * g + 9] = 0;
        v4rhosigmalapl2[18 * g + 10] = 0;
        v4rhosigmalapl2[18 * g + 11] = 0;
        v4rhosigmalapl2[18 * g + 12] = 0;
        v4rhosigmalapl2[18 * g + 13] = 0;
        v4rhosigmalapl2[18 * g + 14] = 0;
        v4rhosigmalapl2[18 * g + 15] = 0;
        v4rhosigmalapl2[18 * g + 16] = 0;
        v4rhosigmalapl2[18 * g + 17] = 0;

        v4rhosigmalapltau[24 * g + 0] = 0;
        v4rhosigmalapltau[24 * g + 1] = 0;
        v4rhosigmalapltau[24 * g + 2] = 0;
        v4rhosigmalapltau[24 * g + 3] = 0;
        v4rhosigmalapltau[24 * g + 4] = 0;
        v4rhosigmalapltau[24 * g + 5] = 0;
        v4rhosigmalapltau[24 * g + 6] = 0;
        v4rhosigmalapltau[24 * g + 7] = 0;
        v4rhosigmalapltau[24 * g + 8] = 0;
        v4rhosigmalapltau[24 * g + 9] = 0;
        v4rhosigmalapltau[24 * g + 10] = 0;
        v4rhosigmalapltau[24 * g + 11] = 0;
        v4rhosigmalapltau[24 * g + 12] = 0;
        v4rhosigmalapltau[24 * g + 13] = 0;
        v4rhosigmalapltau[24 * g + 14] = 0;
        v4rhosigmalapltau[24 * g + 15] = 0;
        v4rhosigmalapltau[24 * g + 16] = 0;
        v4rhosigmalapltau[24 * g + 17] = 0;
        v4rhosigmalapltau[24 * g + 18] = 0;
        v4rhosigmalapltau[24 * g + 19] = 0;
        v4rhosigmalapltau[24 * g + 20] = 0;
        v4rhosigmalapltau[24 * g + 21] = 0;
        v4rhosigmalapltau[24 * g + 22] = 0;
        v4rhosigmalapltau[24 * g + 23] = 0;

        v4rhosigmatau2[36 * g + 0] = 0;
        v4rhosigmatau2[36 * g + 1] = 0;
        v4rhosigmatau2[36 * g + 2] = 0;
        v4rhosigmatau2[36 * g + 3] = 0;
        v4rhosigmatau2[36 * g + 4] = 0;
        v4rhosigmatau2[36 * g + 5] = 0;
        v4rhosigmatau2[36 * g + 6] = 0;
        v4rhosigmatau2[36 * g + 7] = 0;
        v4rhosigmatau2[36 * g + 8] = 0;
        v4rhosigmatau2[36 * g + 9] = 0;
        v4rhosigmatau2[36 * g + 10] = 0;
        v4rhosigmatau2[36 * g + 11] = 0;
        v4rhosigmatau2[36 * g + 12] = 0;
        v4rhosigmatau2[36 * g + 13] = 0;
        v4rhosigmatau2[36 * g + 14] = 0;
        v4rhosigmatau2[36 * g + 15] = 0;
        v4rhosigmatau2[36 * g + 16] = 0;
        v4rhosigmatau2[36 * g + 17] = 0;

        v4rholapl3[8 * g + 0] = 0;
        v4rholapl3[8 * g + 1] = 0;
        v4rholapl3[8 * g + 2] = 0;
        v4rholapl3[8 * g + 3] = 0;
        v4rholapl3[8 * g + 4] = 0;
        v4rholapl3[8 * g + 5] = 0;
        v4rholapl3[8 * g + 6] = 0;
        v4rholapl3[8 * g + 7] = 0;

        v4rholapl2tau[12 * g + 0] = 0;
        v4rholapl2tau[12 * g + 1] = 0;
        v4rholapl2tau[12 * g + 2] = 0;
        v4rholapl2tau[12 * g + 3] = 0;
        v4rholapl2tau[12 * g + 4] = 0;
        v4rholapl2tau[12 * g + 5] = 0;
        v4rholapl2tau[12 * g + 6] = 0;
        v4rholapl2tau[12 * g + 7] = 0;
        v4rholapl2tau[12 * g + 8] = 0;
        v4rholapl2tau[12 * g + 9] = 0;
        v4rholapl2tau[12 * g + 10] = 0;
        v4rholapl2tau[12 * g + 11] = 0;

        v4rholapltau2[12 * g + 0] = 0;
        v4rholapltau2[12 * g + 1] = 0;
        v4rholapltau2[12 * g + 2] = 0;
        v4rholapltau2[12 * g + 3] = 0;
        v4rholapltau2[12 * g + 4] = 0;
        v4rholapltau2[12 * g + 5] = 0;
        v4rholapltau2[12 * g + 6] = 0;
        v4rholapltau2[12 * g + 7] = 0;
        v4rholapltau2[12 * g + 8] = 0;
        v4rholapltau2[12 * g + 9] = 0;
        v4rholapltau2[12 * g + 10] = 0;
        v4rholapltau2[12 * g + 11] = 0;

        v4rhotau3[8 * g + 0] = 0;
        v4rhotau3[8 * g + 1] = 0;
        v4rhotau3[8 * g + 2] = 0;
        v4rhotau3[8 * g + 3] = 0;
        v4rhotau3[8 * g + 4] = 0;
        v4rhotau3[8 * g + 5] = 0;
        v4rhotau3[8 * g + 6] = 0;
        v4rhotau3[8 * g + 7] = 0;

        v4sigma4[15 * g + 0] = 0;
        v4sigma4[15 * g + 1] = 0;
        v4sigma4[15 * g + 2] = 0;
        v4sigma4[15 * g + 3] = 0;
        v4sigma4[15 * g + 4] = 0;
        v4sigma4[15 * g + 5] = 0;
        v4sigma4[15 * g + 6] = 0;
        v4sigma4[15 * g + 7] = 0;
        v4sigma4[15 * g + 8] = 0;
        v4sigma4[15 * g + 9] = 0;
        v4sigma4[15 * g + 10] = 0;
        v4sigma4[15 * g + 11] = 0;
        v4sigma4[15 * g + 12] = 0;
        v4sigma4[15 * g + 13] = 0;
        v4sigma4[15 * g + 14] = 0;

        v4sigma3lapl[20 * g + 0] = 0;
        v4sigma3lapl[20 * g + 1] = 0;
        v4sigma3lapl[20 * g + 2] = 0;
        v4sigma3lapl[20 * g + 3] = 0;
        v4sigma3lapl[20 * g + 4] = 0;
        v4sigma3lapl[20 * g + 5] = 0;
        v4sigma3lapl[20 * g + 6] = 0;
        v4sigma3lapl[20 * g + 7] = 0;
        v4sigma3lapl[20 * g + 8] = 0;
        v4sigma3lapl[20 * g + 9] = 0;
        v4sigma3lapl[20 * g + 10] = 0;
        v4sigma3lapl[20 * g + 11] = 0;
        v4sigma3lapl[20 * g + 12] = 0;
        v4sigma3lapl[20 * g + 13] = 0;
        v4sigma3lapl[20 * g + 14] = 0;
        v4sigma3lapl[20 * g + 15] = 0;
        v4sigma3lapl[20 * g + 16] = 0;
        v4sigma3lapl[20 * g + 17] = 0;
        v4sigma3lapl[20 * g + 18] = 0;
        v4sigma3lapl[20 * g + 19] = 0;

        v4sigma3tau[30 * g + 0] = 0;
        v4sigma3tau[30 * g + 1] = 0;
        v4sigma3tau[30 * g + 2] = 0;
        v4sigma3tau[30 * g + 3] = 0;
        v4sigma3tau[30 * g + 4] = 0;
        v4sigma3tau[30 * g + 5] = 0;
        v4sigma3tau[30 * g + 6] = 0;
        v4sigma3tau[30 * g + 7] = 0;
        v4sigma3tau[30 * g + 8] = 0;
        v4sigma3tau[30 * g + 9] = 0;
        v4sigma3tau[30 * g + 10] = 0;
        v4sigma3tau[30 * g + 11] = 0;
        v4sigma3tau[30 * g + 12] = 0;
        v4sigma3tau[30 * g + 13] = 0;
        v4sigma3tau[30 * g + 14] = 0;
        v4sigma3tau[30 * g + 15] = 0;
        v4sigma3tau[30 * g + 16] = 0;
        v4sigma3tau[30 * g + 17] = 0;
        v4sigma3tau[30 * g + 18] = 0;
        v4sigma3tau[30 * g + 19] = 0;

        v4sigma2lapl2[18 * g + 0] = 0;
        v4sigma2lapl2[18 * g + 1] = 0;
        v4sigma2lapl2[18 * g + 2] = 0;
        v4sigma2lapl2[18 * g + 3] = 0;
        v4sigma2lapl2[18 * g + 4] = 0;
        v4sigma2lapl2[18 * g + 5] = 0;
        v4sigma2lapl2[18 * g + 6] = 0;
        v4sigma2lapl2[18 * g + 7] = 0;
        v4sigma2lapl2[18 * g + 8] = 0;
        v4sigma2lapl2[18 * g + 9] = 0;
        v4sigma2lapl2[18 * g + 10] = 0;
        v4sigma2lapl2[18 * g + 11] = 0;
        v4sigma2lapl2[18 * g + 12] = 0;
        v4sigma2lapl2[18 * g + 13] = 0;
        v4sigma2lapl2[18 * g + 14] = 0;
        v4sigma2lapl2[18 * g + 15] = 0;
        v4sigma2lapl2[18 * g + 16] = 0;
        v4sigma2lapl2[18 * g + 17] = 0;

        v4sigma2lapltau[24 * g + 0] = 0;
        v4sigma2lapltau[24 * g + 1] = 0;
        v4sigma2lapltau[24 * g + 2] = 0;
        v4sigma2lapltau[24 * g + 3] = 0;
        v4sigma2lapltau[24 * g + 4] = 0;
        v4sigma2lapltau[24 * g + 5] = 0;
        v4sigma2lapltau[24 * g + 6] = 0;
        v4sigma2lapltau[24 * g + 7] = 0;
        v4sigma2lapltau[24 * g + 8] = 0;
        v4sigma2lapltau[24 * g + 9] = 0;
        v4sigma2lapltau[24 * g + 10] = 0;
        v4sigma2lapltau[24 * g + 11] = 0;
        v4sigma2lapltau[24 * g + 12] = 0;
        v4sigma2lapltau[24 * g + 13] = 0;
        v4sigma2lapltau[24 * g + 14] = 0;
        v4sigma2lapltau[24 * g + 15] = 0;
        v4sigma2lapltau[24 * g + 16] = 0;
        v4sigma2lapltau[24 * g + 17] = 0;
        v4sigma2lapltau[24 * g + 18] = 0;
        v4sigma2lapltau[24 * g + 19] = 0;
        v4sigma2lapltau[24 * g + 20] = 0;
        v4sigma2lapltau[24 * g + 21] = 0;
        v4sigma2lapltau[24 * g + 22] = 0;
        v4sigma2lapltau[24 * g + 23] = 0;

        v4sigma2tau2[18 * g + 0] = 0;
        v4sigma2tau2[18 * g + 1] = 0;
        v4sigma2tau2[18 * g + 2] = 0;
        v4sigma2tau2[18 * g + 3] = 0;
        v4sigma2tau2[18 * g + 4] = 0;
        v4sigma2tau2[18 * g + 5] = 0;
        v4sigma2tau2[18 * g + 6] = 0;
        v4sigma2tau2[18 * g + 7] = 0;
        v4sigma2tau2[18 * g + 8] = 0;
        v4sigma2tau2[18 * g + 9] = 0;
        v4sigma2tau2[18 * g + 10] = 0;
        v4sigma2tau2[18 * g + 11] = 0;
        v4sigma2tau2[18 * g + 12] = 0;
        v4sigma2tau2[18 * g + 13] = 0;
        v4sigma2tau2[18 * g + 14] = 0;
        v4sigma2tau2[18 * g + 15] = 0;
        v4sigma2tau2[18 * g + 16] = 0;
        v4sigma2tau2[18 * g + 17] = 0;

        v4sigmalapl3[12 * g + 0] = 0;
        v4sigmalapl3[12 * g + 1] = 0;
        v4sigmalapl3[12 * g + 2] = 0;
        v4sigmalapl3[12 * g + 3] = 0;
        v4sigmalapl3[12 * g + 4] = 0;
        v4sigmalapl3[12 * g + 5] = 0;
        v4sigmalapl3[12 * g + 6] = 0;
        v4sigmalapl3[12 * g + 7] = 0;
        v4sigmalapl3[12 * g + 8] = 0;
        v4sigmalapl3[12 * g + 9] = 0;
        v4sigmalapl3[12 * g + 10] = 0;
        v4sigmalapl3[12 * g + 11] = 0;

        v4sigmalapl2tau[18 * g + 0] = 0;
        v4sigmalapl2tau[18 * g + 1] = 0;
        v4sigmalapl2tau[18 * g + 2] = 0;
        v4sigmalapl2tau[18 * g + 3] = 0;
        v4sigmalapl2tau[18 * g + 4] = 0;
        v4sigmalapl2tau[18 * g + 5] = 0;
        v4sigmalapl2tau[18 * g + 6] = 0;
        v4sigmalapl2tau[18 * g + 7] = 0;
        v4sigmalapl2tau[18 * g + 8] = 0;
        v4sigmalapl2tau[18 * g + 9] = 0;
        v4sigmalapl2tau[18 * g + 10] = 0;
        v4sigmalapl2tau[18 * g + 11] = 0;
        v4sigmalapl2tau[18 * g + 12] = 0;
        v4sigmalapl2tau[18 * g + 13] = 0;
        v4sigmalapl2tau[18 * g + 14] = 0;
        v4sigmalapl2tau[18 * g + 15] = 0;
        v4sigmalapl2tau[18 * g + 16] = 0;
        v4sigmalapl2tau[18 * g + 17] = 0;

        v4sigmalapltau2[18 * g + 0] = 0;
        v4sigmalapltau2[18 * g + 1] = 0;
        v4sigmalapltau2[18 * g + 2] = 0;
        v4sigmalapltau2[18 * g + 3] = 0;
        v4sigmalapltau2[18 * g + 4] = 0;
        v4sigmalapltau2[18 * g + 5] = 0;
        v4sigmalapltau2[18 * g + 6] = 0;
        v4sigmalapltau2[18 * g + 7] = 0;
        v4sigmalapltau2[18 * g + 8] = 0;
        v4sigmalapltau2[18 * g + 9] = 0;
        v4sigmalapltau2[18 * g + 10] = 0;
        v4sigmalapltau2[18 * g + 11] = 0;
        v4sigmalapltau2[18 * g + 12] = 0;
        v4sigmalapltau2[18 * g + 13] = 0;
        v4sigmalapltau2[18 * g + 14] = 0;
        v4sigmalapltau2[18 * g + 15] = 0;
        v4sigmalapltau2[18 * g + 16] = 0;
        v4sigmalapltau2[18 * g + 17] = 0;

        v4sigmatau3[12 * g + 0] = 0;
        v4sigmatau3[12 * g + 1] = 0;
        v4sigmatau3[12 * g + 2] = 0;
        v4sigmatau3[12 * g + 3] = 0;
        v4sigmatau3[12 * g + 4] = 0;
        v4sigmatau3[12 * g + 5] = 0;
        v4sigmatau3[12 * g + 6] = 0;
        v4sigmatau3[12 * g + 7] = 0;
        v4sigmatau3[12 * g + 8] = 0;
        v4sigmatau3[12 * g + 9] = 0;
        v4sigmatau3[12 * g + 10] = 0;
        v4sigmatau3[12 * g + 11] = 0;

        v4lapl4[5 * g + 0] = 0;
        v4lapl4[5 * g + 1] = 0;
        v4lapl4[5 * g + 2] = 0;
        v4lapl4[5 * g + 3] = 0;
        v4lapl4[5 * g + 4] = 0;

        v4lapl3tau[8 * g + 0] = 0;
        v4lapl3tau[8 * g + 1] = 0;
        v4lapl3tau[8 * g + 2] = 0;
        v4lapl3tau[8 * g + 3] = 0;
        v4lapl3tau[8 * g + 4] = 0;
        v4lapl3tau[8 * g + 5] = 0;
        v4lapl3tau[8 * g + 6] = 0;
        v4lapl3tau[8 * g + 7] = 0;

        v4lapl2tau2[9 * g + 0] = 0;
        v4lapl2tau2[9 * g + 1] = 0;
        v4lapl2tau2[9 * g + 2] = 0;
        v4lapl2tau2[9 * g + 3] = 0;
        v4lapl2tau2[9 * g + 4] = 0;
        v4lapl2tau2[9 * g + 5] = 0;
        v4lapl2tau2[9 * g + 6] = 0;
        v4lapl2tau2[9 * g + 7] = 0;
        v4lapl2tau2[9 * g + 8] = 0;

        v4lapltau3[8 * g + 0] = 0;
        v4lapltau3[8 * g + 1] = 0;
        v4lapltau3[8 * g + 2] = 0;
        v4lapltau3[8 * g + 3] = 0;
        v4lapltau3[8 * g + 4] = 0;
        v4lapltau3[8 * g + 5] = 0;
        v4lapltau3[8 * g + 6] = 0;
        v4lapltau3[8 * g + 7] = 0;

        v4tau4[5 * g + 0] = 0;
        v4tau4[5 * g + 1] = 0;
        v4tau4[5 * g + 2] = 0;
        v4tau4[5 * g + 3] = 0;
        v4tau4[5 * g + 4] = 0;

    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);
    
    // auto stage_v3rho3 = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[55 * _ldStaging];
    // auto stage_v3rho2sigma = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[59 * _ldStaging];
    // auto stage_v3rho2lapl = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[68 * _ldStaging];
    // auto stage_v3rho2tau = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[74 * _ldStaging];
    // auto stage_v3rhosigma2 = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[80 * _ldStaging];
    // auto stage_v3rhosigmalapl = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[92 * _ldStaging];
    // auto stage_v3rhosigmatau = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[104 * _ldStaging];
    // auto stage_v3rholapl2 = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[116 * _ldStaging];
    // auto stage_v3rholapltau = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[122 * _ldStaging];
    // auto stage_v3rhotau2 = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[130 * _ldStaging];
    // auto stage_v3sigma3 = (alloc) ? mem::malloc<double>(10 * np) : &_stagingBuffer[136 * _ldStaging];
    // auto stage_v3sigma2lapl = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[146 * _ldStaging];
    // auto stage_v3sigma2tau = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[158 * _ldStaging];
    // auto stage_v3sigmalapl2 = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[170 * _ldStaging];
    // auto stage_v3sigmalapltau = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[179 * _ldStaging];
    // auto stage_v3sigmatau2 = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[191 * _ldStaging];
    // auto stage_v3lapl3 = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[200 * _ldStaging];
    // auto stage_v3lapl2tau = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[204 * _ldStaging];
    // auto stage_v3lapltau2 = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[210 * _ldStaging];
    // auto stage_v3tau3 = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[216 * _ldStaging];
    auto stage_v4rho4   = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[220 * _ldStaging];
    auto stage_v4rho3sigma   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[225 * _ldStaging];
    auto stage_v4rho3lapl   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[237 * _ldStaging];
    auto stage_v4rho3tau   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[245 * _ldStaging];
    auto stage_v4rho2sigma2   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[253 * _ldStaging];
    auto stage_v4rho2sigmalapl   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[271 * _ldStaging];
    auto stage_v4rho2sigmatau   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[289 * _ldStaging];
    auto stage_v4rho2lapl2   = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[307 * _ldStaging];
    auto stage_v4rho2lapltau   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[316 * _ldStaging];
    auto stage_v4rho2tau2   = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[328 * _ldStaging];
    auto stage_v4rhosigma3   = (alloc) ? mem::malloc<double>(20 * np) : &_stagingBuffer[337 * _ldStaging];
    auto stage_v4rhosigma2lapl   = (alloc) ? mem::malloc<double>(36 * np) : &_stagingBuffer[357 * _ldStaging];
    auto stage_v4rhosigma2tau   = (alloc) ? mem::malloc<double>(36 * np) : &_stagingBuffer[393 * _ldStaging];
    auto stage_v4rhosigmalapl2   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[429 * _ldStaging];
    auto stage_v4rhosigmalapltau   = (alloc) ? mem::malloc<double>(24 * np) : &_stagingBuffer[447 * _ldStaging];
    auto stage_v4rhosigmatau2   = (alloc) ? mem::malloc<double>(36 * np) : &_stagingBuffer[471 * _ldStaging];
    auto stage_v4rholapl3   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[507 * _ldStaging];
    auto stage_v4rholapl2tau   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[515 * _ldStaging];
    auto stage_v4rholapltau2   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[527 * _ldStaging];
    auto stage_v4rhotau3   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[539 * _ldStaging];
    auto stage_v4sigma4   = (alloc) ? mem::malloc<double>(15 * np) : &_stagingBuffer[547 * _ldStaging];
    auto stage_v4sigma3lapl   = (alloc) ? mem::malloc<double>(20 * np) : &_stagingBuffer[562 * _ldStaging];
    auto stage_v4sigma3tau   = (alloc) ? mem::malloc<double>(30 * np) : &_stagingBuffer[582 * _ldStaging];
    auto stage_v4sigma2lapl2   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[612 * _ldStaging];
    auto stage_v4sigma2lapltau   = (alloc) ? mem::malloc<double>(24 * np) : &_stagingBuffer[630 * _ldStaging];
    auto stage_v4sigma2tau2   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[654 * _ldStaging];
    auto stage_v4sigmalapl3   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[672 * _ldStaging];
    auto stage_v4sigmalapl2tau   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[684 * _ldStaging];
    auto stage_v4sigmalapltau2   = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[702 * _ldStaging];
    auto stage_v4sigmatau3   = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[720 * _ldStaging];
    auto stage_v4lapl4   = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[732 * _ldStaging];
    auto stage_v4lapl3tau   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[737 * _ldStaging];
    auto stage_v4lapl2tau2   = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[745 * _ldStaging];
    auto stage_v4lapltau3   = (alloc) ? mem::malloc<double>(8 * np) : &_stagingBuffer[754 * _ldStaging];
    auto stage_v4tau4   = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[762 * _ldStaging];
    
    
    
    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        auto family = funcptr->info->family;

        if (family == XC_FAMILY_LDA)
        {
            xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

#pragma omp simd aligned(v4rho4, stage_v4rho4 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {

                v4rho4[5 * g + 0] += c * stage_v4rho4[5 * g + 0];
                v4rho4[5 * g + 1] += c * stage_v4rho4[5 * g + 1];
                v4rho4[5 * g + 2] += c * stage_v4rho4[5 * g + 2];
                v4rho4[5 * g + 3] += c * stage_v4rho4[5 * g + 3];
                v4rho4[5 * g + 4] += c * stage_v4rho4[5 * g + 4];
        
            }
        }
        else if (family == XC_FAMILY_GGA)
        {
            xc_gga_lxc(funcptr, np, rho, sigma, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4);

#pragma omp simd aligned(v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4 : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                v4rho4[5 * g + 0] += c * stage_v4rho4[5 * g + 0];
                v4rho4[5 * g + 1] += c * stage_v4rho4[5 * g + 1];
                v4rho4[5 * g + 2] += c * stage_v4rho4[5 * g + 2];
                v4rho4[5 * g + 3] += c * stage_v4rho4[5 * g + 3];
                v4rho4[5 * g + 4] += c * stage_v4rho4[5 * g + 4];

                v4rho3sigma[12 * g + 0]  += c * stage_v4rho3sigma[12 * g + 0] ;
                v4rho3sigma[12 * g + 1]  += c * stage_v4rho3sigma[12 * g + 1] ;
                v4rho3sigma[12 * g + 2]  += c * stage_v4rho3sigma[12 * g + 2] ;
                v4rho3sigma[12 * g + 3]  += c * stage_v4rho3sigma[12 * g + 3] ;
                v4rho3sigma[12 * g + 4]  += c * stage_v4rho3sigma[12 * g + 4] ;
                v4rho3sigma[12 * g + 5]  += c * stage_v4rho3sigma[12 * g + 5] ;
                v4rho3sigma[12 * g + 6]  += c * stage_v4rho3sigma[12 * g + 6] ;
                v4rho3sigma[12 * g + 7]  += c * stage_v4rho3sigma[12 * g + 7] ;
                v4rho3sigma[12 * g + 8]  += c * stage_v4rho3sigma[12 * g + 8] ;
                v4rho3sigma[12 * g + 9]  += c * stage_v4rho3sigma[12 * g + 9] ;
                v4rho3sigma[12 * g + 10] += c * stage_v4rho3sigma[12 * g + 10];
                v4rho3sigma[12 * g + 11] += c * stage_v4rho3sigma[12 * g + 11] ;

                v4rho2sigma2[18 * g + 0]  += c *stage_v4rho2sigma2[18 * g + 0]  ;
                v4rho2sigma2[18 * g + 1]  += c *stage_v4rho2sigma2[18 * g + 1]  ;
                v4rho2sigma2[18 * g + 2]  += c *stage_v4rho2sigma2[18 * g + 2]  ;
                v4rho2sigma2[18 * g + 3]  += c *stage_v4rho2sigma2[18 * g + 3]  ;
                v4rho2sigma2[18 * g + 4]  += c *stage_v4rho2sigma2[18 * g + 4]  ;
                v4rho2sigma2[18 * g + 5]  += c *stage_v4rho2sigma2[18 * g + 5]  ;
                v4rho2sigma2[18 * g + 6]  += c *stage_v4rho2sigma2[18 * g + 6]  ;
                v4rho2sigma2[18 * g + 7]  += c *stage_v4rho2sigma2[18 * g + 7]  ;
                v4rho2sigma2[18 * g + 8]  += c *stage_v4rho2sigma2[18 * g + 8]  ;
                v4rho2sigma2[18 * g + 9]  += c *stage_v4rho2sigma2[18 * g + 9]  ;
                v4rho2sigma2[18 * g + 10] += c *stage_v4rho2sigma2[18 * g + 10] ;
                v4rho2sigma2[18 * g + 11] += c *stage_v4rho2sigma2[18 * g + 11] ;
                v4rho2sigma2[18 * g + 12] += c *stage_v4rho2sigma2[18 * g + 12] ;
                v4rho2sigma2[18 * g + 13] += c *stage_v4rho2sigma2[18 * g + 13] ;
                v4rho2sigma2[18 * g + 14] += c *stage_v4rho2sigma2[18 * g + 14] ;
                v4rho2sigma2[18 * g + 15] += c *stage_v4rho2sigma2[18 * g + 15] ;
                v4rho2sigma2[18 * g + 16] += c *stage_v4rho2sigma2[18 * g + 16] ;
                v4rho2sigma2[18 * g + 17] += c *stage_v4rho2sigma2[18 * g + 17] ;

                v4rhosigma3[20 * g + 0]  += c *stage_v4rhosigma3[20 * g + 0]  ;
                v4rhosigma3[20 * g + 1]  += c *stage_v4rhosigma3[20 * g + 1]  ;
                v4rhosigma3[20 * g + 2]  += c *stage_v4rhosigma3[20 * g + 2]  ;
                v4rhosigma3[20 * g + 3]  += c *stage_v4rhosigma3[20 * g + 3]  ;
                v4rhosigma3[20 * g + 4]  += c *stage_v4rhosigma3[20 * g + 4]  ;
                v4rhosigma3[20 * g + 5]  += c *stage_v4rhosigma3[20 * g + 5]  ;
                v4rhosigma3[20 * g + 6]  += c *stage_v4rhosigma3[20 * g + 6]  ;
                v4rhosigma3[20 * g + 7]  += c *stage_v4rhosigma3[20 * g + 7]  ;
                v4rhosigma3[20 * g + 8]  += c *stage_v4rhosigma3[20 * g + 8]  ;
                v4rhosigma3[20 * g + 9]  += c *stage_v4rhosigma3[20 * g + 9]  ;
                v4rhosigma3[20 * g + 10] += c *stage_v4rhosigma3[20 * g + 10] ;
                v4rhosigma3[20 * g + 11] += c *stage_v4rhosigma3[20 * g + 11] ;
                v4rhosigma3[20 * g + 12] += c *stage_v4rhosigma3[20 * g + 12] ;
                v4rhosigma3[20 * g + 13] += c *stage_v4rhosigma3[20 * g + 13] ;
                v4rhosigma3[20 * g + 14] += c *stage_v4rhosigma3[20 * g + 14] ;
                v4rhosigma3[20 * g + 15] += c *stage_v4rhosigma3[20 * g + 15] ;
                v4rhosigma3[20 * g + 16] += c *stage_v4rhosigma3[20 * g + 16] ;
                v4rhosigma3[20 * g + 17] += c *stage_v4rhosigma3[20 * g + 17] ;
                v4rhosigma3[20 * g + 18] += c *stage_v4rhosigma3[20 * g + 18] ;
                v4rhosigma3[20 * g + 19] += c *stage_v4rhosigma3[20 * g + 19] ;

                v4sigma4[15 * g + 0]  += c *stage_v4sigma4[15 * g + 0]  ;
                v4sigma4[15 * g + 1]  += c *stage_v4sigma4[15 * g + 1]  ;
                v4sigma4[15 * g + 2]  += c *stage_v4sigma4[15 * g + 2]  ;
                v4sigma4[15 * g + 3]  += c *stage_v4sigma4[15 * g + 3]  ;
                v4sigma4[15 * g + 4]  += c *stage_v4sigma4[15 * g + 4]  ;
                v4sigma4[15 * g + 5]  += c *stage_v4sigma4[15 * g + 5]  ;
                v4sigma4[15 * g + 6]  += c *stage_v4sigma4[15 * g + 6]  ;
                v4sigma4[15 * g + 7]  += c *stage_v4sigma4[15 * g + 7]  ;
                v4sigma4[15 * g + 8]  += c *stage_v4sigma4[15 * g + 8]  ;
                v4sigma4[15 * g + 9]  += c *stage_v4sigma4[15 * g + 9]  ;
                v4sigma4[15 * g + 10] += c *stage_v4sigma4[15 * g + 10] ;
                v4sigma4[15 * g + 11] += c *stage_v4sigma4[15 * g + 11] ;
                v4sigma4[15 * g + 12] += c *stage_v4sigma4[15 * g + 12] ;
                v4sigma4[15 * g + 13] += c *stage_v4sigma4[15 * g + 13] ;
                v4sigma4[15 * g + 14] += c *stage_v4sigma4[15 * g + 14] ;

            }
        }
        else if (family == XC_FAMILY_MGGA)
        {
            xc_mgga_lxc(funcptr, np, rho, sigma, lapl, tau,
                         stage_v4rho4, stage_v4rho3sigma, stage_v4rho3lapl, stage_v4rho3tau, stage_v4rho2sigma2,
                         stage_v4rho2sigmalapl, stage_v4rho2sigmatau, stage_v4rho2lapl2, stage_v4rho2lapltau,
                         stage_v4rho2tau2, stage_v4rhosigma3, stage_v4rhosigma2lapl, stage_v4rhosigma2tau,
                         stage_v4rhosigmalapl2, stage_v4rhosigmalapltau, stage_v4rhosigmatau2,
                         stage_v4rholapl3, stage_v4rholapl2tau, stage_v4rholapltau2, stage_v4rhotau3,
                         stage_v4sigma4, stage_v4sigma3lapl, stage_v4sigma3tau, stage_v4sigma2lapl2,
                         stage_v4sigma2lapltau, stage_v4sigma2tau2, stage_v4sigmalapl3, stage_v4sigmalapl2tau,
                         stage_v4sigmalapltau2, stage_v4sigmatau3, stage_v4lapl4, stage_v4lapl3tau,
                         stage_v4lapl2tau2, stage_v4lapltau3, stage_v4tau4);

#pragma omp simd aligned(rho,sigma,tau,lapl,stage_v4rho4,v4rho4,stage_v4rho3sigma,v4rho3sigma,stage_v4rho3lapl,v4rho3lapl,stage_v4rho3tau,v4rho3tau,stage_v4rho2sigma2,v4rho2sigma2,\
                         stage_v4rho2sigmalapl,v4rho2sigmalapl,stage_v4rho2sigmatau,v4rho2sigmatau,stage_v4rho2lapl2,v4rho2lapl2,stage_v4rho2lapltau,v4rho2lapltau,stage_v4rho2tau2,v4rho2tau2,\
                         stage_v4rhosigma3,v4rhosigma3,stage_v4rhosigma2lapl,v4rhosigma2lapl,stage_v4rhosigma2tau,v4rhosigma2tau,stage_v4rhosigmalapl2,v4rhosigmalapl2,stage_v4rhosigmalapltau,v4rhosigmalapltau,\
                         stage_v4rhosigmatau2,v4rhosigmatau2,stage_v4rholapl3,v4rholapl3,stage_v4rholapl2tau,v4rholapl2tau,stage_v4rholapltau2,v4rholapltau2,stage_v4rhotau3,v4rhotau3,\
                         stage_v4sigma4,v4sigma4,stage_v4sigma3lapl,v4sigma3lapl,stage_v4sigma3tau,v4sigma3tau,stage_v4sigma2lapl2,v4sigma2lapl2,stage_v4sigma2lapltau,v4sigma2lapltau,\
                         stage_v4sigma2tau2,v4sigma2tau2,stage_v4sigmalapl3,v4sigmalapl3,stage_v4sigmalapl2tau,v4sigmalapl2tau,stage_v4sigmalapltau2,v4sigmalapltau2,stage_v4sigmatau3,v4sigmatau3,\
                         stage_v4lapl4,v4lapl4,stage_v4lapl3tau,v4lapl3tau,stage_v4lapl2tau2,v4lapl2tau2,stage_v4lapltau3,v4lapltau3,stage_v4tau4,v4tau4 : VLX_ALIGN)
            
            for (auto g = 0; g < np; ++g)
            {

                v4rho4[5 * g + 0] +=  c * stage_v4rho4[5 * g + 0];
                v4rho4[5 * g + 1] +=  c * stage_v4rho4[5 * g + 1];
                v4rho4[5 * g + 2] +=  c * stage_v4rho4[5 * g + 2];
                v4rho4[5 * g + 3] +=  c * stage_v4rho4[5 * g + 3];
                v4rho4[5 * g + 4] +=  c * stage_v4rho4[5 * g + 4];

                v4rho3sigma[12 * g + 0] +=  c * stage_v4rho3sigma[12 * g + 0];
                v4rho3sigma[12 * g + 1] +=  c * stage_v4rho3sigma[12 * g + 1];
                v4rho3sigma[12 * g + 2] +=  c * stage_v4rho3sigma[12 * g + 2];
                v4rho3sigma[12 * g + 3] +=  c * stage_v4rho3sigma[12 * g + 3];
                v4rho3sigma[12 * g + 4] +=  c * stage_v4rho3sigma[12 * g + 4];
                v4rho3sigma[12 * g + 5] +=  c * stage_v4rho3sigma[12 * g + 5];
                v4rho3sigma[12 * g + 6] +=  c * stage_v4rho3sigma[12 * g + 6];
                v4rho3sigma[12 * g + 7] +=  c * stage_v4rho3sigma[12 * g + 7];
                v4rho3sigma[12 * g + 8] +=  c * stage_v4rho3sigma[12 * g + 8];
                v4rho3sigma[12 * g + 9] +=  c * stage_v4rho3sigma[12 * g + 9];
                v4rho3sigma[12 * g + 10] +=  c * stage_v4rho3sigma[12 * g + 10];
                v4rho3sigma[12 * g + 11] +=  c * stage_v4rho3sigma[12 * g + 11];

                v4rho3lapl[8 * g + 0] +=  c * stage_v4rho3lapl[8 * g + 0];
                v4rho3lapl[8 * g + 1] +=  c * stage_v4rho3lapl[8 * g + 1];
                v4rho3lapl[8 * g + 2] +=  c * stage_v4rho3lapl[8 * g + 2];
                v4rho3lapl[8 * g + 3] +=  c * stage_v4rho3lapl[8 * g + 3];
                v4rho3lapl[8 * g + 4] +=  c * stage_v4rho3lapl[8 * g + 4];
                v4rho3lapl[8 * g + 5] +=  c * stage_v4rho3lapl[8 * g + 5];
                v4rho3lapl[8 * g + 6] +=  c * stage_v4rho3lapl[8 * g + 6];
                v4rho3lapl[8 * g + 7] +=  c * stage_v4rho3lapl[8 * g + 7];

                v4rho3tau[8 * g + 0] +=  c * stage_v4rho3tau[8 * g + 0];
                v4rho3tau[8 * g + 1] +=  c * stage_v4rho3tau[8 * g + 1];
                v4rho3tau[8 * g + 2] +=  c * stage_v4rho3tau[8 * g + 2];
                v4rho3tau[8 * g + 3] +=  c * stage_v4rho3tau[8 * g + 3];
                v4rho3tau[8 * g + 4] +=  c * stage_v4rho3tau[8 * g + 4];
                v4rho3tau[8 * g + 5] +=  c * stage_v4rho3tau[8 * g + 5];
                v4rho3tau[8 * g + 6] +=  c * stage_v4rho3tau[8 * g + 6];
                v4rho3tau[8 * g + 7] +=  c * stage_v4rho3tau[8 * g + 7];

                v4rho2sigma2[18 * g + 0] +=  c * stage_v4rho2sigma2[18 * g + 0];
                v4rho2sigma2[18 * g + 1] +=  c * stage_v4rho2sigma2[18 * g + 1];
                v4rho2sigma2[18 * g + 2] +=  c * stage_v4rho2sigma2[18 * g + 2];
                v4rho2sigma2[18 * g + 3] +=  c * stage_v4rho2sigma2[18 * g + 3];
                v4rho2sigma2[18 * g + 4] +=  c * stage_v4rho2sigma2[18 * g + 4];
                v4rho2sigma2[18 * g + 5] +=  c * stage_v4rho2sigma2[18 * g + 5];
                v4rho2sigma2[18 * g + 6] +=  c * stage_v4rho2sigma2[18 * g + 6];
                v4rho2sigma2[18 * g + 7] +=  c * stage_v4rho2sigma2[18 * g + 7];
                v4rho2sigma2[18 * g + 8] +=  c * stage_v4rho2sigma2[18 * g + 8];
                v4rho2sigma2[18 * g + 9] +=  c * stage_v4rho2sigma2[18 * g + 9];
                v4rho2sigma2[18 * g + 10] +=  c * stage_v4rho2sigma2[18 * g + 10];
                v4rho2sigma2[18 * g + 11] +=  c * stage_v4rho2sigma2[18 * g + 11];
                v4rho2sigma2[18 * g + 12] +=  c * stage_v4rho2sigma2[18 * g + 12];
                v4rho2sigma2[18 * g + 13] +=  c * stage_v4rho2sigma2[18 * g + 13];
                v4rho2sigma2[18 * g + 14] +=  c * stage_v4rho2sigma2[18 * g + 14];
                v4rho2sigma2[18 * g + 15] +=  c * stage_v4rho2sigma2[18 * g + 15];
                v4rho2sigma2[18 * g + 16] +=  c * stage_v4rho2sigma2[18 * g + 16];
                v4rho2sigma2[18 * g + 17] +=  c * stage_v4rho2sigma2[18 * g + 17];

                v4rho2sigmalapl[18 * g + 0] +=  c * stage_v4rho2sigmalapl[18 * g + 0];
                v4rho2sigmalapl[18 * g + 1] +=  c * stage_v4rho2sigmalapl[18 * g + 1];
                v4rho2sigmalapl[18 * g + 2] +=  c * stage_v4rho2sigmalapl[18 * g + 2];
                v4rho2sigmalapl[18 * g + 3] +=  c * stage_v4rho2sigmalapl[18 * g + 3];
                v4rho2sigmalapl[18 * g + 4] +=  c * stage_v4rho2sigmalapl[18 * g + 4];
                v4rho2sigmalapl[18 * g + 5] +=  c * stage_v4rho2sigmalapl[18 * g + 5];
                v4rho2sigmalapl[18 * g + 6] +=  c * stage_v4rho2sigmalapl[18 * g + 6];
                v4rho2sigmalapl[18 * g + 7] +=  c * stage_v4rho2sigmalapl[18 * g + 7];
                v4rho2sigmalapl[18 * g + 8] +=  c * stage_v4rho2sigmalapl[18 * g + 8];
                v4rho2sigmalapl[18 * g + 9] +=  c * stage_v4rho2sigmalapl[18 * g + 9];
                v4rho2sigmalapl[18 * g + 10] +=  c * stage_v4rho2sigmalapl[18 * g + 10];
                v4rho2sigmalapl[18 * g + 11] +=  c * stage_v4rho2sigmalapl[18 * g + 11];
                v4rho2sigmalapl[18 * g + 12] +=  c * stage_v4rho2sigmalapl[18 * g + 12];
                v4rho2sigmalapl[18 * g + 13] +=  c * stage_v4rho2sigmalapl[18 * g + 13];
                v4rho2sigmalapl[18 * g + 14] +=  c * stage_v4rho2sigmalapl[18 * g + 14];
                v4rho2sigmalapl[18 * g + 15] +=  c * stage_v4rho2sigmalapl[18 * g + 15];
                v4rho2sigmalapl[18 * g + 16] +=  c * stage_v4rho2sigmalapl[18 * g + 16];
                v4rho2sigmalapl[18 * g + 17] +=  c * stage_v4rho2sigmalapl[18 * g + 17];

                v4rho2sigmatau[18 * g + 0] +=  c * stage_v4rho2sigmatau[18 * g + 0];
                v4rho2sigmatau[18 * g + 1] +=  c * stage_v4rho2sigmatau[18 * g + 1];
                v4rho2sigmatau[18 * g + 2] +=  c * stage_v4rho2sigmatau[18 * g + 2];
                v4rho2sigmatau[18 * g + 3] +=  c * stage_v4rho2sigmatau[18 * g + 3];
                v4rho2sigmatau[18 * g + 4] +=  c * stage_v4rho2sigmatau[18 * g + 4];
                v4rho2sigmatau[18 * g + 5] +=  c * stage_v4rho2sigmatau[18 * g + 5];
                v4rho2sigmatau[18 * g + 6] +=  c * stage_v4rho2sigmatau[18 * g + 6];
                v4rho2sigmatau[18 * g + 7] +=  c * stage_v4rho2sigmatau[18 * g + 7];
                v4rho2sigmatau[18 * g + 8] +=  c * stage_v4rho2sigmatau[18 * g + 8];
                v4rho2sigmatau[18 * g + 9] +=  c * stage_v4rho2sigmatau[18 * g + 9];
                v4rho2sigmatau[18 * g + 10] +=  c * stage_v4rho2sigmatau[18 * g + 10];
                v4rho2sigmatau[18 * g + 11] +=  c * stage_v4rho2sigmatau[18 * g + 11];
                v4rho2sigmatau[18 * g + 12] +=  c * stage_v4rho2sigmatau[18 * g + 12];
                v4rho2sigmatau[18 * g + 13] +=  c * stage_v4rho2sigmatau[18 * g + 13];
                v4rho2sigmatau[18 * g + 14] +=  c * stage_v4rho2sigmatau[18 * g + 14];
                v4rho2sigmatau[18 * g + 15] +=  c * stage_v4rho2sigmatau[18 * g + 15];
                v4rho2sigmatau[18 * g + 16] +=  c * stage_v4rho2sigmatau[18 * g + 16];
                v4rho2sigmatau[18 * g + 17] +=  c * stage_v4rho2sigmatau[18 * g + 17];

                v4rho2lapl2[9 * g + 0] +=  c * stage_v4rho2lapl2[9 * g + 0];
                v4rho2lapl2[9 * g + 1] +=  c * stage_v4rho2lapl2[9 * g + 1];
                v4rho2lapl2[9 * g + 2] +=  c * stage_v4rho2lapl2[9 * g + 2];
                v4rho2lapl2[9 * g + 3] +=  c * stage_v4rho2lapl2[9 * g + 3];
                v4rho2lapl2[9 * g + 4] +=  c * stage_v4rho2lapl2[9 * g + 4];
                v4rho2lapl2[9 * g + 5] +=  c * stage_v4rho2lapl2[9 * g + 5];
                v4rho2lapl2[9 * g + 6] +=  c * stage_v4rho2lapl2[9 * g + 6];
                v4rho2lapl2[9 * g + 7] +=  c * stage_v4rho2lapl2[9 * g + 7];
                v4rho2lapl2[9 * g + 8] +=  c * stage_v4rho2lapl2[9 * g + 8];

                v4rho2lapltau[12 * g + 0] +=  c * stage_v4rho2lapltau[12 * g + 0];
                v4rho2lapltau[12 * g + 1] +=  c * stage_v4rho2lapltau[12 * g + 1];
                v4rho2lapltau[12 * g + 2] +=  c * stage_v4rho2lapltau[12 * g + 2];
                v4rho2lapltau[12 * g + 3] +=  c * stage_v4rho2lapltau[12 * g + 3];
                v4rho2lapltau[12 * g + 4] +=  c * stage_v4rho2lapltau[12 * g + 4];
                v4rho2lapltau[12 * g + 5] +=  c * stage_v4rho2lapltau[12 * g + 5];
                v4rho2lapltau[12 * g + 6] +=  c * stage_v4rho2lapltau[12 * g + 6];
                v4rho2lapltau[12 * g + 7] +=  c * stage_v4rho2lapltau[12 * g + 7];
                v4rho2lapltau[12 * g + 8] +=  c * stage_v4rho2lapltau[12 * g + 8];
                v4rho2lapltau[12 * g + 9] +=  c * stage_v4rho2lapltau[12 * g + 9];
                v4rho2lapltau[12 * g + 10] +=  c * stage_v4rho2lapltau[12 * g + 10];
                v4rho2lapltau[12 * g + 11] +=  c * stage_v4rho2lapltau[12 * g + 11];

                v4rho2tau2[9 * g + 0] +=  c * stage_v4rho2tau2[9 * g + 0];
                v4rho2tau2[9 * g + 1] +=  c * stage_v4rho2tau2[9 * g + 1];
                v4rho2tau2[9 * g + 2] +=  c * stage_v4rho2tau2[9 * g + 2];
                v4rho2tau2[9 * g + 3] +=  c * stage_v4rho2tau2[9 * g + 3];
                v4rho2tau2[9 * g + 4] +=  c * stage_v4rho2tau2[9 * g + 4];
                v4rho2tau2[9 * g + 5] +=  c * stage_v4rho2tau2[9 * g + 5];
                v4rho2tau2[9 * g + 6] +=  c * stage_v4rho2tau2[9 * g + 6];
                v4rho2tau2[9 * g + 7] +=  c * stage_v4rho2tau2[9 * g + 7];
                v4rho2tau2[9 * g + 8] +=  c * stage_v4rho2tau2[9 * g + 8];

                v4rhosigma3[20 * g + 0] +=  c * stage_v4rhosigma3[20 * g + 0];
                v4rhosigma3[20 * g + 1] +=  c * stage_v4rhosigma3[20 * g + 1];
                v4rhosigma3[20 * g + 2] +=  c * stage_v4rhosigma3[20 * g + 2];
                v4rhosigma3[20 * g + 3] +=  c * stage_v4rhosigma3[20 * g + 3];
                v4rhosigma3[20 * g + 4] +=  c * stage_v4rhosigma3[20 * g + 4];
                v4rhosigma3[20 * g + 5] +=  c * stage_v4rhosigma3[20 * g + 5];
                v4rhosigma3[20 * g + 6] +=  c * stage_v4rhosigma3[20 * g + 6];
                v4rhosigma3[20 * g + 7] +=  c * stage_v4rhosigma3[20 * g + 7];
                v4rhosigma3[20 * g + 8] +=  c * stage_v4rhosigma3[20 * g + 8];
                v4rhosigma3[20 * g + 9] +=  c * stage_v4rhosigma3[20 * g + 9];
                v4rhosigma3[20 * g + 10] +=  c * stage_v4rhosigma3[20 * g + 10];
                v4rhosigma3[20 * g + 11] +=  c * stage_v4rhosigma3[20 * g + 11];
                v4rhosigma3[20 * g + 12] +=  c * stage_v4rhosigma3[20 * g + 12];
                v4rhosigma3[20 * g + 13] +=  c * stage_v4rhosigma3[20 * g + 13];
                v4rhosigma3[20 * g + 14] +=  c * stage_v4rhosigma3[20 * g + 14];
                v4rhosigma3[20 * g + 15] +=  c * stage_v4rhosigma3[20 * g + 15];
                v4rhosigma3[20 * g + 16] +=  c * stage_v4rhosigma3[20 * g + 16];
                v4rhosigma3[20 * g + 17] +=  c * stage_v4rhosigma3[20 * g + 17];
                v4rhosigma3[20 * g + 18] +=  c * stage_v4rhosigma3[20 * g + 18];
                v4rhosigma3[20 * g + 19] +=  c * stage_v4rhosigma3[20 * g + 19];

                v4rhosigma2lapl[36 * g + 0] +=  c * stage_v4rhosigma2lapl[36 * g + 0];
                v4rhosigma2lapl[36 * g + 1] +=  c * stage_v4rhosigma2lapl[36 * g + 1];
                v4rhosigma2lapl[36 * g + 2] +=  c * stage_v4rhosigma2lapl[36 * g + 2];
                v4rhosigma2lapl[36 * g + 3] +=  c * stage_v4rhosigma2lapl[36 * g + 3];
                v4rhosigma2lapl[36 * g + 4] +=  c * stage_v4rhosigma2lapl[36 * g + 4];
                v4rhosigma2lapl[36 * g + 5] +=  c * stage_v4rhosigma2lapl[36 * g + 5];
                v4rhosigma2lapl[36 * g + 6] +=  c * stage_v4rhosigma2lapl[36 * g + 6];
                v4rhosigma2lapl[36 * g + 7] +=  c * stage_v4rhosigma2lapl[36 * g + 7];
                v4rhosigma2lapl[36 * g + 8] +=  c * stage_v4rhosigma2lapl[36 * g + 8];
                v4rhosigma2lapl[36 * g + 9] +=  c * stage_v4rhosigma2lapl[36 * g + 9];
                v4rhosigma2lapl[36 * g + 10] +=  c * stage_v4rhosigma2lapl[36 * g + 10];
                v4rhosigma2lapl[36 * g + 11] +=  c * stage_v4rhosigma2lapl[36 * g + 11];
                v4rhosigma2lapl[36 * g + 12] +=  c * stage_v4rhosigma2lapl[36 * g + 12];
                v4rhosigma2lapl[36 * g + 13] +=  c * stage_v4rhosigma2lapl[36 * g + 13];
                v4rhosigma2lapl[36 * g + 14] +=  c * stage_v4rhosigma2lapl[36 * g + 14];
                v4rhosigma2lapl[36 * g + 15] +=  c * stage_v4rhosigma2lapl[36 * g + 15];
                v4rhosigma2lapl[36 * g + 16] +=  c * stage_v4rhosigma2lapl[36 * g + 16];
                v4rhosigma2lapl[36 * g + 17] +=  c * stage_v4rhosigma2lapl[36 * g + 17];
                v4rhosigma2lapl[36 * g + 18] +=  c * stage_v4rhosigma2lapl[36 * g + 18];
                v4rhosigma2lapl[36 * g + 19] +=  c * stage_v4rhosigma2lapl[36 * g + 19];
                v4rhosigma2lapl[36 * g + 20] +=  c * stage_v4rhosigma2lapl[36 * g + 20];
                v4rhosigma2lapl[36 * g + 21] +=  c * stage_v4rhosigma2lapl[36 * g + 21];
                v4rhosigma2lapl[36 * g + 22] +=  c * stage_v4rhosigma2lapl[36 * g + 22];
                v4rhosigma2lapl[36 * g + 23] +=  c * stage_v4rhosigma2lapl[36 * g + 23];

                v4rhosigma2tau[36 * g + 0] +=  c * stage_v4rhosigma2tau[36 * g + 0];
                v4rhosigma2tau[36 * g + 1] +=  c * stage_v4rhosigma2tau[36 * g + 1];
                v4rhosigma2tau[36 * g + 2] +=  c * stage_v4rhosigma2tau[36 * g + 2];
                v4rhosigma2tau[36 * g + 3] +=  c * stage_v4rhosigma2tau[36 * g + 3];
                v4rhosigma2tau[36 * g + 4] +=  c * stage_v4rhosigma2tau[36 * g + 4];
                v4rhosigma2tau[36 * g + 5] +=  c * stage_v4rhosigma2tau[36 * g + 5];
                v4rhosigma2tau[36 * g + 6] +=  c * stage_v4rhosigma2tau[36 * g + 6];
                v4rhosigma2tau[36 * g + 7] +=  c * stage_v4rhosigma2tau[36 * g + 7];
                v4rhosigma2tau[36 * g + 8] +=  c * stage_v4rhosigma2tau[36 * g + 8];
                v4rhosigma2tau[36 * g + 9] +=  c * stage_v4rhosigma2tau[36 * g + 9];
                v4rhosigma2tau[36 * g + 10] +=  c * stage_v4rhosigma2tau[36 * g + 10];
                v4rhosigma2tau[36 * g + 11] +=  c * stage_v4rhosigma2tau[36 * g + 11];
                v4rhosigma2tau[36 * g + 12] +=  c * stage_v4rhosigma2tau[36 * g + 12];
                v4rhosigma2tau[36 * g + 13] +=  c * stage_v4rhosigma2tau[36 * g + 13];
                v4rhosigma2tau[36 * g + 14] +=  c * stage_v4rhosigma2tau[36 * g + 14];
                v4rhosigma2tau[36 * g + 15] +=  c * stage_v4rhosigma2tau[36 * g + 15];
                v4rhosigma2tau[36 * g + 16] +=  c * stage_v4rhosigma2tau[36 * g + 16];
                v4rhosigma2tau[36 * g + 17] +=  c * stage_v4rhosigma2tau[36 * g + 17];
                v4rhosigma2tau[36 * g + 18] +=  c * stage_v4rhosigma2tau[36 * g + 18];
                v4rhosigma2tau[36 * g + 19] +=  c * stage_v4rhosigma2tau[36 * g + 19];
                v4rhosigma2tau[36 * g + 20] +=  c * stage_v4rhosigma2tau[36 * g + 20];
                v4rhosigma2tau[36 * g + 21] +=  c * stage_v4rhosigma2tau[36 * g + 21];
                v4rhosigma2tau[36 * g + 22] +=  c * stage_v4rhosigma2tau[36 * g + 22];
                v4rhosigma2tau[36 * g + 23] +=  c * stage_v4rhosigma2tau[36 * g + 23];

                v4rhosigmalapl2[18 * g + 0] +=  c * stage_v4rhosigmalapl2[18 * g + 0];
                v4rhosigmalapl2[18 * g + 1] +=  c * stage_v4rhosigmalapl2[18 * g + 1];
                v4rhosigmalapl2[18 * g + 2] +=  c * stage_v4rhosigmalapl2[18 * g + 2];
                v4rhosigmalapl2[18 * g + 3] +=  c * stage_v4rhosigmalapl2[18 * g + 3];
                v4rhosigmalapl2[18 * g + 4] +=  c * stage_v4rhosigmalapl2[18 * g + 4];
                v4rhosigmalapl2[18 * g + 5] +=  c * stage_v4rhosigmalapl2[18 * g + 5];
                v4rhosigmalapl2[18 * g + 6] +=  c * stage_v4rhosigmalapl2[18 * g + 6];
                v4rhosigmalapl2[18 * g + 7] +=  c * stage_v4rhosigmalapl2[18 * g + 7];
                v4rhosigmalapl2[18 * g + 8] +=  c * stage_v4rhosigmalapl2[18 * g + 8];
                v4rhosigmalapl2[18 * g + 9] +=  c * stage_v4rhosigmalapl2[18 * g + 9];
                v4rhosigmalapl2[18 * g + 10] +=  c * stage_v4rhosigmalapl2[18 * g + 10];
                v4rhosigmalapl2[18 * g + 11] +=  c * stage_v4rhosigmalapl2[18 * g + 11];
                v4rhosigmalapl2[18 * g + 12] +=  c * stage_v4rhosigmalapl2[18 * g + 12];
                v4rhosigmalapl2[18 * g + 13] +=  c * stage_v4rhosigmalapl2[18 * g + 13];
                v4rhosigmalapl2[18 * g + 14] +=  c * stage_v4rhosigmalapl2[18 * g + 14];
                v4rhosigmalapl2[18 * g + 15] +=  c * stage_v4rhosigmalapl2[18 * g + 15];
                v4rhosigmalapl2[18 * g + 16] +=  c * stage_v4rhosigmalapl2[18 * g + 16];
                v4rhosigmalapl2[18 * g + 17] +=  c * stage_v4rhosigmalapl2[18 * g + 17];

                v4rhosigmalapltau[24 * g + 0] +=  c * stage_v4rhosigmalapltau[24 * g + 0];
                v4rhosigmalapltau[24 * g + 1] +=  c * stage_v4rhosigmalapltau[24 * g + 1];
                v4rhosigmalapltau[24 * g + 2] +=  c * stage_v4rhosigmalapltau[24 * g + 2];
                v4rhosigmalapltau[24 * g + 3] +=  c * stage_v4rhosigmalapltau[24 * g + 3];
                v4rhosigmalapltau[24 * g + 4] +=  c * stage_v4rhosigmalapltau[24 * g + 4];
                v4rhosigmalapltau[24 * g + 5] +=  c * stage_v4rhosigmalapltau[24 * g + 5];
                v4rhosigmalapltau[24 * g + 6] +=  c * stage_v4rhosigmalapltau[24 * g + 6];
                v4rhosigmalapltau[24 * g + 7] +=  c * stage_v4rhosigmalapltau[24 * g + 7];
                v4rhosigmalapltau[24 * g + 8] +=  c * stage_v4rhosigmalapltau[24 * g + 8];
                v4rhosigmalapltau[24 * g + 9] +=  c * stage_v4rhosigmalapltau[24 * g + 9];
                v4rhosigmalapltau[24 * g + 10] +=  c * stage_v4rhosigmalapltau[24 * g + 10];
                v4rhosigmalapltau[24 * g + 11] +=  c * stage_v4rhosigmalapltau[24 * g + 11];
                v4rhosigmalapltau[24 * g + 12] +=  c * stage_v4rhosigmalapltau[24 * g + 12];
                v4rhosigmalapltau[24 * g + 13] +=  c * stage_v4rhosigmalapltau[24 * g + 13];
                v4rhosigmalapltau[24 * g + 14] +=  c * stage_v4rhosigmalapltau[24 * g + 14];
                v4rhosigmalapltau[24 * g + 15] +=  c * stage_v4rhosigmalapltau[24 * g + 15];
                v4rhosigmalapltau[24 * g + 16] +=  c * stage_v4rhosigmalapltau[24 * g + 16];
                v4rhosigmalapltau[24 * g + 17] +=  c * stage_v4rhosigmalapltau[24 * g + 17];
                v4rhosigmalapltau[24 * g + 18] +=  c * stage_v4rhosigmalapltau[24 * g + 18];
                v4rhosigmalapltau[24 * g + 19] +=  c * stage_v4rhosigmalapltau[24 * g + 19];
                v4rhosigmalapltau[24 * g + 20] +=  c * stage_v4rhosigmalapltau[24 * g + 20];
                v4rhosigmalapltau[24 * g + 21] +=  c * stage_v4rhosigmalapltau[24 * g + 21];
                v4rhosigmalapltau[24 * g + 22] +=  c * stage_v4rhosigmalapltau[24 * g + 22];
                v4rhosigmalapltau[24 * g + 23] +=  c * stage_v4rhosigmalapltau[24 * g + 23];

                v4rhosigmatau2[36 * g + 0] +=  c * stage_v4rhosigmatau2[36 * g + 0];
                v4rhosigmatau2[36 * g + 1] +=  c * stage_v4rhosigmatau2[36 * g + 1];
                v4rhosigmatau2[36 * g + 2] +=  c * stage_v4rhosigmatau2[36 * g + 2];
                v4rhosigmatau2[36 * g + 3] +=  c * stage_v4rhosigmatau2[36 * g + 3];
                v4rhosigmatau2[36 * g + 4] +=  c * stage_v4rhosigmatau2[36 * g + 4];
                v4rhosigmatau2[36 * g + 5] +=  c * stage_v4rhosigmatau2[36 * g + 5];
                v4rhosigmatau2[36 * g + 6] +=  c * stage_v4rhosigmatau2[36 * g + 6];
                v4rhosigmatau2[36 * g + 7] +=  c * stage_v4rhosigmatau2[36 * g + 7];
                v4rhosigmatau2[36 * g + 8] +=  c * stage_v4rhosigmatau2[36 * g + 8];
                v4rhosigmatau2[36 * g + 9] +=  c * stage_v4rhosigmatau2[36 * g + 9];
                v4rhosigmatau2[36 * g + 10] +=  c * stage_v4rhosigmatau2[36 * g + 10];
                v4rhosigmatau2[36 * g + 11] +=  c * stage_v4rhosigmatau2[36 * g + 11];
                v4rhosigmatau2[36 * g + 12] +=  c * stage_v4rhosigmatau2[36 * g + 12];
                v4rhosigmatau2[36 * g + 13] +=  c * stage_v4rhosigmatau2[36 * g + 13];
                v4rhosigmatau2[36 * g + 14] +=  c * stage_v4rhosigmatau2[36 * g + 14];
                v4rhosigmatau2[36 * g + 15] +=  c * stage_v4rhosigmatau2[36 * g + 15];
                v4rhosigmatau2[36 * g + 16] +=  c * stage_v4rhosigmatau2[36 * g + 16];
                v4rhosigmatau2[36 * g + 17] +=  c * stage_v4rhosigmatau2[36 * g + 17];

                v4rholapl3[8 * g + 0] +=  c * stage_v4rholapl3[8 * g + 0];
                v4rholapl3[8 * g + 1] +=  c * stage_v4rholapl3[8 * g + 1];
                v4rholapl3[8 * g + 2] +=  c * stage_v4rholapl3[8 * g + 2];
                v4rholapl3[8 * g + 3] +=  c * stage_v4rholapl3[8 * g + 3];
                v4rholapl3[8 * g + 4] +=  c * stage_v4rholapl3[8 * g + 4];
                v4rholapl3[8 * g + 5] +=  c * stage_v4rholapl3[8 * g + 5];
                v4rholapl3[8 * g + 6] +=  c * stage_v4rholapl3[8 * g + 6];
                v4rholapl3[8 * g + 7] +=  c * stage_v4rholapl3[8 * g + 7];

                v4rholapl2tau[12 * g + 0] +=  c * stage_v4rholapl2tau[12 * g + 0];
                v4rholapl2tau[12 * g + 1] +=  c * stage_v4rholapl2tau[12 * g + 1];
                v4rholapl2tau[12 * g + 2] +=  c * stage_v4rholapl2tau[12 * g + 2];
                v4rholapl2tau[12 * g + 3] +=  c * stage_v4rholapl2tau[12 * g + 3];
                v4rholapl2tau[12 * g + 4] +=  c * stage_v4rholapl2tau[12 * g + 4];
                v4rholapl2tau[12 * g + 5] +=  c * stage_v4rholapl2tau[12 * g + 5];
                v4rholapl2tau[12 * g + 6] +=  c * stage_v4rholapl2tau[12 * g + 6];
                v4rholapl2tau[12 * g + 7] +=  c * stage_v4rholapl2tau[12 * g + 7];
                v4rholapl2tau[12 * g + 8] +=  c * stage_v4rholapl2tau[12 * g + 8];
                v4rholapl2tau[12 * g + 9] +=  c * stage_v4rholapl2tau[12 * g + 9];
                v4rholapl2tau[12 * g + 10] +=  c * stage_v4rholapl2tau[12 * g + 10];
                v4rholapl2tau[12 * g + 11] +=  c * stage_v4rholapl2tau[12 * g + 11];

                v4rholapltau2[12 * g + 0] +=  c * stage_v4rholapltau2[12 * g + 0];
                v4rholapltau2[12 * g + 1] +=  c * stage_v4rholapltau2[12 * g + 1];
                v4rholapltau2[12 * g + 2] +=  c * stage_v4rholapltau2[12 * g + 2];
                v4rholapltau2[12 * g + 3] +=  c * stage_v4rholapltau2[12 * g + 3];
                v4rholapltau2[12 * g + 4] +=  c * stage_v4rholapltau2[12 * g + 4];
                v4rholapltau2[12 * g + 5] +=  c * stage_v4rholapltau2[12 * g + 5];
                v4rholapltau2[12 * g + 6] +=  c * stage_v4rholapltau2[12 * g + 6];
                v4rholapltau2[12 * g + 7] +=  c * stage_v4rholapltau2[12 * g + 7];
                v4rholapltau2[12 * g + 8] +=  c * stage_v4rholapltau2[12 * g + 8];
                v4rholapltau2[12 * g + 9] +=  c * stage_v4rholapltau2[12 * g + 9];
                v4rholapltau2[12 * g + 10] +=  c * stage_v4rholapltau2[12 * g + 10];
                v4rholapltau2[12 * g + 11] +=  c * stage_v4rholapltau2[12 * g + 11];

                v4rhotau3[8 * g + 0] +=  c * stage_v4rhotau3[8 * g + 0];
                v4rhotau3[8 * g + 1] +=  c * stage_v4rhotau3[8 * g + 1];
                v4rhotau3[8 * g + 2] +=  c * stage_v4rhotau3[8 * g + 2];
                v4rhotau3[8 * g + 3] +=  c * stage_v4rhotau3[8 * g + 3];
                v4rhotau3[8 * g + 4] +=  c * stage_v4rhotau3[8 * g + 4];
                v4rhotau3[8 * g + 5] +=  c * stage_v4rhotau3[8 * g + 5];
                v4rhotau3[8 * g + 6] +=  c * stage_v4rhotau3[8 * g + 6];
                v4rhotau3[8 * g + 7] +=  c * stage_v4rhotau3[8 * g + 7];

                v4sigma4[15 * g + 0] +=  c * stage_v4sigma4[15 * g + 0];
                v4sigma4[15 * g + 1] +=  c * stage_v4sigma4[15 * g + 1];
                v4sigma4[15 * g + 2] +=  c * stage_v4sigma4[15 * g + 2];
                v4sigma4[15 * g + 3] +=  c * stage_v4sigma4[15 * g + 3];
                v4sigma4[15 * g + 4] +=  c * stage_v4sigma4[15 * g + 4];
                v4sigma4[15 * g + 5] +=  c * stage_v4sigma4[15 * g + 5];
                v4sigma4[15 * g + 6] +=  c * stage_v4sigma4[15 * g + 6];
                v4sigma4[15 * g + 7] +=  c * stage_v4sigma4[15 * g + 7];
                v4sigma4[15 * g + 8] +=  c * stage_v4sigma4[15 * g + 8];
                v4sigma4[15 * g + 9] +=  c * stage_v4sigma4[15 * g + 9];
                v4sigma4[15 * g + 10] +=  c * stage_v4sigma4[15 * g + 10];
                v4sigma4[15 * g + 11] +=  c * stage_v4sigma4[15 * g + 11];
                v4sigma4[15 * g + 12] +=  c * stage_v4sigma4[15 * g + 12];
                v4sigma4[15 * g + 13] +=  c * stage_v4sigma4[15 * g + 13];
                v4sigma4[15 * g + 14] +=  c * stage_v4sigma4[15 * g + 14];

                v4sigma3lapl[20 * g + 0] +=  c * stage_v4sigma3lapl[20 * g + 0];
                v4sigma3lapl[20 * g + 1] +=  c * stage_v4sigma3lapl[20 * g + 1];
                v4sigma3lapl[20 * g + 2] +=  c * stage_v4sigma3lapl[20 * g + 2];
                v4sigma3lapl[20 * g + 3] +=  c * stage_v4sigma3lapl[20 * g + 3];
                v4sigma3lapl[20 * g + 4] +=  c * stage_v4sigma3lapl[20 * g + 4];
                v4sigma3lapl[20 * g + 5] +=  c * stage_v4sigma3lapl[20 * g + 5];
                v4sigma3lapl[20 * g + 6] +=  c * stage_v4sigma3lapl[20 * g + 6];
                v4sigma3lapl[20 * g + 7] +=  c * stage_v4sigma3lapl[20 * g + 7];
                v4sigma3lapl[20 * g + 8] +=  c * stage_v4sigma3lapl[20 * g + 8];
                v4sigma3lapl[20 * g + 9] +=  c * stage_v4sigma3lapl[20 * g + 9];
                v4sigma3lapl[20 * g + 10] +=  c * stage_v4sigma3lapl[20 * g + 10];
                v4sigma3lapl[20 * g + 11] +=  c * stage_v4sigma3lapl[20 * g + 11];
                v4sigma3lapl[20 * g + 12] +=  c * stage_v4sigma3lapl[20 * g + 12];
                v4sigma3lapl[20 * g + 13] +=  c * stage_v4sigma3lapl[20 * g + 13];
                v4sigma3lapl[20 * g + 14] +=  c * stage_v4sigma3lapl[20 * g + 14];
                v4sigma3lapl[20 * g + 15] +=  c * stage_v4sigma3lapl[20 * g + 15];
                v4sigma3lapl[20 * g + 16] +=  c * stage_v4sigma3lapl[20 * g + 16];
                v4sigma3lapl[20 * g + 17] +=  c * stage_v4sigma3lapl[20 * g + 17];
                v4sigma3lapl[20 * g + 18] +=  c * stage_v4sigma3lapl[20 * g + 18];
                v4sigma3lapl[20 * g + 19] +=  c * stage_v4sigma3lapl[20 * g + 19];

                v4sigma3tau[30 * g + 0] +=  c * stage_v4sigma3tau[30 * g + 0];
                v4sigma3tau[30 * g + 1] +=  c * stage_v4sigma3tau[30 * g + 1];
                v4sigma3tau[30 * g + 2] +=  c * stage_v4sigma3tau[30 * g + 2];
                v4sigma3tau[30 * g + 3] +=  c * stage_v4sigma3tau[30 * g + 3];
                v4sigma3tau[30 * g + 4] +=  c * stage_v4sigma3tau[30 * g + 4];
                v4sigma3tau[30 * g + 5] +=  c * stage_v4sigma3tau[30 * g + 5];
                v4sigma3tau[30 * g + 6] +=  c * stage_v4sigma3tau[30 * g + 6];
                v4sigma3tau[30 * g + 7] +=  c * stage_v4sigma3tau[30 * g + 7];
                v4sigma3tau[30 * g + 8] +=  c * stage_v4sigma3tau[30 * g + 8];
                v4sigma3tau[30 * g + 9] +=  c * stage_v4sigma3tau[30 * g + 9];
                v4sigma3tau[30 * g + 10] +=  c * stage_v4sigma3tau[30 * g + 10];
                v4sigma3tau[30 * g + 11] +=  c * stage_v4sigma3tau[30 * g + 11];
                v4sigma3tau[30 * g + 12] +=  c * stage_v4sigma3tau[30 * g + 12];
                v4sigma3tau[30 * g + 13] +=  c * stage_v4sigma3tau[30 * g + 13];
                v4sigma3tau[30 * g + 14] +=  c * stage_v4sigma3tau[30 * g + 14];
                v4sigma3tau[30 * g + 15] +=  c * stage_v4sigma3tau[30 * g + 15];
                v4sigma3tau[30 * g + 16] +=  c * stage_v4sigma3tau[30 * g + 16];
                v4sigma3tau[30 * g + 17] +=  c * stage_v4sigma3tau[30 * g + 17];
                v4sigma3tau[30 * g + 18] +=  c * stage_v4sigma3tau[30 * g + 18];
                v4sigma3tau[30 * g + 19] +=  c * stage_v4sigma3tau[30 * g + 19];

                v4sigma2lapl2[18 * g + 0] +=  c * stage_v4sigma2lapl2[18 * g + 0];
                v4sigma2lapl2[18 * g + 1] +=  c * stage_v4sigma2lapl2[18 * g + 1];
                v4sigma2lapl2[18 * g + 2] +=  c * stage_v4sigma2lapl2[18 * g + 2];
                v4sigma2lapl2[18 * g + 3] +=  c * stage_v4sigma2lapl2[18 * g + 3];
                v4sigma2lapl2[18 * g + 4] +=  c * stage_v4sigma2lapl2[18 * g + 4];
                v4sigma2lapl2[18 * g + 5] +=  c * stage_v4sigma2lapl2[18 * g + 5];
                v4sigma2lapl2[18 * g + 6] +=  c * stage_v4sigma2lapl2[18 * g + 6];
                v4sigma2lapl2[18 * g + 7] +=  c * stage_v4sigma2lapl2[18 * g + 7];
                v4sigma2lapl2[18 * g + 8] +=  c * stage_v4sigma2lapl2[18 * g + 8];
                v4sigma2lapl2[18 * g + 9] +=  c * stage_v4sigma2lapl2[18 * g + 9];
                v4sigma2lapl2[18 * g + 10] +=  c * stage_v4sigma2lapl2[18 * g + 10];
                v4sigma2lapl2[18 * g + 11] +=  c * stage_v4sigma2lapl2[18 * g + 11];
                v4sigma2lapl2[18 * g + 12] +=  c * stage_v4sigma2lapl2[18 * g + 12];
                v4sigma2lapl2[18 * g + 13] +=  c * stage_v4sigma2lapl2[18 * g + 13];
                v4sigma2lapl2[18 * g + 14] +=  c * stage_v4sigma2lapl2[18 * g + 14];
                v4sigma2lapl2[18 * g + 15] +=  c * stage_v4sigma2lapl2[18 * g + 15];
                v4sigma2lapl2[18 * g + 16] +=  c * stage_v4sigma2lapl2[18 * g + 16];
                v4sigma2lapl2[18 * g + 17] +=  c * stage_v4sigma2lapl2[18 * g + 17];

                v4sigma2lapltau[24 * g + 0] +=  c * stage_v4sigma2lapltau[24 * g + 0];
                v4sigma2lapltau[24 * g + 1] +=  c * stage_v4sigma2lapltau[24 * g + 1];
                v4sigma2lapltau[24 * g + 2] +=  c * stage_v4sigma2lapltau[24 * g + 2];
                v4sigma2lapltau[24 * g + 3] +=  c * stage_v4sigma2lapltau[24 * g + 3];
                v4sigma2lapltau[24 * g + 4] +=  c * stage_v4sigma2lapltau[24 * g + 4];
                v4sigma2lapltau[24 * g + 5] +=  c * stage_v4sigma2lapltau[24 * g + 5];
                v4sigma2lapltau[24 * g + 6] +=  c * stage_v4sigma2lapltau[24 * g + 6];
                v4sigma2lapltau[24 * g + 7] +=  c * stage_v4sigma2lapltau[24 * g + 7];
                v4sigma2lapltau[24 * g + 8] +=  c * stage_v4sigma2lapltau[24 * g + 8];
                v4sigma2lapltau[24 * g + 9] +=  c * stage_v4sigma2lapltau[24 * g + 9];
                v4sigma2lapltau[24 * g + 10] +=  c * stage_v4sigma2lapltau[24 * g + 10];
                v4sigma2lapltau[24 * g + 11] +=  c * stage_v4sigma2lapltau[24 * g + 11];
                v4sigma2lapltau[24 * g + 12] +=  c * stage_v4sigma2lapltau[24 * g + 12];
                v4sigma2lapltau[24 * g + 13] +=  c * stage_v4sigma2lapltau[24 * g + 13];
                v4sigma2lapltau[24 * g + 14] +=  c * stage_v4sigma2lapltau[24 * g + 14];
                v4sigma2lapltau[24 * g + 15] +=  c * stage_v4sigma2lapltau[24 * g + 15];
                v4sigma2lapltau[24 * g + 16] +=  c * stage_v4sigma2lapltau[24 * g + 16];
                v4sigma2lapltau[24 * g + 17] +=  c * stage_v4sigma2lapltau[24 * g + 17];
                v4sigma2lapltau[24 * g + 18] +=  c * stage_v4sigma2lapltau[24 * g + 18];
                v4sigma2lapltau[24 * g + 19] +=  c * stage_v4sigma2lapltau[24 * g + 19];
                v4sigma2lapltau[24 * g + 20] +=  c * stage_v4sigma2lapltau[24 * g + 20];
                v4sigma2lapltau[24 * g + 21] +=  c * stage_v4sigma2lapltau[24 * g + 21];
                v4sigma2lapltau[24 * g + 22] +=  c * stage_v4sigma2lapltau[24 * g + 22];
                v4sigma2lapltau[24 * g + 23] +=  c * stage_v4sigma2lapltau[24 * g + 23];

                v4sigma2tau2[18 * g + 0] +=  c * stage_v4sigma2tau2[18 * g + 0];
                v4sigma2tau2[18 * g + 1] +=  c * stage_v4sigma2tau2[18 * g + 1];
                v4sigma2tau2[18 * g + 2] +=  c * stage_v4sigma2tau2[18 * g + 2];
                v4sigma2tau2[18 * g + 3] +=  c * stage_v4sigma2tau2[18 * g + 3];
                v4sigma2tau2[18 * g + 4] +=  c * stage_v4sigma2tau2[18 * g + 4];
                v4sigma2tau2[18 * g + 5] +=  c * stage_v4sigma2tau2[18 * g + 5];
                v4sigma2tau2[18 * g + 6] +=  c * stage_v4sigma2tau2[18 * g + 6];
                v4sigma2tau2[18 * g + 7] +=  c * stage_v4sigma2tau2[18 * g + 7];
                v4sigma2tau2[18 * g + 8] +=  c * stage_v4sigma2tau2[18 * g + 8];
                v4sigma2tau2[18 * g + 9] +=  c * stage_v4sigma2tau2[18 * g + 9];
                v4sigma2tau2[18 * g + 10] +=  c * stage_v4sigma2tau2[18 * g + 10];
                v4sigma2tau2[18 * g + 11] +=  c * stage_v4sigma2tau2[18 * g + 11];
                v4sigma2tau2[18 * g + 12] +=  c * stage_v4sigma2tau2[18 * g + 12];
                v4sigma2tau2[18 * g + 13] +=  c * stage_v4sigma2tau2[18 * g + 13];
                v4sigma2tau2[18 * g + 14] +=  c * stage_v4sigma2tau2[18 * g + 14];
                v4sigma2tau2[18 * g + 15] +=  c * stage_v4sigma2tau2[18 * g + 15];
                v4sigma2tau2[18 * g + 16] +=  c * stage_v4sigma2tau2[18 * g + 16];
                v4sigma2tau2[18 * g + 17] +=  c * stage_v4sigma2tau2[18 * g + 17];

                v4sigmalapl3[12 * g + 0] +=  c * stage_v4sigmalapl3[12 * g + 0];
                v4sigmalapl3[12 * g + 1] +=  c * stage_v4sigmalapl3[12 * g + 1];
                v4sigmalapl3[12 * g + 2] +=  c * stage_v4sigmalapl3[12 * g + 2];
                v4sigmalapl3[12 * g + 3] +=  c * stage_v4sigmalapl3[12 * g + 3];
                v4sigmalapl3[12 * g + 4] +=  c * stage_v4sigmalapl3[12 * g + 4];
                v4sigmalapl3[12 * g + 5] +=  c * stage_v4sigmalapl3[12 * g + 5];
                v4sigmalapl3[12 * g + 6] +=  c * stage_v4sigmalapl3[12 * g + 6];
                v4sigmalapl3[12 * g + 7] +=  c * stage_v4sigmalapl3[12 * g + 7];
                v4sigmalapl3[12 * g + 8] +=  c * stage_v4sigmalapl3[12 * g + 8];
                v4sigmalapl3[12 * g + 9] +=  c * stage_v4sigmalapl3[12 * g + 9];
                v4sigmalapl3[12 * g + 10] +=  c * stage_v4sigmalapl3[12 * g + 10];
                v4sigmalapl3[12 * g + 11] +=  c * stage_v4sigmalapl3[12 * g + 11];

                v4sigmalapl2tau[18 * g + 0] +=  c * stage_v4sigmalapl2tau[18 * g + 0];
                v4sigmalapl2tau[18 * g + 1] +=  c * stage_v4sigmalapl2tau[18 * g + 1];
                v4sigmalapl2tau[18 * g + 2] +=  c * stage_v4sigmalapl2tau[18 * g + 2];
                v4sigmalapl2tau[18 * g + 3] +=  c * stage_v4sigmalapl2tau[18 * g + 3];
                v4sigmalapl2tau[18 * g + 4] +=  c * stage_v4sigmalapl2tau[18 * g + 4];
                v4sigmalapl2tau[18 * g + 5] +=  c * stage_v4sigmalapl2tau[18 * g + 5];
                v4sigmalapl2tau[18 * g + 6] +=  c * stage_v4sigmalapl2tau[18 * g + 6];
                v4sigmalapl2tau[18 * g + 7] +=  c * stage_v4sigmalapl2tau[18 * g + 7];
                v4sigmalapl2tau[18 * g + 8] +=  c * stage_v4sigmalapl2tau[18 * g + 8];
                v4sigmalapl2tau[18 * g + 9] +=  c * stage_v4sigmalapl2tau[18 * g + 9];
                v4sigmalapl2tau[18 * g + 10] +=  c * stage_v4sigmalapl2tau[18 * g + 10];
                v4sigmalapl2tau[18 * g + 11] +=  c * stage_v4sigmalapl2tau[18 * g + 11];
                v4sigmalapl2tau[18 * g + 12] +=  c * stage_v4sigmalapl2tau[18 * g + 12];
                v4sigmalapl2tau[18 * g + 13] +=  c * stage_v4sigmalapl2tau[18 * g + 13];
                v4sigmalapl2tau[18 * g + 14] +=  c * stage_v4sigmalapl2tau[18 * g + 14];
                v4sigmalapl2tau[18 * g + 15] +=  c * stage_v4sigmalapl2tau[18 * g + 15];
                v4sigmalapl2tau[18 * g + 16] +=  c * stage_v4sigmalapl2tau[18 * g + 16];
                v4sigmalapl2tau[18 * g + 17] +=  c * stage_v4sigmalapl2tau[18 * g + 17];

                v4sigmalapltau2[18 * g + 0] +=  c * stage_v4sigmalapltau2[18 * g + 0];
                v4sigmalapltau2[18 * g + 1] +=  c * stage_v4sigmalapltau2[18 * g + 1];
                v4sigmalapltau2[18 * g + 2] +=  c * stage_v4sigmalapltau2[18 * g + 2];
                v4sigmalapltau2[18 * g + 3] +=  c * stage_v4sigmalapltau2[18 * g + 3];
                v4sigmalapltau2[18 * g + 4] +=  c * stage_v4sigmalapltau2[18 * g + 4];
                v4sigmalapltau2[18 * g + 5] +=  c * stage_v4sigmalapltau2[18 * g + 5];
                v4sigmalapltau2[18 * g + 6] +=  c * stage_v4sigmalapltau2[18 * g + 6];
                v4sigmalapltau2[18 * g + 7] +=  c * stage_v4sigmalapltau2[18 * g + 7];
                v4sigmalapltau2[18 * g + 8] +=  c * stage_v4sigmalapltau2[18 * g + 8];
                v4sigmalapltau2[18 * g + 9] +=  c * stage_v4sigmalapltau2[18 * g + 9];
                v4sigmalapltau2[18 * g + 10] +=  c * stage_v4sigmalapltau2[18 * g + 10];
                v4sigmalapltau2[18 * g + 11] +=  c * stage_v4sigmalapltau2[18 * g + 11];
                v4sigmalapltau2[18 * g + 12] +=  c * stage_v4sigmalapltau2[18 * g + 12];
                v4sigmalapltau2[18 * g + 13] +=  c * stage_v4sigmalapltau2[18 * g + 13];
                v4sigmalapltau2[18 * g + 14] +=  c * stage_v4sigmalapltau2[18 * g + 14];
                v4sigmalapltau2[18 * g + 15] +=  c * stage_v4sigmalapltau2[18 * g + 15];
                v4sigmalapltau2[18 * g + 16] +=  c * stage_v4sigmalapltau2[18 * g + 16];
                v4sigmalapltau2[18 * g + 17] +=  c * stage_v4sigmalapltau2[18 * g + 17];

                v4sigmatau3[12 * g + 0] +=  c * stage_v4sigmatau3[12 * g + 0];
                v4sigmatau3[12 * g + 1] +=  c * stage_v4sigmatau3[12 * g + 1];
                v4sigmatau3[12 * g + 2] +=  c * stage_v4sigmatau3[12 * g + 2];
                v4sigmatau3[12 * g + 3] +=  c * stage_v4sigmatau3[12 * g + 3];
                v4sigmatau3[12 * g + 4] +=  c * stage_v4sigmatau3[12 * g + 4];
                v4sigmatau3[12 * g + 5] +=  c * stage_v4sigmatau3[12 * g + 5];
                v4sigmatau3[12 * g + 6] +=  c * stage_v4sigmatau3[12 * g + 6];
                v4sigmatau3[12 * g + 7] +=  c * stage_v4sigmatau3[12 * g + 7];
                v4sigmatau3[12 * g + 8] +=  c * stage_v4sigmatau3[12 * g + 8];
                v4sigmatau3[12 * g + 9] +=  c * stage_v4sigmatau3[12 * g + 9];
                v4sigmatau3[12 * g + 10] +=  c * stage_v4sigmatau3[12 * g + 10];
                v4sigmatau3[12 * g + 11] +=  c * stage_v4sigmatau3[12 * g + 11];

                v4lapl4[5 * g + 0] +=  c * stage_v4lapl4[5 * g + 0];
                v4lapl4[5 * g + 1] +=  c * stage_v4lapl4[5 * g + 1];
                v4lapl4[5 * g + 2] +=  c * stage_v4lapl4[5 * g + 2];
                v4lapl4[5 * g + 3] +=  c * stage_v4lapl4[5 * g + 3];
                v4lapl4[5 * g + 4] +=  c * stage_v4lapl4[5 * g + 4];

                v4lapl3tau[8 * g + 0] +=  c * stage_v4lapl3tau[8 * g + 0];
                v4lapl3tau[8 * g + 1] +=  c * stage_v4lapl3tau[8 * g + 1];
                v4lapl3tau[8 * g + 2] +=  c * stage_v4lapl3tau[8 * g + 2];
                v4lapl3tau[8 * g + 3] +=  c * stage_v4lapl3tau[8 * g + 3];
                v4lapl3tau[8 * g + 4] +=  c * stage_v4lapl3tau[8 * g + 4];
                v4lapl3tau[8 * g + 5] +=  c * stage_v4lapl3tau[8 * g + 5];
                v4lapl3tau[8 * g + 6] +=  c * stage_v4lapl3tau[8 * g + 6];
                v4lapl3tau[8 * g + 7] +=  c * stage_v4lapl3tau[8 * g + 7];

                v4lapl2tau2[9 * g + 0] +=  c * stage_v4lapl2tau2[9 * g + 0];
                v4lapl2tau2[9 * g + 1] +=  c * stage_v4lapl2tau2[9 * g + 1];
                v4lapl2tau2[9 * g + 2] +=  c * stage_v4lapl2tau2[9 * g + 2];
                v4lapl2tau2[9 * g + 3] +=  c * stage_v4lapl2tau2[9 * g + 3];
                v4lapl2tau2[9 * g + 4] +=  c * stage_v4lapl2tau2[9 * g + 4];
                v4lapl2tau2[9 * g + 5] +=  c * stage_v4lapl2tau2[9 * g + 5];
                v4lapl2tau2[9 * g + 6] +=  c * stage_v4lapl2tau2[9 * g + 6];
                v4lapl2tau2[9 * g + 7] +=  c * stage_v4lapl2tau2[9 * g + 7];
                v4lapl2tau2[9 * g + 8] +=  c * stage_v4lapl2tau2[9 * g + 8];

                v4lapltau3[8 * g + 0] +=  c * stage_v4lapltau3[8 * g + 0];
                v4lapltau3[8 * g + 1] +=  c * stage_v4lapltau3[8 * g + 1];
                v4lapltau3[8 * g + 2] +=  c * stage_v4lapltau3[8 * g + 2];
                v4lapltau3[8 * g + 3] +=  c * stage_v4lapltau3[8 * g + 3];
                v4lapltau3[8 * g + 4] +=  c * stage_v4lapltau3[8 * g + 4];
                v4lapltau3[8 * g + 5] +=  c * stage_v4lapltau3[8 * g + 5];
                v4lapltau3[8 * g + 6] +=  c * stage_v4lapltau3[8 * g + 6];
                v4lapltau3[8 * g + 7] +=  c * stage_v4lapltau3[8 * g + 7];

                v4tau4[5 * g + 0] +=  c * stage_v4tau4[5 * g + 0];
                v4tau4[5 * g + 1] +=  c * stage_v4tau4[5 * g + 1];
                v4tau4[5 * g + 2] +=  c * stage_v4tau4[5 * g + 2];
                v4tau4[5 * g + 3] +=  c * stage_v4tau4[5 * g + 3];
                v4tau4[5 * g + 4] +=  c * stage_v4tau4[5 * g + 4];

               
            }

        }        
    }

    if (alloc)
    {
        mem::free(stage_v4rho4);
        mem::free(stage_v4rho3sigma);
        mem::free(stage_v4rho3lapl);
        mem::free(stage_v4rho3tau);
        mem::free(stage_v4rho2sigma2);
        mem::free(stage_v4rho2sigmalapl);
        mem::free(stage_v4rho2sigmatau);
        mem::free(stage_v4rho2lapl2);
        mem::free(stage_v4rho2lapltau);
        mem::free(stage_v4rho2tau2);
        mem::free(stage_v4rhosigma3);
        mem::free(stage_v4rhosigma2lapl);
        mem::free(stage_v4rhosigma2tau);
        mem::free(stage_v4rhosigmalapl2);
        mem::free(stage_v4rhosigmalapltau);
        mem::free(stage_v4rhosigmatau2);
        mem::free(stage_v4rholapl3);
        mem::free(stage_v4rholapl2tau);
        mem::free(stage_v4rholapltau2);
        mem::free(stage_v4rhotau3);
        mem::free(stage_v4sigma4);
        mem::free(stage_v4sigma3lapl);
        mem::free(stage_v4sigma3tau);
        mem::free(stage_v4sigma2lapl2);
        mem::free(stage_v4sigma2lapltau);
        mem::free(stage_v4sigma2tau2);
        mem::free(stage_v4sigmalapl3);
        mem::free(stage_v4sigmalapl2tau);
        mem::free(stage_v4sigmalapltau2);
        mem::free(stage_v4sigmatau3);
        mem::free(stage_v4lapl4);
        mem::free(stage_v4lapl3tau);
        mem::free(stage_v4lapl2tau2);
        mem::free(stage_v4lapltau3);
        mem::free(stage_v4tau4);
    }
}
