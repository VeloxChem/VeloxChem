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
#include <sstream>
#include <string>

#include "ErrorHandler.hpp"
#include "MemAlloc.hpp"

CXCNewFunctional::CXCNewFunctional(const std::vector<std::string>& labels, const std::vector<double>& coeffs, const double fractionOfExactExchange)

    : _fractionOfExactExchange(fractionOfExactExchange)
{
    std::string errmsg(std::string(__func__) + ": Inconsistent sizes of functional labels and coefficients");

    errors::assertMsgCritical(labels.size() == coeffs.size(), errmsg);

    auto is_in_family = [](const auto& options, auto family) {
        return std::any_of(options.begin(), options.end(), [family](auto a) { return a == family; });
    };

    _components.clear();

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
            _components.push_back({coeff, func});
        }
        else if (kind == XC_CORRELATION)
        {
            _components.push_back({coeff, func});
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

    // allocate _stagingBuffer
    _stagingBuffer = mem::malloc<double>(n_xc_outputs * _ldStaging);
}

CXCNewFunctional::~CXCNewFunctional()
{
    // clean up allocated LibXC objects
    std::for_each(std::begin(_components), std::end(_components), [](Component& c) { xc_func_end(&std::get<1>(c)); });

    _components.clear();

    // clean up staging buffer
    mem::free(_stagingBuffer);
}

auto
CXCNewFunctional::compute_exc_for_lda(int32_t np, const double* rho, double* exc) const -> void
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
CXCNewFunctional::compute_vxc_for_lda(int32_t np, const double* rho, double* vrho) const -> void
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
CXCNewFunctional::compute_exc_vxc_for_lda(int32_t np, const double* rho, double* exc, double* vrho) const -> void
{
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc  = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0];
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
CXCNewFunctional::compute_exc_for_gga(int32_t np, const double* rho, const double* sigma, double* exc) const -> void
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
CXCNewFunctional::compute_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void
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
CXCNewFunctional::compute_exc_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0];
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
CXCNewFunctional::compute_exc_for_mgga(int32_t np, const double* rho, const double* sigma, const double* lapl, const double* tau, double* exc) const
    -> void
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
    errors::assertMsgCritical(_hasExc && _hasVxc,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0];
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
