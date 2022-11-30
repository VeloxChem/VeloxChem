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

#include <xc.h>

#include <cstdint>
#include <string>

#include "ErrorHandler.hpp"
#include "MemAlloc.hpp"
#include "PrimitiveFunctional.hpp"

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
