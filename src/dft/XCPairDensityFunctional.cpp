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

#include "XCPairDensityFunctional.hpp"

#include "ErrorHandler.hpp"
#include "MemAlloc.hpp"
#include "PairDensityBecke88.hpp"
#include "PairDensityLYP.hpp"
#include "PairDensityPBE_C.hpp"
#include "PairDensityPBE_X.hpp"
#include "PairDensitySlater.hpp"
#include "PairDensityVWN.hpp"
#include "StringFormat.hpp"

CXCPairDensityFunctional::CXCPairDensityFunctional(const std::string&              nameOfFunctional,
                                                   const std::vector<std::string>& labels,
                                                   const std::vector<double>&      coeffs)

    : _nameOfFunctional(fstr::upcase(nameOfFunctional))
{
    std::string errmsg("XCPairDensityFunctional: Inconsistent sizes of functional labels and coefficients");

    errors::assertMsgCritical(labels.size() == coeffs.size(), errmsg);

    _components.clear();

    bool isPLDA = false, isPGGA = false;

    for (int32_t i = 0; i < static_cast<int32_t>(labels.size()); i++)
    {
        _components.push_back(std::make_tuple(labels[i], coeffs[i]));

        isPGGA = (isPGGA || _isComponentPGGA(labels[i]));
        isPLDA = ((!isPGGA) && (isPLDA || _isComponentPLDA(labels[i])));
    }

    if (isPLDA) _familyOfFunctional = std::string("PLDA");
    if (isPGGA) _familyOfFunctional = std::string("PGGA");

    _allocateStagingBuffer();
}

CXCPairDensityFunctional::CXCPairDensityFunctional(const CXCPairDensityFunctional& source)

    : _nameOfFunctional(source._nameOfFunctional)

    , _familyOfFunctional(source._familyOfFunctional)

    , _ldStaging(source._ldStaging)

    , _components(source._components)
{
    _allocateStagingBuffer();
}

CXCPairDensityFunctional::CXCPairDensityFunctional(CXCPairDensityFunctional&& source) noexcept

    : _nameOfFunctional(std::move(source._nameOfFunctional))

    , _familyOfFunctional(std::move(source._familyOfFunctional))

    , _ldStaging(std::move(source._ldStaging))

    , _components(std::move(source._components))
{
    _allocateStagingBuffer();

    source._freeStagingBuffer();
}

CXCPairDensityFunctional::~CXCPairDensityFunctional()
{
    _components.clear();

    _freeStagingBuffer();
}

void
CXCPairDensityFunctional::_allocateStagingBuffer()
{
    if (_stagingBuffer == nullptr)
    {
        // TODO write function to compute this number based on functional
        // family and order of derivatives available
        int32_t n_xc_outputs = 0;

        if (_familyOfFunctional == std::string("PLDA")) n_xc_outputs = 3;

        if (_familyOfFunctional == std::string("PGGA")) n_xc_outputs = 6;

        if (n_xc_outputs > 0) _stagingBuffer = mem::malloc<double>(n_xc_outputs * _ldStaging);
    }
}

void
CXCPairDensityFunctional::_freeStagingBuffer()
{
    if (_stagingBuffer != nullptr)
    {
        mem::free(_stagingBuffer);

        _stagingBuffer = nullptr;
    }
}

bool
CXCPairDensityFunctional::_isComponentPLDA(const std::string& compName) const
{
    auto upcasename = fstr::upcase(compName);

    if (upcasename == "TSLATER") return true;

    if (upcasename == "TVWN") return true;

    return false;
}

bool
CXCPairDensityFunctional::_isComponentPGGA(const std::string& compName) const
{
    auto upcasename = fstr::upcase(compName);

    if (upcasename == "TPBE_X") return true;

    if (upcasename == "TPBE_C") return true;

    if (upcasename == "TB88") return true;

    if (upcasename == "TLYP") return true;

    return false;
}

void
CXCPairDensityFunctional::_plda_exc_vxc(const std::string& compName, const int32_t np, const double* rho, double* exc, double* vrho) const
{
    auto upcasename = fstr::upcase(compName);

    if (upcasename == "TSLATER")
    {
        pdftslater::compute_exc_vxc(np, rho, exc, vrho);
    }
    else if (upcasename == "TVWN")
    {
        pdftvwn::compute_exc_vxc(np, rho, exc, vrho);
    }
    else
    {
        std::string errmsg("XCPairDensityFunctional._plda_exc_vxc: Invalid functional name");

        errors::assertMsgCritical(false, errmsg);
    }
}

void
CXCPairDensityFunctional::_pgga_exc_vxc(const std::string& compName,
                                        const int32_t      np,
                                        const double*      rho,
                                        const double*      sigma,
                                        double*            exc,
                                        double*            vrho,
                                        double*            vsigma) const
{
    auto upcasename = fstr::upcase(compName);

    if (upcasename == "TPBE_X")
    {
        pdftpbe_x::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else if (upcasename == "TPBE_C")
    {
        pdftpbe_c::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else if (upcasename == "TB88")
    {
        pdftb88::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else if (upcasename == "TLYP")
    {
        pdftlyp::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else
    {
        std::string errmsg("XCPairDensityFunctional._pgga_exc_vxc: Invalid functional name");

        errors::assertMsgCritical(false, errmsg);
    }
}

CXCPairDensityFunctional&
CXCPairDensityFunctional::operator=(const CXCPairDensityFunctional& source)
{
    if (this == &source) return *this;

    _nameOfFunctional = source._nameOfFunctional;

    _familyOfFunctional = source._familyOfFunctional;

    _ldStaging = source._ldStaging;

    _components = source._components;

    _freeStagingBuffer();

    _allocateStagingBuffer();

    return *this;
}

CXCPairDensityFunctional&
CXCPairDensityFunctional::operator=(CXCPairDensityFunctional&& source) noexcept
{
    if (this == &source) return *this;

    _nameOfFunctional = std::move(source._nameOfFunctional);

    _familyOfFunctional = std::move(source._familyOfFunctional);

    _ldStaging = std::move(source._ldStaging);

    _components = std::move(source._components);

    _freeStagingBuffer();

    _allocateStagingBuffer();

    source._freeStagingBuffer();

    return *this;
}

bool
CXCPairDensityFunctional::operator==(const CXCPairDensityFunctional& other) const
{
    if (_nameOfFunctional != other._nameOfFunctional) return false;

    if (_familyOfFunctional != other._familyOfFunctional) return false;

    if (_ldStaging != other._ldStaging) return false;

    if (_components != other._components) return false;

    return true;
}

bool
CXCPairDensityFunctional::operator!=(const CXCPairDensityFunctional& other) const
{
    return !(*this == other);
}

std::string
CXCPairDensityFunctional::getFunctionalLabel() const
{
    return _nameOfFunctional;
}

std::string
CXCPairDensityFunctional::getFunctionalType() const
{
    return _familyOfFunctional;
}

auto
CXCPairDensityFunctional::compute_exc_vxc_for_plda(int32_t np, const double* rho, double* exc, double* vrho) const -> void
{
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

    for (const auto& comp : _components)
    {
        const auto funcname = std::get<0>(comp);

        const auto c = std::get<1>(comp);

        if (_isComponentPLDA(funcname))
        {
            _plda_exc_vxc(funcname, np, rho, stage_exc, stage_vrho);

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                exc[g] += c * stage_exc[g];

                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_exc);
        mem::free(stage_vrho);
    }
}

auto
CXCPairDensityFunctional::compute_exc_vxc_for_pgga(int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma)
    const -> void
{
    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? mem::malloc<double>(1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

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

    for (const auto& comp : _components)
    {
        const auto funcname = std::get<0>(comp);

        const auto c = std::get<1>(comp);

        if (_isComponentPLDA(funcname))
        {
            _plda_exc_vxc(funcname, np, rho, stage_exc, stage_vrho);

#pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                exc[g] += c * stage_exc[g];

                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (_isComponentPGGA(funcname))
        {
            _pgga_exc_vxc(funcname, np, rho, sigma, stage_exc, stage_vrho, stage_vsigma);

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
