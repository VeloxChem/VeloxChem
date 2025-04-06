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

#include "XCPairDensityFunctional.hpp"

#include "ErrorHandler.hpp"
#include "PairDensityBecke88.hpp"
#include "PairDensityB88_erf.hpp"
#include "PairDensityLYP.hpp"
#include "PairDensityLYP_erf.hpp"
#include "PairDensityPBE_C.hpp"
#include "PairDensityPBEC_erf.hpp"
#include "PairDensityPBE_X.hpp"
#include "PairDensityPBEX_erf.hpp"
#include "PairDensitySlater.hpp"
#include "PairDensitySlater_erf.hpp"
#include "PairDensityVWN.hpp"
#include "PairDensityPMGB06.hpp"
#include "PairDensityP86.hpp"
#include "PairDensityHPG20.hpp"
#include "StringFormat.hpp"

CXCPairDensityFunctional::CXCPairDensityFunctional(const std::string&              nameOfFunctional,
                                                   const std::vector<std::string>& labels,
                                                   const std::vector<double>&      coeffs)

    : _nameOfFunctional(format::upper_case(nameOfFunctional))
{
    std::string errmsg("XCPairDensityFunctional: Inconsistent sizes of functional labels and coefficients");

    errors::assertMsgCritical(labels.size() == coeffs.size(), errmsg);

    _components.clear();

    bool isPLDA = false, isPGGA = false;

    for (int i = 0; i < static_cast<int>(labels.size()); i++)
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
        int n_xc_outputs = 0;

        if (_familyOfFunctional == std::string("PLDA")) n_xc_outputs = 3;

        if (_familyOfFunctional == std::string("PGGA")) n_xc_outputs = 6;

        if (n_xc_outputs > 0) _stagingBuffer = (double*)std::malloc(sizeof(double) * n_xc_outputs * _ldStaging);
    }
}

void
CXCPairDensityFunctional::_freeStagingBuffer()
{
    if (_stagingBuffer != nullptr)
    {
        std::free(_stagingBuffer);

        _stagingBuffer = nullptr;
    }
}

bool
CXCPairDensityFunctional::_isComponentPLDA(const std::string& compName) const
{
    auto upcasename = format::upper_case(compName);

    if (upcasename == "TSLATER") return true;

    if (upcasename == "TSLATER_ERF") return true;

    if (upcasename == "TVWN_RPA") return true;

    if (upcasename == "TVWN5") return true;

    if (upcasename == "TPMGB06") return true;

    return false;
}

bool
CXCPairDensityFunctional::_isComponentPGGA(const std::string& compName) const
{
    auto upcasename = format::upper_case(compName);

    if (upcasename == "TPBE_X") return true;

    if (upcasename == "TPBE_C") return true;

    if (upcasename == "TB88") return true;

    if (upcasename == "TLYP") return true;

    if (upcasename == "TP86") return true;

    if (upcasename == "HPG20") return true;

    if (upcasename == "TPBEX_ERF") return true;

    if (upcasename == "TPBEC_ERF") return true;

    if (upcasename == "TB88_ERF") return true;

    if (upcasename == "TLYP_ERF") return true;

    return false;
}

void
CXCPairDensityFunctional::_plda_exc_vxc(const std::string& compName, const int np, const double* rho, double* exc, double* vrho, double rs_omega) const
{
    auto upcasename = format::upper_case(compName);

    if (upcasename == "TSLATER")
    {
        pdftslater::compute_exc_vxc(np, rho, exc, vrho);
    }
    else if (upcasename == "TVWN_RPA")
    {
        pdftvwn_rpa::compute_exc_vxc(np, rho, exc, vrho);
    }
    else if (upcasename == "TVWN5")
    {
        pdftvwn5::compute_exc_vxc(np, rho, exc, vrho);
    }
    else if (upcasename == "TSLATER_ERF")
    {
        pdftslater_erf::compute_exc_vxc(np, rho, exc, vrho, rs_omega);
    }
    else if (upcasename == "TPMGB06")
    {
        pdftpmgb06::compute_exc_vxc(np, rho, exc, vrho, rs_omega);
    }
    else
    {
        std::string errmsg("XCPairDensityFunctional._plda_exc_vxc: Invalid functional name "+upcasename);

        errors::assertMsgCritical(false, errmsg);
    }
}

void
CXCPairDensityFunctional::_pgga_exc_vxc(const std::string& compName,
                                        const int      np,
                                        const double*      rho,
                                        const double*      sigma,
                                        double*            exc,
                                        double*            vrho,
                                        double*            vsigma,
                                        double             rs_omega) const
{
    auto upcasename = format::upper_case(compName);

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
    else if (upcasename == "TP86")
    {
        pdftp86::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else if (upcasename == "HPG20")
    {
        pdfthpg20::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma);
    }
    else if (upcasename == "TPBEX_ERF")
    {
        pdftpbex_erf::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma, rs_omega);
    }
    else if (upcasename == "TPBEC_ERF")
    {
        pdftpbec_erf::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma, rs_omega);
    }
    else if (upcasename == "TB88_ERF")
    {
        pdftb88_erf::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma, rs_omega);
    }
    else if (upcasename == "TLYP_ERF")
    {
        pdftlyp_erf::compute_exc_vxc(np, rho, sigma, exc, vrho, vsigma, rs_omega);
    }
    else
    {
        std::string errmsg("XCPairDensityFunctional._pgga_exc_vxc: Invalid functional name "+upcasename);

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
CXCPairDensityFunctional::compute_exc_vxc_for_plda(int np, const double* rho, double* exc, double* vrho, double rs_omega) const -> void
{
    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc  = (alloc) ? (double*)std::malloc(sizeof(double) * 1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho = (alloc) ? (double*)std::malloc(sizeof(double) * 2 * np) : &_stagingBuffer[1 * _ldStaging];

#pragma omp simd
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
            _plda_exc_vxc(funcname, np, rho, stage_exc, stage_vrho, rs_omega);

#pragma omp simd
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
        std::free(stage_exc);
        std::free(stage_vrho);
    }
}

auto
CXCPairDensityFunctional::compute_exc_vxc_for_pgga(int np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma, double rs_omega)
    const -> void
{
    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_exc    = (alloc) ? (double*)std::malloc(sizeof(double) * 1 * np) : &_stagingBuffer[0 * _ldStaging];
    auto stage_vrho   = (alloc) ? (double*)std::malloc(sizeof(double) * 2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? (double*)std::malloc(sizeof(double) * 3 * np) : &_stagingBuffer[3 * _ldStaging];

#pragma omp simd
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
            _plda_exc_vxc(funcname, np, rho, stage_exc, stage_vrho, rs_omega);

#pragma omp simd
            for (auto g = 0; g < np; ++g)
            {
                exc[g] += c * stage_exc[g];

                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (_isComponentPGGA(funcname))
        {
            _pgga_exc_vxc(funcname, np, rho, sigma, stage_exc, stage_vrho, stage_vsigma, rs_omega);

#pragma omp simd
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
        std::free(stage_exc);
        std::free(stage_vrho);
        std::free(stage_vsigma);
    }
}
