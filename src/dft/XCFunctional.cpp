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

#include "XCFunctional.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "ErrorHandler.hpp"
#include "GridScreener.hpp"
#include "MemAlloc.hpp"
#include "StringFormat.hpp"

CXCFunctional::CXCFunctional(const std::string&              nameOfFunctional,
                             const std::vector<std::string>& labels,
                             const std::vector<double>&      coeffs,
                             const double                    fractionOfExactExchange)

    : _nameOfFunctional(fstr::upcase(nameOfFunctional))

    , _fractionOfExactExchange(fractionOfExactExchange)
{
    std::string errmsg("XCFunctional: Inconsistent sizes of functional labels and coefficients");

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
            errors::assertMsgCritical(labels.size() == 1, std::string("XCFunctional: Cannot mix ") + label + std::string(" with other functionals"));

            _components.push_back(xccomp);
        }
        else
        {
            errors::assertMsgCritical(false, std::string("XCFunctional: Unsupported functional ") + label);
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

        // LDA, GGA, metaGGA
        isMGGA = (isMGGA || xccomp.isMetaGGA());
        isGGA  = ((!isMGGA) && (isGGA || xccomp.isGGA()));
        isLDA  = ((!isMGGA) && (!isGGA) && (isLDA || xccomp.isLDA()));

        // check functional hybrid type
        // TODO use xc_hyb_type when it is available

        auto hyb_exx_coeff = xc_hyb_exx_coef(funcptr);

        if (hyb_exx_coeff != 0.0) _fractionOfExactExchange += coeff * hyb_exx_coeff;

        // TODO process range-separation parameters when range-separted functional is supported

        double omega = 0.0, alpha = 0.0, beta = 0.0;

        xc_hyb_cam_coef(funcptr, &omega, &alpha, &beta);

        errors::assertMsgCritical(std::fabs(beta) < 1.0e-13, std::string("XCFunctional: Range-separated functional is not yet supported"));
    }

    if (hasExc) _maxDerivOrder = 0;
    if (hasVxc) _maxDerivOrder = 1;
    if (hasFxc) _maxDerivOrder = 2;
    if (hasKxc) _maxDerivOrder = 3;
    if (hasLxc) _maxDerivOrder = 4;

    if (isLDA) _familyOfFunctional = std::string("LDA");
    if (isGGA) _familyOfFunctional = std::string("GGA");
    if (isMGGA) _familyOfFunctional = std::string("MGGA");

    errors::assertMsgCritical(!needLaplacian, std::string("XCFunctional: Density Laplacian is not supported"));

    _allocateStagingBuffer();
}

CXCFunctional::CXCFunctional(const CXCFunctional& source)

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

CXCFunctional::CXCFunctional(CXCFunctional&& source) noexcept

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

CXCFunctional::~CXCFunctional()
{
    _components.clear();

    _freeStagingBuffer();
}

void
CXCFunctional::_allocateStagingBuffer()
{
    if (_stagingBuffer == nullptr)
    {
        auto n_xc_outputs = getDimensionOfDerivatives();

        _stagingBuffer = mem::malloc<double>(n_xc_outputs * _ldStaging);
    }
}

void
CXCFunctional::_freeStagingBuffer()
{
    if (_stagingBuffer != nullptr)
    {
        mem::free(_stagingBuffer);

        _stagingBuffer = nullptr;
    }
}

std::unordered_map<std::string, std::array<int32_t, 2>>
CXCFunctional::_getIndicesAndCountsOfDerivatives() const
{
    std::unordered_map<std::string, std::array<int32_t, 2>> indices_and_counts;

    if (_familyOfFunctional == std::string("LDA"))
    {
        auto       ldafunc = getFunctionalPointerToLdaComponent();
        const auto dim     = &(ldafunc->dim);

        int n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int32_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int32_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;

        indices_and_counts["v2rho2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;

        indices_and_counts["v3rho3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;

        indices_and_counts["v4rho4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
    }
    else if (_familyOfFunctional == std::string("GGA"))
    {
        auto       ggafunc = getFunctionalPointerToGgaComponent();
        const auto dim     = &(ggafunc->dim);

        int n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int32_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int32_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;
        indices_and_counts["vsigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->vsigma});
        n_xc_outputs += dim->vsigma;

        indices_and_counts["v2rho2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;
        indices_and_counts["v2rhosigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rhosigma});
        n_xc_outputs += dim->v2rhosigma;
        indices_and_counts["v2sigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2sigma2});
        n_xc_outputs += dim->v2sigma2;

        indices_and_counts["v3rho3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;
        indices_and_counts["v3rho2sigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho2sigma});
        n_xc_outputs += dim->v3rho2sigma;
        indices_and_counts["v3rhosigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rhosigma2});
        n_xc_outputs += dim->v3rhosigma2;
        indices_and_counts["v3sigma3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigma3});
        n_xc_outputs += dim->v3sigma3;

        indices_and_counts["v4rho4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
        indices_and_counts["v4rho3sigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho3sigma});
        n_xc_outputs += dim->v4rho3sigma;
        indices_and_counts["v4rho2sigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2sigma2});
        n_xc_outputs += dim->v4rho2sigma2;
        indices_and_counts["v4rhosigma3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigma3});
        n_xc_outputs += dim->v4rhosigma3;
        indices_and_counts["v4sigma4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma4});
        n_xc_outputs += dim->v4sigma4;
    }
    else if (_familyOfFunctional == std::string("MGGA"))
    {
        auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
        const auto dim      = &(mggafunc->dim);

        int n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int32_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int32_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;
        indices_and_counts["vsigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->vsigma});
        n_xc_outputs += dim->vsigma;
        indices_and_counts["vlapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->vlapl});
        n_xc_outputs += dim->vlapl;
        indices_and_counts["vtau"] = std::array<int32_t, 2>({n_xc_outputs, dim->vtau});
        n_xc_outputs += dim->vtau;

        indices_and_counts["v2rho2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;
        indices_and_counts["v2rhosigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rhosigma});
        n_xc_outputs += dim->v2rhosigma;
        indices_and_counts["v2rholapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rholapl});
        n_xc_outputs += dim->v2rholapl;
        indices_and_counts["v2rhotau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2rhotau});
        n_xc_outputs += dim->v2rhotau;
        indices_and_counts["v2sigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2sigma2});
        n_xc_outputs += dim->v2sigma2;
        indices_and_counts["v2sigmalapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2sigmalapl});
        n_xc_outputs += dim->v2sigmalapl;
        indices_and_counts["v2sigmatau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2sigmatau});
        n_xc_outputs += dim->v2sigmatau;
        indices_and_counts["v2lapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2lapl2});
        n_xc_outputs += dim->v2lapl2;
        indices_and_counts["v2lapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2lapltau});
        n_xc_outputs += dim->v2lapltau;
        indices_and_counts["v2tau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v2tau2});
        n_xc_outputs += dim->v2tau2;

        indices_and_counts["v3rho3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;
        indices_and_counts["v3rho2sigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho2sigma});
        n_xc_outputs += dim->v3rho2sigma;
        indices_and_counts["v3rho2lapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho2lapl});
        n_xc_outputs += dim->v3rho2lapl;
        indices_and_counts["v3rho2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rho2tau});
        n_xc_outputs += dim->v3rho2tau;
        indices_and_counts["v3rhosigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rhosigma2});
        n_xc_outputs += dim->v3rhosigma2;
        indices_and_counts["v3rhosigmalapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rhosigmalapl});
        n_xc_outputs += dim->v3rhosigmalapl;
        indices_and_counts["v3rhosigmatau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rhosigmatau});
        n_xc_outputs += dim->v3rhosigmatau;
        indices_and_counts["v3rholapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rholapl2});
        n_xc_outputs += dim->v3rholapl2;
        indices_and_counts["v3rholapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rholapltau});
        n_xc_outputs += dim->v3rholapltau;
        indices_and_counts["v3rhotau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3rhotau2});
        n_xc_outputs += dim->v3rhotau2;
        indices_and_counts["v3sigma3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigma3});
        n_xc_outputs += dim->v3sigma3;
        indices_and_counts["v3sigma2lapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigma2lapl});
        n_xc_outputs += dim->v3sigma2lapl;
        indices_and_counts["v3sigma2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigma2tau});
        n_xc_outputs += dim->v3sigma2tau;
        indices_and_counts["v3sigmalapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigmalapl2});
        n_xc_outputs += dim->v3sigmalapl2;
        indices_and_counts["v3sigmalapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigmalapltau});
        n_xc_outputs += dim->v3sigmalapltau;
        indices_and_counts["v3sigmatau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3sigmatau2});
        n_xc_outputs += dim->v3sigmatau2;
        indices_and_counts["v3lapl3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3lapl3});
        n_xc_outputs += dim->v3lapl3;
        indices_and_counts["v3lapl2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3lapl2tau});
        n_xc_outputs += dim->v3lapl2tau;
        indices_and_counts["v3lapltau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3lapltau2});
        n_xc_outputs += dim->v3lapltau2;
        indices_and_counts["v3tau3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v3tau3});
        n_xc_outputs += dim->v3tau3;

        indices_and_counts["v4rho4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
        indices_and_counts["v4rho3sigma"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho3sigma});
        n_xc_outputs += dim->v4rho3sigma;
        indices_and_counts["v4rho3lapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho3lapl});
        n_xc_outputs += dim->v4rho3lapl;
        indices_and_counts["v4rho3tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho3tau});
        n_xc_outputs += dim->v4rho3tau;
        indices_and_counts["v4rho2sigma2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2sigma2});
        n_xc_outputs += dim->v4rho2sigma2;
        indices_and_counts["v4rho2sigmalapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2sigmalapl});
        n_xc_outputs += dim->v4rho2sigmalapl;
        indices_and_counts["v4rho2sigmatau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2sigmatau});
        n_xc_outputs += dim->v4rho2sigmatau;
        indices_and_counts["v4rho2lapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2lapl2});
        n_xc_outputs += dim->v4rho2lapl2;
        indices_and_counts["v4rho2lapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2lapltau});
        n_xc_outputs += dim->v4rho2lapltau;
        indices_and_counts["v4rho2tau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rho2tau2});
        n_xc_outputs += dim->v4rho2tau2;
        indices_and_counts["v4rhosigma3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigma3});
        n_xc_outputs += dim->v4rhosigma3;
        indices_and_counts["v4rhosigma2lapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigma2lapl});
        n_xc_outputs += dim->v4rhosigma2lapl;
        indices_and_counts["v4rhosigma2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigma2tau});
        n_xc_outputs += dim->v4rhosigma2tau;
        indices_and_counts["v4rhosigmalapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigmalapl2});
        n_xc_outputs += dim->v4rhosigmalapl2;
        indices_and_counts["v4rhosigmalapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigmalapltau});
        n_xc_outputs += dim->v4rhosigmalapltau;
        indices_and_counts["v4rhosigmatau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhosigmatau2});
        n_xc_outputs += dim->v4rhosigmatau2;
        indices_and_counts["v4rholapl3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rholapl3});
        n_xc_outputs += dim->v4rholapl3;
        indices_and_counts["v4rholapl2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rholapl2tau});
        n_xc_outputs += dim->v4rholapl2tau;
        indices_and_counts["v4rholapltau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rholapltau2});
        n_xc_outputs += dim->v4rholapltau2;
        indices_and_counts["v4rhotau3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4rhotau3});
        n_xc_outputs += dim->v4rhotau3;
        indices_and_counts["v4sigma4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma4});
        n_xc_outputs += dim->v4sigma4;
        indices_and_counts["v4sigma3lapl"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma3lapl});
        n_xc_outputs += dim->v4sigma3lapl;
        indices_and_counts["v4sigma3tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma3tau});
        n_xc_outputs += dim->v4sigma3tau;
        indices_and_counts["v4sigma2lapl2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma2lapl2});
        n_xc_outputs += dim->v4sigma2lapl2;
        indices_and_counts["v4sigma2lapltau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma2lapltau});
        n_xc_outputs += dim->v4sigma2lapltau;
        indices_and_counts["v4sigma2tau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigma2tau2});
        n_xc_outputs += dim->v4sigma2tau2;
        indices_and_counts["v4sigmalapl3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigmalapl3});
        n_xc_outputs += dim->v4sigmalapl3;
        indices_and_counts["v4sigmalapl2tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigmalapl2tau});
        n_xc_outputs += dim->v4sigmalapl2tau;
        indices_and_counts["v4sigmalapltau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigmalapltau2});
        n_xc_outputs += dim->v4sigmalapltau2;
        indices_and_counts["v4sigmatau3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4sigmatau3});
        n_xc_outputs += dim->v4sigmatau3;
        indices_and_counts["v4lapl4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4lapl4});
        n_xc_outputs += dim->v4lapl4;
        indices_and_counts["v4lapl3tau"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4lapl3tau});
        n_xc_outputs += dim->v4lapl3tau;
        indices_and_counts["v4lapl2tau2"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4lapl2tau2});
        n_xc_outputs += dim->v4lapl2tau2;
        indices_and_counts["v4lapltau3"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4lapltau3});
        n_xc_outputs += dim->v4lapltau3;
        indices_and_counts["v4tau4"] = std::array<int32_t, 2>({n_xc_outputs, dim->v4tau4});
        n_xc_outputs += dim->v4tau4;
    }

    return indices_and_counts;
}

CXCFunctional&
CXCFunctional::operator=(const CXCFunctional& source)
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

CXCFunctional&
CXCFunctional::operator=(CXCFunctional&& source) noexcept
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
CXCFunctional::operator==(const CXCFunctional& other) const
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
CXCFunctional::operator!=(const CXCFunctional& other) const
{
    return !(*this == other);
}

std::string
CXCFunctional::getFunctionalLabel() const
{
    return _nameOfFunctional;
}

xcfun
CXCFunctional::getFunctionalType() const
{
    return to_xcfun(_familyOfFunctional);
}

bool
CXCFunctional::isUndefined() const
{
    return (fstr::upcase(_nameOfFunctional) == "UNDEFINED");
}

bool
CXCFunctional::isHybrid() const
{
    return (std::fabs(_fractionOfExactExchange) > 1.0e-13);
}

double
CXCFunctional::getFractionOfExactExchange() const
{
    return _fractionOfExactExchange;
}

auto
CXCFunctional::compute_exc_vxc_for_lda(const int32_t np, const double* rho, double* exc, double* vrho) const -> void
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

    gridscreen::screenExcVxcForLDA(np, rho, exc, vrho);
}

auto
CXCFunctional::compute_vxc_for_lda(const int32_t np, const double* rho, double* vrho) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

    auto stage_vrho = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];

    #pragma omp simd aligned(vrho : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        vrho[2 * g + 0] = 0.0;
        vrho[2 * g + 1] = 0.0;
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        xc_lda_vxc(funcptr, np, rho, stage_vrho);

        const auto c = xccomp.getScalingFactor();

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

    gridscreen::screenVxcForLDA(np, rho, vrho);
}

auto
CXCFunctional::compute_fxc_for_lda(const int32_t np, const double* rho, double* v2rho2) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

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

    gridscreen::screenFxcForLDA(np, rho, v2rho2);
}

auto
CXCFunctional::compute_kxc_for_lda(const int32_t np, const double* rho, double* v3rho3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

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

    gridscreen::screenKxcForLDA(np, rho, v3rho3);
}

auto
CXCFunctional::compute_lxc_for_lda(const int32_t np, const double* rho, double* v4rho4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

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

    gridscreen::screenLxcForLDA(np, rho, v4rho4);
}

auto
CXCFunctional::compute_exc_vxc_for_gga(const int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma) const
    -> void
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

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
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

    gridscreen::screenExcVxcForGGA(np, rho, sigma, exc, vrho, vsigma);
}

auto
CXCFunctional::compute_vxc_for_gga(const int32_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void
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

    auto stage_vrho   = (alloc) ? mem::malloc<double>(2 * np) : &_stagingBuffer[1 * _ldStaging];
    auto stage_vsigma = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[3 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

            #pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (xccomp.isGGA())
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

    gridscreen::screenVxcForGGA(np, rho, sigma, vrho, vsigma);
}

auto
CXCFunctional::compute_fxc_for_gga(const int32_t np, const double* rho, const double* sigma, double* v2rho2, double* v2rhosigma, double* v2sigma2)
    const -> void
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

    auto stage_v2rho2     = (alloc) ? mem::malloc<double>(3 * np) : &_stagingBuffer[6 * _ldStaging];
    auto stage_v2rhosigma = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[9 * _ldStaging];
    auto stage_v2sigma2   = (alloc) ? mem::malloc<double>(6 * np) : &_stagingBuffer[15 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
        {
            xc_gga_fxc(funcptr, np, rho, sigma, stage_v2rho2, stage_v2rhosigma, stage_v2sigma2);

            #pragma omp simd aligned(v2rho2, stage_v2rho2, v2rhosigma, stage_v2rhosigma, \
                                     v2sigma2, stage_v2sigma2 : VLX_ALIGN)
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

    gridscreen::screenFxcForGGA(np, rho, sigma, v2rho2, v2rhosigma, v2sigma2);
}

auto
CXCFunctional::compute_kxc_for_gga(int32_t       np,
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

    auto stage_v3rho3      = (alloc) ? mem::malloc<double>(4 * np) : &_stagingBuffer[21 * _ldStaging];
    auto stage_v3rho2sigma = (alloc) ? mem::malloc<double>(9 * np) : &_stagingBuffer[25 * _ldStaging];
    auto stage_v3rhosigma2 = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[34 * _ldStaging];
    auto stage_v3sigma3    = (alloc) ? mem::malloc<double>(10 * np) : &_stagingBuffer[46 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
        {
            xc_gga_kxc(funcptr, np, rho, sigma, stage_v3rho3, stage_v3rho2sigma, stage_v3rhosigma2, stage_v3sigma3);

            #pragma omp simd aligned(v3rho3, stage_v3rho3, v3rho2sigma, stage_v3rho2sigma, \
                                     v3rhosigma2, stage_v3rhosigma2, v3sigma3, stage_v3sigma3 : VLX_ALIGN)
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

    gridscreen::screenKxcForGGA(np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
}

auto
CXCFunctional::compute_lxc_for_gga(int32_t       np,
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

    auto stage_v4rho4       = (alloc) ? mem::malloc<double>(5 * np) : &_stagingBuffer[56 * _ldStaging];
    auto stage_v4rho3sigma  = (alloc) ? mem::malloc<double>(12 * np) : &_stagingBuffer[61 * _ldStaging];
    auto stage_v4rho2sigma2 = (alloc) ? mem::malloc<double>(18 * np) : &_stagingBuffer[73 * _ldStaging];
    auto stage_v4rhosigma3  = (alloc) ? mem::malloc<double>(20 * np) : &_stagingBuffer[91 * _ldStaging];
    auto stage_v4sigma4     = (alloc) ? mem::malloc<double>(15 * np) : &_stagingBuffer[111 * _ldStaging];

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
        {
            xc_gga_lxc(funcptr, np, rho, sigma, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4);

            #pragma omp simd aligned(v4rho4, stage_v4rho4, v4rho3sigma, stage_v4rho3sigma, \
                                     v4rho2sigma2, stage_v4rho2sigma2, v4rhosigma3, stage_v4rhosigma3, \
                                     v4sigma4, stage_v4sigma4 : VLX_ALIGN)
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

    gridscreen::screenLxcForGGA(np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4);
}

auto
CXCFunctional::compute_exc_vxc_for_mgga(int32_t       np,
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

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
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
        else if (xccomp.isMetaGGA())
        {
            xc_mgga_exc_vxc(funcptr, np, rho, sigma, lapl, tau, stage_exc, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

            #pragma omp simd aligned(exc, stage_exc, vrho, stage_vrho, vsigma, stage_vsigma, \
                                     vlapl, stage_vlapl, vtau, stage_vtau : VLX_ALIGN)
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

    gridscreen::screenExcVxcForMGGA(np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
}

auto
CXCFunctional::compute_vxc_for_mgga(int32_t       np,
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
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

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

        if (xccomp.isLDA())
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

            #pragma omp simd aligned(vrho, stage_vrho : VLX_ALIGN)
            for (auto g = 0; g < np; ++g)
            {
                vrho[2 * g + 0] += c * stage_vrho[2 * g + 0];
                vrho[2 * g + 1] += c * stage_vrho[2 * g + 1];
            }
        }
        else if (xccomp.isGGA())
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
        else if (xccomp.isMetaGGA())
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

    gridscreen::screenVxcForMGGA(np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);
}

auto
CXCFunctional::compute_fxc_for_mgga(int32_t       np,
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

    #pragma omp simd aligned(v2rho2, v2rhosigma, v2rholapl, v2rhotau, \
                             v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, \
                             v2lapltau, v2tau2 : VLX_ALIGN)
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

        v2rholapl[4 * g + 0] = 0.0;
        v2rholapl[4 * g + 1] = 0.0;
        v2rholapl[4 * g + 2] = 0.0;
        v2rholapl[4 * g + 3] = 0.0;

        v2rhotau[4 * g + 0] = 0.0;
        v2rhotau[4 * g + 1] = 0.0;
        v2rhotau[4 * g + 2] = 0.0;
        v2rhotau[4 * g + 3] = 0.0;

        v2sigma2[6 * g + 0] = 0.0;
        v2sigma2[6 * g + 1] = 0.0;
        v2sigma2[6 * g + 2] = 0.0;
        v2sigma2[6 * g + 3] = 0.0;
        v2sigma2[6 * g + 4] = 0.0;
        v2sigma2[6 * g + 5] = 0.0;

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
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

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

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
        {
            xc_gga_fxc(funcptr, np, rho, sigma, stage_v2rho2, stage_v2rhosigma, stage_v2sigma2);

            #pragma omp simd aligned(v2rho2, stage_v2rho2, v2rhosigma, stage_v2rhosigma, \
                                     v2sigma2, stage_v2sigma2 : VLX_ALIGN)
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
        else if (xccomp.isMetaGGA())
        {
            xc_mgga_fxc(funcptr,
                        np,
                        rho,
                        sigma,
                        lapl,
                        tau,
                        stage_v2rho2,
                        stage_v2rhosigma,
                        stage_v2rholapl,
                        stage_v2rhotau,
                        stage_v2sigma2,
                        stage_v2sigmalapl,
                        stage_v2sigmatau,
                        stage_v2lapl2,
                        stage_v2lapltau,
                        stage_v2tau2);

            #pragma omp simd aligned(v2rho2, stage_v2rho2, v2rhosigma, stage_v2rhosigma, \
                                     v2rholapl, stage_v2rholapl, v2rhotau, stage_v2rhotau, \
                                     v2sigma2, stage_v2sigma2, v2sigmalapl, stage_v2sigmalapl, \
                                     v2sigmatau, stage_v2sigmatau, v2lapl2, stage_v2lapl2, \
                                     v2lapltau, stage_v2lapltau, v2tau2, stage_v2tau2 : VLX_ALIGN)
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

                v2rholapl[4 * g + 0] += c * stage_v2rholapl[4 * g + 0];
                v2rholapl[4 * g + 1] += c * stage_v2rholapl[4 * g + 1];
                v2rholapl[4 * g + 2] += c * stage_v2rholapl[4 * g + 2];
                v2rholapl[4 * g + 3] += c * stage_v2rholapl[4 * g + 3];

                v2rhotau[4 * g + 0] += c * stage_v2rhotau[4 * g + 0];
                v2rhotau[4 * g + 1] += c * stage_v2rhotau[4 * g + 1];
                v2rhotau[4 * g + 2] += c * stage_v2rhotau[4 * g + 2];
                v2rhotau[4 * g + 3] += c * stage_v2rhotau[4 * g + 3];

                v2sigma2[6 * g + 0] += c * stage_v2sigma2[6 * g + 0];
                v2sigma2[6 * g + 1] += c * stage_v2sigma2[6 * g + 1];
                v2sigma2[6 * g + 2] += c * stage_v2sigma2[6 * g + 2];
                v2sigma2[6 * g + 3] += c * stage_v2sigma2[6 * g + 3];
                v2sigma2[6 * g + 4] += c * stage_v2sigma2[6 * g + 4];
                v2sigma2[6 * g + 5] += c * stage_v2sigma2[6 * g + 5];

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
            }
        }
    }

    if (alloc)
    {
        mem::free(stage_v2rho2);
        mem::free(stage_v2rhosigma);
        mem::free(stage_v2rholapl);
        mem::free(stage_v2rhotau);
        mem::free(stage_v2sigma2);
        mem::free(stage_v2sigmalapl);
        mem::free(stage_v2sigmatau);
        mem::free(stage_v2lapl2);
        mem::free(stage_v2lapltau);
        mem::free(stage_v2tau2);
    }

    gridscreen::screenFxcForMGGA(
        np, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);
}

auto
CXCFunctional::compute_kxc_for_mgga(int32_t       np,
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

    #pragma omp simd aligned(v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, \
                             v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rholapl2, \
                             v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl, \
                             v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, \
                             v3lapl3, v3lapl2tau, v3lapltau2, v3tau3 : VLX_ALIGN)
    for (auto g = 0; g < np; ++g)
    {
        v3rho3[4 * g + 0] = 0;
        v3rho3[4 * g + 1] = 0;
        v3rho3[4 * g + 2] = 0;
        v3rho3[4 * g + 3] = 0;

        v3rho2sigma[9 * g + 0] = 0;
        v3rho2sigma[9 * g + 1] = 0;
        v3rho2sigma[9 * g + 2] = 0;
        v3rho2sigma[9 * g + 3] = 0;
        v3rho2sigma[9 * g + 4] = 0;
        v3rho2sigma[9 * g + 5] = 0;
        v3rho2sigma[9 * g + 6] = 0;
        v3rho2sigma[9 * g + 7] = 0;
        v3rho2sigma[9 * g + 8] = 0;

        v3rho2lapl[6 * g + 0] = 0;
        v3rho2lapl[6 * g + 1] = 0;
        v3rho2lapl[6 * g + 2] = 0;
        v3rho2lapl[6 * g + 3] = 0;
        v3rho2lapl[6 * g + 4] = 0;
        v3rho2lapl[6 * g + 5] = 0;

        v3rho2tau[6 * g + 0] = 0;
        v3rho2tau[6 * g + 1] = 0;
        v3rho2tau[6 * g + 2] = 0;
        v3rho2tau[6 * g + 3] = 0;
        v3rho2tau[6 * g + 4] = 0;
        v3rho2tau[6 * g + 5] = 0;

        v3rhosigma2[12 * g + 0]  = 0;
        v3rhosigma2[12 * g + 1]  = 0;
        v3rhosigma2[12 * g + 2]  = 0;
        v3rhosigma2[12 * g + 3]  = 0;
        v3rhosigma2[12 * g + 4]  = 0;
        v3rhosigma2[12 * g + 5]  = 0;
        v3rhosigma2[12 * g + 6]  = 0;
        v3rhosigma2[12 * g + 7]  = 0;
        v3rhosigma2[12 * g + 8]  = 0;
        v3rhosigma2[12 * g + 9]  = 0;
        v3rhosigma2[12 * g + 10] = 0;
        v3rhosigma2[12 * g + 11] = 0;

        v3rhosigmalapl[12 * g + 0]  = 0;
        v3rhosigmalapl[12 * g + 1]  = 0;
        v3rhosigmalapl[12 * g + 2]  = 0;
        v3rhosigmalapl[12 * g + 3]  = 0;
        v3rhosigmalapl[12 * g + 4]  = 0;
        v3rhosigmalapl[12 * g + 5]  = 0;
        v3rhosigmalapl[12 * g + 6]  = 0;
        v3rhosigmalapl[12 * g + 7]  = 0;
        v3rhosigmalapl[12 * g + 8]  = 0;
        v3rhosigmalapl[12 * g + 9]  = 0;
        v3rhosigmalapl[12 * g + 10] = 0;
        v3rhosigmalapl[12 * g + 11] = 0;

        v3rhosigmatau[12 * g + 0]  = 0;
        v3rhosigmatau[12 * g + 1]  = 0;
        v3rhosigmatau[12 * g + 2]  = 0;
        v3rhosigmatau[12 * g + 3]  = 0;
        v3rhosigmatau[12 * g + 4]  = 0;
        v3rhosigmatau[12 * g + 5]  = 0;
        v3rhosigmatau[12 * g + 6]  = 0;
        v3rhosigmatau[12 * g + 7]  = 0;
        v3rhosigmatau[12 * g + 8]  = 0;
        v3rhosigmatau[12 * g + 9]  = 0;
        v3rhosigmatau[12 * g + 10] = 0;
        v3rhosigmatau[12 * g + 11] = 0;

        v3rholapl2[6 * g + 0] = 0;
        v3rholapl2[6 * g + 1] = 0;
        v3rholapl2[6 * g + 2] = 0;
        v3rholapl2[6 * g + 3] = 0;
        v3rholapl2[6 * g + 4] = 0;
        v3rholapl2[6 * g + 5] = 0;

        v3rholapltau[8 * g + 0] = 0;
        v3rholapltau[8 * g + 1] = 0;
        v3rholapltau[8 * g + 2] = 0;
        v3rholapltau[8 * g + 3] = 0;
        v3rholapltau[8 * g + 4] = 0;
        v3rholapltau[8 * g + 5] = 0;
        v3rholapltau[8 * g + 6] = 0;
        v3rholapltau[8 * g + 7] = 0;

        v3rhotau2[6 * g + 0] = 0;
        v3rhotau2[6 * g + 1] = 0;
        v3rhotau2[6 * g + 2] = 0;
        v3rhotau2[6 * g + 3] = 0;
        v3rhotau2[6 * g + 4] = 0;
        v3rhotau2[6 * g + 5] = 0;

        v3sigma3[10 * g + 0] = 0;
        v3sigma3[10 * g + 1] = 0;
        v3sigma3[10 * g + 2] = 0;
        v3sigma3[10 * g + 3] = 0;
        v3sigma3[10 * g + 4] = 0;
        v3sigma3[10 * g + 5] = 0;
        v3sigma3[10 * g + 6] = 0;
        v3sigma3[10 * g + 7] = 0;
        v3sigma3[10 * g + 8] = 0;
        v3sigma3[10 * g + 9] = 0;

        v3sigma2lapl[12 * g + 0]  = 0;
        v3sigma2lapl[12 * g + 1]  = 0;
        v3sigma2lapl[12 * g + 2]  = 0;
        v3sigma2lapl[12 * g + 3]  = 0;
        v3sigma2lapl[12 * g + 4]  = 0;
        v3sigma2lapl[12 * g + 5]  = 0;
        v3sigma2lapl[12 * g + 6]  = 0;
        v3sigma2lapl[12 * g + 7]  = 0;
        v3sigma2lapl[12 * g + 8]  = 0;
        v3sigma2lapl[12 * g + 9]  = 0;
        v3sigma2lapl[12 * g + 10] = 0;
        v3sigma2lapl[12 * g + 11] = 0;

        v3sigma2tau[12 * g + 0]  = 0;
        v3sigma2tau[12 * g + 1]  = 0;
        v3sigma2tau[12 * g + 2]  = 0;
        v3sigma2tau[12 * g + 3]  = 0;
        v3sigma2tau[12 * g + 4]  = 0;
        v3sigma2tau[12 * g + 5]  = 0;
        v3sigma2tau[12 * g + 6]  = 0;
        v3sigma2tau[12 * g + 7]  = 0;
        v3sigma2tau[12 * g + 8]  = 0;
        v3sigma2tau[12 * g + 9]  = 0;
        v3sigma2tau[12 * g + 10] = 0;
        v3sigma2tau[12 * g + 11] = 0;

        v3sigmalapl2[9 * g + 0] = 0;
        v3sigmalapl2[9 * g + 1] = 0;
        v3sigmalapl2[9 * g + 2] = 0;
        v3sigmalapl2[9 * g + 3] = 0;
        v3sigmalapl2[9 * g + 4] = 0;
        v3sigmalapl2[9 * g + 5] = 0;
        v3sigmalapl2[9 * g + 6] = 0;
        v3sigmalapl2[9 * g + 7] = 0;
        v3sigmalapl2[9 * g + 8] = 0;

        v3sigmalapltau[12 * g + 0]  = 0;
        v3sigmalapltau[12 * g + 1]  = 0;
        v3sigmalapltau[12 * g + 2]  = 0;
        v3sigmalapltau[12 * g + 3]  = 0;
        v3sigmalapltau[12 * g + 4]  = 0;
        v3sigmalapltau[12 * g + 5]  = 0;
        v3sigmalapltau[12 * g + 6]  = 0;
        v3sigmalapltau[12 * g + 7]  = 0;
        v3sigmalapltau[12 * g + 8]  = 0;
        v3sigmalapltau[12 * g + 9]  = 0;
        v3sigmalapltau[12 * g + 10] = 0;
        v3sigmalapltau[12 * g + 11] = 0;

        v3sigmatau2[9 * g + 0] = 0;
        v3sigmatau2[9 * g + 1] = 0;
        v3sigmatau2[9 * g + 2] = 0;
        v3sigmatau2[9 * g + 3] = 0;
        v3sigmatau2[9 * g + 4] = 0;
        v3sigmatau2[9 * g + 5] = 0;
        v3sigmatau2[9 * g + 6] = 0;
        v3sigmatau2[9 * g + 7] = 0;
        v3sigmatau2[9 * g + 8] = 0;

        v3lapl3[4 * g + 0] = 0;
        v3lapl3[4 * g + 1] = 0;
        v3lapl3[4 * g + 2] = 0;
        v3lapl3[4 * g + 3] = 0;

        v3lapl2tau[6 * g + 0] = 0;
        v3lapl2tau[6 * g + 1] = 0;
        v3lapl2tau[6 * g + 2] = 0;
        v3lapl2tau[6 * g + 3] = 0;
        v3lapl2tau[6 * g + 4] = 0;
        v3lapl2tau[6 * g + 5] = 0;

        v3lapltau2[6 * g + 0] = 0;
        v3lapltau2[6 * g + 1] = 0;
        v3lapltau2[6 * g + 2] = 0;
        v3lapltau2[6 * g + 3] = 0;
        v3lapltau2[6 * g + 4] = 0;
        v3lapltau2[6 * g + 5] = 0;

        v3tau3[4 * g + 0] = 0;
        v3tau3[4 * g + 1] = 0;
        v3tau3[4 * g + 2] = 0;
        v3tau3[4 * g + 3] = 0;
    }

    // should we allocate staging buffers? Or can we use the global one?
    bool alloc = (np > _ldStaging);

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

        if (xccomp.isLDA())
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
        else if (xccomp.isGGA())
        {
            xc_gga_kxc(funcptr, np, rho, sigma, stage_v3rho3, stage_v3rho2sigma, stage_v3rhosigma2, stage_v3sigma3);

            #pragma omp simd aligned(v3rho3, stage_v3rho3, v3rho2sigma, stage_v3rho2sigma, \
                                     v3rhosigma2, stage_v3rhosigma2, v3sigma3, stage_v3sigma3 : VLX_ALIGN)
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
        else if (xccomp.isMetaGGA())
        {
            xc_mgga_kxc(funcptr,
                        np,
                        rho,
                        sigma,
                        lapl,
                        tau,
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

            #pragma omp simd aligned(v3rho3, stage_v3rho3, v3rho2sigma, stage_v3rho2sigma, \
                                     v3rho2lapl, stage_v3rho2lapl, v3rho2tau, stage_v3rho2tau, \
                                     v3rhosigma2, stage_v3rhosigma2, v3rhosigmalapl, stage_v3rhosigmalapl, \
                                     v3rhosigmatau, stage_v3rhosigmatau, v3rholapl2, stage_v3rholapl2, \
                                     v3rholapltau, stage_v3rholapltau, v3rhotau2, stage_v3rhotau2, \
                                     v3sigma3, stage_v3sigma3, v3sigma2lapl, stage_v3sigma2lapl, \
                                     v3sigma2tau, stage_v3sigma2tau, v3sigmalapl2, stage_v3sigmalapl2, \
                                     v3sigmalapltau, stage_v3sigmalapltau, v3sigmatau2, stage_v3sigmatau2, \
                                     v3lapl3, stage_v3lapl3, v3lapl2tau, stage_v3lapl2tau, \
                                     v3lapltau2, stage_v3lapltau2, v3tau3, stage_v3tau3 : VLX_ALIGN)
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

                v3sigmalapl2[9 * g + 0] += c * stage_v3sigmalapl2[9 * g + 0];
                v3sigmalapl2[9 * g + 1] += c * stage_v3sigmalapl2[9 * g + 1];
                v3sigmalapl2[9 * g + 2] += c * stage_v3sigmalapl2[9 * g + 2];
                v3sigmalapl2[9 * g + 3] += c * stage_v3sigmalapl2[9 * g + 3];
                v3sigmalapl2[9 * g + 4] += c * stage_v3sigmalapl2[9 * g + 4];
                v3sigmalapl2[9 * g + 5] += c * stage_v3sigmalapl2[9 * g + 5];
                v3sigmalapl2[9 * g + 6] += c * stage_v3sigmalapl2[9 * g + 6];
                v3sigmalapl2[9 * g + 7] += c * stage_v3sigmalapl2[9 * g + 7];
                v3sigmalapl2[9 * g + 8] += c * stage_v3sigmalapl2[9 * g + 8];

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

    gridscreen::screenKxcForMGGA(np,
                                 rho,
                                 sigma,
                                 lapl,
                                 tau,
                                 v3rho3,
                                 v3rho2sigma,
                                 v3rho2lapl,
                                 v3rho2tau,
                                 v3rhosigma2,
                                 v3rhosigmalapl,
                                 v3rhosigmatau,
                                 v3rholapl2,
                                 v3rholapltau,
                                 v3rhotau2,
                                 v3sigma3,
                                 v3sigma2lapl,
                                 v3sigma2tau,
                                 v3sigmalapl2,
                                 v3sigmalapltau,
                                 v3sigmatau2,
                                 v3lapl3,
                                 v3lapl2tau,
                                 v3lapltau2,
                                 v3tau3);
}

auto
CXCFunctional::compute_lxc_for_mgga(int32_t       np,
                                    const double* rho,
                                    const double* sigma,
                                    const double* lapl,
                                    const double* tau,
                                    double*       v4rho4,
                                    double*       v4rho3sigma,
                                    double*       v4rho3lapl,
                                    double*       v4rho3tau,
                                    double*       v4rho2sigma2,
                                    double*       v4rho2sigmalapl,
                                    double*       v4rho2sigmatau,
                                    double*       v4rho2lapl2,
                                    double*       v4rho2lapltau,
                                    double*       v4rho2tau2,
                                    double*       v4rhosigma3,
                                    double*       v4rhosigma2lapl,
                                    double*       v4rhosigma2tau,
                                    double*       v4rhosigmalapl2,
                                    double*       v4rhosigmalapltau,
                                    double*       v4rhosigmatau2,
                                    double*       v4rholapl3,
                                    double*       v4rholapl2tau,
                                    double*       v4rholapltau2,
                                    double*       v4rhotau3,
                                    double*       v4sigma4,
                                    double*       v4sigma3lapl,
                                    double*       v4sigma3tau,
                                    double*       v4sigma2lapl2,
                                    double*       v4sigma2lapltau,
                                    double*       v4sigma2tau2,
                                    double*       v4sigmalapl3,
                                    double*       v4sigmalapl2tau,
                                    double*       v4sigmalapltau2,
                                    double*       v4sigmatau3,
                                    double*       v4lapl4,
                                    double*       v4lapl3tau,
                                    double*       v4lapl2tau2,
                                    double*       v4lapltau3,
                                    double*       v4tau4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v4rho4            = nullptr;
    double* stage_v4rho3sigma       = nullptr;
    double* stage_v4rho3lapl        = nullptr;
    double* stage_v4rho3tau         = nullptr;
    double* stage_v4rho2sigma2      = nullptr;
    double* stage_v4rho2sigmalapl   = nullptr;
    double* stage_v4rho2sigmatau    = nullptr;
    double* stage_v4rho2lapl2       = nullptr;
    double* stage_v4rho2lapltau     = nullptr;
    double* stage_v4rho2tau2        = nullptr;
    double* stage_v4rhosigma3       = nullptr;
    double* stage_v4rhosigma2lapl   = nullptr;
    double* stage_v4rhosigma2tau    = nullptr;
    double* stage_v4rhosigmalapl2   = nullptr;
    double* stage_v4rhosigmalapltau = nullptr;
    double* stage_v4rhosigmatau2    = nullptr;
    double* stage_v4rholapl3        = nullptr;
    double* stage_v4rholapl2tau     = nullptr;
    double* stage_v4rholapltau2     = nullptr;
    double* stage_v4rhotau3         = nullptr;
    double* stage_v4sigma4          = nullptr;
    double* stage_v4sigma3lapl      = nullptr;
    double* stage_v4sigma3tau       = nullptr;
    double* stage_v4sigma2lapl2     = nullptr;
    double* stage_v4sigma2lapltau   = nullptr;
    double* stage_v4sigma2tau2      = nullptr;
    double* stage_v4sigmalapl3      = nullptr;
    double* stage_v4sigmalapl2tau   = nullptr;
    double* stage_v4sigmalapltau2   = nullptr;
    double* stage_v4sigmatau3       = nullptr;
    double* stage_v4lapl4           = nullptr;
    double* stage_v4lapl3tau        = nullptr;
    double* stage_v4lapl2tau2       = nullptr;
    double* stage_v4lapltau3        = nullptr;
    double* stage_v4tau4            = nullptr;

    if (alloc)
    {
        stage_v4rho4            = mem::malloc<double>(ind_cnt["v4rho4"][1] * np);
        stage_v4rho3sigma       = mem::malloc<double>(ind_cnt["v4rho3sigma"][1] * np);
        stage_v4rho3lapl        = mem::malloc<double>(ind_cnt["v4rho3lapl"][1] * np);
        stage_v4rho3tau         = mem::malloc<double>(ind_cnt["v4rho3tau"][1] * np);
        stage_v4rho2sigma2      = mem::malloc<double>(ind_cnt["v4rho2sigma2"][1] * np);
        stage_v4rho2sigmalapl   = mem::malloc<double>(ind_cnt["v4rho2sigmalapl"][1] * np);
        stage_v4rho2sigmatau    = mem::malloc<double>(ind_cnt["v4rho2sigmatau"][1] * np);
        stage_v4rho2lapl2       = mem::malloc<double>(ind_cnt["v4rho2lapl2"][1] * np);
        stage_v4rho2lapltau     = mem::malloc<double>(ind_cnt["v4rho2lapltau"][1] * np);
        stage_v4rho2tau2        = mem::malloc<double>(ind_cnt["v4rho2tau2"][1] * np);
        stage_v4rhosigma3       = mem::malloc<double>(ind_cnt["v4rhosigma3"][1] * np);
        stage_v4rhosigma2lapl   = mem::malloc<double>(ind_cnt["v4rhosigma2lapl"][1] * np);
        stage_v4rhosigma2tau    = mem::malloc<double>(ind_cnt["v4rhosigma2tau"][1] * np);
        stage_v4rhosigmalapl2   = mem::malloc<double>(ind_cnt["v4rhosigmalapl2"][1] * np);
        stage_v4rhosigmalapltau = mem::malloc<double>(ind_cnt["v4rhosigmalapltau"][1] * np);
        stage_v4rhosigmatau2    = mem::malloc<double>(ind_cnt["v4rhosigmatau2"][1] * np);
        stage_v4rholapl3        = mem::malloc<double>(ind_cnt["v4rholapl3"][1] * np);
        stage_v4rholapl2tau     = mem::malloc<double>(ind_cnt["v4rholapl2tau"][1] * np);
        stage_v4rholapltau2     = mem::malloc<double>(ind_cnt["v4rholapltau2"][1] * np);
        stage_v4rhotau3         = mem::malloc<double>(ind_cnt["v4rhotau3"][1] * np);
        stage_v4sigma4          = mem::malloc<double>(ind_cnt["v4sigma4"][1] * np);
        stage_v4sigma3lapl      = mem::malloc<double>(ind_cnt["v4sigma3lapl"][1] * np);
        stage_v4sigma3tau       = mem::malloc<double>(ind_cnt["v4sigma3tau"][1] * np);
        stage_v4sigma2lapl2     = mem::malloc<double>(ind_cnt["v4sigma2lapl2"][1] * np);
        stage_v4sigma2lapltau   = mem::malloc<double>(ind_cnt["v4sigma2lapltau"][1] * np);
        stage_v4sigma2tau2      = mem::malloc<double>(ind_cnt["v4sigma2tau2"][1] * np);
        stage_v4sigmalapl3      = mem::malloc<double>(ind_cnt["v4sigmalapl3"][1] * np);
        stage_v4sigmalapl2tau   = mem::malloc<double>(ind_cnt["v4sigmalapl2tau"][1] * np);
        stage_v4sigmalapltau2   = mem::malloc<double>(ind_cnt["v4sigmalapltau2"][1] * np);
        stage_v4sigmatau3       = mem::malloc<double>(ind_cnt["v4sigmatau3"][1] * np);
        stage_v4lapl4           = mem::malloc<double>(ind_cnt["v4lapl4"][1] * np);
        stage_v4lapl3tau        = mem::malloc<double>(ind_cnt["v4lapl3tau"][1] * np);
        stage_v4lapl2tau2       = mem::malloc<double>(ind_cnt["v4lapl2tau2"][1] * np);
        stage_v4lapltau3        = mem::malloc<double>(ind_cnt["v4lapltau3"][1] * np);
        stage_v4tau4            = mem::malloc<double>(ind_cnt["v4tau4"][1] * np);
    }
    else
    {
        stage_v4rho4            = &_stagingBuffer[ind_cnt["v4rho4"][0] * _ldStaging];
        stage_v4rho3sigma       = &_stagingBuffer[ind_cnt["v4rho3sigma"][0] * _ldStaging];
        stage_v4rho3lapl        = &_stagingBuffer[ind_cnt["v4rho3lapl"][0] * _ldStaging];
        stage_v4rho3tau         = &_stagingBuffer[ind_cnt["v4rho3tau"][0] * _ldStaging];
        stage_v4rho2sigma2      = &_stagingBuffer[ind_cnt["v4rho2sigma2"][0] * _ldStaging];
        stage_v4rho2sigmalapl   = &_stagingBuffer[ind_cnt["v4rho2sigmalapl"][0] * _ldStaging];
        stage_v4rho2sigmatau    = &_stagingBuffer[ind_cnt["v4rho2sigmatau"][0] * _ldStaging];
        stage_v4rho2lapl2       = &_stagingBuffer[ind_cnt["v4rho2lapl2"][0] * _ldStaging];
        stage_v4rho2lapltau     = &_stagingBuffer[ind_cnt["v4rho2lapltau"][0] * _ldStaging];
        stage_v4rho2tau2        = &_stagingBuffer[ind_cnt["v4rho2tau2"][0] * _ldStaging];
        stage_v4rhosigma3       = &_stagingBuffer[ind_cnt["v4rhosigma3"][0] * _ldStaging];
        stage_v4rhosigma2lapl   = &_stagingBuffer[ind_cnt["v4rhosigma2lapl"][0] * _ldStaging];
        stage_v4rhosigma2tau    = &_stagingBuffer[ind_cnt["v4rhosigma2tau"][0] * _ldStaging];
        stage_v4rhosigmalapl2   = &_stagingBuffer[ind_cnt["v4rhosigmalapl2"][0] * _ldStaging];
        stage_v4rhosigmalapltau = &_stagingBuffer[ind_cnt["v4rhosigmalapltau"][0] * _ldStaging];
        stage_v4rhosigmatau2    = &_stagingBuffer[ind_cnt["v4rhosigmatau2"][0] * _ldStaging];
        stage_v4rholapl3        = &_stagingBuffer[ind_cnt["v4rholapl3"][0] * _ldStaging];
        stage_v4rholapl2tau     = &_stagingBuffer[ind_cnt["v4rholapl2tau"][0] * _ldStaging];
        stage_v4rholapltau2     = &_stagingBuffer[ind_cnt["v4rholapltau2"][0] * _ldStaging];
        stage_v4rhotau3         = &_stagingBuffer[ind_cnt["v4rhotau3"][0] * _ldStaging];
        stage_v4sigma4          = &_stagingBuffer[ind_cnt["v4sigma4"][0] * _ldStaging];
        stage_v4sigma3lapl      = &_stagingBuffer[ind_cnt["v4sigma3lapl"][0] * _ldStaging];
        stage_v4sigma3tau       = &_stagingBuffer[ind_cnt["v4sigma3tau"][0] * _ldStaging];
        stage_v4sigma2lapl2     = &_stagingBuffer[ind_cnt["v4sigma2lapl2"][0] * _ldStaging];
        stage_v4sigma2lapltau   = &_stagingBuffer[ind_cnt["v4sigma2lapltau"][0] * _ldStaging];
        stage_v4sigma2tau2      = &_stagingBuffer[ind_cnt["v4sigma2tau2"][0] * _ldStaging];
        stage_v4sigmalapl3      = &_stagingBuffer[ind_cnt["v4sigmalapl3"][0] * _ldStaging];
        stage_v4sigmalapl2tau   = &_stagingBuffer[ind_cnt["v4sigmalapl2tau"][0] * _ldStaging];
        stage_v4sigmalapltau2   = &_stagingBuffer[ind_cnt["v4sigmalapltau2"][0] * _ldStaging];
        stage_v4sigmatau3       = &_stagingBuffer[ind_cnt["v4sigmatau3"][0] * _ldStaging];
        stage_v4lapl4           = &_stagingBuffer[ind_cnt["v4lapl4"][0] * _ldStaging];
        stage_v4lapl3tau        = &_stagingBuffer[ind_cnt["v4lapl3tau"][0] * _ldStaging];
        stage_v4lapl2tau2       = &_stagingBuffer[ind_cnt["v4lapl2tau2"][0] * _ldStaging];
        stage_v4lapltau3        = &_stagingBuffer[ind_cnt["v4lapltau3"][0] * _ldStaging];
        stage_v4tau4            = &_stagingBuffer[ind_cnt["v4tau4"][0] * _ldStaging];
    }

    // compute derivatives

    auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    auto nthreads = omp_get_max_threads();

    #pragma omp parallel
    {
        auto thread_id = omp_get_thread_num();

        auto grid_batch_size = static_cast<int>(mpi::batch_size(np, thread_id, nthreads));

        auto grid_batch_offset = static_cast<int>(mpi::batch_offset(np, thread_id, nthreads));

        for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; ++g)
        {
            for (int ind = 0; ind < dim->v4rho4; ++ind)
            {
                v4rho4[dim->v4rho4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3sigma; ++ind)
            {
                v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3lapl; ++ind)
            {
                v4rho3lapl[dim->v4rho3lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho3tau; ++ind)
            {
                v4rho3tau[dim->v4rho3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigma2; ++ind)
            {
                v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigmalapl; ++ind)
            {
                v4rho2sigmalapl[dim->v4rho2sigmalapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2sigmatau; ++ind)
            {
                v4rho2sigmatau[dim->v4rho2sigmatau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2lapl2; ++ind)
            {
                v4rho2lapl2[dim->v4rho2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2lapltau; ++ind)
            {
                v4rho2lapltau[dim->v4rho2lapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rho2tau2; ++ind)
            {
                v4rho2tau2[dim->v4rho2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma3; ++ind)
            {
                v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma2lapl; ++ind)
            {
                v4rhosigma2lapl[dim->v4rhosigma2lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigma2tau; ++ind)
            {
                v4rhosigma2tau[dim->v4rhosigma2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmalapl2; ++ind)
            {
                v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmalapltau; ++ind)
            {
                v4rhosigmalapltau[dim->v4rhosigmalapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhosigmatau2; ++ind)
            {
                v4rhosigmatau2[dim->v4rhosigmatau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapl3; ++ind)
            {
                v4rholapl3[dim->v4rholapl3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapl2tau; ++ind)
            {
                v4rholapl2tau[dim->v4rholapl2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rholapltau2; ++ind)
            {
                v4rholapltau2[dim->v4rholapltau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4rhotau3; ++ind)
            {
                v4rhotau3[dim->v4rhotau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma4; ++ind)
            {
                v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma3lapl; ++ind)
            {
                v4sigma3lapl[dim->v4sigma3lapl * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma3tau; ++ind)
            {
                v4sigma3tau[dim->v4sigma3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2lapl2; ++ind)
            {
                v4sigma2lapl2[dim->v4sigma2lapl2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2lapltau; ++ind)
            {
                v4sigma2lapltau[dim->v4sigma2lapltau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigma2tau2; ++ind)
            {
                v4sigma2tau2[dim->v4sigma2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapl3; ++ind)
            {
                v4sigmalapl3[dim->v4sigmalapl3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapl2tau; ++ind)
            {
                v4sigmalapl2tau[dim->v4sigmalapl2tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmalapltau2; ++ind)
            {
                v4sigmalapltau2[dim->v4sigmalapltau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4sigmatau3; ++ind)
            {
                v4sigmatau3[dim->v4sigmatau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl4; ++ind)
            {
                v4lapl4[dim->v4lapl4 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl3tau; ++ind)
            {
                v4lapl3tau[dim->v4lapl3tau * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapl2tau2; ++ind)
            {
                v4lapl2tau2[dim->v4lapl2tau2 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4lapltau3; ++ind)
            {
                v4lapltau3[dim->v4lapltau3 * g + ind] = 0.0;
            }
            for (int ind = 0; ind < dim->v4tau4; ++ind)
            {
                v4tau4[dim->v4tau4 * g + ind] = 0.0;
            }
        }

        for (const auto& xccomp : _components)
        {
            auto funcptr = xccomp.getFunctionalPointer();

            const auto dim = &(funcptr->dim);

            const auto c = xccomp.getScalingFactor();

            if (xccomp.isLDA())
            {
                xc_lda_lxc(funcptr, grid_batch_size, rho + dim->rho * grid_batch_offset, stage_v4rho4 + dim->v4rho4 * grid_batch_offset);

                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; ++g)
                {
                    for (int ind = 0; ind < dim->v4rho4; ++ind)
                    {
                        v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                    }
                }
            }
            else if (xccomp.isGGA())
            {
                xc_gga_lxc(funcptr,
                           grid_batch_size,
                           rho + dim->rho * grid_batch_offset,
                           sigma + dim->sigma * grid_batch_offset,
                           stage_v4rho4 + dim->v4rho4 * grid_batch_offset,
                           stage_v4rho3sigma + dim->v4rho3sigma * grid_batch_offset,
                           stage_v4rho2sigma2 + dim->v4rho2sigma2 * grid_batch_offset,
                           stage_v4rhosigma3 + dim->v4rhosigma3 * grid_batch_offset,
                           stage_v4sigma4 + dim->v4sigma4 * grid_batch_offset);

                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; ++g)
                {
                    for (int ind = 0; ind < dim->v4rho4; ++ind)
                    {
                        v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho3sigma; ++ind)
                    {
                        v4rho3sigma[dim->v4rho3sigma * g + ind] += c * stage_v4rho3sigma[dim->v4rho3sigma * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2sigma2; ++ind)
                    {
                        v4rho2sigma2[dim->v4rho2sigma2 * g + ind] += c * stage_v4rho2sigma2[dim->v4rho2sigma2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigma3; ++ind)
                    {
                        v4rhosigma3[dim->v4rhosigma3 * g + ind] += c * stage_v4rhosigma3[dim->v4rhosigma3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma4; ++ind)
                    {
                        v4sigma4[dim->v4sigma4 * g + ind] += c * stage_v4sigma4[dim->v4sigma4 * g + ind];
                    }
                }
            }
            else if (xccomp.isMetaGGA())
            {
                xc_mgga_lxc(funcptr,
                            grid_batch_size,
                            rho + dim->rho * grid_batch_offset,
                            sigma + dim->sigma * grid_batch_offset,
                            lapl + dim->lapl * grid_batch_offset,
                            tau + dim->tau * grid_batch_offset,
                            stage_v4rho4 + dim->v4rho4 * grid_batch_offset,
                            stage_v4rho3sigma + dim->v4rho3sigma * grid_batch_offset,
                            stage_v4rho3lapl + dim->v4rho3lapl * grid_batch_offset,
                            stage_v4rho3tau + dim->v4rho3tau * grid_batch_offset,
                            stage_v4rho2sigma2 + dim->v4rho2sigma2 * grid_batch_offset,
                            stage_v4rho2sigmalapl + dim->v4rho2sigmalapl * grid_batch_offset,
                            stage_v4rho2sigmatau + dim->v4rho2sigmatau * grid_batch_offset,
                            stage_v4rho2lapl2 + dim->v4rho2lapl2 * grid_batch_offset,
                            stage_v4rho2lapltau + dim->v4rho2lapltau * grid_batch_offset,
                            stage_v4rho2tau2 + dim->v4rho2tau2 * grid_batch_offset,
                            stage_v4rhosigma3 + dim->v4rhosigma3 * grid_batch_offset,
                            stage_v4rhosigma2lapl + dim->v4rhosigma2lapl * grid_batch_offset,
                            stage_v4rhosigma2tau + dim->v4rhosigma2tau * grid_batch_offset,
                            stage_v4rhosigmalapl2 + dim->v4rhosigmalapl2 * grid_batch_offset,
                            stage_v4rhosigmalapltau + dim->v4rhosigmalapltau * grid_batch_offset,
                            stage_v4rhosigmatau2 + dim->v4rhosigmatau2 * grid_batch_offset,
                            stage_v4rholapl3 + dim->v4rholapl3 * grid_batch_offset,
                            stage_v4rholapl2tau + dim->v4rholapl2tau * grid_batch_offset,
                            stage_v4rholapltau2 + dim->v4rholapltau2 * grid_batch_offset,
                            stage_v4rhotau3 + dim->v4rhotau3 * grid_batch_offset,
                            stage_v4sigma4 + dim->v4sigma4 * grid_batch_offset,
                            stage_v4sigma3lapl + dim->v4sigma3lapl * grid_batch_offset,
                            stage_v4sigma3tau + dim->v4sigma3tau * grid_batch_offset,
                            stage_v4sigma2lapl2 + dim->v4sigma2lapl2 * grid_batch_offset,
                            stage_v4sigma2lapltau + dim->v4sigma2lapltau * grid_batch_offset,
                            stage_v4sigma2tau2 + dim->v4sigma2tau2 * grid_batch_offset,
                            stage_v4sigmalapl3 + dim->v4sigmalapl3 * grid_batch_offset,
                            stage_v4sigmalapl2tau + dim->v4sigmalapl2tau * grid_batch_offset,
                            stage_v4sigmalapltau2 + dim->v4sigmalapltau2 * grid_batch_offset,
                            stage_v4sigmatau3 + dim->v4sigmatau3 * grid_batch_offset,
                            stage_v4lapl4 + dim->v4lapl4 * grid_batch_offset,
                            stage_v4lapl3tau + dim->v4lapl3tau * grid_batch_offset,
                            stage_v4lapl2tau2 + dim->v4lapl2tau2 * grid_batch_offset,
                            stage_v4lapltau3 + dim->v4lapltau3 * grid_batch_offset,
                            stage_v4tau4 + dim->v4tau4 * grid_batch_offset);

                for (int g = grid_batch_offset; g < grid_batch_offset + grid_batch_size; ++g)
                {
                    for (int ind = 0; ind < dim->v4rho4; ++ind)
                    {
                        v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho3sigma; ++ind)
                    {
                        v4rho3sigma[dim->v4rho3sigma * g + ind] += c * stage_v4rho3sigma[dim->v4rho3sigma * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho3lapl; ++ind)
                    {
                        v4rho3lapl[dim->v4rho3lapl * g + ind] += c * stage_v4rho3lapl[dim->v4rho3lapl * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho3tau; ++ind)
                    {
                        v4rho3tau[dim->v4rho3tau * g + ind] += c * stage_v4rho3tau[dim->v4rho3tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2sigma2; ++ind)
                    {
                        v4rho2sigma2[dim->v4rho2sigma2 * g + ind] += c * stage_v4rho2sigma2[dim->v4rho2sigma2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2sigmalapl; ++ind)
                    {
                        v4rho2sigmalapl[dim->v4rho2sigmalapl * g + ind] += c * stage_v4rho2sigmalapl[dim->v4rho2sigmalapl * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2sigmatau; ++ind)
                    {
                        v4rho2sigmatau[dim->v4rho2sigmatau * g + ind] += c * stage_v4rho2sigmatau[dim->v4rho2sigmatau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2lapl2; ++ind)
                    {
                        v4rho2lapl2[dim->v4rho2lapl2 * g + ind] += c * stage_v4rho2lapl2[dim->v4rho2lapl2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2lapltau; ++ind)
                    {
                        v4rho2lapltau[dim->v4rho2lapltau * g + ind] += c * stage_v4rho2lapltau[dim->v4rho2lapltau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rho2tau2; ++ind)
                    {
                        v4rho2tau2[dim->v4rho2tau2 * g + ind] += c * stage_v4rho2tau2[dim->v4rho2tau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigma3; ++ind)
                    {
                        v4rhosigma3[dim->v4rhosigma3 * g + ind] += c * stage_v4rhosigma3[dim->v4rhosigma3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigma2lapl; ++ind)
                    {
                        v4rhosigma2lapl[dim->v4rhosigma2lapl * g + ind] += c * stage_v4rhosigma2lapl[dim->v4rhosigma2lapl * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigma2tau; ++ind)
                    {
                        v4rhosigma2tau[dim->v4rhosigma2tau * g + ind] += c * stage_v4rhosigma2tau[dim->v4rhosigma2tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigmalapl2; ++ind)
                    {
                        v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + ind] += c * stage_v4rhosigmalapl2[dim->v4rhosigmalapl2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigmalapltau; ++ind)
                    {
                        v4rhosigmalapltau[dim->v4rhosigmalapltau * g + ind] += c * stage_v4rhosigmalapltau[dim->v4rhosigmalapltau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhosigmatau2; ++ind)
                    {
                        v4rhosigmatau2[dim->v4rhosigmatau2 * g + ind] += c * stage_v4rhosigmatau2[dim->v4rhosigmatau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rholapl3; ++ind)
                    {
                        v4rholapl3[dim->v4rholapl3 * g + ind] += c * stage_v4rholapl3[dim->v4rholapl3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rholapl2tau; ++ind)
                    {
                        v4rholapl2tau[dim->v4rholapl2tau * g + ind] += c * stage_v4rholapl2tau[dim->v4rholapl2tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rholapltau2; ++ind)
                    {
                        v4rholapltau2[dim->v4rholapltau2 * g + ind] += c * stage_v4rholapltau2[dim->v4rholapltau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4rhotau3; ++ind)
                    {
                        v4rhotau3[dim->v4rhotau3 * g + ind] += c * stage_v4rhotau3[dim->v4rhotau3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma4; ++ind)
                    {
                        v4sigma4[dim->v4sigma4 * g + ind] += c * stage_v4sigma4[dim->v4sigma4 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma3lapl; ++ind)
                    {
                        v4sigma3lapl[dim->v4sigma3lapl * g + ind] += c * stage_v4sigma3lapl[dim->v4sigma3lapl * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma3tau; ++ind)
                    {
                        v4sigma3tau[dim->v4sigma3tau * g + ind] += c * stage_v4sigma3tau[dim->v4sigma3tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma2lapl2; ++ind)
                    {
                        v4sigma2lapl2[dim->v4sigma2lapl2 * g + ind] += c * stage_v4sigma2lapl2[dim->v4sigma2lapl2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma2lapltau; ++ind)
                    {
                        v4sigma2lapltau[dim->v4sigma2lapltau * g + ind] += c * stage_v4sigma2lapltau[dim->v4sigma2lapltau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigma2tau2; ++ind)
                    {
                        v4sigma2tau2[dim->v4sigma2tau2 * g + ind] += c * stage_v4sigma2tau2[dim->v4sigma2tau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigmalapl3; ++ind)
                    {
                        v4sigmalapl3[dim->v4sigmalapl3 * g + ind] += c * stage_v4sigmalapl3[dim->v4sigmalapl3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigmalapl2tau; ++ind)
                    {
                        v4sigmalapl2tau[dim->v4sigmalapl2tau * g + ind] += c * stage_v4sigmalapl2tau[dim->v4sigmalapl2tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigmalapltau2; ++ind)
                    {
                        v4sigmalapltau2[dim->v4sigmalapltau2 * g + ind] += c * stage_v4sigmalapltau2[dim->v4sigmalapltau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4sigmatau3; ++ind)
                    {
                        v4sigmatau3[dim->v4sigmatau3 * g + ind] += c * stage_v4sigmatau3[dim->v4sigmatau3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4lapl4; ++ind)
                    {
                        v4lapl4[dim->v4lapl4 * g + ind] += c * stage_v4lapl4[dim->v4lapl4 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4lapl3tau; ++ind)
                    {
                        v4lapl3tau[dim->v4lapl3tau * g + ind] += c * stage_v4lapl3tau[dim->v4lapl3tau * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4lapl2tau2; ++ind)
                    {
                        v4lapl2tau2[dim->v4lapl2tau2 * g + ind] += c * stage_v4lapl2tau2[dim->v4lapl2tau2 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4lapltau3; ++ind)
                    {
                        v4lapltau3[dim->v4lapltau3 * g + ind] += c * stage_v4lapltau3[dim->v4lapltau3 * g + ind];
                    }
                    for (int ind = 0; ind < dim->v4tau4; ++ind)
                    {
                        v4tau4[dim->v4tau4 * g + ind] += c * stage_v4tau4[dim->v4tau4 * g + ind];
                    }
                }
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

    gridscreen::screenLxcForMGGA(this,
                                 np,
                                 rho,
                                 sigma,
                                 lapl,
                                 tau,
                                 v4rho4,
                                 v4rho3sigma,
                                 v4rho3lapl,
                                 v4rho3tau,
                                 v4rho2sigma2,
                                 v4rho2sigmalapl,
                                 v4rho2sigmatau,
                                 v4rho2lapl2,
                                 v4rho2lapltau,
                                 v4rho2tau2,
                                 v4rhosigma3,
                                 v4rhosigma2lapl,
                                 v4rhosigma2tau,
                                 v4rhosigmalapl2,
                                 v4rhosigmalapltau,
                                 v4rhosigmatau2,
                                 v4rholapl3,
                                 v4rholapl2tau,
                                 v4rholapltau2,
                                 v4rhotau3,
                                 v4sigma4,
                                 v4sigma3lapl,
                                 v4sigma3tau,
                                 v4sigma2lapl2,
                                 v4sigma2lapltau,
                                 v4sigma2tau2,
                                 v4sigmalapl3,
                                 v4sigmalapl2tau,
                                 v4sigmalapltau2,
                                 v4sigmatau3,
                                 v4lapl4,
                                 v4lapl3tau,
                                 v4lapl2tau2,
                                 v4lapltau3,
                                 v4tau4);
}

const xc_func_type*
CXCFunctional::getFunctionalPointerToLdaComponent() const
{
    for (const auto& xccomp : _components)
    {
        if (xccomp.isLDA())
        {
            return xccomp.getFunctionalPointer();
        }
    }

    std::string errmsg("XCFunctional.getFunctionalPointerToLdaComponent: Cannot find LDA functional component");

    errors::assertMsgCritical(false, errmsg);

    return nullptr;
}

const xc_func_type*
CXCFunctional::getFunctionalPointerToGgaComponent() const
{
    for (const auto& xccomp : _components)
    {
        if (xccomp.isGGA())
        {
            return xccomp.getFunctionalPointer();
        }
    }

    std::string errmsg("XCFunctional.getFunctionalPointerToGgaComponent: Cannot find GGA functional component");

    errors::assertMsgCritical(false, errmsg);

    return nullptr;
}

const xc_func_type*
CXCFunctional::getFunctionalPointerToMetaGgaComponent() const
{
    for (const auto& xccomp : _components)
    {
        if (xccomp.isMetaGGA())
        {
            return xccomp.getFunctionalPointer();
        }
    }

    std::string errmsg("XCFunctional.getFunctionalPointerToMetaGgaComponent: Cannot find meta-GGA functional component");

    errors::assertMsgCritical(false, errmsg);

    return nullptr;
}

const int32_t
CXCFunctional::getDimensionOfDerivatives() const
{
    auto indices_and_counts = _getIndicesAndCountsOfDerivatives();

    if (_familyOfFunctional == std::string("LDA"))
    {
        auto final_index = indices_and_counts["v4rho4"][0];
        auto final_count = indices_and_counts["v4rho4"][1];

        return static_cast<int32_t>(final_index + final_count);
    }
    else if (_familyOfFunctional == std::string("GGA"))
    {
        auto final_index = indices_and_counts["v4sigma4"][0];
        auto final_count = indices_and_counts["v4sigma4"][1];

        return static_cast<int32_t>(final_index + final_count);
    }
    else if (_familyOfFunctional == std::string("MGGA"))
    {
        auto final_index = indices_and_counts["v4tau4"][0];
        auto final_count = indices_and_counts["v4tau4"][1];

        return static_cast<int32_t>(final_index + final_count);
    }

    return 0;
}
