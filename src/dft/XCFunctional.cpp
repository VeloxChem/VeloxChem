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
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
//#include <iostream>

#include "ErrorHandler.hpp"
#include "GridScreener.hpp"
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

    for (size_t i = 0; i < labels.size(); i++)
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

        xc_hyb_cam_coef(funcptr, &_rangeSeparationParameterOmega, &_rangeSeparationParameterAlpha, &_rangeSeparationParameterBeta);

        errors::assertMsgCritical(std::fabs(beta) < 1.0e-13, std::string("XCFunctional: Range-separated functional is not yet supported"));

        if (isRangeSeparated())
        {
            errors::assertMsgCritical(_components.size() == 1, "XCFunctional: Range-separated functional can only include one component");
        }

        /*
        auto n_ext_param = xc_func_info_get_n_ext_params(funcptr->info);
        std::cout << "n_ext_param: " << n_ext_param << std::endl;
        for (int64_t i_param = 0; i_param < static_cast<int64_t>(n_ext_param); i_param++)
        {
            auto param_name = xc_func_info_get_ext_params_name(funcptr->info, i_param);
            auto param_descr = xc_func_info_get_ext_params_description(funcptr->info, i_param);
            auto param_value = xc_func_info_get_ext_params_default_value(funcptr->info, i_param);
            std::cout << "  param id: " << i_param << std::endl;
            std::cout << "        name: " << param_name << std::endl;
            std::cout << "        descr: " << param_descr << std::endl;
            std::cout << "        value: " << param_value << std::endl;
        }
        */
    }

    if (hasExc) _maxDerivOrder = 0;
    if (hasVxc) _maxDerivOrder = 1;
    if (hasFxc) _maxDerivOrder = 2;
    if (hasKxc) _maxDerivOrder = 3;
    if (hasLxc) _maxDerivOrder = 4;

    if (isLDA) _familyOfFunctional = xcfun::lda;
    if (isGGA) _familyOfFunctional = xcfun::gga;
    if (isMGGA) _familyOfFunctional = xcfun::mgga;

    errors::assertMsgCritical(!needLaplacian, std::string("XCFunctional: Density Laplacian is not supported"));

    _allocateStagingBuffer();
}

CXCFunctional::CXCFunctional(const CXCFunctional& source)

    : _nameOfFunctional(source._nameOfFunctional)

    , _fractionOfExactExchange(source._fractionOfExactExchange)

    , _rangeSeparationParameterAlpha(source._rangeSeparationParameterAlpha)

    , _rangeSeparationParameterBeta(source._rangeSeparationParameterBeta)

    , _rangeSeparationParameterOmega(source._rangeSeparationParameterOmega)

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

    , _rangeSeparationParameterAlpha(std::move(source._rangeSeparationParameterAlpha))

    , _rangeSeparationParameterBeta(std::move(source._rangeSeparationParameterBeta))

    , _rangeSeparationParameterOmega(std::move(source._rangeSeparationParameterOmega))

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

        _stagingBuffer = (double*)std::malloc(sizeof(double) * n_xc_outputs * _ldStaging);
    }
}

void
CXCFunctional::_freeStagingBuffer()
{
    if (_stagingBuffer != nullptr)
    {
        std::free(_stagingBuffer);

        _stagingBuffer = nullptr;
    }
}

std::unordered_map<std::string, std::array<int64_t, 2>>
CXCFunctional::_getIndicesAndCountsOfDerivatives() const
{
    std::unordered_map<std::string, std::array<int64_t, 2>> indices_and_counts;

    if (_familyOfFunctional == xcfun::lda)
    {
        auto       ldafunc = getFunctionalPointerToLdaComponent();
        const auto dim     = &(ldafunc->dim);

        int64_t n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int64_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int64_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;

        indices_and_counts["v2rho2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;

        indices_and_counts["v3rho3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;

        indices_and_counts["v4rho4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
    }
    else if (_familyOfFunctional == xcfun::gga)
    {
        auto       ggafunc = getFunctionalPointerToGgaComponent();
        const auto dim     = &(ggafunc->dim);

        int64_t n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int64_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int64_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;
        indices_and_counts["vsigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->vsigma});
        n_xc_outputs += dim->vsigma;

        indices_and_counts["v2rho2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;
        indices_and_counts["v2rhosigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rhosigma});
        n_xc_outputs += dim->v2rhosigma;
        indices_and_counts["v2sigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2sigma2});
        n_xc_outputs += dim->v2sigma2;

        indices_and_counts["v3rho3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;
        indices_and_counts["v3rho2sigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho2sigma});
        n_xc_outputs += dim->v3rho2sigma;
        indices_and_counts["v3rhosigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rhosigma2});
        n_xc_outputs += dim->v3rhosigma2;
        indices_and_counts["v3sigma3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigma3});
        n_xc_outputs += dim->v3sigma3;

        indices_and_counts["v4rho4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
        indices_and_counts["v4rho3sigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho3sigma});
        n_xc_outputs += dim->v4rho3sigma;
        indices_and_counts["v4rho2sigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2sigma2});
        n_xc_outputs += dim->v4rho2sigma2;
        indices_and_counts["v4rhosigma3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigma3});
        n_xc_outputs += dim->v4rhosigma3;
        indices_and_counts["v4sigma4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma4});
        n_xc_outputs += dim->v4sigma4;
    }
    else if (_familyOfFunctional == xcfun::mgga)
    {
        auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
        const auto dim      = &(mggafunc->dim);

        int64_t n_xc_outputs = 0;

        indices_and_counts["zk"] = std::array<int64_t, 2>({n_xc_outputs, dim->zk});
        n_xc_outputs += dim->zk;

        indices_and_counts["vrho"] = std::array<int64_t, 2>({n_xc_outputs, dim->vrho});
        n_xc_outputs += dim->vrho;
        indices_and_counts["vsigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->vsigma});
        n_xc_outputs += dim->vsigma;
        indices_and_counts["vlapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->vlapl});
        n_xc_outputs += dim->vlapl;
        indices_and_counts["vtau"] = std::array<int64_t, 2>({n_xc_outputs, dim->vtau});
        n_xc_outputs += dim->vtau;

        indices_and_counts["v2rho2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rho2});
        n_xc_outputs += dim->v2rho2;
        indices_and_counts["v2rhosigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rhosigma});
        n_xc_outputs += dim->v2rhosigma;
        indices_and_counts["v2rholapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rholapl});
        n_xc_outputs += dim->v2rholapl;
        indices_and_counts["v2rhotau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2rhotau});
        n_xc_outputs += dim->v2rhotau;
        indices_and_counts["v2sigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2sigma2});
        n_xc_outputs += dim->v2sigma2;
        indices_and_counts["v2sigmalapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2sigmalapl});
        n_xc_outputs += dim->v2sigmalapl;
        indices_and_counts["v2sigmatau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2sigmatau});
        n_xc_outputs += dim->v2sigmatau;
        indices_and_counts["v2lapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2lapl2});
        n_xc_outputs += dim->v2lapl2;
        indices_and_counts["v2lapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2lapltau});
        n_xc_outputs += dim->v2lapltau;
        indices_and_counts["v2tau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v2tau2});
        n_xc_outputs += dim->v2tau2;

        indices_and_counts["v3rho3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho3});
        n_xc_outputs += dim->v3rho3;
        indices_and_counts["v3rho2sigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho2sigma});
        n_xc_outputs += dim->v3rho2sigma;
        indices_and_counts["v3rho2lapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho2lapl});
        n_xc_outputs += dim->v3rho2lapl;
        indices_and_counts["v3rho2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rho2tau});
        n_xc_outputs += dim->v3rho2tau;
        indices_and_counts["v3rhosigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rhosigma2});
        n_xc_outputs += dim->v3rhosigma2;
        indices_and_counts["v3rhosigmalapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rhosigmalapl});
        n_xc_outputs += dim->v3rhosigmalapl;
        indices_and_counts["v3rhosigmatau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rhosigmatau});
        n_xc_outputs += dim->v3rhosigmatau;
        indices_and_counts["v3rholapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rholapl2});
        n_xc_outputs += dim->v3rholapl2;
        indices_and_counts["v3rholapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rholapltau});
        n_xc_outputs += dim->v3rholapltau;
        indices_and_counts["v3rhotau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3rhotau2});
        n_xc_outputs += dim->v3rhotau2;
        indices_and_counts["v3sigma3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigma3});
        n_xc_outputs += dim->v3sigma3;
        indices_and_counts["v3sigma2lapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigma2lapl});
        n_xc_outputs += dim->v3sigma2lapl;
        indices_and_counts["v3sigma2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigma2tau});
        n_xc_outputs += dim->v3sigma2tau;
        indices_and_counts["v3sigmalapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigmalapl2});
        n_xc_outputs += dim->v3sigmalapl2;
        indices_and_counts["v3sigmalapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigmalapltau});
        n_xc_outputs += dim->v3sigmalapltau;
        indices_and_counts["v3sigmatau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3sigmatau2});
        n_xc_outputs += dim->v3sigmatau2;
        indices_and_counts["v3lapl3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3lapl3});
        n_xc_outputs += dim->v3lapl3;
        indices_and_counts["v3lapl2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3lapl2tau});
        n_xc_outputs += dim->v3lapl2tau;
        indices_and_counts["v3lapltau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3lapltau2});
        n_xc_outputs += dim->v3lapltau2;
        indices_and_counts["v3tau3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v3tau3});
        n_xc_outputs += dim->v3tau3;

        indices_and_counts["v4rho4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho4});
        n_xc_outputs += dim->v4rho4;
        indices_and_counts["v4rho3sigma"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho3sigma});
        n_xc_outputs += dim->v4rho3sigma;
        indices_and_counts["v4rho3lapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho3lapl});
        n_xc_outputs += dim->v4rho3lapl;
        indices_and_counts["v4rho3tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho3tau});
        n_xc_outputs += dim->v4rho3tau;
        indices_and_counts["v4rho2sigma2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2sigma2});
        n_xc_outputs += dim->v4rho2sigma2;
        indices_and_counts["v4rho2sigmalapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2sigmalapl});
        n_xc_outputs += dim->v4rho2sigmalapl;
        indices_and_counts["v4rho2sigmatau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2sigmatau});
        n_xc_outputs += dim->v4rho2sigmatau;
        indices_and_counts["v4rho2lapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2lapl2});
        n_xc_outputs += dim->v4rho2lapl2;
        indices_and_counts["v4rho2lapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2lapltau});
        n_xc_outputs += dim->v4rho2lapltau;
        indices_and_counts["v4rho2tau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rho2tau2});
        n_xc_outputs += dim->v4rho2tau2;
        indices_and_counts["v4rhosigma3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigma3});
        n_xc_outputs += dim->v4rhosigma3;
        indices_and_counts["v4rhosigma2lapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigma2lapl});
        n_xc_outputs += dim->v4rhosigma2lapl;
        indices_and_counts["v4rhosigma2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigma2tau});
        n_xc_outputs += dim->v4rhosigma2tau;
        indices_and_counts["v4rhosigmalapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigmalapl2});
        n_xc_outputs += dim->v4rhosigmalapl2;
        indices_and_counts["v4rhosigmalapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigmalapltau});
        n_xc_outputs += dim->v4rhosigmalapltau;
        indices_and_counts["v4rhosigmatau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhosigmatau2});
        n_xc_outputs += dim->v4rhosigmatau2;
        indices_and_counts["v4rholapl3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rholapl3});
        n_xc_outputs += dim->v4rholapl3;
        indices_and_counts["v4rholapl2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rholapl2tau});
        n_xc_outputs += dim->v4rholapl2tau;
        indices_and_counts["v4rholapltau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rholapltau2});
        n_xc_outputs += dim->v4rholapltau2;
        indices_and_counts["v4rhotau3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4rhotau3});
        n_xc_outputs += dim->v4rhotau3;
        indices_and_counts["v4sigma4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma4});
        n_xc_outputs += dim->v4sigma4;
        indices_and_counts["v4sigma3lapl"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma3lapl});
        n_xc_outputs += dim->v4sigma3lapl;
        indices_and_counts["v4sigma3tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma3tau});
        n_xc_outputs += dim->v4sigma3tau;
        indices_and_counts["v4sigma2lapl2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma2lapl2});
        n_xc_outputs += dim->v4sigma2lapl2;
        indices_and_counts["v4sigma2lapltau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma2lapltau});
        n_xc_outputs += dim->v4sigma2lapltau;
        indices_and_counts["v4sigma2tau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigma2tau2});
        n_xc_outputs += dim->v4sigma2tau2;
        indices_and_counts["v4sigmalapl3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigmalapl3});
        n_xc_outputs += dim->v4sigmalapl3;
        indices_and_counts["v4sigmalapl2tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigmalapl2tau});
        n_xc_outputs += dim->v4sigmalapl2tau;
        indices_and_counts["v4sigmalapltau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigmalapltau2});
        n_xc_outputs += dim->v4sigmalapltau2;
        indices_and_counts["v4sigmatau3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4sigmatau3});
        n_xc_outputs += dim->v4sigmatau3;
        indices_and_counts["v4lapl4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4lapl4});
        n_xc_outputs += dim->v4lapl4;
        indices_and_counts["v4lapl3tau"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4lapl3tau});
        n_xc_outputs += dim->v4lapl3tau;
        indices_and_counts["v4lapl2tau2"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4lapl2tau2});
        n_xc_outputs += dim->v4lapl2tau2;
        indices_and_counts["v4lapltau3"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4lapltau3});
        n_xc_outputs += dim->v4lapltau3;
        indices_and_counts["v4tau4"] = std::array<int64_t, 2>({n_xc_outputs, dim->v4tau4});
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

    _rangeSeparationParameterAlpha = source._rangeSeparationParameterAlpha;

    _rangeSeparationParameterBeta = source._rangeSeparationParameterBeta;

    _rangeSeparationParameterOmega = source._rangeSeparationParameterOmega;

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

    _rangeSeparationParameterAlpha = std::move(source._rangeSeparationParameterAlpha);

    _rangeSeparationParameterBeta = std::move(source._rangeSeparationParameterBeta);

    _rangeSeparationParameterOmega = std::move(source._rangeSeparationParameterOmega);

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

    if (_rangeSeparationParameterAlpha != other._rangeSeparationParameterAlpha) return false;

    if (_rangeSeparationParameterBeta != other._rangeSeparationParameterBeta) return false;

    if (_rangeSeparationParameterOmega != other._rangeSeparationParameterOmega) return false;

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
    return _familyOfFunctional;
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

bool
CXCFunctional::isRangeSeparated() const
{
    return (std::fabs(_rangeSeparationParameterBeta) > 1.0e-13);
}

double
CXCFunctional::getRangeSeparationParameterAlpha() const
{
    return _rangeSeparationParameterAlpha;
}

double
CXCFunctional::getRangeSeparationParameterBeta() const
{
    return _rangeSeparationParameterBeta;
}

double
CXCFunctional::getRangeSeparationParameterOmega() const
{
    return _rangeSeparationParameterOmega;
}

auto
CXCFunctional::compute_exc_vxc_for_lda(const int64_t np, const double* rho, double* zk, double* vrho) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_zk   = nullptr;
    double* stage_vrho = nullptr;

    if (alloc)
    {
        stage_zk   = (double*)std::malloc(sizeof(double) * ind_cnt["zk"][1] * np);
        stage_vrho = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
    }
    else
    {
        stage_zk   = &_stagingBuffer[ind_cnt["zk"][0] * _ldStaging];
        stage_vrho = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ldafunc = getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->zk; ++ind)
        {
            zk[dim->zk * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_exc_vxc(funcptr, np, rho, stage_zk, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
    }

    gridscreen::screenExcVxcForLDA(this, np, rho, zk, vrho);

    if (alloc)
    {
        std::free(stage_zk);
        std::free(stage_vrho);
    }
}

auto
CXCFunctional::compute_vxc_for_lda(const int64_t np, const double* rho, double* vrho) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_vrho = nullptr;

    if (alloc)
    {
        stage_vrho = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
    }
    else
    {
        stage_vrho = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ldafunc = getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
    }

    gridscreen::screenVxcForLDA(this, np, rho, vrho);

    if (alloc)
    {
        std::free(stage_vrho);
    }
}

auto
CXCFunctional::compute_fxc_for_lda(const int64_t np, const double* rho, double* v2rho2) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v2rho2 = nullptr;

    if (alloc)
    {
        stage_v2rho2 = (double*)std::malloc(sizeof(double) * ind_cnt["v2rho2"][1] * np);
    }
    else
    {
        stage_v2rho2 = &_stagingBuffer[ind_cnt["v2rho2"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ldafunc = getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v2rho2; ++ind)
        {
            v2rho2[dim->v2rho2 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_fxc(funcptr, np, rho, stage_v2rho2);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
            }
        }
    }

    gridscreen::screenFxcForLDA(this, np, rho, v2rho2);

    if (alloc)
    {
        std::free(stage_v2rho2);
    }
}

auto
CXCFunctional::compute_kxc_for_lda(const int64_t np, const double* rho, double* v3rho3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v3rho3 = nullptr;

    if (alloc)
    {
        stage_v3rho3 = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho3"][1] * np);
    }
    else
    {
        stage_v3rho3 = &_stagingBuffer[ind_cnt["v3rho3"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ldafunc = getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v3rho3; ++ind)
        {
            v3rho3[dim->v3rho3 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_kxc(funcptr, np, rho, stage_v3rho3);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
            }
        }
    }

    gridscreen::screenKxcForLDA(this, np, rho, v3rho3);

    if (alloc)
    {
        std::free(stage_v3rho3);
    }
}

auto
CXCFunctional::compute_lxc_for_lda(const int64_t np, const double* rho, double* v4rho4) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 4,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Lxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v4rho4 = nullptr;

    if (alloc)
    {
        stage_v4rho4 = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho4"][1] * np);
    }
    else
    {
        stage_v4rho4 = &_stagingBuffer[ind_cnt["v4rho4"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ldafunc = getFunctionalPointerToLdaComponent();
    const auto dim     = &(ldafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v4rho4; ++ind)
        {
            v4rho4[dim->v4rho4 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v4rho4; ++ind)
                {
                    v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                }
            }
        }
    }

    gridscreen::screenLxcForLDA(this, np, rho, v4rho4);

    if (alloc)
    {
        std::free(stage_v4rho4);
    }
}

auto
CXCFunctional::compute_exc_vxc_for_gga(const int64_t np, const double* rho, const double* sigma, double* zk, double* vrho, double* vsigma) const
    -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_zk     = nullptr;
    double* stage_vrho   = nullptr;
    double* stage_vsigma = nullptr;

    if (alloc)
    {
        stage_zk     = (double*)std::malloc(sizeof(double) * ind_cnt["zk"][1] * np);
        stage_vrho   = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
        stage_vsigma = (double*)std::malloc(sizeof(double) * ind_cnt["vsigma"][1] * np);
    }
    else
    {
        stage_zk     = &_stagingBuffer[ind_cnt["zk"][0] * _ldStaging];
        stage_vrho   = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
        stage_vsigma = &_stagingBuffer[ind_cnt["vsigma"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ggafunc = getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->zk; ++ind)
        {
            zk[dim->zk * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vsigma; ++ind)
        {
            vsigma[dim->vsigma * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_exc_vxc(funcptr, np, rho, stage_zk, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_exc_vxc(funcptr, np, rho, sigma, stage_zk, stage_vrho, stage_vsigma);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
            }
        }
    }

    gridscreen::screenExcVxcForGGA(this, np, rho, sigma, zk, vrho, vsigma);

    if (alloc)
    {
        std::free(stage_zk);
        std::free(stage_vrho);
        std::free(stage_vsigma);
    }
}

auto
CXCFunctional::compute_vxc_for_gga(const int64_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Vxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_vrho   = nullptr;
    double* stage_vsigma = nullptr;

    if (alloc)
    {
        stage_vrho   = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
        stage_vsigma = (double*)std::malloc(sizeof(double) * ind_cnt["vsigma"][1] * np);
    }
    else
    {
        stage_vrho   = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
        stage_vsigma = &_stagingBuffer[ind_cnt["vsigma"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ggafunc = getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vsigma; ++ind)
        {
            vsigma[dim->vsigma * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_vxc(funcptr, np, rho, sigma, stage_vrho, stage_vsigma);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
            }
        }
    }

    gridscreen::screenVxcForGGA(this, np, rho, sigma, vrho, vsigma);

    if (alloc)
    {
        std::free(stage_vrho);
        std::free(stage_vsigma);
    }
}

auto
CXCFunctional::compute_fxc_for_gga(const int64_t np, const double* rho, const double* sigma, double* v2rho2, double* v2rhosigma, double* v2sigma2)
    const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 2,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Fxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v2rho2     = nullptr;
    double* stage_v2rhosigma = nullptr;
    double* stage_v2sigma2   = nullptr;

    if (alloc)
    {
        stage_v2rho2     = (double*)std::malloc(sizeof(double) * ind_cnt["v2rho2"][1] * np);
        stage_v2rhosigma = (double*)std::malloc(sizeof(double) * ind_cnt["v2rhosigma"][1] * np);
        stage_v2sigma2   = (double*)std::malloc(sizeof(double) * ind_cnt["v2sigma2"][1] * np);
    }
    else
    {
        stage_v2rho2     = &_stagingBuffer[ind_cnt["v2rho2"][0] * _ldStaging];
        stage_v2rhosigma = &_stagingBuffer[ind_cnt["v2rhosigma"][0] * _ldStaging];
        stage_v2sigma2   = &_stagingBuffer[ind_cnt["v2sigma2"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ggafunc = getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v2rho2; ++ind)
        {
            v2rho2[dim->v2rho2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2rhosigma; ++ind)
        {
            v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2sigma2; ++ind)
        {
            v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_fxc(funcptr, np, rho, stage_v2rho2);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_fxc(funcptr, np, rho, sigma, stage_v2rho2, stage_v2rhosigma, stage_v2sigma2);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
                for (int ind = 0; ind < dim->v2rhosigma; ++ind)
                {
                    v2rhosigma[dim->v2rhosigma * g + ind] += c * stage_v2rhosigma[dim->v2rhosigma * g + ind];
                }
                for (int ind = 0; ind < dim->v2sigma2; ++ind)
                {
                    v2sigma2[dim->v2sigma2 * g + ind] += c * stage_v2sigma2[dim->v2sigma2 * g + ind];
                }
            }
        }
    }

    gridscreen::screenFxcForGGA(this, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2);

    if (alloc)
    {
        std::free(stage_v2rho2);
        std::free(stage_v2rhosigma);
        std::free(stage_v2sigma2);
    }
}

auto
CXCFunctional::compute_kxc_for_gga(const int64_t np,
                                   const double* rho,
                                   const double* sigma,
                                   double*       v3rho3,
                                   double*       v3rho2sigma,
                                   double*       v3rhosigma2,
                                   double*       v3sigma3) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 3,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Kxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v3rho3      = nullptr;
    double* stage_v3rho2sigma = nullptr;
    double* stage_v3rhosigma2 = nullptr;
    double* stage_v3sigma3    = nullptr;

    if (alloc)
    {
        stage_v3rho3      = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho3"][1] * np);
        stage_v3rho2sigma = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho2sigma"][1] * np);
        stage_v3rhosigma2 = (double*)std::malloc(sizeof(double) * ind_cnt["v3rhosigma2"][1] * np);
        stage_v3sigma3    = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigma3"][1] * np);
    }
    else
    {
        stage_v3rho3      = &_stagingBuffer[ind_cnt["v3rho3"][0] * _ldStaging];
        stage_v3rho2sigma = &_stagingBuffer[ind_cnt["v3rho2sigma"][0] * _ldStaging];
        stage_v3rhosigma2 = &_stagingBuffer[ind_cnt["v3rhosigma2"][0] * _ldStaging];
        stage_v3sigma3    = &_stagingBuffer[ind_cnt["v3sigma3"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ggafunc = getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v3rho3; ++ind)
        {
            v3rho3[dim->v3rho3 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rho2sigma; ++ind)
        {
            v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rhosigma2; ++ind)
        {
            v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigma3; ++ind)
        {
            v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_kxc(funcptr, np, rho, stage_v3rho3);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_kxc(funcptr, np, rho, sigma, stage_v3rho3, stage_v3rho2sigma, stage_v3rhosigma2, stage_v3sigma3);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
                for (int ind = 0; ind < dim->v3rho2sigma; ++ind)
                {
                    v3rho2sigma[dim->v3rho2sigma * g + ind] += c * stage_v3rho2sigma[dim->v3rho2sigma * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhosigma2; ++ind)
                {
                    v3rhosigma2[dim->v3rhosigma2 * g + ind] += c * stage_v3rhosigma2[dim->v3rhosigma2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigma3; ++ind)
                {
                    v3sigma3[dim->v3sigma3 * g + ind] += c * stage_v3sigma3[dim->v3sigma3 * g + ind];
                }
            }
        }
    }

    gridscreen::screenKxcForGGA(this, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

    if (alloc)
    {
        std::free(stage_v3rho3);
        std::free(stage_v3rho2sigma);
        std::free(stage_v3rhosigma2);
        std::free(stage_v3sigma3);
    }
}

auto
CXCFunctional::compute_lxc_for_gga(const int64_t np,
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

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v4rho4       = nullptr;
    double* stage_v4rho3sigma  = nullptr;
    double* stage_v4rho2sigma2 = nullptr;
    double* stage_v4rhosigma3  = nullptr;
    double* stage_v4sigma4     = nullptr;

    if (alloc)
    {
        stage_v4rho4       = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho4"][1] * np);
        stage_v4rho3sigma  = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho3sigma"][1] * np);
        stage_v4rho2sigma2 = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2sigma2"][1] * np);
        stage_v4rhosigma3  = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigma3"][1] * np);
        stage_v4sigma4     = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma4"][1] * np);
    }
    else
    {
        stage_v4rho4       = &_stagingBuffer[ind_cnt["v4rho4"][0] * _ldStaging];
        stage_v4rho3sigma  = &_stagingBuffer[ind_cnt["v4rho3sigma"][0] * _ldStaging];
        stage_v4rho2sigma2 = &_stagingBuffer[ind_cnt["v4rho2sigma2"][0] * _ldStaging];
        stage_v4rhosigma3  = &_stagingBuffer[ind_cnt["v4rhosigma3"][0] * _ldStaging];
        stage_v4sigma4     = &_stagingBuffer[ind_cnt["v4sigma4"][0] * _ldStaging];
    }

    // compute derivatives

    auto       ggafunc = getFunctionalPointerToGgaComponent();
    const auto dim     = &(ggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v4rho4; ++ind)
        {
            v4rho4[dim->v4rho4 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v4rho3sigma; ++ind)
        {
            v4rho3sigma[dim->v4rho3sigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v4rho2sigma2; ++ind)
        {
            v4rho2sigma2[dim->v4rho2sigma2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v4rhosigma3; ++ind)
        {
            v4rhosigma3[dim->v4rhosigma3 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v4sigma4; ++ind)
        {
            v4sigma4[dim->v4sigma4 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v4rho4; ++ind)
                {
                    v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_lxc(funcptr, np, rho, sigma, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4);

            for (int64_t g = 0; g < np; ++g)
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
    }

    gridscreen::screenLxcForGGA(this, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4);

    if (alloc)
    {
        std::free(stage_v4rho4);
        std::free(stage_v4rho3sigma);
        std::free(stage_v4rho2sigma2);
        std::free(stage_v4rhosigma3);
        std::free(stage_v4sigma4);
    }
}

auto
CXCFunctional::compute_exc_vxc_for_mgga(const int64_t np,
                                        const double* rho,
                                        const double* sigma,
                                        const double* lapl,
                                        const double* tau,
                                        double*       zk,
                                        double*       vrho,
                                        double*       vsigma,
                                        double*       vlapl,
                                        double*       vtau) const -> void
{
    errors::assertMsgCritical(_maxDerivOrder >= 1,
                              std::string(__func__) + ": exchange-correlation functional does not provide evaluators for Exc and Vxc on grid");

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_zk     = nullptr;
    double* stage_vrho   = nullptr;
    double* stage_vsigma = nullptr;
    double* stage_vlapl  = nullptr;
    double* stage_vtau   = nullptr;

    if (alloc)
    {
        stage_zk     = (double*)std::malloc(sizeof(double) * ind_cnt["zk"][1] * np);
        stage_vrho   = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
        stage_vsigma = (double*)std::malloc(sizeof(double) * ind_cnt["vsigma"][1] * np);
        stage_vlapl  = (double*)std::malloc(sizeof(double) * ind_cnt["vlapl"][1] * np);
        stage_vtau   = (double*)std::malloc(sizeof(double) * ind_cnt["vtau"][1] * np);
    }
    else
    {
        stage_zk     = &_stagingBuffer[ind_cnt["zk"][0] * _ldStaging];
        stage_vrho   = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
        stage_vsigma = &_stagingBuffer[ind_cnt["vsigma"][0] * _ldStaging];
        stage_vlapl  = &_stagingBuffer[ind_cnt["vlapl"][0] * _ldStaging];
        stage_vtau   = &_stagingBuffer[ind_cnt["vtau"][0] * _ldStaging];
    }

    // compute derivatives

    auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->zk; ++ind)
        {
            zk[dim->zk * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vsigma; ++ind)
        {
            vsigma[dim->vsigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vlapl; ++ind)
        {
            vlapl[dim->vlapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vtau; ++ind)
        {
            vtau[dim->vtau * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_exc_vxc(funcptr, np, rho, stage_zk, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_exc_vxc(funcptr, np, rho, sigma, stage_zk, stage_vrho, stage_vsigma);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
            }
        }
        else if (xccomp.isMetaGGA())
        {
            xc_mgga_exc_vxc(funcptr, np, rho, sigma, lapl, tau, stage_zk, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->zk; ++ind)
                {
                    zk[dim->zk * g + ind] += c * stage_zk[dim->zk * g + ind];
                }
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
                for (int ind = 0; ind < dim->vlapl; ++ind)
                {
                    vlapl[dim->vlapl * g + ind] += c * stage_vlapl[dim->vlapl * g + ind];
                }
                for (int ind = 0; ind < dim->vtau; ++ind)
                {
                    vtau[dim->vtau * g + ind] += c * stage_vtau[dim->vtau * g + ind];
                }
            }
        }
    }

    gridscreen::screenExcVxcForMGGA(this, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau);

    if (alloc)
    {
        std::free(stage_zk);
        std::free(stage_vrho);
        std::free(stage_vsigma);
        std::free(stage_vlapl);
        std::free(stage_vtau);
    }
}

auto
CXCFunctional::compute_vxc_for_mgga(const int64_t np,
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

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_vrho   = nullptr;
    double* stage_vsigma = nullptr;
    double* stage_vlapl  = nullptr;
    double* stage_vtau   = nullptr;

    if (alloc)
    {
        stage_vrho   = (double*)std::malloc(sizeof(double) * ind_cnt["vrho"][1] * np);
        stage_vsigma = (double*)std::malloc(sizeof(double) * ind_cnt["vsigma"][1] * np);
        stage_vlapl  = (double*)std::malloc(sizeof(double) * ind_cnt["vlapl"][1] * np);
        stage_vtau   = (double*)std::malloc(sizeof(double) * ind_cnt["vtau"][1] * np);
    }
    else
    {
        stage_vrho   = &_stagingBuffer[ind_cnt["vrho"][0] * _ldStaging];
        stage_vsigma = &_stagingBuffer[ind_cnt["vsigma"][0] * _ldStaging];
        stage_vlapl  = &_stagingBuffer[ind_cnt["vlapl"][0] * _ldStaging];
        stage_vtau   = &_stagingBuffer[ind_cnt["vtau"][0] * _ldStaging];
    }

    // compute derivatives

    auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->vrho; ++ind)
        {
            vrho[dim->vrho * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vsigma; ++ind)
        {
            vsigma[dim->vsigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vlapl; ++ind)
        {
            vlapl[dim->vlapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->vtau; ++ind)
        {
            vtau[dim->vtau * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_vxc(funcptr, np, rho, stage_vrho);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_vxc(funcptr, np, rho, sigma, stage_vrho, stage_vsigma);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
            }
        }
        else if (xccomp.isMetaGGA())
        {
            xc_mgga_vxc(funcptr, np, rho, sigma, lapl, tau, stage_vrho, stage_vsigma, stage_vlapl, stage_vtau);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->vrho; ++ind)
                {
                    vrho[dim->vrho * g + ind] += c * stage_vrho[dim->vrho * g + ind];
                }
                for (int ind = 0; ind < dim->vsigma; ++ind)
                {
                    vsigma[dim->vsigma * g + ind] += c * stage_vsigma[dim->vsigma * g + ind];
                }
                for (int ind = 0; ind < dim->vlapl; ++ind)
                {
                    vlapl[dim->vlapl * g + ind] += c * stage_vlapl[dim->vlapl * g + ind];
                }
                for (int ind = 0; ind < dim->vtau; ++ind)
                {
                    vtau[dim->vtau * g + ind] += c * stage_vtau[dim->vtau * g + ind];
                }
            }
        }
    }

    gridscreen::screenVxcForMGGA(this, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau);

    if (alloc)
    {
        std::free(stage_vrho);
        std::free(stage_vsigma);
        std::free(stage_vlapl);
        std::free(stage_vtau);
    }
}

auto
CXCFunctional::compute_fxc_for_mgga(const int64_t np,
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

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v2rho2      = nullptr;
    double* stage_v2rhosigma  = nullptr;
    double* stage_v2rholapl   = nullptr;
    double* stage_v2rhotau    = nullptr;
    double* stage_v2sigma2    = nullptr;
    double* stage_v2sigmalapl = nullptr;
    double* stage_v2sigmatau  = nullptr;
    double* stage_v2lapl2     = nullptr;
    double* stage_v2lapltau   = nullptr;
    double* stage_v2tau2      = nullptr;

    if (alloc)
    {
        stage_v2rho2      = (double*)std::malloc(sizeof(double) * ind_cnt["v2rho2"][1] * np);
        stage_v2rhosigma  = (double*)std::malloc(sizeof(double) * ind_cnt["v2rhosigma"][1] * np);
        stage_v2rholapl   = (double*)std::malloc(sizeof(double) * ind_cnt["v2rholapl"][1] * np);
        stage_v2rhotau    = (double*)std::malloc(sizeof(double) * ind_cnt["v2rhotau"][1] * np);
        stage_v2sigma2    = (double*)std::malloc(sizeof(double) * ind_cnt["v2sigma2"][1] * np);
        stage_v2sigmalapl = (double*)std::malloc(sizeof(double) * ind_cnt["v2sigmalapl"][1] * np);
        stage_v2sigmatau  = (double*)std::malloc(sizeof(double) * ind_cnt["v2sigmatau"][1] * np);
        stage_v2lapl2     = (double*)std::malloc(sizeof(double) * ind_cnt["v2lapl2"][1] * np);
        stage_v2lapltau   = (double*)std::malloc(sizeof(double) * ind_cnt["v2lapltau"][1] * np);
        stage_v2tau2      = (double*)std::malloc(sizeof(double) * ind_cnt["v2tau2"][1] * np);
    }
    else
    {
        stage_v2rho2      = &_stagingBuffer[ind_cnt["v2rho2"][0] * _ldStaging];
        stage_v2rhosigma  = &_stagingBuffer[ind_cnt["v2rhosigma"][0] * _ldStaging];
        stage_v2rholapl   = &_stagingBuffer[ind_cnt["v2rholapl"][0] * _ldStaging];
        stage_v2rhotau    = &_stagingBuffer[ind_cnt["v2rhotau"][0] * _ldStaging];
        stage_v2sigma2    = &_stagingBuffer[ind_cnt["v2sigma2"][0] * _ldStaging];
        stage_v2sigmalapl = &_stagingBuffer[ind_cnt["v2sigmalapl"][0] * _ldStaging];
        stage_v2sigmatau  = &_stagingBuffer[ind_cnt["v2sigmatau"][0] * _ldStaging];
        stage_v2lapl2     = &_stagingBuffer[ind_cnt["v2lapl2"][0] * _ldStaging];
        stage_v2lapltau   = &_stagingBuffer[ind_cnt["v2lapltau"][0] * _ldStaging];
        stage_v2tau2      = &_stagingBuffer[ind_cnt["v2tau2"][0] * _ldStaging];
    }

    // compute derivatives

    auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v2rho2; ++ind)
        {
            v2rho2[dim->v2rho2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2rhosigma; ++ind)
        {
            v2rhosigma[dim->v2rhosigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2rholapl; ++ind)
        {
            v2rholapl[dim->v2rholapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2rhotau; ++ind)
        {
            v2rhotau[dim->v2rhotau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2sigma2; ++ind)
        {
            v2sigma2[dim->v2sigma2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2sigmalapl; ++ind)
        {
            v2sigmalapl[dim->v2sigmalapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2sigmatau; ++ind)
        {
            v2sigmatau[dim->v2sigmatau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2lapl2; ++ind)
        {
            v2lapl2[dim->v2lapl2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2lapltau; ++ind)
        {
            v2lapltau[dim->v2lapltau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v2tau2; ++ind)
        {
            v2tau2[dim->v2tau2 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_fxc(funcptr, np, rho, stage_v2rho2);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_fxc(funcptr, np, rho, sigma, stage_v2rho2, stage_v2rhosigma, stage_v2sigma2);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
                for (int ind = 0; ind < dim->v2rhosigma; ++ind)
                {
                    v2rhosigma[dim->v2rhosigma * g + ind] += c * stage_v2rhosigma[dim->v2rhosigma * g + ind];
                }
                for (int ind = 0; ind < dim->v2sigma2; ++ind)
                {
                    v2sigma2[dim->v2sigma2 * g + ind] += c * stage_v2sigma2[dim->v2sigma2 * g + ind];
                }
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

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v2rho2; ++ind)
                {
                    v2rho2[dim->v2rho2 * g + ind] += c * stage_v2rho2[dim->v2rho2 * g + ind];
                }
                for (int ind = 0; ind < dim->v2rhosigma; ++ind)
                {
                    v2rhosigma[dim->v2rhosigma * g + ind] += c * stage_v2rhosigma[dim->v2rhosigma * g + ind];
                }
                for (int ind = 0; ind < dim->v2rholapl; ++ind)
                {
                    v2rholapl[dim->v2rholapl * g + ind] += c * stage_v2rholapl[dim->v2rholapl * g + ind];
                }
                for (int ind = 0; ind < dim->v2rhotau; ++ind)
                {
                    v2rhotau[dim->v2rhotau * g + ind] += c * stage_v2rhotau[dim->v2rhotau * g + ind];
                }
                for (int ind = 0; ind < dim->v2sigma2; ++ind)
                {
                    v2sigma2[dim->v2sigma2 * g + ind] += c * stage_v2sigma2[dim->v2sigma2 * g + ind];
                }
                for (int ind = 0; ind < dim->v2sigmalapl; ++ind)
                {
                    v2sigmalapl[dim->v2sigmalapl * g + ind] += c * stage_v2sigmalapl[dim->v2sigmalapl * g + ind];
                }
                for (int ind = 0; ind < dim->v2sigmatau; ++ind)
                {
                    v2sigmatau[dim->v2sigmatau * g + ind] += c * stage_v2sigmatau[dim->v2sigmatau * g + ind];
                }
                for (int ind = 0; ind < dim->v2lapl2; ++ind)
                {
                    v2lapl2[dim->v2lapl2 * g + ind] += c * stage_v2lapl2[dim->v2lapl2 * g + ind];
                }
                for (int ind = 0; ind < dim->v2lapltau; ++ind)
                {
                    v2lapltau[dim->v2lapltau * g + ind] += c * stage_v2lapltau[dim->v2lapltau * g + ind];
                }
                for (int ind = 0; ind < dim->v2tau2; ++ind)
                {
                    v2tau2[dim->v2tau2 * g + ind] += c * stage_v2tau2[dim->v2tau2 * g + ind];
                }
            }
        }
    }

    gridscreen::screenFxcForMGGA(
        this, np, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2);

    if (alloc)
    {
        std::free(stage_v2rho2);
        std::free(stage_v2rhosigma);
        std::free(stage_v2rholapl);
        std::free(stage_v2rhotau);
        std::free(stage_v2sigma2);
        std::free(stage_v2sigmalapl);
        std::free(stage_v2sigmatau);
        std::free(stage_v2lapl2);
        std::free(stage_v2lapltau);
        std::free(stage_v2tau2);
    }
}

auto
CXCFunctional::compute_kxc_for_mgga(const int64_t np,
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

    // set up staging buffers

    bool alloc = (np > _ldStaging);

    auto ind_cnt = _getIndicesAndCountsOfDerivatives();

    double* stage_v3rho3         = nullptr;
    double* stage_v3rho2sigma    = nullptr;
    double* stage_v3rho2lapl     = nullptr;
    double* stage_v3rho2tau      = nullptr;
    double* stage_v3rhosigma2    = nullptr;
    double* stage_v3rhosigmalapl = nullptr;
    double* stage_v3rhosigmatau  = nullptr;
    double* stage_v3rholapl2     = nullptr;
    double* stage_v3rholapltau   = nullptr;
    double* stage_v3rhotau2      = nullptr;
    double* stage_v3sigma3       = nullptr;
    double* stage_v3sigma2lapl   = nullptr;
    double* stage_v3sigma2tau    = nullptr;
    double* stage_v3sigmalapl2   = nullptr;
    double* stage_v3sigmalapltau = nullptr;
    double* stage_v3sigmatau2    = nullptr;
    double* stage_v3lapl3        = nullptr;
    double* stage_v3lapl2tau     = nullptr;
    double* stage_v3lapltau2     = nullptr;
    double* stage_v3tau3         = nullptr;

    if (alloc)
    {
        stage_v3rho3         = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho3"][1] * np);
        stage_v3rho2sigma    = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho2sigma"][1] * np);
        stage_v3rho2lapl     = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho2lapl"][1] * np);
        stage_v3rho2tau      = (double*)std::malloc(sizeof(double) * ind_cnt["v3rho2tau"][1] * np);
        stage_v3rhosigma2    = (double*)std::malloc(sizeof(double) * ind_cnt["v3rhosigma2"][1] * np);
        stage_v3rhosigmalapl = (double*)std::malloc(sizeof(double) * ind_cnt["v3rhosigmalapl"][1] * np);
        stage_v3rhosigmatau  = (double*)std::malloc(sizeof(double) * ind_cnt["v3rhosigmatau"][1] * np);
        stage_v3rholapl2     = (double*)std::malloc(sizeof(double) * ind_cnt["v3rholapl2"][1] * np);
        stage_v3rholapltau   = (double*)std::malloc(sizeof(double) * ind_cnt["v3rholapltau"][1] * np);
        stage_v3rhotau2      = (double*)std::malloc(sizeof(double) * ind_cnt["v3rhotau2"][1] * np);
        stage_v3sigma3       = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigma3"][1] * np);
        stage_v3sigma2lapl   = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigma2lapl"][1] * np);
        stage_v3sigma2tau    = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigma2tau"][1] * np);
        stage_v3sigmalapl2   = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigmalapl2"][1] * np);
        stage_v3sigmalapltau = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigmalapltau"][1] * np);
        stage_v3sigmatau2    = (double*)std::malloc(sizeof(double) * ind_cnt["v3sigmatau2"][1] * np);
        stage_v3lapl3        = (double*)std::malloc(sizeof(double) * ind_cnt["v3lapl3"][1] * np);
        stage_v3lapl2tau     = (double*)std::malloc(sizeof(double) * ind_cnt["v3lapl2tau"][1] * np);
        stage_v3lapltau2     = (double*)std::malloc(sizeof(double) * ind_cnt["v3lapltau2"][1] * np);
        stage_v3tau3         = (double*)std::malloc(sizeof(double) * ind_cnt["v3tau3"][1] * np);
    }
    else
    {
        stage_v3rho3         = &_stagingBuffer[ind_cnt["v3rho3"][0] * _ldStaging];
        stage_v3rho2sigma    = &_stagingBuffer[ind_cnt["v3rho2sigma"][0] * _ldStaging];
        stage_v3rho2lapl     = &_stagingBuffer[ind_cnt["v3rho2lapl"][0] * _ldStaging];
        stage_v3rho2tau      = &_stagingBuffer[ind_cnt["v3rho2tau"][0] * _ldStaging];
        stage_v3rhosigma2    = &_stagingBuffer[ind_cnt["v3rhosigma2"][0] * _ldStaging];
        stage_v3rhosigmalapl = &_stagingBuffer[ind_cnt["v3rhosigmalapl"][0] * _ldStaging];
        stage_v3rhosigmatau  = &_stagingBuffer[ind_cnt["v3rhosigmatau"][0] * _ldStaging];
        stage_v3rholapl2     = &_stagingBuffer[ind_cnt["v3rholapl2"][0] * _ldStaging];
        stage_v3rholapltau   = &_stagingBuffer[ind_cnt["v3rholapltau"][0] * _ldStaging];
        stage_v3rhotau2      = &_stagingBuffer[ind_cnt["v3rhotau2"][0] * _ldStaging];
        stage_v3sigma3       = &_stagingBuffer[ind_cnt["v3sigma3"][0] * _ldStaging];
        stage_v3sigma2lapl   = &_stagingBuffer[ind_cnt["v3sigma2lapl"][0] * _ldStaging];
        stage_v3sigma2tau    = &_stagingBuffer[ind_cnt["v3sigma2tau"][0] * _ldStaging];
        stage_v3sigmalapl2   = &_stagingBuffer[ind_cnt["v3sigmalapl2"][0] * _ldStaging];
        stage_v3sigmalapltau = &_stagingBuffer[ind_cnt["v3sigmalapltau"][0] * _ldStaging];
        stage_v3sigmatau2    = &_stagingBuffer[ind_cnt["v3sigmatau2"][0] * _ldStaging];
        stage_v3lapl3        = &_stagingBuffer[ind_cnt["v3lapl3"][0] * _ldStaging];
        stage_v3lapl2tau     = &_stagingBuffer[ind_cnt["v3lapl2tau"][0] * _ldStaging];
        stage_v3lapltau2     = &_stagingBuffer[ind_cnt["v3lapltau2"][0] * _ldStaging];
        stage_v3tau3         = &_stagingBuffer[ind_cnt["v3tau3"][0] * _ldStaging];
    }

    // compute derivatives

    auto       mggafunc = getFunctionalPointerToMetaGgaComponent();
    const auto dim      = &(mggafunc->dim);

    for (int64_t g = 0; g < np; ++g)
    {
        for (int ind = 0; ind < dim->v3rho3; ++ind)
        {
            v3rho3[dim->v3rho3 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rho2sigma; ++ind)
        {
            v3rho2sigma[dim->v3rho2sigma * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rho2lapl; ++ind)
        {
            v3rho2lapl[dim->v3rho2lapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rho2tau; ++ind)
        {
            v3rho2tau[dim->v3rho2tau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rhosigma2; ++ind)
        {
            v3rhosigma2[dim->v3rhosigma2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rhosigmalapl; ++ind)
        {
            v3rhosigmalapl[dim->v3rhosigmalapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rhosigmatau; ++ind)
        {
            v3rhosigmatau[dim->v3rhosigmatau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rholapl2; ++ind)
        {
            v3rholapl2[dim->v3rholapl2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rholapltau; ++ind)
        {
            v3rholapltau[dim->v3rholapltau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3rhotau2; ++ind)
        {
            v3rhotau2[dim->v3rhotau2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigma3; ++ind)
        {
            v3sigma3[dim->v3sigma3 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigma2lapl; ++ind)
        {
            v3sigma2lapl[dim->v3sigma2lapl * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigma2tau; ++ind)
        {
            v3sigma2tau[dim->v3sigma2tau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigmalapl2; ++ind)
        {
            v3sigmalapl2[dim->v3sigmalapl2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigmalapltau; ++ind)
        {
            v3sigmalapltau[dim->v3sigmalapltau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3sigmatau2; ++ind)
        {
            v3sigmatau2[dim->v3sigmatau2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3lapl3; ++ind)
        {
            v3lapl3[dim->v3lapl3 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3lapl2tau; ++ind)
        {
            v3lapl2tau[dim->v3lapl2tau * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3lapltau2; ++ind)
        {
            v3lapltau2[dim->v3lapltau2 * g + ind] = 0.0;
        }
        for (int ind = 0; ind < dim->v3tau3; ++ind)
        {
            v3tau3[dim->v3tau3 * g + ind] = 0.0;
        }
    }

    for (const auto& xccomp : _components)
    {
        auto funcptr = xccomp.getFunctionalPointer();

        const auto dim = &(funcptr->dim);

        const auto c = xccomp.getScalingFactor();

        if (xccomp.isLDA())
        {
            xc_lda_kxc(funcptr, np, rho, stage_v3rho3);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_kxc(funcptr, np, rho, sigma, stage_v3rho3, stage_v3rho2sigma, stage_v3rhosigma2, stage_v3sigma3);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
                for (int ind = 0; ind < dim->v3rho2sigma; ++ind)
                {
                    v3rho2sigma[dim->v3rho2sigma * g + ind] += c * stage_v3rho2sigma[dim->v3rho2sigma * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhosigma2; ++ind)
                {
                    v3rhosigma2[dim->v3rhosigma2 * g + ind] += c * stage_v3rhosigma2[dim->v3rhosigma2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigma3; ++ind)
                {
                    v3sigma3[dim->v3sigma3 * g + ind] += c * stage_v3sigma3[dim->v3sigma3 * g + ind];
                }
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

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v3rho3; ++ind)
                {
                    v3rho3[dim->v3rho3 * g + ind] += c * stage_v3rho3[dim->v3rho3 * g + ind];
                }
                for (int ind = 0; ind < dim->v3rho2sigma; ++ind)
                {
                    v3rho2sigma[dim->v3rho2sigma * g + ind] += c * stage_v3rho2sigma[dim->v3rho2sigma * g + ind];
                }
                for (int ind = 0; ind < dim->v3rho2lapl; ++ind)
                {
                    v3rho2lapl[dim->v3rho2lapl * g + ind] += c * stage_v3rho2lapl[dim->v3rho2lapl * g + ind];
                }
                for (int ind = 0; ind < dim->v3rho2tau; ++ind)
                {
                    v3rho2tau[dim->v3rho2tau * g + ind] += c * stage_v3rho2tau[dim->v3rho2tau * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhosigma2; ++ind)
                {
                    v3rhosigma2[dim->v3rhosigma2 * g + ind] += c * stage_v3rhosigma2[dim->v3rhosigma2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhosigmalapl; ++ind)
                {
                    v3rhosigmalapl[dim->v3rhosigmalapl * g + ind] += c * stage_v3rhosigmalapl[dim->v3rhosigmalapl * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhosigmatau; ++ind)
                {
                    v3rhosigmatau[dim->v3rhosigmatau * g + ind] += c * stage_v3rhosigmatau[dim->v3rhosigmatau * g + ind];
                }
                for (int ind = 0; ind < dim->v3rholapl2; ++ind)
                {
                    v3rholapl2[dim->v3rholapl2 * g + ind] += c * stage_v3rholapl2[dim->v3rholapl2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3rholapltau; ++ind)
                {
                    v3rholapltau[dim->v3rholapltau * g + ind] += c * stage_v3rholapltau[dim->v3rholapltau * g + ind];
                }
                for (int ind = 0; ind < dim->v3rhotau2; ++ind)
                {
                    v3rhotau2[dim->v3rhotau2 * g + ind] += c * stage_v3rhotau2[dim->v3rhotau2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigma3; ++ind)
                {
                    v3sigma3[dim->v3sigma3 * g + ind] += c * stage_v3sigma3[dim->v3sigma3 * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigma2lapl; ++ind)
                {
                    v3sigma2lapl[dim->v3sigma2lapl * g + ind] += c * stage_v3sigma2lapl[dim->v3sigma2lapl * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigma2tau; ++ind)
                {
                    v3sigma2tau[dim->v3sigma2tau * g + ind] += c * stage_v3sigma2tau[dim->v3sigma2tau * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigmalapl2; ++ind)
                {
                    v3sigmalapl2[dim->v3sigmalapl2 * g + ind] += c * stage_v3sigmalapl2[dim->v3sigmalapl2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigmalapltau; ++ind)
                {
                    v3sigmalapltau[dim->v3sigmalapltau * g + ind] += c * stage_v3sigmalapltau[dim->v3sigmalapltau * g + ind];
                }
                for (int ind = 0; ind < dim->v3sigmatau2; ++ind)
                {
                    v3sigmatau2[dim->v3sigmatau2 * g + ind] += c * stage_v3sigmatau2[dim->v3sigmatau2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3lapl3; ++ind)
                {
                    v3lapl3[dim->v3lapl3 * g + ind] += c * stage_v3lapl3[dim->v3lapl3 * g + ind];
                }
                for (int ind = 0; ind < dim->v3lapl2tau; ++ind)
                {
                    v3lapl2tau[dim->v3lapl2tau * g + ind] += c * stage_v3lapl2tau[dim->v3lapl2tau * g + ind];
                }
                for (int ind = 0; ind < dim->v3lapltau2; ++ind)
                {
                    v3lapltau2[dim->v3lapltau2 * g + ind] += c * stage_v3lapltau2[dim->v3lapltau2 * g + ind];
                }
                for (int ind = 0; ind < dim->v3tau3; ++ind)
                {
                    v3tau3[dim->v3tau3 * g + ind] += c * stage_v3tau3[dim->v3tau3 * g + ind];
                }
            }
        }
    }

    gridscreen::screenKxcForMGGA(this,
                                 np,
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

    if (alloc)
    {
        std::free(stage_v3rho3);
        std::free(stage_v3rho2sigma);
        std::free(stage_v3rho2lapl);
        std::free(stage_v3rho2tau);
        std::free(stage_v3rhosigma2);
        std::free(stage_v3rhosigmalapl);
        std::free(stage_v3rhosigmatau);
        std::free(stage_v3rholapl2);
        std::free(stage_v3rholapltau);
        std::free(stage_v3rhotau2);
        std::free(stage_v3sigma3);
        std::free(stage_v3sigma2lapl);
        std::free(stage_v3sigma2tau);
        std::free(stage_v3sigmalapl2);
        std::free(stage_v3sigmalapltau);
        std::free(stage_v3sigmatau2);
        std::free(stage_v3lapl3);
        std::free(stage_v3lapl2tau);
        std::free(stage_v3lapltau2);
        std::free(stage_v3tau3);
    }
}

auto
CXCFunctional::compute_lxc_for_mgga(const int64_t np,
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
        stage_v4rho4            = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho4"][1] * np);
        stage_v4rho3sigma       = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho3sigma"][1] * np);
        stage_v4rho3lapl        = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho3lapl"][1] * np);
        stage_v4rho3tau         = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho3tau"][1] * np);
        stage_v4rho2sigma2      = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2sigma2"][1] * np);
        stage_v4rho2sigmalapl   = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2sigmalapl"][1] * np);
        stage_v4rho2sigmatau    = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2sigmatau"][1] * np);
        stage_v4rho2lapl2       = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2lapl2"][1] * np);
        stage_v4rho2lapltau     = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2lapltau"][1] * np);
        stage_v4rho2tau2        = (double*)std::malloc(sizeof(double) * ind_cnt["v4rho2tau2"][1] * np);
        stage_v4rhosigma3       = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigma3"][1] * np);
        stage_v4rhosigma2lapl   = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigma2lapl"][1] * np);
        stage_v4rhosigma2tau    = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigma2tau"][1] * np);
        stage_v4rhosigmalapl2   = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigmalapl2"][1] * np);
        stage_v4rhosigmalapltau = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigmalapltau"][1] * np);
        stage_v4rhosigmatau2    = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhosigmatau2"][1] * np);
        stage_v4rholapl3        = (double*)std::malloc(sizeof(double) * ind_cnt["v4rholapl3"][1] * np);
        stage_v4rholapl2tau     = (double*)std::malloc(sizeof(double) * ind_cnt["v4rholapl2tau"][1] * np);
        stage_v4rholapltau2     = (double*)std::malloc(sizeof(double) * ind_cnt["v4rholapltau2"][1] * np);
        stage_v4rhotau3         = (double*)std::malloc(sizeof(double) * ind_cnt["v4rhotau3"][1] * np);
        stage_v4sigma4          = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma4"][1] * np);
        stage_v4sigma3lapl      = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma3lapl"][1] * np);
        stage_v4sigma3tau       = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma3tau"][1] * np);
        stage_v4sigma2lapl2     = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma2lapl2"][1] * np);
        stage_v4sigma2lapltau   = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma2lapltau"][1] * np);
        stage_v4sigma2tau2      = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigma2tau2"][1] * np);
        stage_v4sigmalapl3      = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigmalapl3"][1] * np);
        stage_v4sigmalapl2tau   = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigmalapl2tau"][1] * np);
        stage_v4sigmalapltau2   = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigmalapltau2"][1] * np);
        stage_v4sigmatau3       = (double*)std::malloc(sizeof(double) * ind_cnt["v4sigmatau3"][1] * np);
        stage_v4lapl4           = (double*)std::malloc(sizeof(double) * ind_cnt["v4lapl4"][1] * np);
        stage_v4lapl3tau        = (double*)std::malloc(sizeof(double) * ind_cnt["v4lapl3tau"][1] * np);
        stage_v4lapl2tau2       = (double*)std::malloc(sizeof(double) * ind_cnt["v4lapl2tau2"][1] * np);
        stage_v4lapltau3        = (double*)std::malloc(sizeof(double) * ind_cnt["v4lapltau3"][1] * np);
        stage_v4tau4            = (double*)std::malloc(sizeof(double) * ind_cnt["v4tau4"][1] * np);
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

    for (int64_t g = 0; g < np; ++g)
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
            xc_lda_lxc(funcptr, np, rho, stage_v4rho4);

            for (int64_t g = 0; g < np; ++g)
            {
                for (int ind = 0; ind < dim->v4rho4; ++ind)
                {
                    v4rho4[dim->v4rho4 * g + ind] += c * stage_v4rho4[dim->v4rho4 * g + ind];
                }
            }
        }
        else if (xccomp.isGGA())
        {
            xc_gga_lxc(funcptr, np, rho, sigma, stage_v4rho4, stage_v4rho3sigma, stage_v4rho2sigma2, stage_v4rhosigma3, stage_v4sigma4);

            for (int64_t g = 0; g < np; ++g)
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
                        np,
                        rho,
                        sigma,
                        lapl,
                        tau,
                        stage_v4rho4,
                        stage_v4rho3sigma,
                        stage_v4rho3lapl,
                        stage_v4rho3tau,
                        stage_v4rho2sigma2,
                        stage_v4rho2sigmalapl,
                        stage_v4rho2sigmatau,
                        stage_v4rho2lapl2,
                        stage_v4rho2lapltau,
                        stage_v4rho2tau2,
                        stage_v4rhosigma3,
                        stage_v4rhosigma2lapl,
                        stage_v4rhosigma2tau,
                        stage_v4rhosigmalapl2,
                        stage_v4rhosigmalapltau,
                        stage_v4rhosigmatau2,
                        stage_v4rholapl3,
                        stage_v4rholapl2tau,
                        stage_v4rholapltau2,
                        stage_v4rhotau3,
                        stage_v4sigma4,
                        stage_v4sigma3lapl,
                        stage_v4sigma3tau,
                        stage_v4sigma2lapl2,
                        stage_v4sigma2lapltau,
                        stage_v4sigma2tau2,
                        stage_v4sigmalapl3,
                        stage_v4sigmalapl2tau,
                        stage_v4sigmalapltau2,
                        stage_v4sigmatau3,
                        stage_v4lapl4,
                        stage_v4lapl3tau,
                        stage_v4lapl2tau2,
                        stage_v4lapltau3,
                        stage_v4tau4);

            for (int64_t g = 0; g < np; ++g)
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

    if (alloc)
    {
        std::free(stage_v4rho4);
        std::free(stage_v4rho3sigma);
        std::free(stage_v4rho3lapl);
        std::free(stage_v4rho3tau);
        std::free(stage_v4rho2sigma2);
        std::free(stage_v4rho2sigmalapl);
        std::free(stage_v4rho2sigmatau);
        std::free(stage_v4rho2lapl2);
        std::free(stage_v4rho2lapltau);
        std::free(stage_v4rho2tau2);
        std::free(stage_v4rhosigma3);
        std::free(stage_v4rhosigma2lapl);
        std::free(stage_v4rhosigma2tau);
        std::free(stage_v4rhosigmalapl2);
        std::free(stage_v4rhosigmalapltau);
        std::free(stage_v4rhosigmatau2);
        std::free(stage_v4rholapl3);
        std::free(stage_v4rholapl2tau);
        std::free(stage_v4rholapltau2);
        std::free(stage_v4rhotau3);
        std::free(stage_v4sigma4);
        std::free(stage_v4sigma3lapl);
        std::free(stage_v4sigma3tau);
        std::free(stage_v4sigma2lapl2);
        std::free(stage_v4sigma2lapltau);
        std::free(stage_v4sigma2tau2);
        std::free(stage_v4sigmalapl3);
        std::free(stage_v4sigmalapl2tau);
        std::free(stage_v4sigmalapltau2);
        std::free(stage_v4sigmatau3);
        std::free(stage_v4lapl4);
        std::free(stage_v4lapl3tau);
        std::free(stage_v4lapl2tau2);
        std::free(stage_v4lapltau3);
        std::free(stage_v4tau4);
    }
}

xc_func_type*
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

xc_func_type*
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

xc_func_type*
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

int64_t
CXCFunctional::getDimensionOfDerivatives() const
{
    auto indices_and_counts = _getIndicesAndCountsOfDerivatives();

    if (_familyOfFunctional == xcfun::lda)
    {
        auto final_index = indices_and_counts["v4rho4"][0];
        auto final_count = indices_and_counts["v4rho4"][1];

        return final_index + final_count;
    }
    else if (_familyOfFunctional == xcfun::gga)
    {
        auto final_index = indices_and_counts["v4sigma4"][0];
        auto final_count = indices_and_counts["v4sigma4"][1];

        return final_index + final_count;
    }
    else if (_familyOfFunctional == xcfun::mgga)
    {
        auto final_index = indices_and_counts["v4tau4"][0];
        auto final_count = indices_and_counts["v4tau4"][1];

        return final_index + final_count;
    }

    return 0;
}

auto
CXCFunctional::setRangeSeparatedParameterOmega(const double omega) -> void
{
    errors::assertMsgCritical(isRangeSeparated(), "XCFunctional.setRangeSeparatedParameterOmega: Only applicable to range-separated functional");

    errors::assertMsgCritical(_components.size() == 1, "XCFunctional.setRangeSeparatedParameterOmega: Only applicable to single-component functional");

    std::string param_name("_omega");

    auto funcptr = _components[0].getFunctionalPointer();

    xc_func_set_ext_params_name(funcptr, param_name.c_str(), omega);

    xc_hyb_cam_coef(funcptr, &_rangeSeparationParameterOmega, &_rangeSeparationParameterAlpha, &_rangeSeparationParameterBeta);
}
