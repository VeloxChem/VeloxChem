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

#ifndef XCNewFunctional_hpp
#define XCNewFunctional_hpp

#include <xc.h>

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "Buffer.hpp"
#include "DensityGrid.hpp"
#include "XCComponent.hpp"
#include "XCFuncType.hpp"

/**
 * Class CXCNewFunctional is a wrapper to the C functions and structs provided by LibXC.
 *
 * @author R. Di Remigio Eikås, Z. Rinkevicius, X. Li
 */
class CXCNewFunctional
{
   private:
    /** Name of functional. */
    std::string _nameOfFunctional{std::string("Undefined")};

    /** The fraction of exact Hartree-Fock exchange in functional. */
    double _fractionOfExactExchange{0.0};

    /** The range-separation parameters in functional. */
    std::array<double, 3> _rangeSeparationParameters{0.0, 0.0, 0.0};

    /** Highest order of available derivatives. */
    int32_t _maxDerivOrder{0};

    /** Family of functional. */
    std::string _familyOfFunctional{std::string("LDA")};

    /** Leading dimension for initial allocation of staging buffer. */
    int32_t _ldStaging{1024};

    /** Buffer to stage output results from LibXC invocations. */
    double* _stagingBuffer{nullptr};

    /** The functional components and their coefficients. */
    std::vector<CXCComponent> _components;

    /** Allocates the staging buffer. */
    void _allocateStagingBuffer();

    /** Frees the staging buffer. */
    void _freeStagingBuffer();

   public:
    /** Creates an exchange-correlation functional object.
     *
     * @param[in] nameOfFunctional name of functional.
     * @param[in] labels list of labels of functional components.
     * @param[in] coeffs list of coefficients for functional components.
     * @param[in] fractionOfExactExchange fraction of exact exchange.
     */
    CXCNewFunctional(const std::string&              nameOfFunctional,
                     const std::vector<std::string>& labels,
                     const std::vector<double>&      coeffs,
                     const double                    fractionOfExactExchange = 0.0);

    /**
     Creates an XC functional object by copying other XC functional object.

     @param source the XC functional object.
     */
    CXCNewFunctional(const CXCNewFunctional& source);

    /**
     Creates an XC functional object by moving other XC functional object.

     @param source the XC functional object.
     */
    CXCNewFunctional(CXCNewFunctional&& source) noexcept;

    /**
     * Destroys an exchange-correlation functional object.
     */
    ~CXCNewFunctional();

    /**
     Assigns an XC functional object by copying other XC functional object.

     @param source the XC functional object.
     */
    CXCNewFunctional& operator=(const CXCNewFunctional& source);

    /**
     Assigns an XC functional object by moving other XC functional object.

     @param source the XC functional object.
     */
    CXCNewFunctional& operator=(CXCNewFunctional&& source) noexcept;

    /**
     Compares XC functional object with other XC functional object.

     @param other the XC functional object.
     @return true if XC functional objects are equal, false otherwise.
     */
    bool operator==(const CXCNewFunctional& other) const;

    /**
     Compares XC functional object with other XC functional object.

     @param other the XC functional object.
     @return true if XC functional objects are not equal, false otherwise.
     */
    bool operator!=(const CXCNewFunctional& other) const;

    /**
     Gets XC functional name.
     */
    std::string getFunctionalLabel() const;

    /**
     Gets XC functional type.
     */
    xcfun getFunctionalType() const;

    /**
     Determines whether the XC functional is undefined.
     */
    bool isUndefined() const;

    /**
     Determines whether the XC functional is hybrid.
     */
    bool isHybrid() const;

    /**
     Gets the fraction of exact exchange.
     */
    double getFractionOfExactExchange() const;

    /** Computes values and first derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     */
    auto compute_exc_vxc_for_lda(int32_t np, const double* rho, double* exc, double* vrho) const -> void;

    /** Computes first derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     */
    auto compute_vxc_for_lda(int32_t np, const double* rho, double* vrho) const -> void;

    /** Computes second derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] v2rho2 values of the second derivative of the
     * exchange-correlation kernel wrt density. Size: 3*np, order:
     * [(0, 0), (0, 1), (1, 1)].
     */
    auto compute_fxc_for_lda(int32_t np, const double* rho, double* v2rho2) const -> void;

    /** Computes third derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] v3rho3 values of the third derivative of the
     * exchange-correlation kernel wrt density. Size: 4*np, order:
     * [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)].
     */
    auto compute_kxc_for_lda(int32_t np, const double* rho, double* v3rho3) const -> void;

    /** Computes fourth derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] v4rho4 values of the fourth derivative of the
     * exchange-correlation kernel wrt density. Size: 5*np, order:
     * [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 1, 1), (1, 1, 1, 1)].
     */
    auto compute_lxc_for_lda(int32_t np, const double* rho, double* v4rho4) const -> void;

    /** Computes values and first derivative of GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     * @param[in,out] vsigma values of the first derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 3*np, order: [(0), (1), (2)].
     */
    auto compute_exc_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma) const -> void;

    /** Computes first derivative of GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     * @param[in,out] vsigma values of the first derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 3*np, order: [(0), (1), (2)].
     */
    auto compute_vxc_for_gga(int32_t np, const double* rho, const double* sigma, double* vrho, double* vsigma) const -> void;

    /** Computes second derivative of GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] v2rho2 values of the second derivative of the
     * exchange-correlation kernel wrt density. Size: 3*np, order:
     * [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] v2rhosigma values of the second derivative of the
     * exchange-correlation kernel wrt density and contracted gradients. Size: 6*np, order:
     * [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)].
     * @param[in,out] v2sigma2 values of the second derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 6*np, order:
     * [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)].
     */
    auto compute_fxc_for_gga(int32_t np, const double* rho, const double* sigma, double* v2rho2, double* v2rhosigma, double* v2sigma2) const -> void;

    /** Computes third derivative of GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] v3rho3 values of the third derivative of the
     * exchange-correlation kernel wrt density. Size: 4*np, order:
     * [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)].
     * @param[in,out] v3rho2sigma values of the third derivative of the
     * exchange-correlation kernel wrt density and contracted gradients. Size: 9*np, order:
     * [(0, 0, 0), (0, 0, 1), (0, 0, 2),
     *  (0, 1, 0), (0, 1, 1), (0, 1, 2),
     *  (1, 1, 0), (1, 1, 1), (1, 1, 2)]
     * @param[in,out] v3rhosigma2 values of the third derivative of the
     * exchange-correlation kernel wrt density and contracted gradients. Size: 12*np, order:
     * [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2),
     *  (1, 0, 0), (1, 0, 1), (1, 0, 2), (1, 1, 1), (1, 1, 2), (1, 2, 2)]
     * @param[in,out] v3sigma3 values of the third derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 10*np, order:
     * [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2), (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
     */
    auto compute_kxc_for_gga(int32_t       np,
                             const double* rho,
                             const double* sigma,
                             double*       v3rho3,
                             double*       v3rho2sigma,
                             double*       v3rhosigma2,
                             double*       v3sigma3) const -> void;

    /** Computes fourth derivative of GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] v4rho4 values of the fourth derivative of the
     * exchange-correlation kernel wrt density. Size: 5*np, order:
     * [(0,0,0,0), (0,0,0,1), (0,0,1,1), (0,1,1,1), (1,1,1,1)].
     * @param[in,out] v4rho3sigma values of the fourth derivative of the
     * exchange-correlation kernel wrt density and contracted gradients. Size: 12*np, order:
     * [(0,0,0,0), (0,0,0,1), (0,0,0,2),
     *  (0,0,1,0), (0,0,1,1), (0,0,1,2),
     *  (0,1,1,0), (0,1,1,1), (0,1,1,2),
     *  (1,1,1,0), (1,1,1,1), (1,1,1,2)]
     * @param[in,out] v4rho2sigma2 values of the fourth derivative of the
     * exchange-correlation kernel wrt density and contracted gradients. Size: 18*np, order:
     * [(0,0,0,0), (0,0,0,1), (0,0,0,2), (0,0,1,1), (0,0,1,2), (0,0,2,2),
     *  (0,1,0,0), (0,1,0,1), (0,1,0,2), (0,1,1,1), (0,1,1,2), (0,1,2,2),
     *  (1,1,0,0), (1,1,0,1), (1,1,0,2), (1,1,1,1), (1,1,1,2), (1,1,2,2)]
     * @param[in,out] v4rhosigma3 values of the fourth derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 20*np, order:
     * [(0,0,0,0),(0,0,0,1),(0,0,0,2),(0,0,1,1),(0,0,1,2),(0,0,2,2),(0,1,1,1),(0,1,1,2),(0,1,2,2),(0,2,2,2),
     *  (1,0,0,0),(1,0,0,1),(1,0,0,2),(1,0,1,1),(1,0,1,2),(1,0,2,2),(1,1,1,1),(1,1,1,2),(1,1,2,2),(1,2,2,2)]
     * @param[in,out] v4sigma4 values of the fourth derivative of the
     * exchange-correlation kernel wrt contracted gradients. Size: 15*np, order:
     * [(0,0,0,0),(0,0,0,1),(0,0,0,2),(0,0,1,1),(0,0,1,2),(0,0,2,2),(0,1,1,1),(0,1,1,2),
     *  (0,1,2,2),(0,2,2,2),(1,1,1,1),(1,1,1,2),(1,1,2,2),(1,2,2,2),(2,2,2,2)]
     */
    auto compute_lxc_for_gga(int32_t       np,
                             const double* rho,
                             const double* sigma,
                             double*       v4rho4,
                             double*       v4rho3sigma,
                             double*       v4rho2sigma2,
                             double*       v4rhosigma3,
                             double*       v4sigma4) const -> void;

    /** Computes values and first derivative of metaGGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in] lapl values of the density Laplacian at grid points. Order: [(0), (1)].
     * @param[in] tau values of the kinetic energy density at grid points. Order: [(0), (1)].
     * @param[in,out] exc values of the exchange-correlation kernel.
     * @param[in,out] vrho values of the first derivative
     * @param[in,out] vsigma values of the first derivative
     * @param[in,out] vlapl values of the first derivative
     * @param[in,out] vtau alues of the first derivative
     */
    auto compute_exc_vxc_for_mgga(int32_t       np,
                                  const double* rho,
                                  const double* sigma,
                                  const double* lapl,
                                  const double* tau,
                                  double*       exc,
                                  double*       vrho,
                                  double*       vsigma,
                                  double*       vlapl,
                                  double*       vtau) const -> void;

    /** Computes first derivative of metaGGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in] lapl values of the density Laplacian at grid points. Order: [(0), (1)].
     * @param[in] tau values of the kinetic energy density at grid points. Order: [(0), (1)].
     * @param[in,out] vrho values of the first derivative
     * @param[in,out] vsigma values of the first derivative
     * @param[in,out] vlapl values of the first derivative
     * @param[in,out] vtau alues of the first derivative
     */
    auto compute_vxc_for_mgga(int32_t       np,
                              const double* rho,
                              const double* sigma,
                              const double* lapl,
                              const double* tau,
                              double*       vrho,
                              double*       vsigma,
                              double*       vlapl,
                              double*       vtau) const -> void;

    /** Computes second derivative of metaGGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in] lapl values of the density Laplacian at grid points. Order: [(0), (1)].
     * @param[in] tau values of the kinetic energy density at grid points. Order: [(0), (1)].
     * @param[in,out] v2rho2 values of the second derivative
     * @param[in,out] v2rhosigma values of the second derivative
     * @param[in,out] v2rholapl values of the second derivative
     * @param[in,out] v2rhotau values of the second derivative
     * @param[in,out] v2sigma2 values of the second derivative
     * @param[in,out] v2sigmalapl values of the second derivative
     * @param[in,out] v2sigmatau values of the second derivative
     * @param[in,out] v2lapl2 values of the second derivative
     * @param[in,out] v2lapltau values of the second derivative
     * @param[in,out] v2tau2 values of the second derivative
     */
    auto compute_fxc_for_mgga(int32_t       np,
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
                              double*       v2tau2) const -> void;

    /** Computes third derivative of metaGGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in] lapl values of the density Laplacian at grid points. Order: [(0), (1)].
     * @param[in] tau values of the kinetic energy density at grid points. Order: [(0), (1)].
     * @param[in,out] v3rho3 values of the third derivative
     * @param[in,out] v3rho2sigma values of the third derivative
     * @param[in,out] v3rho2lapl values of the third derivative
     * @param[in,out] v3rho2tau values of the third derivative
     * @param[in,out] v3rhosigma2 values of the third derivative
     * @param[in,out] v3rhosigmalapl values of the third derivative
     * @param[in,out] v3rhosigmatau values of the third derivative
     * @param[in,out] v3rholapl2 values of the third derivative
     * @param[in,out] v3rholapltau values of the third derivative
     * @param[in,out] v3rhotau2 values of the third derivative
     * @param[in,out] v3sigma3 values of the third derivative
     * @param[in,out] v3sigma2lapl values of the third derivative
     * @param[in,out] v3sigma2tau values of the third derivative
     * @param[in,out] v3sigmalapl2 values of the third derivative
     * @param[in,out] v3sigmalapltau values of the third derivative
     * @param[in,out] v3sigmatau2 values of the third derivative
     * @param[in,out] v3lapl3 values of the third derivative
     * @param[in,out] v3lapl2tau values of the third derivative
     * @param[in,out] v3lapltau2 values of the third derivative
     * @param[in,out] v3tau3 values of the third derivative
     */
    auto compute_kxc_for_mgga(int32_t       np,
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
                              double*       v3tau3) const -> void;

    /** Computes fourth derivative of metaGGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in] lapl values of the density Laplacian at grid points. Order: [(0), (1)].
     * @param[in] tau values of the kinetic energy density at grid points. Order: [(0), (1)].
     * @param[in,out] v4rho4 values of the fourth derivative
     * @param[in,out] v4rho3sigma values of the fourth derivative
     * @param[in,out] v4rho3lapl values of the fourth derivative
     * @param[in,out] v4rho3tau values of the fourth derivative
     * @param[in,out] v4rho2sigma2 values of the fourth derivative
     * @param[in,out] v4rho2sigmalapl values of the fourth derivative
     * @param[in,out] v4rho2sigmatau values of the fourth derivative
     * @param[in,out] v4rho2lapl2 values of the fourth derivative
     * @param[in,out] v4rho2lapltau values of the fourth derivative
     * @param[in,out] v4rho2tau2 values of the fourth derivative
     * @param[in,out] v4rhosigma3 values of the fourth derivative
     * @param[in,out] v4rhosigma2lapl values of the fourth derivative
     * @param[in,out] v4rhosigma2tau values of the fourth derivative
     * @param[in,out] v4rhosigmalapl2 values of the fourth derivative
     * @param[in,out] v4rhosigmalapltau values of the fourth derivative
     * @param[in,out] v4rhosigmatau2 values of the fourth derivative
     * @param[in,out] v4rholapl3 values of the fourth derivative
     * @param[in,out] v4rholapl2tau values of the fourth derivative
     * @param[in,out] v4rholapltau2 values of the fourth derivative
     * @param[in,out] v4rhotau3 values of the fourth derivative
     * @param[in,out] v4sigma4 values of the fourth derivative
     * @param[in,out] v4sigma3lapl values of the fourth derivative
     * @param[in,out] v4sigma3tau values of the fourth derivative
     * @param[in,out] v4sigma2lapl2 values of the fourth derivative
     * @param[in,out] v4sigma2lapltau values of the fourth derivative
     * @param[in,out] v4sigma2tau2 values of the fourth derivative
     * @param[in,out] v4sigmalapl3 values of the fourth derivative
     * @param[in,out] v4sigmalapl2tau values of the fourth derivative
     * @param[in,out] v4sigmalapltau2 values of the fourth derivative
     * @param[in,out] v4sigmatau3 values of the fourth derivative
     * @param[in,out] v4lapl4 values of the fourth derivative
     * @param[in,out] v4lapl3tau values of the fourth derivative
     * @param[in,out] v4lapl2tau2 values of the fourth derivative
     * @param[in,out] v4lapltau3 values of the fourth derivative
     * @param[in,out] v4tau4 values of the fourth derivative
     */
    auto compute_lxc_for_mgga(int32_t       np,
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
                              double*       v4tau4) const -> void;
};

#endif /* XCNewFunctional_hpp */
