//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef XCPairDensityFunctional_hpp
#define XCPairDensityFunctional_hpp

#include <xc.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

/**
 * Class CXCPairDensityFunctional implements pair density functionals.
 */
class CXCPairDensityFunctional
{
   private:
    /** Name of functional. */
    std::string _nameOfFunctional{std::string("Undefined")};

    /** Family of functional. */
    std::string _familyOfFunctional{std::string("PLDA")};

    /** Leading dimension for initial allocation of staging buffer. */
    int _ldStaging{1024};

    /** Buffer to stage output results. */
    double* _stagingBuffer{nullptr};

    /** The functional components and their coefficients. */
    std::vector<std::tuple<std::string, double>> _components;

    /** Allocates the staging buffer. */
    void _allocateStagingBuffer();

    /** Frees the staging buffer. */
    void _freeStagingBuffer();

    /** Checks if a component is PLDA.
     *
     * @param[in] compName name of the functional component.
     * @return whether the component is PLDA.
     */
    bool _isComponentPLDA(const std::string& compName) const;

    /** Checks if a component is PGGA.
     *
     * @param[in] compName name of the functional component.
     * @return whether the component is PGGA.
     */
    bool _isComponentPGGA(const std::string& compName) const;

    /** Computes values and first derivative of pair-LDA exchange-correlation
     * functional component on grid.
     *
     * @param[in] compName name of the functional component.
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * @param[in] rs_omega the range-separation parameter.
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     */
    void _plda_exc_vxc(const std::string& compName, const int np, const double* rho, double* exc, double* vrho, double rs_omega) const;

    /** Computes values and first derivative of pair-GGA exchange-correlation
     * functional component on grid.
     *
     * @param[in] compName name of the functional component.
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     * @param[in,out] vsigma values of the first derivative of the
     * @param[in] rs_omega the range-separation parameter.
     * exchange-correlation kernel wrt contracted gradients. Size: 3*np, order: [(0), (1), (2)].
     */
    void _pgga_exc_vxc(const std::string& compName,
                       const int      np,
                       const double*      rho,
                       const double*      sigma,
                       double*            exc,
                       double*            vrho,
                       double*            vsigma,
                       double             rs_omega) const;

   public:
    /** Creates an exchange-correlation functional object.
     *
     * @param[in] name of functional.
     * @param[in] labels list of labels of functional components.
     * @param[in] coeffs list of coefficients for functional components.
     */
    CXCPairDensityFunctional(const std::string&, const std::vector<std::string>& labels, const std::vector<double>& coeffs);

    /**
     Creates an XC functional object by copying other XC functional object.

     @param source the XC functional object.
     */
    CXCPairDensityFunctional(const CXCPairDensityFunctional& source);

    /**
     Creates an XC functional object by moving other XC functional object.

     @param source the XC functional object.
     */
    CXCPairDensityFunctional(CXCPairDensityFunctional&& source) noexcept;

    /**
     * Destroys an exchange-correlation functional object.
     */
    ~CXCPairDensityFunctional();

    /**
     Assigns an XC functional object by copying other XC functional object.

     @param source the XC functional object.
     */
    CXCPairDensityFunctional& operator=(const CXCPairDensityFunctional& source);

    /**
     Assigns an XC functional object by moving other XC functional object.

     @param source the XC functional object.
     */
    CXCPairDensityFunctional& operator=(CXCPairDensityFunctional&& source) noexcept;

    /**
     Compares XC functional object with other XC functional object.

     @param other the XC functional object.
     @return true if XC functional objects are equal, false otherwise.
     */
    bool operator==(const CXCPairDensityFunctional& other) const;

    /**
     Compares XC functional object with other XC functional object.

     @param other the XC functional object.
     @return true if XC functional objects are not equal, false otherwise.
     */
    bool operator!=(const CXCPairDensityFunctional& other) const;

    /**
     Gets functional name.
     */
    std::string getFunctionalLabel() const;

    /**
     Gets XC functional type.
     */
    std::string getFunctionalType() const;

    /** Computes values and first derivative of pair-LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * @param[in] rs_omega the range-separation parameter.
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     */
    auto compute_exc_vxc_for_plda(int np, const double* rho, double* exc, double* vrho, double rs_omega) const -> void;

    /** Computes values and first derivative of pair-GGA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in] sigma values of the contracted gradient of density at grid points. Order: [(0, 0), (0, 1), (1, 1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     * @param[in,out] vsigma values of the first derivative of the
     * @param[in] rs_omega the range-separation parameter.
     * exchange-correlation kernel wrt contracted gradients. Size: 3*np, order: [(0), (1), (2)].
     */
    auto compute_exc_vxc_for_pgga(int np, const double* rho, const double* sigma, double* exc, double* vrho, double* vsigma, double rs_omega) const -> void;
};

#endif /* XCPairDensityFunctional_hpp */
