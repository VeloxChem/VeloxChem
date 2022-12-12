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

#ifndef XCPairDensityFunctional_hpp
#define XCPairDensityFunctional_hpp

#include <xc.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

/**
 * Class CXCPairDensityFunctional is a wrapper to the C functions and structs provided by LibXC.
 *
 * @author X. Li
 */
class CXCPairDensityFunctional
{
   private:
    /** Name of functional. */
    std::string _nameOfFunctional{std::string("Undefined")};

    /** Family of functional. */
    std::string _familyOfFunctional{std::string("PLDA")};

    /** Leading dimension for initial allocation of staging buffer. */
    int32_t _ldStaging{1024};

    /** Buffer to stage output results from LibXC invocations. */
    double* _stagingBuffer{nullptr};

    /** The functional components and their coefficients. */
    std::vector<std::tuple<std::string, double>> _components;

    /** Allocates the staging buffer. */
    void _allocateStagingBuffer();

    /** Frees the staging buffer. */
    void _freeStagingBuffer();

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

    /** Computes values and first derivative of LDA exchange-correlation functional on grid.
     *
     * @param[in] np number of grid points.
     * @param[in] rho values of the density at grid points. Order: [(0), (1)].
     * @param[in,out] exc values of the exchange-correlation kernel. Size: np.
     * @param[in,out] vrho values of the first derivative of the
     * exchange-correlation kernel wrt density. Size: 2*np, order: [(0), (1)].
     */
    auto compute_exc_vxc_for_plda(int32_t np, const double* rho, double* exc, double* vrho) const -> void;
};

#endif /* XCPairDensityFunctional_hpp */
