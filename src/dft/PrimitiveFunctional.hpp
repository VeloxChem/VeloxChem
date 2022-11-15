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

#ifndef PrimitiveFunctional_hpp
#define PrimitiveFunctional_hpp

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "DensityGrid.hpp"
#include "XCCubicHessianGrid.hpp"
#include "XCFuncType.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"

using def_vxc_func_typ = void(CXCGradientGrid&, const double factor, const CDensityGrid&);

using def_vxc2_func_typ = void(CXCHessianGrid&, const double factor, const CDensityGrid&);

using def_vxc3_func_typ = void(CXCCubicHessianGrid&, const double factor, const CDensityGrid&);

// using def_vxc2_func_typ = void(CXCGradientGrid&, const CDensityGrid&);

/**
 Class CPrimitiveFunctional stores information about primitive exchange-correlation functional
 and provides methods to perform actions with stored exchange-correlation functional.

 @author Z. Rinkevicius
 */
class CPrimitiveFunctional
{
    /**
     The label of primitive exchange-correlation functional.
     */
    std::string _label;

    /**
     The type of primitive exchange-correlation functional.
     */
    xcfun _xcFuncType;

    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc_func_typ> _abFirstOrderFunction;

    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::lima type.
     */
    std::function<def_vxc_func_typ> _aFirstOrderFunction;

    /**
     The functions for computing first derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::limb type.
     */
    std::function<def_vxc_func_typ> _bFirstOrderFunction;

    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc2_func_typ> _abSecondOrderFunction;

    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::lima type.
     */
    std::function<def_vxc2_func_typ> _aSecondOrderFunction;

    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::limb type.
     */
    std::function<def_vxc2_func_typ> _bSecondOrderFunction;

    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc3_func_typ> _abThirdOrderFunction;

    /**
    The functions for computing second derrivatives of primitive exchange-correlation functional for density
    grid of dengrid::ab type.
    */
    std::function<def_vxc3_func_typ> _aThirdOrderFunction;

    /**
     The functions for computing second derrivatives of primitive exchange-correlation functional for density
     grid of dengrid::ab type.
     */
    std::function<def_vxc3_func_typ> _bThirdOrderFunction;

   public:
    /**
     Creates an empty primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional();

    /**
     Creates a primitive exchange-correlation functional object.

     @param label the label of primitive exchange-correlation functional.
     @param xcFuncType the type of primitive exchange-correlation functional.
     @param abFirstOrderFunction the first-order derivative function (dengrid::ab).
     @param aFirstOrderFunction the first-order derivative function (dengrid::lima).
     @param bFirstOrderFunction the first-order derivative function (dengrid::limb).
     @param abSecondOrderFunction the second-order derivative function (dengrid::ab).
     @param aSecondOrderFunction the second-order derivative function (dengrid::lima).
     @param bSecondOrderFunction the second-order derivative function (dengrid::limb).
     */
    CPrimitiveFunctional(const std::string&                      label,
                         const xcfun                             xcFuncType,
                         const std::function<def_vxc_func_typ>&  abFirstOrderFunction,
                         const std::function<def_vxc_func_typ>&  aFirstOrderFunction,
                         const std::function<def_vxc_func_typ>&  bFirstOrderFunction,
                         const std::function<def_vxc2_func_typ>& abSecondOrderFunction,
                         const std::function<def_vxc2_func_typ>& aSecondOrderFunction,
                         const std::function<def_vxc2_func_typ>& bSecondOrderFunction,
                         const std::function<def_vxc3_func_typ>& abThirdOrderFunction,
                         const std::function<def_vxc3_func_typ>& aThirdOrderFunction,
                         const std::function<def_vxc3_func_typ>& bThirdOrderFunction);

    /**
     Creates a primitive exchange-correlation functional object by copying other primitive exchange-correlation functional
     object.

     @param source the recursion term object.
     */
    CPrimitiveFunctional(const CPrimitiveFunctional& source);

    /**
     Creates a primitive exchange-correlation functional object by moving other primitive exchange-correlation functional
     object.

     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional(CPrimitiveFunctional&& source) noexcept;

    /**
     Destroys a primitive exchange-correlation functional object.
     */
    ~CPrimitiveFunctional();

    /**
     Assigns a primitive exchange-correlation functional object by copying other primitive exchange-correlation functional
     object.

     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional& operator=(const CPrimitiveFunctional& source);

    /**
     Assigns a primitive exchange-correlation functional object by moving other primitive exchange-correlation functional
     object.

     @param source the primitive exchange-correlation functional object.
     */
    CPrimitiveFunctional& operator=(CPrimitiveFunctional&& source) noexcept;

    /**
     Compares primitive exchange-correlation functional object with other primitive exchange-correlation functional object.

     @param other the primitive exchange-correlation functional object.
     @return true if primitive exchange-correlation functional objects are equal, false otherwise.
     */
    bool operator==(const CPrimitiveFunctional& other) const;

    /**
     Compares primitive exchange-correlation functional object with other primitive exchange-correlation functional object.

     @param other the primitive exchange-correlation functional object.
     @return true if primitive exchange-correlation functional objects are not equal, false otherwise.
     */
    bool operator!=(const CPrimitiveFunctional& other) const;

    /**
     Gets label of primitive exchange-correlation functional.

     @return the label.
     */
    std::string getLabel() const;

    /**
     Gets type of primitive exchange-correlation functional.

     @return the type of exchange-correlation functional.
     */
    xcfun getFunctionalType() const;

    /**
     Computes first derivative of exchange-correlation functional for given density grid.

     @param xcGradientGrid the exchange-correlation gradient grid object.
     @param factor the scaling factor of functional contribution.
     @param densityGrid the density grid object.
     */
    void compute(CXCGradientGrid& xcGradientGrid, const double factor, const CDensityGrid& densityGrid) const;

    /**
     Computes second derivative of exchange-correlation functional for given density grid.

     @param xcHessianGrid the exchange-correlation hessian grid object.
     @param factor the scaling factor of functional contribution.
     @param densityGrid the density grid object.
     */
    void compute(CXCHessianGrid& xcHessianGrid, const double factor, const CDensityGrid& densityGrid) const;

    /**
     Computes second derivative of exchange-correlation functional for given density grid.

     @param xcCubicHessianGrid the exchange-correlation hessian grid object.
     @param factor the scaling factor of functional contribution.
     @param densityGrid the density grid object.
     */
    void compute(CXCCubicHessianGrid& xcCubicHessianGrid, const double factor, const CDensityGrid& densityGrid) const;

    /**
     Converts primitive exchange-correlation functional object to text format and insert it into output text stream.

     @param output the output text stream.
     @param source the primitive exchange-correlation functional object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CPrimitiveFunctional& source);
};

// Forward-declaration from LibXC
struct xc_func_type;

using Component = std::tuple<double, xc_func_type>;

/**
 * Class CPrimitiveFunctional is a wrapper to the C functions and structs provided by LibXC.
 *
 * @author R. Di Remigio Eikås, Z. Rinkevicius
 */
class Functional
{
   private:
    /** Whether the functional is spin-polarized or not. */
    bool _polarized{false};

    /** @{ Exchange-correlation family. */
    /** Whether the functional is of generalized-gradient approximation type. */
    bool _isGGA{false};

    /** Whether the functional is of meta-generalized-gradient approximation type. */
    bool _isMetaGGA{false};

    /** Whether the functional is of global hybrid type. */
    bool _isGlobalHybrid{false};

    /** Whether the functional is of range-separated hybrid type. */
    bool _isRangeSeparatedHybrid{false};
    /** }@ */

    /** @{ Available derivative orders. */
    /** Whether has zeroth-order derivatives available. */
    bool _hasExc{true};

    /** Whether has first-order derivatives available. */
    bool _hasVxc{true};

    /** Whether has second-order derivatives available. */
    bool _hasFxc{true};

    /** Whether has third-order derivatives available. */
    bool _hasKxc{true};

    /** Whether has fourth-order derivatives available. */
    bool _hasLxc{true};
    /** }@ */

    int32_t _numberOfExchangeFunctionals{0};
    int32_t _numberOfCorrelationFunctionals{0};

    std::string _xcName{""};

    std::string _citation{""};

    /** The primitive exchange-correlation functionals and their coefficients. */
    std::vector<Component> _components;

   public:
    /** Creates an exchange-correlation functional object.
     *
     * @param[in] labels list of labels of component exchange and correlation functionals.
     * @param[in] coeffs list of coefficients for the components of the functional.
     * @param[in] polarized whether the functional is spin-polarized.
     */
    Functional(const std::vector<std::string>& labels, const std::vector<double>& coeffs, bool polarized);

    /** Copy-constructor.
     *
     * @param[in] src the exchange-correlation functional.
     */
    Functional(const Functional& src);

    /** Move-constructor.
     *
     *  @param[in] src the exchange-correlation functional object.
     */
    Functional(Functional&& src) noexcept;

    /** Destroys an exchange-correlation functional object. */
    ~Functional();

    /** Copy-assignment operator.
     *
     * @param[in] src the exchange-correlation functional object.
     */
    auto operator=(const Functional& src) -> Functional&;

    /** Move-assignment operator.
     *
     * @param[in] src the exchange-correlation functional object.
     */
    auto operator=(Functional&& src) noexcept -> Functional&;

    /** Equality operator.
     *
     * @param[in] other the exchange-correlation functional object.
     * @return true if exchange-correlation functional objects are equal, false otherwise.
     */
    auto operator==(const Functional& other) const -> bool;

    /** Inequality operator.
     *
     * @param[in] other the exchange-correlation functional object.
     * @return true if exchange-correlation functional objects are different, false otherwise.
     */
    auto operator!=(const Functional& other) const -> bool;

    /** Gets label of primitive exchange-correlation functional.
     *
     * @return the label.
     */
    auto getXCName() const -> std::string { return _xcName; }

    /** Get bibliographic reference of primitive functional.
     *
     * @return the reference.
     */
    auto getCitation() const -> std::string { return _citation; }

    auto getNumberOfExchangeFunctionals() const -> int32_t { return _numberOfExchangeFunctionals; }

    auto getNumberOfCorrelationFunctionals() const -> int32_t { return _numberOfCorrelationFunctionals; }

    /** String representation of primitive functional. */
    auto repr() const -> std::string;

    /** Computes first derivative of exchange-correlation functional for given density grid.
     *
     * @param[in,out] xcGradientGrid the exchange-correlation gradient grid object.
     * @param[in] densityGrid the density grid object.
     */
    auto compute(CXCGradientGrid& xcGradientGrid, const CDensityGrid& densityGrid) const -> void;

    /** Computes second derivative of exchange-correlation functional for given density grid.
     *
     * @param[in,out] xcHessianGrid the exchange-correlation hessian grid object.
     * @param[in] densityGrid the density grid object.
     */
    auto compute(CXCHessianGrid& xcHessianGrid, const CDensityGrid& densityGrid) const -> void;

    /** Computes third derivative of exchange-correlation functional for given density grid.
     *
     * @param[in,out] xcCubicHessianGrid the exchange-correlation third-derivative grid object.
     * @param[in] densityGrid the density grid object.
     */
    auto compute(CXCCubicHessianGrid& xcCubicHessianGrid, const CDensityGrid& densityGrid) const -> void;
};

std::ostream& operator<<(std::ostream& output, const Functional& source);

/** Get descriptiton, with bibliographic reference, LibXC in use.
 *
 * @return the description.
 */
auto getLibXCDescription() -> std::string;

namespace detail {
auto getXCName(const std::vector<Component>& x_funcs, const std::vector<Component>& c_funcs) -> std::string;

auto getFunctionalCitation(const xc_func_type& func) -> std::string;

auto getCitation(const std::vector<Component>& x_funcs, const std::vector<Component>& c_funcs) -> std::string;
}  // namespace detail

#endif /* PrimitiveFunctional_hpp */
