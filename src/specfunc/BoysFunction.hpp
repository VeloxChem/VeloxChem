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

#ifndef BoysFunction_hpp
#define BoysFunction_hpp

#include <array>
#include <cstdint>

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"

/**
 Defines alias to STL 2D array used for storing Boys function recursion
 parameters.
 */
using CBFTable = std::array<std::array<double, 7>, 121>;

/**
 Class CBoysFunction implements computation of Boys function.

 @author Z. Rinkevicius
 */
class CBoysFunction
{
    /**
     The order of Boys function.
     */
    int32_t _order;

    /**
     The vector of Boys function recursion parameters for arguments between
     0.0-12.0.
     */
    CMemBlock2D<double> _table;

    /**
     Sets vector of recursion parameters required for Boys function computation.
     */
    void _setTable();

    /**
     Computes vector of 0-order Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF00(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..1 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF01(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..2 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF02(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..3 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF03(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..4 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF04(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..5 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF05(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..6 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF06(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..7 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF07(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..8 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF08(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..9 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF09(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..10 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF10(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..11 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF11(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..12 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF12(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..13 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF13(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..14 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF14(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..15 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF15(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..16 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF16(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..17 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF17(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..18 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF18(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..19 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF19(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..20 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF20(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..21 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF21(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..22 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF22(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..23 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF23(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..24 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF24(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..25 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF25(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..26 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF26(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..27 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF27(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Computes vector of 0..28 orders Boys function values for given vector of
     arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     */
    void _computeBF28(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments) const;

    /**
     Loads recursion paramaters table into specified part of recursion parameters
     vector.

     @param bfData the recursion parameters table.
     @param identifier the identifier of part of recursion parameters vector.
     */
    void _loadTable(const CBFTable& bfData, const int32_t identifier);

    /**
     Generates and loads recursion parameters for Boys function of 0-th order.
     */
    void _generateTable00();

    /**
     Generates and loads recursion parameters for Boys function of 1-th order.
     */
    void _generateTable01();

    /**
     Generates and loads recursion parameters for Boys function of 2-th order.
     */
    void _generateTable02();

    /**
     Generates and loads recursion parameters for Boys function of 3-th order.
     */
    void _generateTable03();

    /**
     Generates and loads recursion parameters for Boys function of 4-th order.
     */
    void _generateTable04();

    /**
     Generates and loads recursion parameters for Boys function of 5-th order.
     */
    void _generateTable05();

    /**
     Generates and loads recursion parameters for Boys function of 6-th order.
     */
    void _generateTable06();

    /**
     Generates and loads recursion parameters for Boys function of 7-th order.
     */
    void _generateTable07();

    /**
     Generates and loads recursion parameters for Boys function of 8-th order.
     */
    void _generateTable08();

    /**
     Generates and loads recursion parameters for Boys function of 9-th order.
     */
    void _generateTable09();

    /**
     Generates and loads recursion parameters for Boys function of 10-th order.
     */
    void _generateTable10();

    /**
     Generates and loads recursion parameters for Boys function of 11-th order.
     */
    void _generateTable11();

    /**
     Generates and loads recursion parameters for Boys function of 12-th order.
     */
    void _generateTable12();

    /**
     Generates and loads recursion parameters for Boys function of 13-th order.
     */
    void _generateTable13();

    /**
     Generates and loads recursion parameters for Boys function of 14-th order.
     */
    void _generateTable14();

    /**
     Generates and loads recursion parameters for Boys function of 15-th order.
     */
    void _generateTable15();

    /**
     Generates and loads recursion parameters for Boys function of 16-th order.
     */
    void _generateTable16();

    /**
     Generates and loads recursion parameters for Boys function of 17-th order.
     */
    void _generateTable17();

    /**
     Generates and loads recursion parameters for Boys function of 18-th order.
     */
    void _generateTable18();

    /**
     Generates and loads recursion parameters for Boys function of 19-th order.
     */
    void _generateTable19();

    /**
     Generates and loads recursion parameters for Boys function of 20-th order.
     */
    void _generateTable20();

    /**
     Generates and loads recursion parameters for Boys function of 21-th order.
     */
    void _generateTable21();

    /**
     Generates and loads recursion parameters for Boys function of 22-th order.
     */
    void _generateTable22();

    /**
     Generates and loads recursion parameters for Boys function of 23-th order.
     */
    void _generateTable23();

    /**
     Generates and loads recursion parameters for Boys function of 24-th order.
     */
    void _generateTable24();

    /**
     Generates and loads recursion parameters for Boys function of 25-th order.
     */
    void _generateTable25();

    /**
     Generates and loads recursion parameters for Boys function of 26-th order.
     */
    void _generateTable26();

    /**
     Generates and loads recursion parameters for Boys function of 27-th order.
     */
    void _generateTable27();

    /**
     Generates and loads recursion parameters for Boys function of 28-th order.
     */
    void _generateTable28();

   public:
    /**
     Creates an unitialized Boys function evaluator object.
     */
    CBoysFunction();

    /**
     Creates a Boys function evaluator object.

     @param order the max. order of Boys function.
     */
    CBoysFunction(const int32_t order);

    /**
     Destroys a Boys function evaluator object.
     */
    ~CBoysFunction();

    /**
     Computes Boys function values up to specified order (inclusively) for given
     vector of arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param iOrder the order of Boys function.
     */
    void compute(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t iOrder) const;

    /**
     Computes Boys function values up to specified order (inclusively) for given
     vector of arguments.

     @param values the vector of Boys function values.
     @param arguments the vector of Boys function arguments.
     @param nArguments the number of Boys function arguments used in
            computations.
     @param iOrder the order of Boys function.
     */
    void compute(CMemBlock2D<double>& values, const CMemBlock<double>& arguments, const int32_t nArguments, const int32_t iOrder) const;
};

#endif /* BoysFunction_hpp */
