//
//                           VELOXCHEM 1.0-RC
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

#ifndef RecursionFunction_hpp
#define RecursionFunction_hpp

#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <vector>

#include "RecursionTerm.hpp"

using def_rec_func_typ = std::vector<CRecursionTerm>(const CRecursionTerm&);

/**
 Class CRecursionFunction stores information about recursion function and
 provides methods to perform actions with stored recursion function.

 @author Z. Rinkevicius
 */
class CRecursionFunction
{
    /**
     The label of recursion function.
     */
    std::string _label;

    /**
     The action of recursion function.
     */
    std::function<def_rec_func_typ> _funcAction;

   public:
    /**
     Creates an empty recursion function object.
     */
    CRecursionFunction();

    /**
     Creates a recursion function object.

     @param label the label of recursion function.
     @param funcAction the action of recursion function.
     */
    CRecursionFunction(const std::string& label, const std::function<def_rec_func_typ>& funcAction);

    /**
     Creates a recursion function object by copying other recursion function
     object.

     @param source the recursion term object.
     */
    CRecursionFunction(const CRecursionFunction& source);

    /**
     Creates a recursion function object by moving other recursion function
     object.

     @param source the recursion function object.
     */
    CRecursionFunction(CRecursionFunction&& source) noexcept;

    /**
     Destroys a recursion function object.
     */
    ~CRecursionFunction();

    /**
     Assigns a recursion function object by copying other recursion function
     object.

     @param source the recursion function object.
     */
    CRecursionFunction& operator=(const CRecursionFunction& source);

    /**
     Assigns a recursion function object by moving other recursion function
     object.

     @param source the recursion function object.
     */
    CRecursionFunction& operator=(CRecursionFunction&& source) noexcept;

    /**
     Compares recursion function object with other recursion function object.

     @param other the recursion function object.
     @return true if recursion function objects are equal, false otherwise.
     */
    bool operator==(const CRecursionFunction& other) const;

    /**
     Compares recursion function object with other recursion function object.

     @param other the recursion function object.
     @return true if recursion function objects are not equal, false otherwise.
     */
    bool operator!=(const CRecursionFunction& other) const;

    /**
     Gets label of recursion function.

     @return the label.
     */
    std::string getLabel() const;

    /**
     Applies function defined by function action to recursion term object.

     @param recursionTerm the recursion term object.
     @return the vector of recursion term objects.
     */
    std::vector<CRecursionTerm> compute(const CRecursionTerm& recursionTerm) const;

    /**
     Checks if recursion function object has matching label.

     @param label the label of function.
     @return true if labels match, false otherwise.
     */
    bool isMatch(const std::string label) const;

    /**
     Converts recursion function object to text format and insert it into output text stream.

     @param output the output text stream.
     @param source the recursion function object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CRecursionFunction& source);
};

#endif /* RecursionFunction_hpp */
