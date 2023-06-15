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

#ifndef XCComponent_hpp
#define XCComponent_hpp

#include <xc.h>

#include <string>

/**
 Class CXCComponent implements XC functional component.

 @author X. Li
 */
class CXCComponent
{
    /**
     Libxc name of the XC functional component.
     */
    std::string _label{std::string("Undefined")};

    /**
     Scaling factor of the functional component.
     */
    double _scalingFactor{1.0};

    /**
     The Libxc functional component.
     */
    xc_func_type _func;

    /**
     Whether Libxc functional component has been initalized.
     */
    bool _initialized{false};

    /**
     Initializes libxc functional component.
     */
    void _init_libxc_func();

    /**
     Finalizes libxc functional component.
     */
    void _end_libxc_func();

    /**
     Resets libxc functional component.
     */
    void _reset_libxc_func();

   public:
    /**
     Creates an XC component object.
     */
    CXCComponent(const std::string& label, const double scalingFactor);

    /**
     Creates an XC component object by copying other XC component object.

     @param source the XC component object.
     */
    CXCComponent(const CXCComponent& source);

    /**
     Creates an XC component object by moving other XC component object.

     @param source the XC component object.
     */
    CXCComponent(CXCComponent&& source) noexcept;

    /**
     Destroys an XC component object.
     */
    ~CXCComponent();

    /**
     Assigns an XC component object by copying other XC component object.

     @param source the XC component object.
     */
    CXCComponent& operator=(const CXCComponent& source);

    /**
     Assigns an XC component object by moving other XC component object.

     @param source the XC component object.
     */
    CXCComponent& operator=(CXCComponent&& source) noexcept;

    /**
     Compares XC component object with other XC component object.

     @param other the XC component object.
     @return true if XC component objects are equal, false otherwise.
     */
    bool operator==(const CXCComponent& other) const;

    /**
     Compares XC component object with other XC component object.

     @param other the XC component object.
     @return true if XC component objects are not equal, false otherwise.
     */
    bool operator!=(const CXCComponent& other) const;

    /**
     Gets the Libxc name.

     @return the Libxc name.
     */
    std::string getLabel() const;

    /**
     Gets the scaling factor.

     @return the scaling factor.
     */
    double getScalingFactor() const;

    /**
     Gets pointer to the Libxc functional component.

     @return the pointer.
     */
    const xc_func_type* getFunctionalPointer() const;

    /**
     Checks if the functional component is LDA.

     @return true if the funcitonal is LDA, otherwise false.
     */
    bool isLDA() const;

    /**
     Checks if the functional component is GGA.

     @return true if the funcitonal is GGA, otherwise false.
     */
    bool isGGA() const;

    /**
     Checks if the functional component is meta-GGA.

     @return true if the funcitonal is meta-GGA, otherwise false.
     */
    bool isMetaGGA() const;
};

#endif /* XCComponent_hpp */
