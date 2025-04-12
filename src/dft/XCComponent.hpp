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
    xc_func_type* _func;

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
    xc_func_type* getFunctionalPointer() const;

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
