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

#ifndef GpuDevices_hpp
#define GpuDevices_hpp

#include <cstdint>
#include <string>
#include <vector>

/**
 Class GpuDevices stores data about available GPU devices and provides these
 devices managment functions.

 @author Z. Rinkevicius
 */
class CGpuDevices
{
    /**
     The vector of names of GPU devices.
     */
    std::vector<std::string> _namesOfDevices;

    /**
     The vector of available global memory on GPU devices.
     */
    std::vector<int64_t> _globalMemoryOfDevices;

    /**
     The vector of compute major capabilities on GPU devices.
     */
    std::vector<int64_t> _computeMajorCapabilityOfDevices;

    /**
     The vector of compute minor capabilities on GPU devices.
     */
    std::vector<int64_t> _computeMinorCapabilityOfDevices;

   public:
    /**
     Creates an GPU devices object.
     */
    CGpuDevices();

    /**
     Gets number of GPU devices.

     @return the number of devices.
     */
    int64_t getNumberOfDevices() const;

    /**
     Gets string representation of GPU devices object.

     @return the representation string.
     */
    std::string getString() const;
};

#endif /* GpuDevices_hpp */
