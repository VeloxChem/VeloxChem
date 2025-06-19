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

#include "GpuRuntime.hpp"


#include "DeviceProp.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"

namespace gpu {  // gpu namespace

void
getDevicesProperty(std::vector<std::string>& namesOfDevices,
                   std::vector<int64_t>&     globalMemoryOfDevices,
                   std::vector<int64_t>&     computeMajorCapabilityOfDevices,
                   std::vector<int64_t>&     computeMinorCapabilityOfDevices)
{
    int device_count = 0;

    gpuSafe(gpuGetDeviceCount(&device_count));

    for (int i = 0; i < device_count; i++)
    {
        gpuDeviceProp prop;

        gpuSafe(gpuGetDeviceProperties(&prop, i));

        namesOfDevices.push_back(std::string(prop.name));

        globalMemoryOfDevices.push_back(static_cast<int64_t>(prop.totalGlobalMem));

        computeMajorCapabilityOfDevices.push_back(static_cast<int64_t>(prop.major));

        computeMinorCapabilityOfDevices.push_back(static_cast<int64_t>(prop.minor));
    }
}

}  // namespace gpu
