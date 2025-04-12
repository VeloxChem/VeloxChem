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

#include "GpuDevices.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

#include "DeviceProp.hpp"
#include "StringFormat.hpp"

CGpuDevices::CGpuDevices()

    : _namesOfDevices(std::vector<std::string>())

    , _globalMemoryOfDevices(std::vector<int64_t>())

    , _computeMajorCapabilityOfDevices(std::vector<int64_t>())

    , _computeMinorCapabilityOfDevices(std::vector<int64_t>())
{
    gpu::getDevicesProperty(_namesOfDevices, _globalMemoryOfDevices, _computeMajorCapabilityOfDevices, _computeMinorCapabilityOfDevices);
}

int64_t
CGpuDevices::getNumberOfDevices() const
{
    return static_cast<int64_t>(_namesOfDevices.size());
}

std::string
CGpuDevices::getString() const
{
    std::stringstream ss;

    std::string str("GPU Devices");

    ss << str << "\n" << std::string(str.size(), '=') << "\n\n";

    for (size_t i = 0; i < _namesOfDevices.size(); i++)
    {
        ss << "GPU device ID: " << std::to_string(i) << "\n";

        ss << "  Device name:             " << _namesOfDevices[i] << "\n";

        ss << "  Compute capability:      ";

        ss << std::to_string(_computeMajorCapabilityOfDevices[i]) << "." << std::to_string(_computeMinorCapabilityOfDevices[i]) << "\n";

        ss << "  Global memory:           ";

        ss << std::fixed << std::setprecision(0) << static_cast<double>(_globalMemoryOfDevices[i]) / std::pow(1024.0, 3.0) << " GB\n";
    }

    return ss.str();
}
