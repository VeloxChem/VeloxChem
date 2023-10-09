//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#include "CudaDevices.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

#include "DeviceProp.hpp"
#include "StringFormat.hpp"

CCudaDevices::CCudaDevices()

    : _namesOfDevices(std::vector<std::string>())

    , _globalMemoryOfDevices(std::vector<int64_t>())

    , _computeMajorCapabilityOfDevices(std::vector<int64_t>())

    , _computeMinorCapabilityOfDevices(std::vector<int64_t>())
{
#ifdef ENABLE_GPU
    gpu::getDevicesProperty(_namesOfDevices, _globalMemoryOfDevices, _computeMajorCapabilityOfDevices, _computeMinorCapabilityOfDevices);
#endif
}

int64_t
CCudaDevices::getNumberOfDevices() const
{
    return static_cast<int64_t>(_namesOfDevices.size());
}

std::string
CCudaDevices::getString() const
{
    std::stringstream ss;

    std::string str("CUDA Compute Capable Devices");

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
