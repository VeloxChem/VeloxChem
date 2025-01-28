//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include <pybind11/pybind11.h>

#include "ExportDft.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ExportMoldata.hpp"
#include "ExportOneElecInts.hpp"
#include "ExportOrbdata.hpp"
#include "ExportVisualization.hpp"
#include "ExportT2CIntegrals.hpp"
#include "ExportT4CIntegrals.hpp"

PYBIND11_MODULE(veloxchemlib, m)
{
    vlx_general::export_general(m);

    vlx_math::export_math(m);

    vlx_moldata::export_moldata(m);

    vlx_orbdata::export_orbdata(m);

    vlx_dft::export_dft(m);

    vlx_oneeints::export_oneeints(m);

    vlx_visualization::export_visualization(m);

    vlx_t2cintegrals::export_t2cintegrals(m);

    vlx_t4cintegrals::export_t4cintegrals(m);
}
