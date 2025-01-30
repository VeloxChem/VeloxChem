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

#include "NewDispersionModel.hpp"

#include <vector>

#include "dftd4.h"

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CNewDispersionModel::CNewDispersionModel()
{
}

CNewDispersionModel::~CNewDispersionModel()
{
}

void
CNewDispersionModel::compute(const CMolecule& molecule, const std::string& xcLabel)
{
    auto net_charge = molecule.get_charge();

    auto natoms = molecule.number_of_atoms();

    _gradient = CDenseMatrix(natoms, 3);

    std::vector<double> atom_coords(natoms * 3);

    auto xyzcoord = molecule.coordinates("bohr");

    for (int i = 0; i < natoms; i++)
    {
        const auto rixyz = xyzcoord[i].coordinates();

        for (int d = 0; d < 3; d++)
        {
            atom_coords[i * 3 + d] = rixyz[d];
        }
    }

    auto identifiers = molecule.identifiers();

    std::string input_xc_label = std::string(xcLabel);

    if (format::lower_case(xcLabel) == std::string("lrc-wpbeh")) input_xc_label = std::string("lc-wpbeh");

    if (format::lower_case(xcLabel) == std::string("m06-l")) input_xc_label = std::string("m06l");

    std::vector<char> xc_label_char_vec(input_xc_label.begin(), input_xc_label.end());
    xc_label_char_vec.push_back('\0');

    dftd4_error error_handler = dftd4_new_error();

    dftd4_structure disp_mol = dftd4_new_structure(error_handler, natoms, identifiers.data(), atom_coords.data(), &net_charge, nullptr, nullptr);
    check_error_code(dftd4_check_error(error_handler), std::string("creating molecule"));

    dftd4_model disp_model = dftd4_new_d4_model(error_handler, disp_mol);
    check_error_code(dftd4_check_error(error_handler), std::string("creating d4 model"));

    dftd4_param disp_param = dftd4_load_rational_damping(error_handler, xc_label_char_vec.data(), true);
    check_error_code(dftd4_check_error(error_handler), std::string("loading parameters for " + xcLabel));

    dftd4_get_dispersion(error_handler, disp_mol, disp_model, disp_param, &_energy, _gradient.values(), nullptr);
    check_error_code(dftd4_check_error(error_handler), std::string("computing dispersion"));

    dftd4_delete_param(&disp_param);
    dftd4_delete_model(&disp_model);
    dftd4_delete_structure(&disp_mol);
    dftd4_delete_error(&error_handler);
}

void
CNewDispersionModel::check_error_code(const int error_code, const std::string& msg)
{
    errors::assertMsgCritical(error_code == 0, std::string("DFTD4 dispersion model error in ") + msg);
}

double
CNewDispersionModel::getEnergy() const
{
    return _energy;
}

CDenseMatrix
CNewDispersionModel::getGradient() const
{
    return _gradient;
}
