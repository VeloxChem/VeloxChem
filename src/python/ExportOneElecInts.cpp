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

#include "ExportOneElecInts.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AngularMomentumIntegrals.hpp"
#include "ElectricFieldIntegrals.hpp"
#include "ElectricFieldValues.hpp"
#include "ExportGeneral.hpp"
#include "ErrorHandler.hpp"
#include "LinearMomentumIntegrals.hpp"
#include "NuclearPotentialValues.hpp"
#include "QuadrupoleIntegrals.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_oneeints {

// Exports classes/functions in src/onee_ints to python

void
export_oneeints(py::module& m)
{
    m.def("compute_linear_momentum_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis) -> py::list {
                auto linmom = onee::computeLinearMomentumIntegrals(molecule, basis);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(linmom[0].values(), {linmom[0].getNumberOfRows(), linmom[0].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(linmom[1].values(), {linmom[1].getNumberOfRows(), linmom[1].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(linmom[2].values(), {linmom[2].getNumberOfRows(), linmom[2].getNumberOfColumns()}));
                return ret;
            },
            "Computes linear momentum integrals.");

    m.def("compute_angular_momentum_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<double>& origin) -> py::list {
                auto angmom = onee::computeAngularMomentumIntegrals(molecule, basis, origin);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(angmom[0].values(), {angmom[0].getNumberOfRows(), angmom[0].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(angmom[1].values(), {angmom[1].getNumberOfRows(), angmom[1].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(angmom[2].values(), {angmom[2].getNumberOfRows(), angmom[2].getNumberOfColumns()}));
                return ret;
            },
            "Computes angular momentum integrals.");

    m.def("compute_electric_field_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& dipole_moments) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_integrals: Expecting contiguous numpy arrays");
                auto        c_style_coords  = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                auto        c_style_moments = py::detail::check_flags(dipole_moments.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_coords && c_style_moments), errstyle);
                std::string errsize("compute_electric_field_integrals: Inconsistent number of dipole sites");
                errors::assertMsgCritical(
                        ((dipole_coords.shape(0) == dipole_moments.shape(0)) &&
                         (dipole_coords.shape(1) == 3) &&
                         (dipole_moments.shape(1) == 3)),
                        errsize);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto efield = onee::computeElectricFieldIntegrals(molecule, basis, dipole_coords.data(), dipole_moments.data(), ndipoles);
                return vlx_general::pointer_to_numpy(efield.values(), {efield.getNumberOfRows(), efield.getNumberOfColumns()});
            },
            "Computes electric field integrals.");

    m.def("compute_electric_field_values",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_values: Expecting contiguous numpy arrays");
                auto        c_style = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errsize("compute_electric_field_values: Inconsistent dimension of dipole coordinates");
                errors::assertMsgCritical(dipole_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_dipole_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto ef_vals = onee::computeElectricFieldValues(molecule, basis, dipole_coords.data(), ndipoles, D.data(), naos);
                return vlx_general::pointer_to_numpy(ef_vals.values(), {ef_vals.getNumberOfRows(), ef_vals.getNumberOfColumns()});
            },
            "Computes electric field values.");

    m.def("compute_nuclear_potential_values",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_values: Expecting contiguous numpy arrays");
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_electric_point_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto npot_vals = onee::computeNuclearPotentialValues(molecule, basis, point_coords.data(), npoints, D.data(), naos);
                return vlx_general::pointer_to_numpy(npot_vals.data(), {static_cast<int>(npot_vals.size())});
            },
            "Computes nuclear potential values.");

    m.def("compute_quadrupole_integrals",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const std::vector<double>& origin) -> py::list {
                auto mu = onee::computeQuadrupoleIntegrals(molecule, basis, origin);
                py::list ret;
                ret.append(vlx_general::pointer_to_numpy(mu[0].values(), {mu[0].getNumberOfRows(), mu[0].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(mu[1].values(), {mu[1].getNumberOfRows(), mu[1].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(mu[2].values(), {mu[2].getNumberOfRows(), mu[2].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(mu[3].values(), {mu[3].getNumberOfRows(), mu[3].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(mu[4].values(), {mu[4].getNumberOfRows(), mu[4].getNumberOfColumns()}));
                ret.append(vlx_general::pointer_to_numpy(mu[5].values(), {mu[5].getNumberOfRows(), mu[5].getNumberOfColumns()}));
                return ret;
            },
            "Computes quadrupole integrals.");
}

}  // namespace vlx_oneeints