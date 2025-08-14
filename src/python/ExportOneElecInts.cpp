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

#include "ExportOneElecInts.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AngularMomentumIntegrals.hpp"
#include "ElectricFieldFockGradient.hpp"
#include "ElectricFieldIntegrals.hpp"
#include "ElectricFieldPotentialGradient.hpp"
#include "ElectricFieldPotentialGradientAtMMSites.hpp"
#include "ElectricFieldPotentialHessian.hpp"
#include "ElectricFieldValues.hpp"
#include "ExportGeneral.hpp"
#include "ErrorHandler.hpp"
#include "LinearMomentumIntegrals.hpp"
#include "NuclearPotentialErfGradient.hpp"
#include "NuclearPotentialErfValues.hpp"
#include "NuclearPotentialHessian110.hpp"
#include "NuclearPotentialHessian101.hpp"
#include "QuadrupoleIntegrals.hpp"
#include "OldOneElecIntsDrivers.hpp"
#include "OverlapHessian101.hpp"
#include "KineticEnergyHessian101.hpp"
#include "TCOSFock1GOST.hpp"
#include "TCOSValues.hpp"
#include "TCOPFock2GOST.hpp"
#include "TCOPValues.hpp"
#include "TCODValues.hpp"
#include "TCOFValues.hpp"
#include "TCOSGradient.hpp"
#include "TCOPGradient.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_oneeints {

// COldElectricDipoleIntegralsDriver constructor (for backward compatibility only)
static std::shared_ptr<COldElectricDipoleIntegralsDriver>
COldElectricDipoleIntegralsDriver_from_obj(py::object py_comm)
{
    return std::make_shared<COldElectricDipoleIntegralsDriver>();
}

// COldLinearMomentumIntegralsDriver constructor (for backward compatibility only)
static std::shared_ptr<COldLinearMomentumIntegralsDriver>
COldLinearMomentumIntegralsDriver_from_obj(py::object py_comm)
{
    return std::make_shared<COldLinearMomentumIntegralsDriver>();
}

// COldAngularMomentumIntegralsDriver constructor (for backward compatibility only)
static std::shared_ptr<COldAngularMomentumIntegralsDriver>
COldAngularMomentumIntegralsDriver_from_obj(py::object py_comm)
{
    return std::make_shared<COldAngularMomentumIntegralsDriver>();
}

// Exports classes/functions in src/onee_ints to python

void
export_oneeints(py::module& m)
{
    // COldOneElecIntsMatrix class (for backward compatibility only)
    PyClass<COldOneElecIntsMatrix>(m, "OldOneElecIntsMatrix")
        .def(py::init<const std::vector<CDenseMatrix>&>())
        .def(
            "x_to_numpy",
            [](const COldOneElecIntsMatrix& self) -> py::array_t<double> {
                const auto mat = self.get_matrix(0);
                return vlx_general::pointer_to_numpy(mat.values(), {mat.getNumberOfRows(), mat.getNumberOfColumns()});
            },
            "Gets x matrix as numpy array.")
        .def(
            "y_to_numpy",
            [](const COldOneElecIntsMatrix& self) -> py::array_t<double> {
                const auto mat = self.get_matrix(1);
                return vlx_general::pointer_to_numpy(mat.values(), {mat.getNumberOfRows(), mat.getNumberOfColumns()});
            },
            "Gets y matrix as numpy array.")
        .def(
            "z_to_numpy",
            [](const COldOneElecIntsMatrix& self) -> py::array_t<double> {
                const auto mat = self.get_matrix(2);
                return vlx_general::pointer_to_numpy(mat.values(), {mat.getNumberOfRows(), mat.getNumberOfColumns()});
            },
            "Gets z matrix as numpy array.");

    // COldElectricDipoleIntegralsDriver class (for backward compatibility only)
    PyClass<COldElectricDipoleIntegralsDriver>(m, "ElectricDipoleIntegralsDriver")
        .def(py::init<>())
        .def(py::init(&COldElectricDipoleIntegralsDriver_from_obj), "comm"_a)
        .def("compute", &COldElectricDipoleIntegralsDriver::compute, "Computes electric dipole integrals.");

    // COldLinearMomentumIntegralsDriver class (for backward compatibility only)
    PyClass<COldLinearMomentumIntegralsDriver>(m, "LinearMomentumIntegralsDriver")
        .def(py::init<>())
        .def(py::init(&COldLinearMomentumIntegralsDriver_from_obj), "comm"_a)
        .def("compute", &COldLinearMomentumIntegralsDriver::compute, "Computes electric dipole integrals.");

    // COldAngularMomentumIntegralsDriver class (for backward compatibility only)
    PyClass<COldAngularMomentumIntegralsDriver>(m, "AngularMomentumIntegralsDriver")
        .def(py::init<>())
        .def(py::init(&COldAngularMomentumIntegralsDriver_from_obj), "comm"_a)
        .def("compute", &COldAngularMomentumIntegralsDriver::compute, "Computes electric dipole integrals.");

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
            "Computes linear momentum integrals.",
             "molecule"_a,
             "basis"_a);

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
            "Computes angular momentum integrals.",
             "molecule"_a,
             "basis"_a,
             "origin"_a = std::vector<double>({0.0, 0.0, 0.0}));

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
            "Computes electric field integrals.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "dipole_moments"_a);

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
            "Computes electric field values.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "density"_a);

    m.def("compute_electric_field_potential_gradient",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& dipole_moments,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_potential_gradient: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(dipole_moments.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2), errstyle);
                std::string errsize("compute_electric_field_potential_gradient: Inconsistent dimension of dipole coordinates/moments");
                errors::assertMsgCritical(dipole_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(dipole_moments.shape(1) == 3, errsize);
                std::string errshape("compute_electric_field_potential_gradient: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto ef_grad = onee::computeElectricFieldPotentialGradient(molecule, basis, dipole_coords.data(), dipole_moments.data(), ndipoles, D.data(), naos);
                return vlx_general::pointer_to_numpy(ef_grad.values(), {ef_grad.getNumberOfRows(), ef_grad.getNumberOfColumns()});
            },
            "Computes electric field integrals contribution to molecular gradient.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "dipole_moments"_a,
             "density"_a);

    m.def("compute_electric_field_fock_gradient",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& dipole_moments,
               const int                  qm_atom_index) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_fock_gradient: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(dipole_moments.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2), errstyle);
                std::string errsize("compute_electric_field_fock_gradient: Inconsistent dimension of dipole coordinates/moments");
                errors::assertMsgCritical(dipole_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(dipole_moments.shape(1) == 3, errsize);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto ef_fock_grad = onee::computeElectricFieldFockGradient(molecule, basis, dipole_coords.data(), dipole_moments.data(), ndipoles, qm_atom_index);
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                CDenseMatrix ret(3, naos * naos);
                for (int n = 0; n < 3; n++) std::memcpy(ret.row(n), ef_fock_grad[n].values(), naos * naos * sizeof(double));
                return vlx_general::pointer_to_numpy(ret.values(), {3, naos, naos});
            },
            "Computes electric field integrals contribution to fock gradient.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "dipole_moments"_a,
             "qm_atom_index"_a);

    m.def("compute_electric_field_potential_gradient_for_mm",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& D,
               const int                  atom_idx) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_potential_gradient_for_mm: Expecting contiguous numpy arrays");
                auto        c_style = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errsize("compute_electric_field_potential_gradient_for_mm: Inconsistent dimension of dipole coordinates/moments");
                errors::assertMsgCritical(dipole_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_field_potential_gradient_for_mm: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto ef_grad_for_mm = onee::computeElectricFieldPotentialGradientAtMMSites(molecule, basis, dipole_coords.data(), ndipoles, D.data(), naos, atom_idx);
                return vlx_general::pointer_to_numpy(ef_grad_for_mm.values(), {3, ndipoles, 3});
            },
            "Computes electric field potential gradient for solving MM induced dipoles.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "density"_a,
             "qm_atom_index"_a);

    m.def("compute_electric_field_potential_hessian",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& dipole_coords,
               const py::array_t<double>& dipole_moments,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_electric_field_potential_hessian: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(dipole_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(dipole_moments.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2), errstyle);
                std::string errsize("compute_electric_field_potential_hessian: Inconsistent dimension of dipole coordinates/moments");
                errors::assertMsgCritical(dipole_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(dipole_moments.shape(1) == 3, errsize);
                std::string errshape("compute_electric_field_potential_hessian: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto ndipoles = static_cast<int>(dipole_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto ef_hess = onee::computeElectricFieldPotentialHessian(molecule, basis, dipole_coords.data(), dipole_moments.data(), ndipoles, D.data(), naos);
                return vlx_general::pointer_to_numpy(ef_hess.values(), {ef_hess.getNumberOfRows(), ef_hess.getNumberOfColumns()});
            },
            "Computes electric field integrals contribution to molecular Hessian.",
             "molecule"_a,
             "basis"_a,
             "dipole_coords"_a,
             "dipole_moments"_a,
             "density"_a);

    m.def("compute_nuclear_potential_values",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_values: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(D.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2), errstyle);
                std::string errshape("compute_electric_point_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto npot_vals = onee::computeNuclearPotentialErfValues(molecule, basis, point_coords.data(), npoints, D.data(), naos);
                return vlx_general::pointer_to_numpy(npot_vals.data(), {static_cast<int>(npot_vals.size())});
            },
            "Computes nuclear potential values.",
             "molecule"_a,
             "basis"_a,
             "point_coords"_a,
             "density"_a);

    m.def("compute_nuclear_potential_erf_values",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& D,
               const py::array_t<double>& omega) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_erf_values: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(D.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(omega.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_electric_point_erf_values: Inconsistent sizes");
                errors::assertMsgCritical(omega.shape(0) == point_coords.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_point_erf_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto npot_vals = onee::computeNuclearPotentialErfValues(molecule, basis, point_coords.data(), npoints, D.data(), naos, omega.data());
                return vlx_general::pointer_to_numpy(npot_vals.data(), {static_cast<int>(npot_vals.size())});
            },
            "Computes nuclear potential values.",
             "molecule"_a,
             "basis"_a,
             "point_coords"_a,
             "density"_a,
             "omega"_a);

    m.def("compute_nuclear_potential_erf_gradient",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& point_charges,
               const py::array_t<double>& D,
               const py::array_t<double>& omega,
               const py::array_t<int>&    atom_indices) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_erf_gradient: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_charges.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(D.ptr(), py::array::c_style);
                auto        c_style_4 = py::detail::check_flags(omega.ptr(), py::array::c_style);
                auto        c_style_5 = py::detail::check_flags(atom_indices.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4 && c_style_5), errstyle);
                std::string errsize("compute_electric_point_erf_gradient: Inconsistent sizes");
                errors::assertMsgCritical(omega.shape(0) == point_coords.shape(0), errsize);
                errors::assertMsgCritical(point_charges.shape(0) == point_coords.shape(0), errsize);
                errors::assertMsgCritical(point_charges.shape(0) == atom_indices.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_point_erf_gradient: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto grad = onee::computeNuclearPotentialErfGradient(
                    molecule, basis, point_coords.data(), point_charges.data(), npoints, D.data(), naos, omega.data(), atom_indices.data());
                return vlx_general::pointer_to_numpy(grad.values(), {grad.getNumberOfRows(), grad.getNumberOfColumns()});
            },
            "Computes nuclear potential gradient on charges.",
             "molecule"_a,
             "basis"_a,
             "point_coords"_a,
             "point_charges"_a,
             "density"_a,
             "omega"_a,
             "atom_indices"_a);

    m.def("compute_nuclear_potential_hessian_110",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& point_charges,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_hessian_110: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_charges.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(D.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_electric_point_hessian_110: Inconsistent sizes");
                errors::assertMsgCritical(point_charges.shape(0) == point_coords.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_point_hessian_110: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto hess = onee::computeNuclearPotentialHessian110(
                    molecule, basis, point_coords.data(), point_charges.data(), npoints, D.data(), naos);
                return vlx_general::pointer_to_numpy(hess.values(), {hess.getNumberOfRows(), hess.getNumberOfColumns()});
            },
            "Computes nuclear potential Hessian contribution.",
             "molecule"_a,
             "basis"_a,
             "point_coords"_a,
             "point_charges"_a,
             "density"_a);

    m.def("compute_nuclear_potential_hessian_101",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& point_coords,
               const py::array_t<double>& point_charges,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_nuclear_potential_hessian_101: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_charges.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(D.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_electric_point_hessian_101: Inconsistent sizes");
                errors::assertMsgCritical(point_charges.shape(0) == point_coords.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                std::string errshape("compute_electric_point_hessian_101: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(D.shape(0));
                auto hess = onee::computeNuclearPotentialHessian101(
                    molecule, basis, point_coords.data(), point_charges.data(), npoints, D.data(), naos);
                return vlx_general::pointer_to_numpy(hess.values(), {hess.getNumberOfRows(), hess.getNumberOfColumns()});
            },
            "Computes nuclear potential Hessian contribution.",
             "molecule"_a,
             "basis"_a,
             "point_coords"_a,
             "point_charges"_a,
             "density"_a);

    m.def("compute_overlap_hessian_101",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& W) -> py::array_t<double> {
                std::string errstyle("compute_overlap_hessian_101: Expecting contiguous numpy arrays");
                auto        c_style = py::detail::check_flags(W.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_electric_point_hessian_101: Expecting square matrix W");
                errors::assertMsgCritical(W.shape(0) == W.shape(1), errshape);
                auto naos = static_cast<int>(W.shape(0));
                auto hess = onee::computeOverlapHessian101(molecule, basis, W.data(), naos);
                return vlx_general::pointer_to_numpy(hess.values(), {hess.getNumberOfRows(), hess.getNumberOfColumns()});
            },
            "Computes overlap Hessian contribution.",
             "molecule"_a,
             "basis"_a,
             "omega_ao"_a);

    m.def("compute_kinetic_energy_hessian_101",
            [](const CMolecule&           molecule,
               const CMolecularBasis&     basis,
               const py::array_t<double>& D) -> py::array_t<double> {
                std::string errstyle("compute_kinetic_energy_hessian_101: Expecting contiguous numpy arrays");
                auto        c_style = py::detail::check_flags(D.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_electric_point_hessian_101: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto naos = static_cast<int>(D.shape(0));
                auto hess = onee::computeKineticEnergyHessian101(molecule, basis, D.data(), naos);
                return vlx_general::pointer_to_numpy(hess.values(), {hess.getNumberOfRows(), hess.getNumberOfColumns()});
            },
            "Computes kinetic energy Hessian contribution.",
             "molecule"_a,
             "basis"_a,
             "density"_a);

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
            "Computes quadrupole integrals.",
             "molecule"_a,
             "basis"_a,
             "origin"_a = std::vector<double>({0.0, 0.0, 0.0}));

    m.def("compute_tco_s_fock",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp) -> py::array_t<double> {
                std::string errstyle("compute_tco_s_fock: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_tco_s_fock: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_s = onee::computeFock1GOSTcontrib(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), naos);
                return vlx_general::pointer_to_numpy(tco_s.values(), {tco_s.getNumberOfRows(), tco_s.getNumberOfColumns()});
            },
            "Computes three-center s-type overlap integrals.",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a);

    m.def("compute_tco_s_values",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& D
                ) -> py::array_t<double> {
                std::string errstyle("compute_tco_s_values: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_tco_s_values: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_tco_s_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_s_vals = onee::computeTCOSValues(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), D.data(), naos);
                return vlx_general::pointer_to_numpy(tco_s_vals.data(), {static_cast<int>(tco_s_vals.size())});
            },
            "Computes three-center s-type overlap integrals and contracts with density matrix",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "density"_a);

    m.def("compute_tco_p_fock",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& point_norms) -> py::array_t<double> {
                std::string errstyle("compute_tco_p_fock: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                auto        c_style_4 = py::detail::check_flags(point_norms.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
                std::string errsize("compute_tco_p_fock: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_norms.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(point_norms.shape(1) == 3, errsize);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_p = onee::computeFock2GOSTcontrib(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), point_norms.data(), naos);
                return vlx_general::pointer_to_numpy(tco_p.values(), {tco_p.getNumberOfRows(), tco_p.getNumberOfColumns()});
            },
            "Computes three-center p-type overlap integrals and takes dot product with grid point normal vector.",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "point_norms"_a);

    m.def("compute_tco_p_values",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& point_norms,
                const py::array_t<double>& D
                ) -> py::array_t<double> {
                std::string errstyle("compute_tco_p_values: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                auto        c_style_4 = py::detail::check_flags(point_norms.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
                std::string errsize("compute_tco_p_values: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_norms.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(point_norms.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_tco_p_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_p_vals = onee::computeTCOPValues(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), point_norms.data(), D.data(), naos);
                return vlx_general::pointer_to_numpy(tco_p_vals.data(), {static_cast<int>(tco_p_vals.size())});
            },
            "Computes three-center p-type overlap integrals.",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "point_norms"_a,
                "density"_a);

    m.def("compute_tco_d_values",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& D
                ) -> py::array_t<double> {
                std::string errstyle("compute_tco_d_fock: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_tco_d_fock: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_electric_point_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_d_vals = onee::computeTCODValues(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), D.data(), naos);
                return vlx_general::pointer_to_numpy(tco_d_vals.data(), {static_cast<int>(tco_d_vals.size())});
            },
            "Computes three-center s-type d-function overlap integrals.",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "density"_a);

    m.def("compute_tco_f_values",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& point_norms,
                const py::array_t<double>& D
                ) -> py::array_t<double> {
                std::string errstyle("compute_tco_f_fock: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                auto        c_style_4 = py::detail::check_flags(point_norms.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
                std::string errsize("compute_tco_f_fock: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_norms.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(point_norms.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_electric_point_values: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto e_tilde_vals = onee::computeTCOFValues(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), point_norms.data(), D.data(), naos);
                return vlx_general::pointer_to_numpy(e_tilde_vals.data(), {static_cast<int>(e_tilde_vals.size())});
            },
            "Computes linear combination of three-center p-type f-function overlap integrals.",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "point_norms"_a,
                "density"_a);

    m.def("compute_tco_s_gradient",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& D
                ) -> py::list {
                std::string errstyle("compute_tco_p_gradient: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
                std::string errsize("compute_tco_p_gradient: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_tco_p_gradient: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_s_grad = onee::computeTCOSGradient(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), D.data(), naos);
                py::list ret;
                for (int c = 0; c < npoints; c++) 
                ret.append(vlx_general::pointer_to_numpy(tco_s_grad[c].values(), {tco_s_grad[c].getNumberOfRows(), tco_s_grad[c].getNumberOfColumns()}));
                return ret;
            },
            "Computes three-center s-type overlap integrals gradient contribution and contracts with density matrix",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "density"_a);

    m.def("compute_tco_p_gradient",
            [](const CMolecule&           molecule,
                const CMolecularBasis&     basis,
                const py::array_t<double>& point_coords,
                const py::array_t<double>& point_exp,
                const py::array_t<double>& point_amp,
                const py::array_t<double>& point_norms,
                const py::array_t<double>& D
                ) -> py::list {
                std::string errstyle("compute_tco_p_gradient: Expecting contiguous numpy arrays");
                auto        c_style_1 = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                auto        c_style_2 = py::detail::check_flags(point_exp.ptr(), py::array::c_style);
                auto        c_style_3 = py::detail::check_flags(point_amp.ptr(), py::array::c_style);
                auto        c_style_4 = py::detail::check_flags(point_norms.ptr(), py::array::c_style);
                errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
                std::string errsize("compute_tco_p_gradient: Inconsistent sizes");
                errors::assertMsgCritical(point_coords.shape(0) == point_exp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(0) == point_amp.shape(0), errsize);
                errors::assertMsgCritical(point_coords.shape(1) == 3, errsize);
                errors::assertMsgCritical(point_norms.shape(1) == 3, errsize);
                auto        c_style = py::detail::check_flags(point_coords.ptr(), py::array::c_style);
                errors::assertMsgCritical(c_style, errstyle);
                std::string errshape("compute_tco_p_gradient: Expecting square matrix D");
                errors::assertMsgCritical(D.shape(0) == D.shape(1), errshape);
                auto npoints = static_cast<int>(point_coords.shape(0));
                auto naos = static_cast<int>(basis.dimensions_of_basis());
                auto tco_p_grad = onee::computeTCOPGradient(molecule, basis, point_coords.data(), npoints, point_exp.data(), point_amp.data(), point_norms.data(), D.data(), naos);
                py::list ret;
                for (int c = 0; c < npoints; c++) 
                ret.append(vlx_general::pointer_to_numpy(tco_p_grad[c].values(), {tco_p_grad[c].getNumberOfRows(), tco_p_grad[c].getNumberOfColumns()}));
                return ret;
            },
            "Computes three-center p-type overlap integrals gradient contribution and contracts with density matrix",
                "molecule"_a,
                "basis"_a,
                "point_coords"_a,
                "point_exp"_a,
                "point_amp"_a,
                "point_norms"_a,
                "density"_a);

}

}  // namespace vlx_oneeints
