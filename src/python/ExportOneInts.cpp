//
//                           VELOXCHEM 1.0-RC
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

#include "ExportOneInts.hpp"

#include <mpi.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <memory>
#include <string>

#include "AngularMomentumIntegralsDriver.hpp"
#include "AngularMomentumMatrix.hpp"
#include "CartesianComponents.hpp"
#include "DenseMatrix.hpp"
#include "ElectricDipoleIntegralsDriver.hpp"
#include "ElectricDipoleMatrix.hpp"
#include "ElectricFieldIntegralsDriver.hpp"
#include "ElectricFieldMatrix.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "KineticEnergyIntegralsDriver.hpp"
#include "KineticEnergyMatrix.hpp"
#include "LinearMomentumIntegralsDriver.hpp"
#include "LinearMomentumMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "NuclearPotentialIntegralsDriver.hpp"
#include "NuclearPotentialMatrix.hpp"
#include "OverlapIntegralsDriver.hpp"
#include "OverlapMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_oneints {  // vlx_oneints namespace
// Exports classes/functions in src/oneints to python

void
export_oneints(py::module& m)
{
    // COverlapMatrix class

    auto py_overlap_matrix = bind_operator_matrix<COverlapMatrix>(m, "OverlapMatrix");
    py_overlap_matrix.def("get_ortho_matrix", &COverlapMatrix::getOrthogonalizationMatrix);

    // COverlapIntegralsDriver class

    auto py_ovl_driver = bind_integrals_driver<COverlapIntegralsDriver>(m, "OverlapIntegralsDriver", "overlap");

    // CKineticEnergyMatrix class

    auto py_kinetic_matrix = bind_operator_matrix<CKineticEnergyMatrix>(m, "KineticEnergyMatrix");
    py_kinetic_matrix.def("get_energy", &CKineticEnergyMatrix::getKineticEnergy);

    // CKineticEnergyIntegralsDriver class

    auto py_kin_driver = bind_integrals_driver<CKineticEnergyIntegralsDriver>(m, "KineticEnergyIntegralsDriver", "kinetic energy");

    // CNuclearPotentialMatrix class

    auto py_npot_matrix = bind_operator_matrix<CNuclearPotentialMatrix>(m, "NuclearPotentialMatrix");
    py_npot_matrix
        .def(
            "reduce_sum",
            [](CNuclearPotentialMatrix& obj, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                obj.reduce_sum(rank, nodes, comm);
            },
            "Sum-reduce nuclear potential matrix object from all MPI ranks within the communicator into nuclear potential matrix object on "
            "master node.",
            "rank"_a,
            "nodes"_a,
            "comm"_a)
        .def("get_energy", &CNuclearPotentialMatrix::getNuclearPotentialEnergy);

    // CNuclearPotentialIntegralsDriver class

    PyClass<CNuclearPotentialIntegralsDriver>(m, "NuclearPotentialIntegralsDriver")
        .def(py::init(&vlx_general::create<CNuclearPotentialIntegralsDriver>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&>()(&CNuclearPotentialIntegralsDriver::compute, py::const_),
             "Compute AO-basis nuclear potential integrals for given molecule and basis",
             "molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecule&>()(&CNuclearPotentialIntegralsDriver::compute,
                                                                                                       py::const_),
             "Compute AO-basis nuclear potential integrals for given molecule and basis, using charges from the second molecule object as sources",
             "molecule"_a,
             "basis"_a,
             "charges"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&, const CMolecule&>()(
                 &CNuclearPotentialIntegralsDriver::compute, py::const_),
             "Compute AO-basis nuclear potential integrals for given molecule and bases, using charges from the second molecule object as sources",
             "molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a,
             "charges"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecule&>()(
                 &CNuclearPotentialIntegralsDriver::compute, py::const_),
             "Compute AO-basis nuclear potential integrals for two molecules in given basis, using charges from the third molecule object as sources",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a,
             "charges"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&, const CMolecule&>()(
                 &CNuclearPotentialIntegralsDriver::compute, py::const_),
             "Compute AO-basis nuclear potential integrals for two molecules, each with its own basis. Use the charges from the third molecule "
             "object as sources",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a,
             "charges"_a)
        .def("compute",
             [](const CNuclearPotentialIntegralsDriver&        obj,
                const CMolecule&                               molecule,
                const CMolecularBasis&                         basis,
                const py::array_t<double>&                     charges,
                const py::array_t<double, py::array::f_style>& coordinates) -> CNuclearPotentialMatrix {
                 auto chg = vlx_general::numpy_to_memblock(charges);
                 auto xyz = vlx_general::numpy_to_memblock2d(coordinates);

                 return obj.compute(molecule, basis, chg, xyz);
             });

    // CElectricDipoleMatrix class

    PyClass<CElectricDipoleMatrix>(m, "ElectricDipoleMatrix")
        .def(py::init<>())
        .def(py::init<const std::array<CDenseMatrix, 3>&, const std::array<double, 3>&>(), "matrices"_a, "origin"_a)
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&, const double, const double, const double>(),
             "xMatrix"_a,
             "yMatrix"_a,
             "zMatrix"_a,
             "x_0"_a,
             "y_0"_a,
             "z_0"_a)
        .def(py::init<const CElectricDipoleMatrix&>())
        .def_property("origin",
                      &CElectricDipoleMatrix::getOriginCoordinates,
                      vlx_general::overload_cast_<const std::array<double, 3>&>()(&CElectricDipoleMatrix::setOriginCoordinates))
        .def("__str__", &CElectricDipoleMatrix::getString)
        .def(
            "to_numpy",
            [](const CElectricDipoleMatrix& obj) {
                return std::array<py::array_t<double>, 3>{
                    {matrix_to_numpy(obj, cartesians::X), matrix_to_numpy(obj, cartesians::Y), matrix_to_numpy(obj, cartesians::Z)}};
            },
            "Return a list of NumPy arrays holding the components of the electric dipole moment integrals")
        .def("x_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::X>)
        .def("y_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::Y>)
        .def("z_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::Z>)
        .def(py::self == py::self);

    // CElectricDipoleIntegralsDriver class

    auto py_dip_driver = bind_integrals_driver<CElectricDipoleIntegralsDriver>(m, "ElectricDipoleIntegralsDriver", "electric dipole");
    py_dip_driver.def_property("origin",
                               &CElectricDipoleIntegralsDriver::getElectricDipoleOrigin,
                               vlx_general::overload_cast_<const std::array<double, 3>&>()(&CElectricDipoleIntegralsDriver::setElectricDipoleOrigin));

    // CLinearMomentumMatrix class

    PyClass<CLinearMomentumMatrix>(m, "LinearMomentumMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&>(), "xMatrix"_a, "yMatrix"_a, "zMatrix"_a)
        .def(py::init<const std::array<CDenseMatrix, 3>&>(), "matrices"_a)
        .def(py::init<const CLinearMomentumMatrix&>())
        .def("__str__", &CLinearMomentumMatrix::getString)
        .def(
            "to_numpy",
            [](const CLinearMomentumMatrix& obj) {
                return std::array<py::array_t<double>, 3>{
                    {matrix_to_numpy(obj, cartesians::X), matrix_to_numpy(obj, cartesians::Y), matrix_to_numpy(obj, cartesians::Z)}};
            },
            "Return a list of NumPy arrays holding the components of the electric dipole moment integrals")
        .def("x_to_numpy", &matrix_to_numpy<CLinearMomentumMatrix, cartesians::X>)
        .def("y_to_numpy", &matrix_to_numpy<CLinearMomentumMatrix, cartesians::Y>)
        .def("z_to_numpy", &matrix_to_numpy<CLinearMomentumMatrix, cartesians::Z>)
        .def(py::self == py::self);

    // CLinearMomentumIntegralsDriver class

    auto py_lmom_driver = bind_integrals_driver<CLinearMomentumIntegralsDriver>(m, "LinearMomentumIntegralsDriver", "linear momentum");

    // CAngularMomentumMatrix class

    PyClass<CAngularMomentumMatrix>(m, "AngularMomentumMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&, const double, const double, const double>(),
             "xMatrix"_a,
             "yMatrix"_a,
             "zMatrix"_a,
             "x_0"_a,
             "y_0"_a,
             "z_0"_a)
        .def(py::init<const std::array<CDenseMatrix, 3>&, const std::array<double, 3>&>(), "matrices"_a, "origin"_a)
        .def(py::init<const CAngularMomentumMatrix&>())
        .def_property("origin",
                      &CAngularMomentumMatrix::getOriginCoordinates,
                      vlx_general::overload_cast_<const std::array<double, 3>&>()(&CAngularMomentumMatrix::setOriginCoordinates))
        .def("__str__", &CAngularMomentumMatrix::getString)
        .def(
            "to_numpy",
            [](const CAngularMomentumMatrix& obj) {
                return std::array<py::array_t<double>, 3>{
                    {matrix_to_numpy(obj, cartesians::X), matrix_to_numpy(obj, cartesians::Y), matrix_to_numpy(obj, cartesians::Z)}};
            },
            "Return a list of NumPy arrays holding the components of the electric dipole moment integrals")
        .def("x_to_numpy", &matrix_to_numpy<CAngularMomentumMatrix, cartesians::X>)
        .def("y_to_numpy", &matrix_to_numpy<CAngularMomentumMatrix, cartesians::Y>)
        .def("z_to_numpy", &matrix_to_numpy<CAngularMomentumMatrix, cartesians::Z>)
        .def(py::self == py::self);

    // CAngularMomentumIntegralsDriver class
    auto py_amom_driver = bind_integrals_driver<CAngularMomentumIntegralsDriver>(m, "AngularMomentumIntegralsDriver", "angular momentum");
    py_amom_driver.def_property(
        "origin",
        &CAngularMomentumIntegralsDriver::getAngularMomentumOrigin,
        vlx_general::overload_cast_<const std::array<double, 3>&>()(&CAngularMomentumIntegralsDriver::setAngularMomentumOrigin));

    // CElectricFieldMatrix class

    PyClass<CElectricFieldMatrix>(m, "ElectricFieldMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&>(), "xMatrix"_a, "yMatrix"_a, "zMatrix"_a)
        .def(py::init<const std::array<CDenseMatrix, 3>&>(), "matrices"_a)
        .def(py::init<const CElectricFieldMatrix&>())
        .def("__str__", &CElectricFieldMatrix::getString)
        .def(
            "to_numpy",
            [](const CElectricFieldMatrix& obj) {
                return std::array<py::array_t<double>, 3>{
                    {matrix_to_numpy(obj, cartesians::X), matrix_to_numpy(obj, cartesians::Y), matrix_to_numpy(obj, cartesians::Z)}};
            },
            "Return a list of NumPy arrays holding the components of the electric dipole moment integrals")
        .def("x_to_numpy", &matrix_to_numpy<CElectricFieldMatrix, cartesians::X>)
        .def("y_to_numpy", &matrix_to_numpy<CElectricFieldMatrix, cartesians::Y>)
        .def("z_to_numpy", &matrix_to_numpy<CElectricFieldMatrix, cartesians::Z>)
        .def(py::self == py::self);

    // CElectricFieldIntegralsDriver class

    PyClass<CElectricFieldIntegralsDriver>(m, "ElectricFieldIntegralsDriver")
        .def(py::init(&vlx_general::create<CElectricFieldIntegralsDriver>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const double, const double, const double>()(
                 &CElectricFieldIntegralsDriver::compute, py::const_),
             "Compute AO-basis electric field integrals for given molecule and basis, at given point",
             "molecule"_a,
             "basis"_a,
             "X"_a,
             "Y"_a,
             "Z"_a)
        .def(
            "compute",
            vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&, const double, const double, const double>()(
                &CElectricFieldIntegralsDriver::compute, py::const_),
            "Compute AO-basis electric field integrals for given molecule and bases, at given point",
            "molecule"_a,
            "bra_basis"_a,
            "ket_basis"_a,
            "X"_a,
            "Y"_a,
            "Z"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const double, const double, const double>()(
                 &CElectricFieldIntegralsDriver::compute, py::const_),
             "Compute AO-basis electric field integrals for two molecules in given basis, at given point",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a,
             "X"_a,
             "Y"_a,
             "Z"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&,
                                         const CMolecule&,
                                         const CMolecularBasis&,
                                         const CMolecularBasis&,
                                         const double,
                                         const double,
                                         const double>()(&CElectricFieldIntegralsDriver::compute, py::const_),
             "Compute AO-basis electric field integrals for two molecules, each with its own basis, at given point",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a,
             "X"_a,
             "Y"_a,
             "Z"_a)
        .def(
            "compute",
            [](const CElectricFieldIntegralsDriver&           obj,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const py::array_t<double, py::array::f_style>& dipoles,
               const py::array_t<double, py::array::f_style>& coordinates) {
                std::string errdims("ElectricFieldIntegralsDirver.compute: Inconsistent size of dipoles and/or their coordinates");

                errors::assertMsgCritical(coordinates.shape(0) == dipoles.shape(0), errdims);

                errors::assertMsgCritical(coordinates.shape(0) > 0, errdims);

                errors::assertMsgCritical(coordinates.shape(1) == 3, errdims);

                errors::assertMsgCritical(dipoles.shape(0) > 0, errdims);

                errors::assertMsgCritical(dipoles.shape(1) == 3, errdims);

                auto ds  = vlx_general::numpy_to_memblock2d(dipoles);
                auto xyz = vlx_general::numpy_to_memblock2d(coordinates);

                return obj.compute(molecule, basis, ds, xyz);
            },
            // the wonky format of the raw string literal is to get help(...) in Python to look nice
            R"pbdoc(
Compute electric field integrals for a collection of point dipoles.

:param molecule:
    The molecule.
:param basis:
    The basis set.
:param dipoles:
    Values of the point dipoles. This is a NumpPy array of shape (N,3), i.e. [..., [dx_i, dy_i, dz_i], ...]
:param coordinates:
    Coordinates of the point dipoles. This is a NumpPy array of shape (N,3), i.e. [..., [x_i, y_i, z_i], ...]
)pbdoc",
            "molecule"_a,
            "basis"_a,
            "dipoles"_a,
            "coordinates"_a)
        /* Overloads of compute accepting a 3-vector */
        .def(
            "compute",
            [](const CElectricFieldIntegralsDriver& obj, const CMolecule& mol, const CMolecularBasis& bas, const std::array<double, 3>& point) {
                return obj.compute(mol, bas, point[0], point[1], point[2]);
            },
            "Compute AO-basis electric field integrals for given molecule and basis, at given point",
            "molecule"_a,
            "basis"_a,
            "point"_a)
        .def(
            "compute",
            [](const CElectricFieldIntegralsDriver& obj,
               const CMolecule&                     mol,
               const CMolecularBasis&               bra,
               const CMolecularBasis&               ket,
               const std::array<double, 3>&         point) { return obj.compute(mol, bra, ket, point[0], point[1], point[2]); },
            "Compute AO-basis electric field integrals for given molecule and bases, at given point",
            "molecule"_a,
            "bra_basis"_a,
            "ket_basis"_a,
            "point"_a)
        .def(
            "compute",
            [](const CElectricFieldIntegralsDriver& obj,
               const CMolecule&                     bra_mol,
               const CMolecule&                     ket_mol,
               const CMolecularBasis&               bas,
               const std::array<double, 3>&         point) { return obj.compute(bra_mol, ket_mol, bas, point[0], point[1], point[2]); },
            "Compute AO-basis electric field integrals for two molecules in given basis, at given point",
            "bra_molecule"_a,
            "ket_molecule"_a,
            "basis"_a,
            "point"_a)
        .def(
            "compute",
            [](const CElectricFieldIntegralsDriver& obj,
               const CMolecule&                     bra_mol,
               const CMolecule&                     ket_mol,
               const CMolecularBasis&               bra,
               const CMolecularBasis&               ket,
               const std::array<double, 3>&         point) { return obj.compute(bra_mol, ket_mol, bra, ket, point[0], point[1], point[2]); },
            "Compute AO-basis electric field integrals for two molecules, each with its own basis, at given point",
            "bra_molecule"_a,
            "ket_molecule"_a,
            "bra_basis"_a,
            "ket_basis"_a,
            "point"_a);
}
}  // namespace vlx_oneints
