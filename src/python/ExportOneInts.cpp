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
// Helper function for CLinearMomentumIntegralsDriver constructor

static std::shared_ptr<CLinearMomentumIntegralsDriver>
CLinearMomentumIntegralsDriver_create(py::object py_comm)
{
    if (py_comm.is_none())
    {
        return std::make_shared<CLinearMomentumIntegralsDriver>(MPI_COMM_WORLD);
    }
    else
    {
        auto comm = vlx_general::get_mpi_comm(py_comm);

        return std::make_shared<CLinearMomentumIntegralsDriver>(comm);
    }
}

// Helper function for printing CLinearMomentumMatrix

static std::string
CLinearMomentumMatrix_str(const CLinearMomentumMatrix& self)
{
    return self.getStringForComponentX() + self.getStringForComponentY() + self.getStringForComponentZ();
}

// Helper function for converting CLinearMomentumMatrix to numpy array

static py::array_t<double>
CLinearMomentumMatrix_x_to_numpy(const CLinearMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.xvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CLinearMomentumMatrix_y_to_numpy(const CLinearMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.yvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CLinearMomentumMatrix_z_to_numpy(const CLinearMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.zvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

// Helper function for CAngularMomentumIntegralsDriver constructor

static std::shared_ptr<CAngularMomentumIntegralsDriver>
CAngularMomentumIntegralsDriver_create(py::object py_comm)
{
    if (py_comm.is_none())
    {
        return std::make_shared<CAngularMomentumIntegralsDriver>(MPI_COMM_WORLD);
    }
    else
    {
        auto comm = vlx_general::get_mpi_comm(py_comm);

        return std::make_shared<CAngularMomentumIntegralsDriver>(comm);
    }
}

// Helper function for printing CAngularMomentumMatrix

static std::string
CAngularMomentumMatrix_str(const CAngularMomentumMatrix& self)
{
    return self.getStringForComponentX() + self.getStringForComponentY() + self.getStringForComponentZ();
}

// Helper function for converting CAngularMomentumMatrix to numpy array

static py::array_t<double>
CAngularMomentumMatrix_x_to_numpy(const CAngularMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.xvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CAngularMomentumMatrix_y_to_numpy(const CAngularMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.yvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CAngularMomentumMatrix_z_to_numpy(const CAngularMomentumMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.zvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

// Helper function for CElectricFieldIntegralsDriver constructor

static std::shared_ptr<CElectricFieldIntegralsDriver>
CElectricFieldIntegralsDriver_create(py::object py_comm)
{
    if (py_comm.is_none())
    {
        return std::make_shared<CElectricFieldIntegralsDriver>(MPI_COMM_WORLD);
    }
    else
    {
        auto comm = vlx_general::get_mpi_comm(py_comm);

        return std::make_shared<CElectricFieldIntegralsDriver>(comm);
    }
}

// Helper function for printing CElectricDipoleMatrix

static std::string
CElectricFieldMatrix_str(const CElectricFieldMatrix& self)
{
    return self.getStringForComponentX() + self.getStringForComponentY() + self.getStringForComponentZ();
}

// Helper function for converting CElectricDipoleMatrix to numpy array

static py::array_t<double>
CElectricFieldMatrix_x_to_numpy(const CElectricFieldMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.xvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CElectricFieldMatrix_y_to_numpy(const CElectricFieldMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.yvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CElectricFieldMatrix_z_to_numpy(const CElectricFieldMatrix& self)
{
    return vlx_general::pointer_to_numpy(self.zvalues(), self.getNumberOfRows(), self.getNumberOfColumns());
}

CElectricFieldMatrix
CElectricFieldIntegralsDirver_compute(const CElectricFieldIntegralsDriver&           self,
                                      const CMolecule&                               molecule,
                                      const CMolecularBasis&                         basis,
                                      const py::array_t<double, py::array::f_style>& py_dipoles,
                                      const py::array_t<double, py::array::f_style>& py_coords)
{
    // NOTE: Dipoles data order
    // Dipole values:
    // namely np.array([[dx1, dy1, dz1], [dx2, dy2, dz2], [dx3, dy3, dz3], [dx4, dy4, dz4], ...])
    // Positions of dipoles:
    // namely np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], [d4, y4, z4], ...])

    // sanity check

    std::string errdims("ElectricFieldIntegralsDirver.compute: Inconsistent size of dipoles and/or their coordinates");

    errors::assertMsgCritical(py_coords.shape(0) == py_dipoles.shape(0), errdims);

    errors::assertMsgCritical(py_coords.shape(0) > 0, errdims);

    errors::assertMsgCritical(py_coords.shape(1) == 3, errdims);

    errors::assertMsgCritical(py_dipoles.shape(0) > 0, errdims);

    errors::assertMsgCritical(py_dipoles.shape(1) == 3, errdims);

    // form coordinate vector

    std::vector<double> coords(py_coords.size());

    std::vector<double> dipoles(py_dipoles.size());

    std::memcpy(coords.data(), py_coords.data(), py_coords.size() * sizeof(double));

    std::memcpy(dipoles.data(), py_dipoles.data(), py_dipoles.size() * sizeof(double));

    CMemBlock2D<double> dipdat(dipoles, static_cast<int32_t>(py_dipoles.shape(0)), 3);

    CMemBlock2D<double> crddat(coords, static_cast<int32_t>(py_dipoles.shape(0)), 3);

    return self.compute(molecule, basis, &dipdat, &crddat);
}

// Exports classes/functions in src/oneints to python

void
export_oneints(py::module& m)
{
    // COverlapMatrix class

    py::class_<COverlapMatrix, std::shared_ptr<COverlapMatrix>>(m, "OverlapMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init<const COverlapMatrix&>())
        .def(py::init(&matrix_from_numpy<COverlapMatrix>))
        .def("__str__", &COverlapMatrix::getString)
        .def("to_numpy", &matrix_to_numpy<COverlapMatrix>)
        .def("get_ortho_matrix", &COverlapMatrix::getOrthogonalizationMatrix)
        .def(py::self == py::self);

    // COverlapIntegralsDriver class

    py::class_<COverlapIntegralsDriver, std::shared_ptr<COverlapIntegralsDriver>>(m, "OverlapIntegralsDriver")
        .def(py::init(&vlx_general::create<COverlapIntegralsDriver>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&>()(&COverlapIntegralsDriver::compute, py::const_),
             "Compute AO-basis overlap integrals for given molecule and basis",
             "molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(&COverlapIntegralsDriver::compute,
                                                                                                             py::const_),
             "Compute mixed AO-basis overlap integrals for given molecule",
             "molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&>()(&COverlapIntegralsDriver::compute, py::const_),
             "Compute AO-basis overlap integrals for two molecules in given basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(
                 &COverlapIntegralsDriver::compute, py::const_),
             "Compute AO-basis overlap integrals for two molecules, each with its own basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a);

    // CKineticEnergyMatrix class

    py::class_<CKineticEnergyMatrix, std::shared_ptr<CKineticEnergyMatrix>>(m, "KineticEnergyMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init<const CKineticEnergyMatrix&>())
        .def(py::init(&matrix_from_numpy<CKineticEnergyMatrix>))
        .def("__str__", &CKineticEnergyMatrix::getString)
        .def("to_numpy", &matrix_to_numpy<CKineticEnergyMatrix>)
        .def("get_energy", &CKineticEnergyMatrix::getKineticEnergy)
        .def(py::self == py::self);

    // CKineticEnergyIntegralsDriver class

    py::class_<CKineticEnergyIntegralsDriver, std::shared_ptr<CKineticEnergyIntegralsDriver>>(m, "KineticEnergyIntegralsDriver")
        .def(py::init(&vlx_general::create<CKineticEnergyIntegralsDriver>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&>()(&CKineticEnergyIntegralsDriver::compute, py::const_),
             "Compute AO-basis kinetic energy integrals for given molecule and basis",
             "molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(&CKineticEnergyIntegralsDriver::compute,
                                                                                                             py::const_),
             "Compute mixed AO-basis kinetic energy integrals for given molecule",
             "molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&>()(&CKineticEnergyIntegralsDriver::compute,
                                                                                                       py::const_),
             "Compute AO-basis kinetic energy integrals for two molecules in given basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(
                 &CKineticEnergyIntegralsDriver::compute, py::const_),
             "Compute AO-basis kinetic energy integrals for two molecules, each with its own basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a);

    // CNuclearPotentialMatrix class

    py::class_<CNuclearPotentialMatrix, std::shared_ptr<CNuclearPotentialMatrix>>(m, "NuclearPotentialMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init<const CNuclearPotentialMatrix&>())
        .def(py::init(&matrix_from_numpy<CNuclearPotentialMatrix>))
        .def("__str__", &CNuclearPotentialMatrix::getString)
        .def("to_numpy", &matrix_to_numpy<CNuclearPotentialMatrix>)
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
        .def("get_energy", &CNuclearPotentialMatrix::getNuclearPotentialEnergy)
        .def(py::self == py::self);

    // CNuclearPotentialIntegralsDriver class

    py::class_<CNuclearPotentialIntegralsDriver, std::shared_ptr<CNuclearPotentialIntegralsDriver>>(m, "NuclearPotentialIntegralsDriver")
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

    // cartesians enum class

    py::enum_<cartesians>(m, "cartesians")
        .value("X", cartesians::X)
        .value("Y", cartesians::Y)
        .value("Z", cartesians::Z);

    // CElectricDipoleMatrix class

    py::class_<CElectricDipoleMatrix, std::shared_ptr<CElectricDipoleMatrix>>(m, "ElectricDipoleMatrix")
        .def(py::init<>())
        .def(py::init<const std::array<CDenseMatrix, 3>&, const std::array<double, 3>&>(), "matrices"_a, "origin"_a)
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&, const double, const double, const double>(), "xMatrix"_a, "yMatrix"_a, "zMatrix"_a, "x_0"_a, "y_0"_a, "z_0"_a)
        .def(py::init<const CElectricDipoleMatrix&>())
        .def_property("origin", &CElectricDipoleMatrix::getOriginCoordinates, &CElectricDipoleMatrix::setOriginCoordinates)
        .def("__str__", &CElectricDipoleMatrix::getString)
        .def("to_numpy", [](const CElectricDipoleMatrix& obj, cartesians cart) { return matrix_to_numpy(obj, cart); }, "component"_a)
        .def("to_numpy", [](const CElectricDipoleMatrix& obj, int32_t cart) { auto c = static_cast<cartesians>(cart); return matrix_to_numpy(obj, c); }, "component"_a)
        .def("x_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::X>)
        .def("y_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::Y>)
        .def("z_to_numpy", &matrix_to_numpy<CElectricDipoleMatrix, cartesians::Z>)
        .def(py::self == py::self);

    // CElectricDipoleIntegralsDriver class

    py::class_<CElectricDipoleIntegralsDriver, std::shared_ptr<CElectricDipoleIntegralsDriver>>(m, "ElectricDipoleIntegralsDriver")
        .def(py::init(&vlx_general::create<CElectricDipoleIntegralsDriver>), "comm"_a = py::none())
        .def_property("origin", &CElectricDipoleIntegralsDriver::getElectricDipoleOrigin, &CElectricDipoleIntegralsDriver::setElectricDipoleOrigin)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&>()(&CElectricDipoleIntegralsDriver::compute, py::const_),
             "Compute AO-basis electric dipole integrals for given molecule and basis",
             "molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(&CElectricDipoleIntegralsDriver::compute,
                                                                                                             py::const_),
             "Compute mixed AO-basis electric dipole integrals for given molecule",
             "molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&>()(&CElectricDipoleIntegralsDriver::compute,
                                                                                                       py::const_),
             "Compute AO-basis electric dipole integrals for two molecules in given basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(
                 &CElectricDipoleIntegralsDriver::compute, py::const_),
             "Compute AO-basis electric dipole integrals for two molecules, each with its own basis",
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a);

    // CLinearMomentumMatrix class

    py::class_<CLinearMomentumMatrix, std::shared_ptr<CLinearMomentumMatrix>>(m, "LinearMomentumMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&>())
        .def(py::init<const CLinearMomentumMatrix&>())
        .def("__str__", &CLinearMomentumMatrix_str)
        .def("x_to_numpy", &CLinearMomentumMatrix_x_to_numpy)
        .def("y_to_numpy", &CLinearMomentumMatrix_y_to_numpy)
        .def("z_to_numpy", &CLinearMomentumMatrix_z_to_numpy)
        .def(py::self == py::self);

    // CLinearMomentumIntegralsDriver class

    py::class_<CLinearMomentumIntegralsDriver, std::shared_ptr<CLinearMomentumIntegralsDriver>>(m, "LinearMomentumIntegralsDriver")
        .def(py::init(&CLinearMomentumIntegralsDriver_create), py::arg("py_comm") = py::none())
        .def("compute",
             (CLinearMomentumMatrix(CLinearMomentumIntegralsDriver::*)(const CMolecule&, const CMolecularBasis&) const) &
                 CLinearMomentumIntegralsDriver::compute)
        .def("compute",
             (CLinearMomentumMatrix(CLinearMomentumIntegralsDriver::*)(const CMolecule&, const CMolecularBasis&, const CMolecularBasis&) const) &
                 CLinearMomentumIntegralsDriver::compute)
        .def("compute",
             (CLinearMomentumMatrix(CLinearMomentumIntegralsDriver::*)(const CMolecule&, const CMolecule&, const CMolecularBasis&) const) &
                 CLinearMomentumIntegralsDriver::compute)
        .def("compute",
             (CLinearMomentumMatrix(CLinearMomentumIntegralsDriver::*)(
                 const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&) const) &
                 CLinearMomentumIntegralsDriver::compute);

    // CAngularMomentumMatrix class

    py::class_<CAngularMomentumMatrix, std::shared_ptr<CAngularMomentumMatrix>>(m, "AngularMomentumMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&, const double, const double, const double>())
        .def(py::init<const CAngularMomentumMatrix&>())
        .def("__str__", &CAngularMomentumMatrix_str)
        .def("x_to_numpy", &CAngularMomentumMatrix_x_to_numpy)
        .def("y_to_numpy", &CAngularMomentumMatrix_y_to_numpy)
        .def("z_to_numpy", &CAngularMomentumMatrix_z_to_numpy)
        .def(py::self == py::self);

    // CAngularMomentumIntegralsDriver class

    py::class_<CAngularMomentumIntegralsDriver, std::shared_ptr<CAngularMomentumIntegralsDriver>>(m, "AngularMomentumIntegralsDriver")
        .def(py::init(&CAngularMomentumIntegralsDriver_create), py::arg("py_comm") = py::none())
        .def("set_origin", &CAngularMomentumIntegralsDriver::setAngularMomentumOrigin)
        .def("compute",
             (CAngularMomentumMatrix(CAngularMomentumIntegralsDriver::*)(const CMolecule&, const CMolecularBasis&) const) &
                 CAngularMomentumIntegralsDriver::compute)
        .def("compute",
             (CAngularMomentumMatrix(CAngularMomentumIntegralsDriver::*)(const CMolecule&, const CMolecularBasis&, const CMolecularBasis&) const) &
                 CAngularMomentumIntegralsDriver::compute)
        .def("compute",
             (CAngularMomentumMatrix(CAngularMomentumIntegralsDriver::*)(const CMolecule&, const CMolecule&, const CMolecularBasis&) const) &
                 CAngularMomentumIntegralsDriver::compute)
        .def("compute",
             (CAngularMomentumMatrix(CAngularMomentumIntegralsDriver::*)(
                 const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&) const) &
                 CAngularMomentumIntegralsDriver::compute);

    // CElectricFieldMatrix class

    py::class_<CElectricFieldMatrix, std::shared_ptr<CElectricFieldMatrix>>(m, "ElectricFieldMatrix")
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&, const CDenseMatrix&, const CDenseMatrix&>())
        .def(py::init<const CElectricFieldMatrix&>())
        .def("__str__", &CElectricFieldMatrix_str)
        .def("x_to_numpy", &CElectricFieldMatrix_x_to_numpy)
        .def("y_to_numpy", &CElectricFieldMatrix_y_to_numpy)
        .def("z_to_numpy", &CElectricFieldMatrix_z_to_numpy)
        .def(py::self == py::self);

    // CElectricFieldIntegralsDriver class

    py::class_<CElectricFieldIntegralsDriver, std::shared_ptr<CElectricFieldIntegralsDriver>>(m, "ElectricFieldIntegralsDriver")
        .def(py::init(&CElectricFieldIntegralsDriver_create), py::arg("py_comm") = py::none())
        .def("compute",
             (CElectricFieldMatrix(CElectricFieldIntegralsDriver::*)(
                 const CMolecule&, const CMolecularBasis&, const double, const double, const double) const) &
                 CElectricFieldIntegralsDriver::compute)
        .def("compute",
             (CElectricFieldMatrix(CElectricFieldIntegralsDriver::*)(
                 const CMolecule&, const CMolecularBasis&, const CMolecularBasis&, const double, const double, const double) const) &
                 CElectricFieldIntegralsDriver::compute)
        .def("compute",
             (CElectricFieldMatrix(CElectricFieldIntegralsDriver::*)(
                 const CMolecule&, const CMolecule&, const CMolecularBasis&, const double, const double, const double) const) &
                 CElectricFieldIntegralsDriver::compute)
        .def(
            "compute",
            (CElectricFieldMatrix(CElectricFieldIntegralsDriver::*)(
                const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&, const double, const double, const double) const) &
                CElectricFieldIntegralsDriver::compute)
        .def("compute", &CElectricFieldIntegralsDirver_compute);
}

}  // namespace vlx_oneints
