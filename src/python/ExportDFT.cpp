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

#include "ExportDFT.hpp"

#include <mpi.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <memory>

#include "DensityGridDriver.hpp"
#include "ExportGeneral.hpp"
#include "FunctionalParser.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"
#include "XCFunctional.hpp"
#include "XCIntegrator.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

// Exports classes/functions in src/dft to python

// Helper function for CAOKohnShamMatrix constructor

static std::shared_ptr<CAOKohnShamMatrix>
CAOKohnShamMatrix_from_dimensions(const int32_t nrows, const int32_t ncols, const bool is_rest)
{
    return std::make_shared<CAOKohnShamMatrix>(nrows, ncols, is_rest);
}

// Helper function for reduce_sum CAOKohnShamMatrix object

static void
CAOKohnShamMatrix_reduce_sum(CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.reduce_sum(rank, nodes, *comm);
}

// Helper function for collect CAOKohnShamMatrix object

static void
CAOKohnShamMatrix_collect(CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm, int32_t source)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.collect(rank, nodes, *comm, source);
}

// Helper function for getting grid coordinates and weigths as numpy array

static py::array_t<double>
CMolecularGrid_x_to_numpy(const CMolecularGrid& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesX(), self.getNumberOfGridPoints());
}

static py::array_t<double>
CMolecularGrid_y_to_numpy(const CMolecularGrid& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesY(), self.getNumberOfGridPoints());
}

static py::array_t<double>
CMolecularGrid_z_to_numpy(const CMolecularGrid& self)
{
    return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), self.getNumberOfGridPoints());
}

static py::array_t<double>
CMolecularGrid_w_to_numpy(const CMolecularGrid& self)
{
    return vlx_general::pointer_to_numpy(self.getWeights(), self.getNumberOfGridPoints());
}

// Helper function for distributing CMolecularGrid object

static void
CMolecularGrid_distribute(CMolecularGrid& self, int32_t rank, int32_t nodes, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.distribute(rank, nodes, *comm);
}

// Helper function for broadcasting CMolecularGrid object

static void
CMolecularGrid_broadcast(CMolecularGrid& self, int32_t rank, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm);
}

void
export_dft(py::module& m)
{
    // xcfun enum class

    // clang-format off
    py::enum_<xcfun>(m, "xcfun")
        .value("lda", xcfun::lda)
        .value("gga", xcfun::gga)
        .value("mgga", xcfun::mgga);
    // clang-format on

    // CAOKohnShamMatrix class

    PyClass<CAOKohnShamMatrix>(m, "AOKohnShamMatrix")
        .def(py::init<>())
        .def(py::init(&CAOKohnShamMatrix_from_dimensions))
        .def("__str__", &CAOKohnShamMatrix::getString)
        .def("get_matrix", &CAOKohnShamMatrix::getReferenceToKohnSham, "beta"_a = false)
        .def("reduce_sum", &CAOKohnShamMatrix_reduce_sum)
        .def("collect", &CAOKohnShamMatrix_collect)
        .def("get_electrons", &CAOKohnShamMatrix::getNumberOfElectrons)
        .def("get_energy", &CAOKohnShamMatrix::getExchangeCorrelationEnergy)
        .def(py::self == py::self);

    // CXCFunctional class

    PyClass<CXCFunctional>(m, "XCFunctional")
        .def(py::init<>())
        .def("get_frac_exact_exchange", &CXCFunctional::getFractionOfExactExchange)
        .def("get_func_type", &CXCFunctional::getFunctionalType)
        .def("get_func_label", &CXCFunctional::getLabel)
        .def("is_hybrid", &CXCFunctional::isHybridFunctional)
        .def("is_undefined", &CXCFunctional::isUndefined)
        .def(py::self == py::self);

    // CMolecularGrid class

    PyClass<CMolecularGrid>(m, "MolecularGrid")
        .def(py::init<>())
        .def(py::init<const CMolecularGrid&>())
        .def("number_of_points", &CMolecularGrid::getNumberOfGridPoints)
        .def("x_to_numpy", &CMolecularGrid_x_to_numpy)
        .def("y_to_numpy", &CMolecularGrid_y_to_numpy)
        .def("z_to_numpy", &CMolecularGrid_z_to_numpy)
        .def("w_to_numpy", &CMolecularGrid_w_to_numpy)
        .def("distribute", &CMolecularGrid_distribute)
        .def("broadcast", &CMolecularGrid_broadcast)
        .def(py::self == py::self);

    // CGridDriver class

    PyClass<CGridDriver>(m, "GridDriver")
        .def(py::init(&vlx_general::create<CGridDriver>), "comm"_a = py::none())
        .def("generate", &CGridDriver::generate)
        .def("set_level", &CGridDriver::setLevel);

    // CDensityGridDriver class

    PyClass<CDensityGridDriver>(m, "DensityGridDriver")
        .def(py::init(&vlx_general::create<CDensityGridDriver>), "comm"_a)
        .def("generate", &CDensityGridDriver::generate);

    // CXCIntegrator class

    PyClass<CXCIntegrator>(m, "XCIntegrator")
        .def(py::init(&vlx_general::create<CXCIntegrator>), "comm"_a = py::none())
        .def("integrate",
             vlx_general::
                 overload_cast_<const CAODensityMatrix&, const CMolecule&, const CMolecularBasis&, const CMolecularGrid&, const std::string&>()(
                     &CXCIntegrator::integrate, py::const_),
             "Integrate exchange-correlation functional contribution to zero order Kohn-Sham matrix.",
             "ao_density"_a,
             "molecule"_a,
             "ao_basis"_a,
             "grid"_a,
             "xcfun"_a)
        .def("integrate",
             vlx_general::overload_cast_<CAOFockMatrix&,
                                         const CAODensityMatrix&,
                                         const CAODensityMatrix&,
                                         const CMolecule&,
                                         const CMolecularBasis&,
                                         const CMolecularGrid&,
                                         const std::string&>()(&CXCIntegrator::integrate, py::const_),
             "Integrate exchange-correlation functional contribution to first order Fock matrices and adds it to AO Fock matrix.",
             "fock_matrix"_a,
             "rw_density"_a,
             "gs_density"_a,
             "molecule"_a,
             "ao_basis"_a,
             "grid"_a,
             "xcfun"_a);

    // exposing functions

    m.def("to_xcfun", &to_xcfun);

    m.def("parse_xc_func", &vxcfuncs::getExchangeCorrelationFunctional);
}

}  // namespace vlx_dft
