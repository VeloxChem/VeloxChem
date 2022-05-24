//
//                           VELOXCHEM 1.0-RC2
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

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "DensityGrid.hpp"
#include "DensityGridDriver.hpp"
#include "ExportGeneral.hpp"
#include "FunctionalParser.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"
#include "XCFunctional.hpp"
#include "XCIntegrator.hpp"
#include "XCMolecularGradient.hpp"
#include "DenseMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_dft {  // vlx_dft namespace

// Exports classes/functions in src/dft to python

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
        .def(py::init<int32_t, int32_t, bool>(), "nrows"_a, "ncols"_a, "is_rest"_a)
        .def("__str__", &CAOKohnShamMatrix::getString)
        .def("get_matrix", &CAOKohnShamMatrix::getReferenceToKohnSham, "Gets constant reference to specific Kohn-Sham matrix.", "beta"_a = false)
        .def(
            "reduce_sum",
            [](CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.reduce_sum(rank, nodes, *comm);
            },
            "Performs reduce_sum for AOKohnShamMatrix object",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a)
        .def(
            "collect",
            [](CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm, int32_t source) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.collect(rank, nodes, *comm, source);
            },
            "Collects AOKohnShamMatrix object.",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a,
            "source"_a)
        .def("get_electrons", &CAOKohnShamMatrix::getNumberOfElectrons, "Gets number of electrons obtained by integrating Kohn-Sham matrix.")
        .def("get_energy", &CAOKohnShamMatrix::getExchangeCorrelationEnergy, "Gets exchange-correlation energy associated with Kohn-Sham matrix.")
        .def(py::self == py::self);

    // CXCFunctional class

    PyClass<CXCFunctional>(m, "XCFunctional")
        .def(py::init<>())
        .def("get_frac_exact_exchange",
             &CXCFunctional::getFractionOfExactExchange,
             "Gets fraction of exact Hatree-Fock exchange in exchange-correlation functional.")
        .def("get_func_type", &CXCFunctional::getFunctionalType, "Gets type of exchange-correlation functional.")
        .def("get_func_label", &CXCFunctional::getLabel, "Gets label of exchange-correlation functional.")
        .def("is_hybrid",
             &CXCFunctional::isHybridFunctional,
             "Determines if exchange-correlation functional is of hybrid type i.e. non-zero fraction of exact Hatree-Fock exchange.")
        .def("is_undefined", &CXCFunctional::isUndefined, "Determines if exchange-correlation function is undefined.")
        .def(py::self == py::self);

    // CMolecularGrid class

    PyClass<CMolecularGrid>(m, "MolecularGrid")
        .def(py::init<>())
        .def(py::init<const CMolecularGrid&>())
        .def("number_of_points", &CMolecularGrid::getNumberOfGridPoints)
        .def(
            "x_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesX(), self.getNumberOfGridPoints());
            },
            "Gets X coordinates of grid as numpy array.")
        .def(
            "y_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesY(), self.getNumberOfGridPoints());
            },
            "Gets Y coordinates of grid as numpy array.")
        .def(
            "z_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoordinatesZ(), self.getNumberOfGridPoints());
            },
            "Gets Z coordinates of grid as numpy array.")
        .def(
            "w_to_numpy",
            [](const CMolecularGrid& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getWeights(), self.getNumberOfGridPoints());
            },
            "Gets weights of grid as numpy array.")
        .def(
            "distribute",
            [](CMolecularGrid& self, int32_t rank, int32_t nodes, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.distribute(rank, nodes, *comm);
            },
            "Distributes MolecularGrid object.",
            "rank"_a,
            "nodes"_a,
            "py_comm"_a)
        .def(
            "broadcast",
            [](CMolecularGrid& self, int32_t rank, py::object py_comm) -> void {
                auto comm = vlx_general::get_mpi_comm(py_comm);
                self.broadcast(rank, *comm);
            },
            "Broadcasts MolecularGrid object.",
            "rank"_a,
            "py_comm"_a)
        .def(py::self == py::self);

    // CGridDriver class

    PyClass<CGridDriver>(m, "GridDriver")
        .def(py::init(&vlx_general::create<CGridDriver>), "comm"_a = py::none())
        .def("generate",
             &CGridDriver::generate,
             "Generates molecular grid for molecule. Errors are printed to output stream. Grid generation is distributed within domain of MPI "
             "communicator.",
             "molecule"_a)
        .def("set_level",
             &CGridDriver::setLevel,
             "Sets accuracy level for grid generation. Level: 1-6, where 1 is coarse grid, 5 is ultrafine grid, 6 special benchmarking grid.",
             "gridLevel"_a);

    // CDensityGridDriver class

    PyClass<CDensityGridDriver>(m, "DensityGridDriver")
        .def(py::init(&vlx_general::create<CDensityGridDriver>), "comm"_a = py::none())
        .def(
            "generate",
            [](CDensityGridDriver&     self,
               const CAODensityMatrix& aoDensityMatrix,
               const CMolecule&        molecule,
               const CMolecularBasis&  basis,
               const CMolecularGrid&   molecularGrid,
               const std::string&      xcFunctionalType) -> CDensityGrid {
                return self.generate(aoDensityMatrix, molecule, basis, molecularGrid, to_xcfun(xcFunctionalType));
            },
            "Generates partitioned density grid for given molecule and type of exchange-correlation functional. Density grid generation is "
            "distributed within domain of MPI communicator.",
            "aoDensityMatrix"_a,
            "molecule"_a,
            "basis"_a,
            "molecularGrid"_a,
            "xcFunctionalType"_a);

    // CXCIntegrator class

    PyClass<CXCIntegrator>(m, "XCIntegrator")
        .def(py::init(&vlx_general::create<CXCIntegrator>), "comm"_a = py::none())
        .def("integrate",
             py::overload_cast<const CAODensityMatrix&, const CMolecule&, const CMolecularBasis&, const CMolecularGrid&, const std::string&>(
                 &CXCIntegrator::integrate, py::const_),
             "Integrate exchange-correlation functional contribution to zero order Kohn-Sham matrix.",
             "aoDensityMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate",
             py::overload_cast<CAOFockMatrix&,
                               const CAODensityMatrix&,
                               const CAODensityMatrix&,
                               const CMolecule&,
                               const CMolecularBasis&,
                               const CMolecularGrid&,
                               const std::string&>(&CXCIntegrator::integrate, py::const_),
             "Integrate exchange-correlation functional contribution to first order Fock matrices and adds it to AO Fock matrix.",
             "aoFockMatrix"_a,
             "rwDensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a)
        .def("integrate",
             py::overload_cast<CAOFockMatrix&,
                               const CAODensityMatrix&,
                               const CAODensityMatrix&,
                               const CAODensityMatrix&,
                               const CMolecule&,
                               const CMolecularBasis&,
                               const CMolecularGrid&,
                               const std::string&,
                               const std::string&>(&CXCIntegrator::integrate, py::const_),
             "Integrate exchange-correlation functional contribution to first order Fock matrices and adds it to AO Fock matrix.",
             "aoFockMatrix"_a,
             "rwDensityMatrix"_a,
             "rw12DensityMatrix"_a,
             "gsDensityMatrix"_a,
             "molecule"_a,
             "basis"_a,
             "molecularGrid"_a,
             "xcFuncLabel"_a,
             "quadMode"_a);
    
    // CXCMolecularGradient class

    PyClass<CXCMolecularGradient>(m, "XCMolecularGradient")
        .def(py::init(&vlx_general::create<CXCMolecularGradient>), "comm"_a = py::none())
        .def(
            "integrate",
            [](      CXCMolecularGradient& self,
               const std::vector<int32_t>& idsAtomic,
               const CAODensityMatrix&     aoDensityMatrix,
               const CMolecule&            molecule,
               const CMolecularBasis&      basis,
               const CMolecularGrid&       molecularGrid,
               const std::string&          xcFuncLabel) -> py::array_t<double> {
                    auto molgrad = self.integrate(idsAtomic, aoDensityMatrix,
                                                  molecule, basis,
                                                  molecularGrid, xcFuncLabel);
                    return vlx_general::pointer_to_numpy(molgrad.values(),
                                                         molgrad.getNumberOfRows(),
                                                         molgrad.getNumberOfColumns());
             },
            "Integrates exchange-correlation contribution to moleculargradient.",
            "idsAtomic"_a,
            "aoDensityMatrix"_a,
            "molecule"_a,
            "basis"_a, 
            "molecularGrid"_a,
            "xcFuncLabel"_a);

    // CDensityGrid class

    PyClass<CDensityGrid>(m, "DensityGrid")
        .def(py::init<>())
        .def(py::init<const CDensityGrid&>())
        .def("number_of_points", &CDensityGrid::getNumberOfGridPoints)
        .def("number_of_density_matrices", &CDensityGrid::getNumberOfDensityMatrices)
        .def(
            "density_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensity(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "density_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensity(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_norm_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient norm on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_norm_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient norm on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_product_ab_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.mixedDensityGradient(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets mixed density gradient product on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_x_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientX(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient X component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_y_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientY(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient Y component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_z_aa_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.alphaDensityGradientZ(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets alpha density gradient Z component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_x_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientX(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient X component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_y_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientY(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient Y component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(
            "gradient_z_bb_to_numpy",
            [](const CDensityGrid& self, int32_t iDensityMatrix) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.betaDensityGradientZ(iDensityMatrix), self.getNumberOfGridPoints());
            },
            "Gets beta density gradient Z component on grid as numpy array.",
            "iDensityMatrix"_a)
        .def(py::self == py::self);

    // exposing functions

    m.def("to_xcfun", &to_xcfun, "Converts string label to its enumerate class value.", "label"_a);

    m.def("available_functionals", &vxcfuncs::getAvailableFunctionals, "Gets a list of available exchange-correlation functionals.");

    m.def("parse_xc_func",
          &vxcfuncs::getExchangeCorrelationFunctional,
          "Converts exchange-correlation functional label to exchange-correlation functional object.",
          "xcLabel"_a);
}

}  // namespace vlx_dft
