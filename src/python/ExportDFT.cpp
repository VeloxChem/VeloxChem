//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include "ExportDFT.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <mpi.h>
#include <memory>

#include "DensityGridDriver.hpp"
#include "ExportGeneral.hpp"
#include "GridDriver.hpp"
#include "MolecularGrid.hpp"
#include "XCFuncType.hpp"
#include "XCIntegrator.hpp"
#include "XCFunctional.hpp"
#include "FunctionalParser.hpp"

namespace py = pybind11;

namespace vlx_dft {  // vlx_dft namespace

// Exports classes/functions in src/dft to python

// Helper function for CAOKohnShamMatrix constructor

static std::shared_ptr<CAOKohnShamMatrix>
CAOKohnShamMatrix_from_dimensions(const int32_t nrows, const int32_t ncols, const bool is_rest)
{
    return std::shared_ptr<CAOKohnShamMatrix>(new CAOKohnShamMatrix(nrows, ncols, is_rest));
}
    
// Helper function for printing CAOKohnShamkMatrix

static std::string
CAOKohnShamMatrix_str(const CAOKohnShamMatrix& self)
{
    return self.getString();
}
    
// Helper function for reduce_sum CAOKohnShamMatrix object

static void
CAOKohnShamMatrix_reduce_sum(CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.reduce_sum(rank, nodes, *comm_ptr);
}

// Helper function for collect CAOKohnShamMatrix object

static void
CAOKohnShamMatrix_collect(CAOKohnShamMatrix& self, int32_t rank, int32_t nodes, py::object py_comm, int32_t source)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.collect(rank, nodes, *comm_ptr, source);
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
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.distribute(rank, nodes, *comm_ptr);
}

// Helper function for broadcasting CMolecularGrid object

static void
CMolecularGrid_broadcast(CMolecularGrid& self, int32_t rank, py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Helper function for CGridDriver constructor

static std::shared_ptr<CGridDriver>
CGridDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CGridDriver>(new CGridDriver(*comm_ptr));
}

// Helper function for CDensityGridDriver constructor

static std::shared_ptr<CDensityGridDriver>
CDensityGridDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CDensityGridDriver>(new CDensityGridDriver(*comm_ptr));
}
    
// Helper function for CXCIntegrator constructor

static std::shared_ptr<CXCIntegrator>
CXCIntegrator_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CXCIntegrator>(new CXCIntegrator(*comm_ptr));
}

void
export_dft(py::module& m)
{
    // xcfun enum class

    py::enum_<xcfun>(m, "xcfun").value("lda", xcfun::lda).value("gga", xcfun::gga).value("mgga", xcfun::mgga);
    
    // CAOKohnShamMatrix class

    py::class_<CAOKohnShamMatrix, std::shared_ptr<CAOKohnShamMatrix>>(m, "AOKohnShamMatrix")
        .def(py::init<>())
        .def(py::init(&CAOKohnShamMatrix_from_dimensions))
        .def("__str__", &CAOKohnShamMatrix_str)
        .def("get_matrix", &CAOKohnShamMatrix::getReferenceToKohnSham, py::arg("beta") = false)
        .def("reduce_sum", &CAOKohnShamMatrix_reduce_sum)
        .def("collect", &CAOKohnShamMatrix_collect)
        .def("get_electrons", &CAOKohnShamMatrix::getNumberOfElectrons)
        .def("get_energy", &CAOKohnShamMatrix::getExchangeCorrelationEnergy)
        .def(py::self == py::self);

    // CXCFunctional class
    
     py::class_<CXCFunctional, std::shared_ptr<CXCFunctional>>(m, "XCFunctional")
        .def(py::init<>())
        .def("get_frac_exact_exchange", &CXCFunctional::getFractionOfExactExchange)
        .def("get_func_type", &CXCFunctional::getFunctionalType)
        .def("get_func_label", &CXCFunctional::getLabel)
        .def("is_hybrid", &CXCFunctional::isHybridFunctional)
        .def("is_undefined", &CXCFunctional::isUndefined)
        .def(py::self == py::self);
    
    // CMolecularGrid class

    py::class_<CMolecularGrid, std::shared_ptr<CMolecularGrid>>(m, "MolecularGrid")
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

    py::class_<CGridDriver, std::shared_ptr<CGridDriver>>(m, "GridDriver")
        .def(py::init(&CGridDriver_create))
        .def("generate", &CGridDriver::generate)
        .def("set_level", &CGridDriver::setLevel);

    // CDensityGridDriver class

    py::class_<CDensityGridDriver, std::shared_ptr<CDensityGridDriver>>(m, "DensityGridDriver")
        .def(py::init(&CDensityGridDriver_create))
        .def("generate", &CDensityGridDriver::generate);

    // CXCIntegrator class

    py::class_<CXCIntegrator, std::shared_ptr<CXCIntegrator>>(m, "XCIntegrator")
        .def(py::init(&CXCIntegrator_create))
        .def("integrate", (CAOKohnShamMatrix (CXCIntegrator::*)(const CAODensityMatrix&,
                                                                const COverlapMatrix& 
                                                                const CMolecule&,
                                                                const CMolecularBasis&,
                                                                const CMolecularGrid&,
                                                                const std::string&) const) & CXCIntegrator::integrate)
        .def("integrate", (void (CXCIntegrator::*)(      CAOFockMatrix&,
                                                   const CAODensityMatrix&,
                                                   const CAODensityMatrix&,
                                                   const CMolecule&,
                                                   const CMolecularBasis&,
                                                   const CMolecularGrid&,
                                                   const std::string&) const) &  CXCIntegrator::integrate);
    
    // exposing functions

    m.def("to_xcfun", &to_xcfun);
    
    m.def("parse_xc_func", &vxcfuncs::getExchangeCorrelationFunctional);
}

}  // namespace vlx_dft
