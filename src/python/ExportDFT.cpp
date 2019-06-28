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

#include "XCFuncType.hpp"
#include "MolecularGrid.hpp"
#include "GridDriver.hpp"
#include "DensityGridDriver.hpp"
#include "ExportGeneral.hpp"


namespace py = pybind11;

namespace vlx_dft {  // vlx_dft namespace
    
// Exports classes/functions in src/dft to python
    
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

void
export_dft(py::module& m)
{
    // xcfun enum class
    
    py::enum_<xcfun>(m, "xcfun")
        .value("lda", xcfun::lda)
        .value("gga", xcfun::gga)
        .value("mgga", xcfun::mgga);
    
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

}
    
}  // namespace vlx_dft
