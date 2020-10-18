//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ExportXTB.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <mpi.h>
#include <memory>

#include "XTBDriver.hpp"
#include "ExportGeneral.hpp"

namespace py = pybind11;

namespace vlx_xtb {  // vlx_xtb namespace

// Exports classes/functions in src/xtb to python

// Helper function for CXTBDriver constructor

static std::shared_ptr<CXTBDriver>
CXTBDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CXTBDriver>(new CXTBDriver(*comm_ptr));
}

void
export_xtb(py::module& m)
{
    // CXTBDriver class

    py::class_<CXTBDriver, std::shared_ptr<CXTBDriver>>(m, "XTBDriver")
        .def(py::init(&CXTBDriver_create))
	.def("is_master_node", &CXTBDriver::isMasterNode)
	.def("set_max_iter", &CXTBDriver::setMaxIterations)
	.def("set_elec_temp", &CXTBDriver::setElectronicTemp)
        .def("compute", &CXTBDriver::compute);
}

}  // namespace vlx_xtb
