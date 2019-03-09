//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <mpi.h>
#include <mpi4py/mpi4py.h>
#include <string>

#include "MpiFunc.hpp"
#include "ErrorHandler.hpp"
#include "Codata.hpp"
#include "StringFormat.hpp"
#include "SpinBlock.hpp"
#include "ExportGeneral.hpp"

namespace py = pybind11;

namespace vlx_general { // vlx_general namespace

// Gets MPI_Comm pointer from a mpi4py communicator object
// Not a static function; used in other files

MPI_Comm*
get_mpi_comm(py::object py_comm)
{
    PyObject* py_obj = py_comm.ptr();

    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);

    if (comm_ptr == NULL) throw py::error_already_set();

    return comm_ptr;
}

// Gets numpy array from double pointer and dimensions
// Not a static function; used in other files

py::array_t<double>
pointer_to_numpy(const double*               ptr,
                 const std::vector<int32_t>& dimension)
{
    std::vector<ssize_t> shape, strides;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        shape.push_back(static_cast<ssize_t>(dimension[i]));

        size_t strd = 1;

        for (size_t j = i + 1; j < dimension.size(); j++)
        {
            strd *= dimension[j];
        }

        strides.push_back(static_cast<ssize_t>(strd * sizeof(double)));
    }

    return py::array_t<double>(shape, strides, ptr);
}

py::array_t<double>
pointer_to_numpy(const double* ptr, int32_t nElements)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nElements}));
}

py::array_t<double>
pointer_to_numpy(const double* ptr, int32_t nRows, int32_t nColumns)
{
    return pointer_to_numpy(ptr, std::vector<int32_t>({nRows, nColumns}));
}

// Helper function for converting angular momentum

static std::string
string_to_angular_momentum(const int32_t angl)
{
    return fstr::to_AngularMomentum(angl);
}

static int32_t
integer_to_angular_momentum(const std::string& label)
{
    return fstr::to_AngularMomentum(label);
}

// Exports classes/functions in src/general to python

void export_general(py::module& m)
{
    // initialize mpi4py's C-API

    auto err = import_mpi4py();

    std::string errmpi4py("mpi4py: failed to import mpi4py");

    errors::assertMsgCritical(err == 0, errmpi4py);
    
    // szblock enum class
    
    py::enum_<szblock> (m, "szblock")
        .value("aa", szblock::aa)
        .value("ab", szblock::ab)
        .value("ba", szblock::ba)
        .value("bb", szblock::bb)
    ;

    // exposing functions

    m.def("mpi_master", &mpi::master);

    m.def("mpi_initialized", &mpi::initialized);

    m.def("assert_msg_critical", &errors::assertMsgCritical);

    m.def("bohr_in_angstroms", &units::getBohrValueInAngstroms);

    m.def("hartree_in_ev", &units::getHatreeValueInElectronVolts);

    m.def("to_angular_momentum", &string_to_angular_momentum);

    m.def("to_angular_momentum", &integer_to_angular_momentum);
}

} // vlx_general namespace
