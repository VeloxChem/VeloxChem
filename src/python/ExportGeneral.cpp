//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/tuple.hpp>
#include <mpi.h>
#include <mpi4py/mpi4py.h>

#include "MpiFunc.hpp"
#include "ErrorHandler.hpp"

#include "ExportGeneral.hpp"

namespace bp = boost::python;

namespace np = boost::python::numpy;

namespace bp_general { // bp_general namespace

// Gets MPI_Comm pointer from a boost python communicator object

MPI_Comm*
get_mpi_comm(bp::object py_comm)
{
    PyObject* py_obj = py_comm.ptr();

    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);

    if (comm_ptr == NULL) bp::throw_error_already_set();

    return comm_ptr;
}

// Gets numpy array from double pointer and dimensions

np::ndarray
pointer_to_numpy(const double* ptr, int32_t nRows, int32_t nColumns)
{
    np::dtype dt_double = np::dtype::get_builtin<double>();

    if (ptr == nullptr || nRows == 0 || nColumns == 0)
    {
        bp::tuple shape = bp::make_tuple(0, 0);

        return np::empty(shape, dt_double);
    }

    bp::tuple shape = bp::make_tuple(nRows, nColumns);

    bp::tuple stride = bp::make_tuple(sizeof(double) * nColumns,
                                      sizeof(double) * 1);

    return np::from_data(ptr, dt_double, shape, stride, bp::object());
}

// Exports classes/functions in src/general to python

void export_general()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // initialize numpy

    Py_Initialize();

    np::initialize();

    // exposing functions

    bp::def("mpi_master", &mpi::master);

    bp::def("mpi_initialized", &mpi::initialized);

    bp::def("assert_msg_critical", &errors::assertMsgCritical);
}

} // bp_general namespace
