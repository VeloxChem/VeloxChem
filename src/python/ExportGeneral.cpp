//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <mpi.h>
#include <mpi4py/mpi4py.h>

#include "MpiFunc.hpp"

#include "ExportGeneral.hpp"

namespace bp = boost::python;

namespace bp_general { // bp_general namespace

// Gets MPI_Comm pointer from a boost python communicator object

MPI_Comm* get_mpi_comm(bp::object py_comm)
{
    PyObject* py_obj = py_comm.ptr();

    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);

    if (comm_ptr == NULL) bp::throw_error_already_set();

    return comm_ptr;
}

// Exports classes/functions in src/general to python

void export_general()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // exposing functions

    bp::def("mpi_master", &mpi::master);
}

} // bp_general namespace
