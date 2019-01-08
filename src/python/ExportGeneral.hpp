//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <mpi.h>

namespace py = pybind11;

namespace bp_general { // bp_general namespace

/**
 Gets MPI_Comm pointer from a mpi4py communicator object.

 @param py_comm mpi4py communicator object.
 @return the MPI_Comm pointer.
 */
MPI_Comm* get_mpi_comm(py::object py_comm);

/**
 Gets numpy array from double pointer and dimension.

 @param ptr the double pointer.
 @param nRows number of rows.
 @param nColumns number of columns.
 @return numpy array.
 */
py::array_t<double> pointer_to_numpy(const double* ptr,
                                     const int32_t nRows,
                                     const int32_t nColumns);

/**
 Exports classes/functions in src/general to python.
 */
void export_general(py::module& m);

} // bp_general namespace

#endif /* ExportGeneral_hpp */
