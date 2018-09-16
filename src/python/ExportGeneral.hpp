//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <boost/python.hpp>
#include <mpi.h>

namespace bp = boost::python;

namespace bp_general { // bp_general namespace

/**
 Gets MPI_Comm pointer from a boost python communicator object.

 @param py_comm boost python communicator object.
 @return the MPI_Comm pointer.
 */
MPI_Comm* get_mpi_comm(bp::object py_comm);

/**
 Exports classes/functions in src/general to python.
 */
void export_general();

} // bp_general namespace

#endif /* ExportGeneral_hpp */
