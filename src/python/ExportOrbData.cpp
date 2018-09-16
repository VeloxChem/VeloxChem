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

#include "MolecularBasis.hpp"

#include "ExportGeneral.hpp"
#include "ExportOrbData.hpp"

namespace bp = boost::python;

namespace bp_orbdata { // bp_orbdata namespace

// Helper function for broadcasting a CMolecularBasis object

void
CMolecularBasis_broadcast(CMolecularBasis& self,
                          int32_t          rank,
                          bp::object       py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Exports classes/functions in src/orbdata to python

void export_orbdata()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // CMolecularBasis class

    bp::class_< CMolecularBasis >
        (
            "CMolecularBasis",
            bp::init<>()
        )
        .def("get_label", &CMolecularBasis::getLabel)
        .def("broadcast", &bp_orbdata::CMolecularBasis_broadcast)
    ;
}

} // bp_orbdata namespace
