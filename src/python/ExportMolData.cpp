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

#include <memory>
#include <vector>
#include <string>

#include "Molecule.hpp"

#include "ExportGeneral.hpp"
#include "ExportMolData.hpp"

namespace bp = boost::python;

namespace bp_moldata { // bp_moldata namespace

// Helper function for broadcasting a CMolecule object

void
CMolecule_broadcast(CMolecule& self,
                    int32_t    rank,
                    bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}

// Exports classes/functions in src/moldata to python

void export_moldata()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // CMolecule class
    // Note: CMolecule has several constructors

    bp::class_< CMolecule, std::shared_ptr<CMolecule> >
        (
            "Molecule",
            bp::init<
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<std::string>&,
                const std::vector<int32_t>&
                >()
        )
        .def(bp::init<>())
        .def(bp::init<const CMolecule&>())
        .def(bp::init<const CMolecule&, const CMolecule&>())
        .def("print_geometry", &CMolecule::printGeometry)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
        .def("broadcast", &CMolecule_broadcast)
    ;
}

} // bp_moldata namespace
