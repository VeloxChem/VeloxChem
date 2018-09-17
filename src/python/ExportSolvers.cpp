//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "SADGuessDriver.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"

#include "ExportGeneral.hpp"
#include "ExportSolvers.hpp"

namespace bp = boost::python;

namespace bp_solvers { // bp_solvers namespace

// Helper function for creating a CSADGuessDriver object

static std::shared_ptr<CSADGuessDriver>
CSADGuessDriver_create(int32_t    globRank,
                       int32_t    globNodes,
                       bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CSADGuessDriver>(
        new CSADGuessDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading CSADGuessDriver::compute

CDenseMatrix
CSADGuessDriver_compute(
          CSADGuessDriver& self,
    const CMolecule&       molecule,
    const CMolecularBasis& basis_1,
    const CMolecularBasis& basis_2,
    const COverlapMatrix&  S12,
    const COverlapMatrix&  S22,
          COutputStream&   oStream,
          bp::object       py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis_1, basis_2, S12, S22, oStream, *comm_ptr);
}

// Exports classes/functions in src/solvers to python

void export_solvers()
{
    // CSADGuessDriver class

    bp::class_< CSADGuessDriver, std::shared_ptr<CSADGuessDriver> >
        (
            "SADGuessDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create", &CSADGuessDriver_create)
        .staticmethod("create")
        .def("compute", &CSADGuessDriver_compute)
    ;
}

} // bp_solvers namespace
