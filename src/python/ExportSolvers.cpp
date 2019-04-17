//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>

#include <mpi.h>

#include "SADGuessDriver.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "OverlapMatrix.hpp"
#include "ExportGeneral.hpp"
#include "ExportSolvers.hpp"

namespace py = pybind11;

namespace vlx_solvers { // vlx_solvers namespace

// Helper function for CSADGuessDriver constructor

static std::shared_ptr<CSADGuessDriver>
CSADGuessDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CSADGuessDriver>(
        new CSADGuessDriver(*comm_ptr)
        );
}

// Exports classes/functions in src/solvers to python

void export_solvers(py::module& m)
{
    // CSADGuessDriver class

    py::class_< CSADGuessDriver, std::shared_ptr<CSADGuessDriver> >
        (
            m, "SADGuessDriver"
        )
        .def(py::init(&CSADGuessDriver_create))
        .def("compute", &CSADGuessDriver::compute)
    ;
}

} // vlx_solvers namespace
