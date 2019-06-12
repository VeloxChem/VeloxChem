//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "ExportGeneral.hpp"
#include "ExportSolvers.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OverlapMatrix.hpp"
#include "SADGuessDriver.hpp"

namespace py = pybind11;

namespace vlx_solvers {  // vlx_solvers namespace

// Helper function for CSADGuessDriver constructor

static std::shared_ptr<CSADGuessDriver>
CSADGuessDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CSADGuessDriver>(new CSADGuessDriver(*comm_ptr));
}

// Exports classes/functions in src/solvers to python

void
export_solvers(py::module& m)
{
    // CSADGuessDriver class

    py::class_<CSADGuessDriver, std::shared_ptr<CSADGuessDriver>>(m, "SADGuessDriver")
        .def(py::init(&CSADGuessDriver_create))
        .def("compute", &CSADGuessDriver::compute);
}

}  // namespace vlx_solvers
