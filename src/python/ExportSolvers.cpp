//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
CSADGuessDriver_create(int32_t    globRank,
                       int32_t    globNodes,
                       py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CSADGuessDriver>(
        new CSADGuessDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading CSADGuessDriver::compute

static CAODensityMatrix
CSADGuessDriver_compute(
          CSADGuessDriver& self,
    const CMolecule&       molecule,
    const CMolecularBasis& basis_1,
    const CMolecularBasis& basis_2,
    const COverlapMatrix&  S12,
    const COverlapMatrix&  S22,
          py::object       py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis_1, basis_2, S12, S22, *comm_ptr);
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
        .def("compute", &CSADGuessDriver_compute)
    ;
}

} // vlx_solvers namespace
