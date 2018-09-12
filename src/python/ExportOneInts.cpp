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
#include <string>

#include "MpiFunc.hpp"
#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "OverlapIntegralsDriver.hpp"

namespace bp = boost::python;

// ==> boost python helper function <==
// for creating COverlapIntegralsDriver object

static std::shared_ptr<COverlapIntegralsDriver>
COverlapIntegralsDriver_create(int32_t    globRank,
                               int32_t    globNodes,
                               bp::object py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return std::shared_ptr<COverlapIntegralsDriver>(
        new COverlapIntegralsDriver (globRank, globNodes, *comm_ptr)
        );
}

// ==> boost python helper function <==
// for overloading COverlapIntegralsDriver::compute

COverlapMatrix
COverlapIntegralsDriver_compute_1(
          COverlapIntegralsDriver& self,
    const CMolecule&               molecule,
    const CMolecularBasis&         basis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return self.compute(molecule, basis, oStream, *comm_ptr);
}

COverlapMatrix
COverlapIntegralsDriver_compute_2(
          COverlapIntegralsDriver& self,
    const CMolecule&               molecule,
    const CMolecularBasis&         braBasis,
    const CMolecularBasis&         ketBasis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return self.compute(molecule, braBasis, ketBasis, oStream, *comm_ptr);
}

COverlapMatrix
COverlapIntegralsDriver_compute_3(
          COverlapIntegralsDriver& self,
    const CMolecule&               braMolecule,
    const CMolecule&               ketMolecule,
    const CMolecularBasis&         basis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return self.compute(braMolecule, ketMolecule, basis, oStream, *comm_ptr);
}

COverlapMatrix
COverlapIntegralsDriver_compute_4(
          COverlapIntegralsDriver& self,
    const CMolecule&               braMolecule,
    const CMolecule&               ketMolecule,
    const CMolecularBasis&         braBasis,
    const CMolecularBasis&         ketBasis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis,
                        oStream, *comm_ptr);
}

// ==> boost python helper function <==
// for printing COverlapMatrix in python

std::string
COverlapMatrix_str (const COverlapMatrix& self)
{
    return self.getString();
}

// ==> boost python <==
// functions and classes

void export_oneints()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // COverlapMatrix class

    bp::class_< COverlapMatrix, std::shared_ptr<COverlapMatrix> >
        (
            "COverlapMatrix",
            bp::init<const CDenseMatrix&>()
        )
        .def(bp::init<>())
        .def("__str__", &COverlapMatrix_str)
        .def(bp::self == bp::other<COverlapMatrix>())
    ;

    // COverlapIntegralsDriver class

    bp::class_< COverlapIntegralsDriver, std::shared_ptr<COverlapIntegralsDriver> >
        (
            "COverlapIntegralsDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create",  &COverlapIntegralsDriver_create)
        .def("compute", &COverlapIntegralsDriver_compute_1)
        .def("compute", &COverlapIntegralsDriver_compute_2)
        .def("compute", &COverlapIntegralsDriver_compute_3)
        .def("compute", &COverlapIntegralsDriver_compute_4)
        .staticmethod("create")
    ;
}
