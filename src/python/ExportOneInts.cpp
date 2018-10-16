//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
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
#include "KineticEnergyMatrix.hpp"
#include "KineticEnergyIntegralsDriver.hpp"
#include "NuclearPotentialMatrix.hpp"
#include "NuclearPotentialIntegralsDriver.hpp"

#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ExportOneInts.hpp"

namespace bp = boost::python;

namespace np = boost::python::numpy;

namespace bp_oneints { // bp_oneints namespace

// Helper function for creating a COverlapIntegralsDriver object

static std::shared_ptr<COverlapIntegralsDriver>
COverlapIntegralsDriver_create(int32_t    globRank,
                               int32_t    globNodes,
                               bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return std::shared_ptr<COverlapIntegralsDriver>(
        new COverlapIntegralsDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading COverlapIntegralsDriver::compute

COverlapMatrix
COverlapIntegralsDriver_compute_1(
          COverlapIntegralsDriver& self,
    const CMolecule&               molecule,
    const CMolecularBasis&         basis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

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
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

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
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

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
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis,
                        oStream, *comm_ptr);
}

// Helper function for printing COverlapMatrix

std::string
COverlapMatrix_str (const COverlapMatrix& self)
{
    return self.getString();
}

// Helper function for converting COverlapMatrix to numpy array

np::ndarray
COverlapMatrix_to_numpy(const COverlapMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for converting numpy array to COverlapMatrix

COverlapMatrix
COverlapMatrix_from_numpy(const np::ndarray& arr)
{
    CDenseMatrix m = bp_math::CDenseMatrix_from_numpy(arr);

    return COverlapMatrix(m);
}

// Helper function for creating a CKineticEnergyIntegralsDriver object

static std::shared_ptr<CKineticEnergyIntegralsDriver>
CKineticEnergyIntegralsDriver_create(int32_t    globRank,
                                     int32_t    globNodes,
                                     bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CKineticEnergyIntegralsDriver>(
        new CKineticEnergyIntegralsDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading CKineticEnergyIntegralsDriver::compute

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_1(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     molecule,
    const CMolecularBasis&               basis,
          COutputStream&                 oStream, 
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, oStream, *comm_ptr);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_2(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     molecule,
    const CMolecularBasis&               braBasis,
    const CMolecularBasis&               ketBasis,
          COutputStream&                 oStream, 
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, braBasis, ketBasis, oStream, *comm_ptr);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_3(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     braMolecule,
    const CMolecule&                     ketMolecule,
    const CMolecularBasis&               basis,
          COutputStream&                 oStream, 
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, basis, oStream, *comm_ptr);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_4(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     braMolecule,
    const CMolecule&                     ketMolecule,
    const CMolecularBasis&               braBasis,
    const CMolecularBasis&               ketBasis,
          COutputStream&                 oStream, 
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis,
                        oStream, *comm_ptr);
}

// Helper function for printing CKineticEnergyMatrix

std::string
CKineticEnergyMatrix_str (const CKineticEnergyMatrix& self)
{
    return self.getString();
}

// Helper function for converting CKineticEnergyMatrix to numpy array

np::ndarray
CKineticEnergyMatrix_to_numpy(const CKineticEnergyMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for converting numpy array to CKineticEnergyMatrix

CKineticEnergyMatrix
CKineticEnergyMatrix_from_numpy(const np::ndarray& arr)
{
    CDenseMatrix m = bp_math::CDenseMatrix_from_numpy(arr);

    return CKineticEnergyMatrix(m);
}

// Helper function for creating a CKineticEnergyIntegralsDriver object

static std::shared_ptr<CNuclearPotentialIntegralsDriver>
CNuclearPotentialIntegralsDriver_create(int32_t    globRank,
                                        int32_t    globNodes,
                                        bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CNuclearPotentialIntegralsDriver>(
        new CNuclearPotentialIntegralsDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading CNuclearPotentialIntegralsDriver::compute

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_0(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  basis,
          COutputStream&                    oStream, 
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, oStream, *comm_ptr);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_1(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  basis,
    const CMolecule&                        pchgMolecule,
          COutputStream&                    oStream, 
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, pchgMolecule, oStream, *comm_ptr);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_2(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  braBasis,
    const CMolecularBasis&                  ketBasis,
    const CMolecule&                        pchgMolecule,
          COutputStream&                    oStream, 
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, braBasis, ketBasis, pchgMolecule,
                        oStream, *comm_ptr);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_3(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        braMolecule,
    const CMolecule&                        ketMolecule,
    const CMolecularBasis&                  basis,
    const CMolecule&                        pchgMolecule,
          COutputStream&                    oStream, 
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, basis, pchgMolecule,
                        oStream, *comm_ptr);
}

CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_4(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        braMolecule,
    const CMolecule&                        ketMolecule,
    const CMolecularBasis&                  braBasis,
    const CMolecularBasis&                  ketBasis,
    const CMolecule&                        pchgMolecule,
          COutputStream&                    oStream, 
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis, pchgMolecule,
                        oStream, *comm_ptr);
}

// Helper function for printing CNuclearPotentialMatrix

std::string
CNuclearPotentialMatrix_str (const CNuclearPotentialMatrix& self)
{
    return self.getString();
}

// Helper function for converting CNuclearPotentialMatrix to numpy array

np::ndarray
CNuclearPotentialMatrix_to_numpy(const CNuclearPotentialMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for converting numpy array to CNuclearPotentialMatrix

CNuclearPotentialMatrix
CNuclearPotentialMatrix_from_numpy(const np::ndarray& arr)
{
    CDenseMatrix m = bp_math::CDenseMatrix_from_numpy(arr);

    return CNuclearPotentialMatrix(m);
}

// Exports classes/functions in src/oneints to python

void export_oneints()
{
    // initialize mpi4py's C-API

    if (import_mpi4py() < 0) return;

    // initialize numpy

    Py_Initialize();

    np::initialize();

    // COverlapMatrix class

    bp::class_< COverlapMatrix, std::shared_ptr<COverlapMatrix> >
        (
            "OverlapMatrix",
            bp::init<const CDenseMatrix&>()
        )
        .def(bp::init<>())
        .def("__str__", &COverlapMatrix_str)
        .def("to_numpy", &COverlapMatrix_to_numpy)
        .def("from_numpy", &COverlapMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<COverlapMatrix>())
        .def("get_ortho_matrix", &COverlapMatrix::getOrthogonalizationMatrix)
    ;

    // COverlapIntegralsDriver class

    bp::class_< COverlapIntegralsDriver, std::shared_ptr<COverlapIntegralsDriver> >
        (
            "OverlapIntegralsDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create", &COverlapIntegralsDriver_create)
        .staticmethod("create")
        .def("compute", &COverlapIntegralsDriver_compute_1)
        .def("compute", &COverlapIntegralsDriver_compute_2)
        .def("compute", &COverlapIntegralsDriver_compute_3)
        .def("compute", &COverlapIntegralsDriver_compute_4)
    ;

    // CKineticEnergyMatrix class

    bp::class_< CKineticEnergyMatrix, std::shared_ptr<CKineticEnergyMatrix> >
        (
            "KineticEnergyMatrix",
            bp::init<const CDenseMatrix&>()
        )
        .def(bp::init<>())
        .def("__str__", &CKineticEnergyMatrix_str)
        .def("to_numpy", &CKineticEnergyMatrix_to_numpy)
        .def("from_numpy", &CKineticEnergyMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<CKineticEnergyMatrix>())
    ;

    // CKineticEnergyIntegralsDriver class

    bp::class_< CKineticEnergyIntegralsDriver, std::shared_ptr<CKineticEnergyIntegralsDriver> >
        (
            "KineticEnergyIntegralsDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create", &CKineticEnergyIntegralsDriver_create)
        .staticmethod("create")
        .def("compute", &CKineticEnergyIntegralsDriver_compute_1)
        .def("compute", &CKineticEnergyIntegralsDriver_compute_2)
        .def("compute", &CKineticEnergyIntegralsDriver_compute_3)
        .def("compute", &CKineticEnergyIntegralsDriver_compute_4)
    ;

    // CNuclearPotentialMatrix class

    bp::class_< CNuclearPotentialMatrix, std::shared_ptr<CNuclearPotentialMatrix> >
        (
            "NuclearPotentialMatrix",
            bp::init<const CDenseMatrix&>()
        )
        .def(bp::init<>())
        .def("__str__", &CNuclearPotentialMatrix_str)
        .def("to_numpy", &CNuclearPotentialMatrix_to_numpy)
        .def("from_numpy", &CNuclearPotentialMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<CNuclearPotentialMatrix>())
    ;

    // CNuclearPotentialIntegralsDriver class

    bp::class_< CNuclearPotentialIntegralsDriver, std::shared_ptr<CNuclearPotentialIntegralsDriver> >
        (
            "NuclearPotentialIntegralsDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create", &CNuclearPotentialIntegralsDriver_create)
        .staticmethod("create")
        .def("compute", &CNuclearPotentialIntegralsDriver_compute_0)
        .def("compute", &CNuclearPotentialIntegralsDriver_compute_1)
        .def("compute", &CNuclearPotentialIntegralsDriver_compute_2)
        .def("compute", &CNuclearPotentialIntegralsDriver_compute_3)
        .def("compute", &CNuclearPotentialIntegralsDriver_compute_4)
    ;
}

} // bp_oneints namespace
