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
#include <memory>
#include <string>

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

// Helper function for COverlapIntegralsDriver constructor

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

static COverlapMatrix
COverlapIntegralsDriver_compute_1(
          COverlapIntegralsDriver& self,
    const CMolecule&               molecule,
    const CMolecularBasis&         basis,
          bp::object               py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, *comm_ptr);
}

static COverlapMatrix
COverlapIntegralsDriver_compute_2(
          COverlapIntegralsDriver& self,
    const CMolecule&               molecule,
    const CMolecularBasis&         braBasis,
    const CMolecularBasis&         ketBasis,
          bp::object               py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, braBasis, ketBasis, *comm_ptr);
}

static COverlapMatrix
COverlapIntegralsDriver_compute_3(
          COverlapIntegralsDriver& self,
    const CMolecule&               braMolecule,
    const CMolecule&               ketMolecule,
    const CMolecularBasis&         basis,
          bp::object               py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, basis, *comm_ptr);
}

static COverlapMatrix
COverlapIntegralsDriver_compute_4(
          COverlapIntegralsDriver& self,
    const CMolecule&               braMolecule,
    const CMolecule&               ketMolecule,
    const CMolecularBasis&         braBasis,
    const CMolecularBasis&         ketBasis,
          bp::object               py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis,
                        *comm_ptr);
}

// Helper function for printing COverlapMatrix

static std::string
COverlapMatrix_str (const COverlapMatrix& self)
{
    return self.getString();
}

// Helper function for converting COverlapMatrix to numpy array

static np::ndarray
COverlapMatrix_to_numpy(const COverlapMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for COverlapMatrix constructor

static std::shared_ptr<COverlapMatrix>
COverlapMatrix_from_numpy(const np::ndarray& arr)
{
    std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);

    return std::shared_ptr<COverlapMatrix>(new COverlapMatrix(*mp));
}

// Helper function for CKineticEnergyIntegralsDriver constructor

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

static CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_1(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     molecule,
    const CMolecularBasis&               basis,
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, *comm_ptr);
}

static CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_2(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     molecule,
    const CMolecularBasis&               braBasis,
    const CMolecularBasis&               ketBasis,
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, braBasis, ketBasis, *comm_ptr);
}

static CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_3(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     braMolecule,
    const CMolecule&                     ketMolecule,
    const CMolecularBasis&               basis,
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, basis, *comm_ptr);
}

static CKineticEnergyMatrix
CKineticEnergyIntegralsDriver_compute_4(
          CKineticEnergyIntegralsDriver& self,
    const CMolecule&                     braMolecule,
    const CMolecule&                     ketMolecule,
    const CMolecularBasis&               braBasis,
    const CMolecularBasis&               ketBasis,
          bp::object                     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis,
                        *comm_ptr);
}

// Helper function for printing CKineticEnergyMatrix

static std::string
CKineticEnergyMatrix_str (const CKineticEnergyMatrix& self)
{
    return self.getString();
}

// Helper function for converting CKineticEnergyMatrix to numpy array

static np::ndarray
CKineticEnergyMatrix_to_numpy(const CKineticEnergyMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for CKineticEnergyMatrix constructor

static std::shared_ptr<CKineticEnergyMatrix>
CKineticEnergyMatrix_from_numpy(const np::ndarray& arr)
{
    std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);

    return std::shared_ptr<CKineticEnergyMatrix>(new CKineticEnergyMatrix(*mp));
}

// Helper function for CKineticEnergyIntegralsDriver constructor

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

static CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_0(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  basis,
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, *comm_ptr);
}

static CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_1(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  basis,
    const CMolecule&                        pchgMolecule,
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, basis, pchgMolecule, *comm_ptr);
}

static CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_2(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        molecule,
    const CMolecularBasis&                  braBasis,
    const CMolecularBasis&                  ketBasis,
    const CMolecule&                        pchgMolecule,
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(molecule, braBasis, ketBasis, pchgMolecule,
                        *comm_ptr);
}

static CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_3(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        braMolecule,
    const CMolecule&                        ketMolecule,
    const CMolecularBasis&                  basis,
    const CMolecule&                        pchgMolecule,
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, basis, pchgMolecule,
                        *comm_ptr);
}

static CNuclearPotentialMatrix
CNuclearPotentialIntegralsDriver_compute_4(
          CNuclearPotentialIntegralsDriver& self,
    const CMolecule&                        braMolecule,
    const CMolecule&                        ketMolecule,
    const CMolecularBasis&                  braBasis,
    const CMolecularBasis&                  ketBasis,
    const CMolecule&                        pchgMolecule,
          bp::object                        py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis, pchgMolecule,
                        *comm_ptr);
}

// Helper function for printing CNuclearPotentialMatrix

static std::string
CNuclearPotentialMatrix_str (const CNuclearPotentialMatrix& self)
{
    return self.getString();
}

// Helper function for converting CNuclearPotentialMatrix to numpy array

static np::ndarray
CNuclearPotentialMatrix_to_numpy(const CNuclearPotentialMatrix& self)
{
    return bp_general::pointer_to_numpy(self.values(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

// Helper function for CNuclearPotentialMatrix constructor

static std::shared_ptr<CNuclearPotentialMatrix>
CNuclearPotentialMatrix_from_numpy(const np::ndarray& arr)
{
    std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);

    return std::shared_ptr<CNuclearPotentialMatrix>(new CNuclearPotentialMatrix(*mp));
}

// Exports classes/functions in src/oneints to python

void export_oneints()
{
    // initialize numpy

    Py_Initialize();

    np::initialize();

    // COverlapMatrix class

    bp::class_< COverlapMatrix,
                std::shared_ptr<COverlapMatrix> >
        (
            "OverlapMatrix",
             bp::init<>()
        )
        .def(bp::init<const CDenseMatrix&>())
        .def("__str__", &COverlapMatrix_str)
        .def("to_numpy", &COverlapMatrix_to_numpy)
        .def("from_numpy", &COverlapMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<COverlapMatrix>())
        .def("get_ortho_matrix", &COverlapMatrix::getOrthogonalizationMatrix)
    ;

    // COverlapIntegralsDriver class

    bp::class_< COverlapIntegralsDriver,
                std::shared_ptr<COverlapIntegralsDriver> >
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

    bp::class_< CKineticEnergyMatrix,
                std::shared_ptr<CKineticEnergyMatrix> >
        (
            "KineticEnergyMatrix",
            bp::init<>()
        )
        .def(bp::init<const CDenseMatrix&>())
        .def("__str__", &CKineticEnergyMatrix_str)
        .def("to_numpy", &CKineticEnergyMatrix_to_numpy)
        .def("from_numpy", &CKineticEnergyMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<CKineticEnergyMatrix>())
    ;

    // CKineticEnergyIntegralsDriver class

    bp::class_< CKineticEnergyIntegralsDriver,
                std::shared_ptr<CKineticEnergyIntegralsDriver> >
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

    bp::class_< CNuclearPotentialMatrix,
                std::shared_ptr<CNuclearPotentialMatrix> >
        (
            "NuclearPotentialMatrix",
            bp::init<>()
        )
        .def(bp::init<const CDenseMatrix&>())
        .def("__str__", &CNuclearPotentialMatrix_str)
        .def("to_numpy", &CNuclearPotentialMatrix_to_numpy)
        .def("from_numpy", &CNuclearPotentialMatrix_from_numpy)
        .staticmethod("from_numpy")
        .def(bp::self == bp::other<CNuclearPotentialMatrix>())
    ;

    // CNuclearPotentialIntegralsDriver class

    bp::class_< CNuclearPotentialIntegralsDriver,
                std::shared_ptr<CNuclearPotentialIntegralsDriver> >
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
