//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include <mpi.h>
#include <mpi4py/mpi4py.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "MpiFunc.hpp"
#include "AppManager.hpp"
#include "InputStream.hpp"
#include "OutputStream.hpp"
#include "MolXYZReader.hpp"
#include "EnvironmentReader.hpp"
#include "BasisReader.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "OverlapIntegralsDriver.hpp"

namespace bp = boost::python;

namespace PyVLX { // PyVLX namespace

    int
    vlx_run(int    argc,
            char** argv)
    {
        // initialize global MPI communicator

        if (!mpi::init(argc, argv)) return EXIT_FAILURE;

        // initialize application manager

        CAppManager appManager(argc, argv);

        if (!appManager.getState())
        {
            mpi::abort(MPI_ERR_OTHER, "main()");
        
            return EXIT_FAILURE;
        }

        // execute application manager

        appManager.execute();

        if (!appManager.getState())
        {
            mpi::abort(MPI_ERR_OTHER, "main()");

            return EXIT_FAILURE;
        }

        // finalize global MPI communicator 

        if (!mpi::finalize()) return EXIT_FAILURE;

        return EXIT_SUCCESS;
    }

    int
    py_vlx_run(bp::list inputs)
    {
        // prepare argc and argv from boost python list

        int argc = len(inputs) + 1;

        if (argc < 3) return EXIT_FAILURE;

        char** argv;

        argv = (char**)malloc(sizeof(char*) * argc);

        argv[0] = (char*)malloc(sizeof(char) * (strlen("exe")+1));

        memset(argv[0], '\0', sizeof(char) * (strlen("exe")+1));

        strcpy(argv[0], "exe");

        for (int iarg = 1; iarg < argc; iarg++)
        {
            std::string line = bp::extract<std::string>(inputs[iarg - 1]);

            const char* text = line.c_str();

            argv[iarg] = (char*)malloc(sizeof(char) * (strlen(text) + 1));

            memset(argv[iarg], '\0', sizeof(char) * (strlen(text) + 1));

            memcpy(argv[iarg], text, sizeof(char) * strlen(text));
        }

        // run VeloxChem

        return vlx_run(argc, argv);
    }
}

// ==> boost python helper function <==
// for creating CAppManager object

static std::shared_ptr<CAppManager>
CAppManager_create(std::string input_string, std::string output_string)
{
    int argc = 3;

    char* argv[argc];

    std::vector<std::string> inputs ({std::string("exe"), input_string, output_string});

    for (int i = 0; i < argc; i++)
    {
        const char* text = inputs[i].c_str();

        argv[i] = (char*)malloc(sizeof(char) * (strlen(text) + 1));

        memset(argv[i], '\0', sizeof(char) * (strlen(text) + 1));

        memcpy(argv[i], text, sizeof(char) * strlen(text));
    }

    mpi::init(argc, argv);

    return std::shared_ptr<CAppManager>(new CAppManager (argc, argv));
}

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
    const CMolecularBasis&         braBasis,
    const CMolecularBasis&         ketBasis,
          COutputStream&           oStream,
          bp::object               py_comm)
{
    PyObject* py_obj = py_comm.ptr();
    MPI_Comm* comm_ptr = PyMPIComm_Get(py_obj);
    if (comm_ptr == NULL) bp::throw_error_already_set();

    return self.compute(braMolecule, ketMolecule, braBasis, ketBasis, oStream, *comm_ptr);
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

BOOST_PYTHON_MODULE(VeloxChemMP)
{
    // initialize mpi4py's C-API
    if (import_mpi4py() < 0) return;

    // -----------------
    // ==> Functions <==
    // -----------------

    // run VeloxChem
    bp::def("run", PyVLX::py_vlx_run);

    // -----------------------
    // ==> Manager classes <==
    // -----------------------

    // CAppManager class
    // Note: "create" is a static method that returns a CAppManager object

    bp::class_< CAppManager, std::shared_ptr<CAppManager> >
        (
            "CAppManager",
            bp::init<int, char**>()
        )
        .def("create",    &CAppManager_create)
        .def("execute",   &CAppManager::execute)
        .def("get_state", &CAppManager::getState)
        .staticmethod("create")
    ;

    // ----------------------
    // ==> Stream classes <==
    // ----------------------

    // COutputStream class

    bp::class_< COutputStream, std::shared_ptr<COutputStream> >
        (
            "COutputStream",
            bp::init<const std::string&>()
        )
        .def("get_state", &COutputStream::getState)
        .def("flush",     &COutputStream::flush)
    ;

    // CInputStream class

    bp::class_< CInputStream, std::shared_ptr<CInputStream> >
        (
            "CInputStream",
            bp::init<const std::string&, COutputStream&>()
        )
        .def("read",      &CInputStream::read)
        .def("get_state", &CInputStream::getState)
    ;

    // CInputData class

    bp::class_< CInputData, std::shared_ptr<CInputData> >
        (
            "CInputData", bp::init<>()
        )
    ;

    // ----------------------
    // ==> Reader classes <==
    // ----------------------

    // CMolXYZReader class
    // Note: Need member function pointers for proper overloading

    void (CMolXYZReader::*parse_1)(      CMolecule&     molecule,
                                   const CInputData&    inputData,
                                         COutputStream& oStream)
        = &CMolXYZReader::parse;

    void (CMolXYZReader::*parse_2)(      CMolecule&     molecule,
                                   const CInputData&    inputData,
                                   const int32_t        iGroup,
                                         COutputStream& oStream)
        = &CMolXYZReader::parse;

    bp::class_< CMolXYZReader >
        (
            "CMolXYZReader", bp::init<>()
        )
        .def("parse", parse_1)
        .def("parse", parse_2)
        .def("get_state", &CMolXYZReader::getState)
    ;

    // CEnvironmentReader class

    bp::class_< CEnvironmentReader >
        (
            "CEnvironmentReader", bp::init<>()
        )
        .def("parse", &CEnvironmentReader::parse)
        .def("get_state", &CEnvironmentReader::getState)
        .def("get_path_to_basis_sets", &CEnvironmentReader::getPathToBasisSets)
    ;

    // CBasisReader class

    bp::class_< CBasisReader >
        (
            "CBasisReader", bp::init<>()
        )
        .def("parse", &CBasisReader::parse)
        .def("get_state", &CBasisReader::getState)
        .def("get_ao_basis", &CBasisReader::getAOBasis)
        .def("get_rij_basis", &CBasisReader::getRIJBasis)
        .def("get_min_basis", &CBasisReader::getMinBasis)
    ;

    // ------------------------------
    // ==> Molecule/Basis classes <==
    // ------------------------------

    // CMolecule class
    // Note: CMolecule has two constructors

    bp::class_< CMolecule, std::shared_ptr<CMolecule> >
        (
            "CMolecule", bp::init<
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<std::string>&,
                const std::vector<int32_t>&
                >()
        )
        .def(bp::init<>())
        .def("print_geometry", &CMolecule::printGeometry)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
    ;

    // CMolecularBasis class

    bp::class_< CMolecularBasis >
        (
            "CMolecularBasis", bp::init<>()
        )
        .def("get_label", &CMolecularBasis::getLabel)
    ;

    // ----------------------
    // ==> Matrix classes <==
    // ----------------------

    // CDenseMatrix class

    bp::class_< CDenseMatrix, std::shared_ptr<CDenseMatrix> >
        (
            "CDenseMatrix", bp::init<
                const std::vector<double>&,
                const int32_t,
                const int32_t
                >()
        )
        .def(bp::init<>())
        .def(bp::init<const int32_t, const int32_t>())
        .def(bp::init<const int32_t>())
        .def("get_number_of_rows", &CDenseMatrix::getNumberOfRows)
        .def("get_number_of_columns", &CDenseMatrix::getNumberOfColumns)
    ;

    // COverlapMatrix class

    bp::class_< COverlapMatrix, std::shared_ptr<COverlapMatrix> >
        (
            "COverlapMatrix",
            bp::init<const CDenseMatrix&>()
        )
        .def(bp::init<>())
        .def("__str__", &COverlapMatrix_str);
    ;

    // -------------------------------
    // ==> Integral driver classes <==
    // -------------------------------

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
        .staticmethod("create")
    ;
}
