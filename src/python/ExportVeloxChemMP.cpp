//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

//#include <mpi4py/mpi4py.h>
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>

#include "MpiFunc.hpp"
#include "AppManager.hpp"

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

        return vlx_run(argc, argv);
    }

}

BOOST_PYTHON_MODULE(VeloxChemMP)
{
    // initialize mpi4py's C-API
    //if (import_mpi4py() < 0) return;

    // ==> initialize module <==

    // VeloxChem
    bp::def("run", PyVLX::py_vlx_run);

    // classes
    bp::class_< CAppManager, std::shared_ptr<CAppManager> > ("CAppManager", bp::init<int, char**>())
        .def("create", &CAppManager::create)
        .staticmethod("create")
        .def("execute", &CAppManager::execute)
    ;
}
