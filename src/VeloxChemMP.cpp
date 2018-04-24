//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cstdlib>

#include "MpiFunc.hpp"
#include "AppManager.hpp"

/**
 Executes a VeloxChemMP program.

 @param argc the number of command line arguments.
 @param argv the array of command line arguments.
 @return the program state: EXIT_SUCCESS for normal state and EXIT_FAILURE
         for abnormal state
 */
int
main(int    argc,
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
