//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "SADGuess.hpp"

#include "ExportSolvers.hpp"

namespace bp = boost::python;

namespace bp_solvers { // bp_solvers namespace

// Exports classes/functions in src/solvers to python

void export_solvers()
{
    // exposing functions

    bp::def("get_sad_initial_guess",
            &sad_guess::getSADInitialGuess);
}

} // bp_solvers namespace
