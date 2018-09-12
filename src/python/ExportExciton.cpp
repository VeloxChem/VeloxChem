//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "AssembleMatrices.hpp"

namespace bp = boost::python;

// ==> boost python <==
// functions and classes

void export_exciton()
{
    bp::def("assemble_overlap_matrices", &assembleOverlapMatrices);
}
