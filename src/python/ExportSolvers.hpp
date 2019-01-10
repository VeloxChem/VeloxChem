//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ExportSolvers_hpp
#define ExportSolvers_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_solvers { // vlx_solvers namespace

/**
 Exports classes/functions in src/solvers to python.
 */
void export_solvers(py::module& m);

} // vlx_solvers namespace

#endif /* ExportSolvers_hpp */
