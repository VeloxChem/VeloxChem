//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportSolvers_hpp
#define ExportSolvers_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace vlx_solvers {  // vlx_solvers namespace

/**
 Exports classes/functions in src/solvers to python.
 */
void export_solvers(py::module& m);

}  // namespace vlx_solvers

#endif /* ExportSolvers_hpp */
