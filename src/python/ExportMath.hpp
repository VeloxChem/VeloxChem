//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>

#include "DenseMatrix.hpp"

namespace py = pybind11;

namespace vlx_math { // vlx_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
std::shared_ptr<CDenseMatrix>
CDenseMatrix_from_numpy(const py::array_t<double>& arr);

/**
 Exports classes/functions in src/math to python.
 */
void export_math(py::module& m);

} // vlx_math namespace

#endif /* ExportMath_hpp */
