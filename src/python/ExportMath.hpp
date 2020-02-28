//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>

#include "DenseMatrix.hpp"

namespace py = pybind11;

namespace vlx_math {  // vlx_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
std::shared_ptr<CDenseMatrix> CDenseMatrix_from_numpy(const py::array_t<double>& arr);

/**
 Exports classes/functions in src/math to python.
 */
void export_math(py::module& m);

}  // namespace vlx_math

#endif /* ExportMath_hpp */
