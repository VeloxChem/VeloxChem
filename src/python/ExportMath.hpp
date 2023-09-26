#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>

#include "DenseMatrix.hpp"
#include "ExportHelpers.hpp"

namespace vlx_math {  // vlx_math namespace

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
auto CDenseMatrix_from_numpy(const py::array_t<double>& arr) -> std::shared_ptr<CDenseMatrix>;

/**
 Exports classes/functions in src/math to python.
 */
auto export_math(py::module& m) -> void;

}  // namespace vlx_math

#endif /* ExportMath_hpp */
