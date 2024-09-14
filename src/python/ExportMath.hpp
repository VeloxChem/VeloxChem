#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "DenseMatrix.hpp"
#include "ExportHelpers.hpp"

namespace vlx_math {

/**
 Converts numpy array to CDenseMatrix

 @param arr the numpy array.
 @return a CDenseMatrix object.
 */
auto CDenseMatrix_from_numpy(const py::array_t<double>& arr) -> std::shared_ptr<CDenseMatrix>;

/// @brief Exports classes/functions in src/math to Python module.
/// @param m The Python module.
auto export_math(py::module &m) -> void;

}  // namespace vlx_math

#endif /* ExportMath_hpp */
