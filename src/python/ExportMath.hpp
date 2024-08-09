#ifndef ExportMath_hpp
#define ExportMath_hpp

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "ExportHelpers.hpp"

namespace vlx_math {

/// @brief Exports classes/functions in src/math to Python module.
/// @param m The Python module.
auto export_math(py::module &m) -> void;
}  // namespace vlx_math

#endif /* ExportMath_hpp */
