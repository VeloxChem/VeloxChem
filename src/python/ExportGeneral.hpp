#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ExportHelpers.hpp"

namespace py = pybind11;

namespace vlx_general {

/// @brief Exports classes/functions in src/orbdata to Python module.
/// @param m The Python module.
auto export_general(py::module &m) -> void;
}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
