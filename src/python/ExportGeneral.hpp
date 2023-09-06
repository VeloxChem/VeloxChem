#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <pybind11/numpy.h>

#include "ExportHelpers.hpp"

namespace py = pybind11;

namespace vlx_general {  // vlx_general namespace

/** Gets numpy array from pointer and shape.
 *
 * @tparam T scalar type of array.
 * @param ptr pointer to data.
 * @param dimension the shape of numpy array.
 * @return numpy array.
 */
auto pointer_to_numpy(const double* ptr, const std::vector<int64_t>& dimension) -> py::array_t<double>;

/**
 Exports classes/functions in src/general to python.
 */
auto export_general(py::module& m) -> void;

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
