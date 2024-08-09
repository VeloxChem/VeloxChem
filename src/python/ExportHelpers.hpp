#ifndef ExportHelpers_hpp
#define ExportHelpers_hpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename T>
using PyClass = py::class_<T, std::shared_ptr<T>>;

#endif /* ExportHelpers_hpp */
