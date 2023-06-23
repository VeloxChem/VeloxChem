#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include "ExportHelpers.hpp"

namespace vlx_general {  // vlx_general namespace

/**
 Exports classes/functions in src/general to python.
 */
auto export_general(py::module& m) -> void;

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
