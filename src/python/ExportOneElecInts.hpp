#ifndef ExportOneElecInts_hpp
#define ExportOneElecInts_hpp

#include "ExportHelpers.hpp"

namespace vlx_oneeints {

/**
 Exports classes/functions in src/t2c_*  to python.
 */
auto export_oneeints(py::module& m) -> void;

}  // namespace vlx_oneeints

#endif /* ExportOneElecInts_hpp */
