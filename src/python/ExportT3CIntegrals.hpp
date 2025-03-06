#ifndef ExportT3CIntegrals_hpp
#define ExportT3CIntegrals_hpp

#include "ExportHelpers.hpp"

namespace vlx_t3cintegrals {

/**
 Exports classes/functions in src/t3c_*  to python.
 */
auto export_t3cintegrals(py::module& m) -> void;

}  // namespace vlx_t3cintegrals

#endif /* ExportT3CIntegrals_hpp */
