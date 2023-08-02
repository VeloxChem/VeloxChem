#ifndef ExportT4CIntegrals_hpp
#define ExportT4CIntegrals_hpp

#include "ExportHelpers.hpp"

namespace vlx_t4cintegrals {  // vlx_orbdata namespace

/**
 Exports classes/functions in src/t4c_*  to python.
 */
auto export_t4cintegrals(py::module& m) -> void;

}  // namespace vlx_t4cintegrals

#endif /* ExportT4CIntegrals_hpp */
