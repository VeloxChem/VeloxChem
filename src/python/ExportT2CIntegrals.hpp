#ifndef ExportT2CIntegrals_hpp
#define ExportT2CIntegrals_hpp

#include "ExportHelpers.hpp"

namespace vlx_t2cintegrals {  // vlx_orbdata namespace

/**
 Exports classes/functions in src/t2c_*  to python.
 */
auto export_t2cintegrals(py::module& m) -> void;

}  // namespace vlx_t2cintegrals

#endif /* ExportT2CIntegrals_hpp */
