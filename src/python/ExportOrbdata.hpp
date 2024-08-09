#ifndef ExportOrbdata_hpp
#define ExportOrbdata_hpp

#include "ExportHelpers.hpp"

namespace vlx_orbdata {  // vlx_orbdata namespace

/// @brief Exports classes/functions in src/orbdata to Python module.
/// @param m The Python module.
auto export_orbdata(py::module &m) -> void;

}  // namespace vlx_orbdata

#endif /* ExportOrbdata_hpp */
