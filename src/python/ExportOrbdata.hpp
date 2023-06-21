#ifndef ExportOrbdata_hpp
#define ExportOrbdata_hpp

#include "ExportHelpers.hpp"

namespace vlx_orbdata {  // vlx_orbdata namespace

/**
 Exports classes/functions in src/orbdata to python.
 */
auto
export_orbdata(py::module& m) -> void;

}  // namespace vlx_orbdata

#endif /* ExportOrbdata_hpp */
