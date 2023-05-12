#ifndef ExportMoldata_hpp
#define ExportMoldata_hpp

#include "ExportHelpers.hpp"

namespace vlx_moldata {  // vlx_moldata namespace

/**
 Exports classes/functions in src/moldata to python.
 */
auto
export_moldata(py::module& m) -> void;

}  // namespace vlx_moldata

#endif /* ExportMoldata_hpp */
