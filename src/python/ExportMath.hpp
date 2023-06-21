#ifndef ExportMath_hpp
#define ExportMath_hpp

#include "ExportHelpers.hpp"

namespace vlx_math {  // vlx_math namespace

/**
 Exports classes/functions in src/math to python.
 */
auto
export_math(py::module& m) -> void;

}  // namespace vlx_math

#endif /* ExportMath_hpp */
