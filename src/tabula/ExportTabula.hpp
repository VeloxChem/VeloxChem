//
//  Tabula — custom-recursion molecular-integral machinery.
//  Python (pybind11) bindings.
//

#ifndef ExportTabula_hpp
#define ExportTabula_hpp

#include <pybind11/pybind11.h>

namespace tabula {  // tabula namespace

/// @brief Exports Tabula classes/functions to the Python module.
/// @param m The Python module.
auto export_tabula(pybind11::module& m) -> void;

}  // namespace tabula

#endif /* ExportTabula_hpp */
