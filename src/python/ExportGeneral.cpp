#include "ExportGeneral.hpp"

#include "Codata.hpp"
#include "StringFormat.hpp"

namespace vlx_general {  // vlx_general namespace

// Exports classes/functions in src/general to python

auto
export_general(py::module& m) -> void
{
    // exposing enum from FmtType.hpp
    
    py::enum_<fmt>(m, "fmt")
        .value("center", fmt::center)
        .value("left",   fmt::left)
        .value("right",  fmt::right);

    // exposing functions from Codata.hpp

    m.def("bohr_in_angstroms",
          &units::getBohrValueInAngstroms,
          "Gets Bohr value in Angstroms.");
    
    m.def("hartree_in_ev",
          &units::getHartreeValueInElectronVolts,
          "Gets Hartree value in electronvolts.");
    
    // exposing functions from StringFormat.hpp
    
    m.def("upcase",
          &fstr::upcase,
          "Convers string to upper case string.");
    
    m.def("format",
          &fstr::format,
          "Formats string to string with specified alignment.");
    
    m.def("to_string",
          [](const double source,
             const size_t presicion,
             const size_t width,
             const fmt    aligment) -> std::string
          { return fstr::to_string(source, presicion, width, aligment); },
          "Formats double precision number to string with specified alignment.");
    
    m.def("to_string",
          [](const double source,
             const size_t presicion) -> std::string
          { return fstr::to_string(source, presicion); },
          "Formats double precision number to string with specified alignment.");
    
    m.def("to_string",
          [](const int64_t source,
             const size_t  width,
             const fmt     aligment) -> std::string
          { return fstr::to_string(source, width, aligment); },
          "Formats integer number to string with specified alignment.");
    
    m.def("to_string",
          [](const bool source) -> std::string
          { return fstr::to_string(source); },
          "Formats bool to string.");
    
    m.def("to_angular_momentum",
          [](const int64_t angl) -> std::string
          { return fstr::to_AngularMomentum(angl); },
          "Converts angular momentum integer to string.");
    
    m.def("to_angular_momentum",
          [](const std::string& label) -> int64_t
          { return fstr::to_AngularMomentum(label); },
          "Converts angular momentum string to integer.");
}

}  // namespace vlx_general
