#include "ExportGeneral.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "Codata.hpp"
#include "StringFormat.hpp"
#include "BatchFunc.hpp"
#include "OpenMPFunc.hpp"
#include "GtoBlock.hpp"

namespace vlx_general {  // vlx_general namespace

// Exports classes/functions in src/general to python

auto
export_general(py::module& m) -> void
{
    // exposing enum from FmtType.hpp
    
    py::enum_<fmt_t>(m, "fmt_t")
        .value("center", fmt_t::center)
        .value("left",   fmt_t::left)
        .value("right",  fmt_t::right);

    // exposing functions from Codata.hpp

    m.def("bohr_in_angstroms",
          &units::getBohrValueInAngstroms,
          "Gets Bohr value in Angstroms.");
    
    m.def("hartree_in_ev",
          &units::getHartreeValueInElectronVolts,
          "Gets Hartree value in electronvolts.");
    
    // exposing functions from BatchFunc.hpp

    m.def("get_batch_index",
          &batch::getBatchIndex,
          "Gets starting index of batch.");
    
    m.def("number_of_batches",
          &batch::getNumberOfBatches,
          "Gets number of batches for given vector size.");
    
    m.def("get_batch_range",
          &batch::getBatchRange,
          "Gets range of specific batch.");
    
    // exposing functions from OpenMPFunc.hpp

    m.def("set_number_of_threads",
          &omp::setNumberOfThreads,
          "Sets number of OMP threads to requested value.");
    
    m.def("get_number_of_threads",
          &omp::getNumberOfThreads,
          "Gets number of OMP threads available.");
    
    m.def("make_workgroup",
         py::overload_cast<const std::vector<CGtoBlock>&>(&omp::makeWorkGroup),
         "Gets work group for given vector of basis function blocks.");
    
    m.def("make_workgroup",
         py::overload_cast<const std::vector<CGtoBlock>&,
                           const std::vector<CGtoBlock>&>(&omp::makeWorkGroup),
         "Gets work group for given two vectors of basis function blocks.");
    
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
             const fmt_t  aligment) -> std::string
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
             const fmt_t   aligment) -> std::string
          { return fstr::to_string(source, width, aligment); },
          "Formats integer number to string with specified alignment.");
    
    m.def("to_string",
          [](const bool source) -> std::string
          { return fstr::to_string(source); },
          "Formats bool to string.");
    
    m.def("to_angular_momentum",
          [](const int64_t angmom) -> std::string
          { return fstr::to_AngularMomentum(angmom); },
          "Converts angular momentum integer to string.");
    
    m.def("to_angular_momentum",
          [](const std::string& label) -> int64_t
          { return fstr::to_AngularMomentum(label); },
          "Converts angular momentum string to integer.");
}

}  // namespace vlx_general
