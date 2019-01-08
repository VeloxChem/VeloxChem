//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>

#include <memory>
#include <string>

#include "OutputStream.hpp"
#include "StringFormat.hpp"
#include "ExportStreams.hpp"

namespace py = pybind11;

namespace bp_streams { // bp_streams namespace

// Helper function for writing to output stream

static void
COutputStream_print_line(      COutputStream& self,
                         const std::string&   source)
{
    self << source.c_str() << fmt::end;
}

static void
COutputStream_print_info(      COutputStream& self,
                         const std::string&   source)
{
    self << fmt::info << source.c_str() << fmt::end;
}
    
static void
COutputStream_print_title(      COutputStream& self,
                          const std::string&   source)
{
    self << fmt::title << source.c_str() << fmt::end;
}
    
static void
COutputStream_print_header(      COutputStream& self,
                           const std::string&   source)
{
    self << fmt::header << source.c_str() << fmt::end;
}

static void
COutputStream_print_separator(COutputStream& self)
{
    self << fmt::title << fmt::tsep;
}
    
static void
COutputStream_print_blank(COutputStream& self)
{
    self << fmt::blank;
}

// Helper function for converting angular momentum

static std::string
string_to_angular_momentum(const int32_t angl)
{
    return fstr::to_AngularMomentum(angl);
}

static int32_t
integer_to_angular_momentum(const std::string& label)
{
    return fstr::to_AngularMomentum(label);
}

// Exports classes/functions in src/streams to python

void export_streams(py::module& m)
{
    // COutputStream class

    py::class_< COutputStream, std::shared_ptr<COutputStream> >
        (
            m, "OutputStream"
        )
        .def(py::init<const std::string&>())
        .def("get_state", &COutputStream::getState)
        .def("flush", &COutputStream::flush)
        .def("print_line", &COutputStream_print_line)
        .def("print_info", &COutputStream_print_info)
        .def("print_title", &COutputStream_print_title)
        .def("print_header", &COutputStream_print_header)
        .def("print_separator", &COutputStream_print_separator)
        .def("print_blank", &COutputStream_print_blank)
    ;

    // to_angular_momentum methods

    m.def("to_angular_momentum", &string_to_angular_momentum);

    m.def("to_angular_momentum", &integer_to_angular_momentum);
}

} // bp_streams namespace
