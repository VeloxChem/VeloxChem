//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <memory>
#include <string>
#include <iostream>

#include "OutputStream.hpp"
#include "StringFormat.hpp"
#include "ExportStreams.hpp"

namespace bp = boost::python;

namespace bp_streams { // bp_streams namespace

// Helper function for writing to output stream

static void
COutputStream_put_info(      COutputStream& self,
                       const std::string&   source)
{
    self << fmt::info << source.c_str() << fmt::end;
}
    
static void
COutputStream_put_title(      COutputStream& self,
                        const std::string&   source)
{
    self << fmt::title << source.c_str() << fmt::end;
}
    
static void
COutputStream_put_header(      COutputStream& self,
                         const std::string&   source)
{
    self << fmt::header << source.c_str() << fmt::end;
}

static void
COutputStream_put_separator(COutputStream& self)
{
    self << fmt::title << fmt::tsep;
}
    
static void
COutputStream_new_line(COutputStream& self)
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

void export_streams()
{
    // COutputStream class

    bp::class_< COutputStream, std::shared_ptr<COutputStream> >
        (
            "OutputStream",
            bp::init<const std::string&>()
        )
        .def("get_state", &COutputStream::getState)
        .def("flush", &COutputStream::flush)
        .def("put_info", &COutputStream_put_info)
        .def("put_title", &COutputStream_put_title)
        .def("put_header", &COutputStream_put_header)
        .def("put_separator", &COutputStream_put_separator)
        .def("new_line", &COutputStream_new_line)
    ;

    // to_angular_momentum methods

    bp::def("to_angular_momentum", &string_to_angular_momentum);

    bp::def("to_angular_momentum", &integer_to_angular_momentum);
}

} // bp_streams namespace
