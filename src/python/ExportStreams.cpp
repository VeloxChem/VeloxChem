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

#include "InputStream.hpp"
#include "OutputStream.hpp"

#include "ExportStreams.hpp"

namespace bp = boost::python;

namespace bp_streams { // bp_streams namespace

// Helper function for writing to output stream

void COutputStream_put_info(      COutputStream& self,
                            const std::string&   source)
{
    self << fmt::info << source.c_str() << fmt::end;
}

void COutputStream_new_line(COutputStream& self)
{
    self << fmt::blank;

    self.flush();
}

// Exports classes/functions in src/streams to python

void export_streams()
{
    // COutputStream class

    bp::class_< COutputStream, std::shared_ptr<COutputStream> >
        (
            "COutputStream",
            bp::init<const std::string&>()
        )
        .def("get_state", &COutputStream::getState)
        .def("flush", &COutputStream::flush)
        .def("put_info", &COutputStream_put_info)
        .def("new_line", &COutputStream_new_line)
    ;

    // CInputStream class

    bp::class_< CInputStream, std::shared_ptr<CInputStream> >
        (
            "CInputStream",
            bp::init<const std::string&, COutputStream&>()
        )
        .def("read",      &CInputStream::read)
        .def("get_state", &CInputStream::getState)
    ;

    // CInputData class

    bp::class_< CInputData, std::shared_ptr<CInputData> >
        (
            "CInputData",
            bp::init<>()
        )
    ;
}

} // bp_streams namespace
