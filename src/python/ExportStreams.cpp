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

#include "InputStream.hpp"
#include "OutputStream.hpp"

namespace bp = boost::python;

// ==> boost python <==
// functions and classes

void export_streams()
{
    // COutputStream class

    bp::class_< COutputStream, std::shared_ptr<COutputStream> >
        (
            "COutputStream",
            bp::init<const std::string&>()
        )
        .def("get_state", &COutputStream::getState)
        .def("flush",     &COutputStream::flush)
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
