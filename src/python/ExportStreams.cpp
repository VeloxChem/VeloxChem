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

#include "StringFormat.hpp"
#include "ExportStreams.hpp"

namespace py = pybind11;

namespace vlx_streams { // vlx_streams namespace

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
    // to_angular_momentum methods

    m.def("to_angular_momentum", &string_to_angular_momentum);

    m.def("to_angular_momentum", &integer_to_angular_momentum);
}

} // vlx_streams namespace
