//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "InputData.hpp"
#include "OutputStream.hpp"
#include "BasisReader.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "ExportReaders.hpp"

namespace bp = boost::python;

namespace bp_readers { // bp_readers namespace

// Exports classes/functions in src/readers to python

void export_readers()
{
    // CBasisReader class

    bp::class_< CBasisReader, std::shared_ptr<CBasisReader> >
        (
            "BasisReader",
            bp::init<>()
        )
        .def("parse", &CBasisReader::parse)
        .def("get_state", &CBasisReader::getState)
        .def("get_ao_basis", &CBasisReader::getAOBasis)
        .def("get_rij_basis", &CBasisReader::getRIJBasis)
        .def("get_min_basis", &CBasisReader::getMinBasis)
    ;
}

} // bp_readers namespace
