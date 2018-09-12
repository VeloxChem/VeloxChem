//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "MolecularBasis.hpp"

namespace bp = boost::python;

// ==> boost python <==
// functions and classes

void export_orbdata()
{
    // CMolecularBasis class

    bp::class_< CMolecularBasis >
        (
            "CMolecularBasis",
            bp::init<>()
        )
        .def("get_label", &CMolecularBasis::getLabel)
    ;
}
