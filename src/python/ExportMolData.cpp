//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include <memory>
#include <vector>
#include <string>

#include "Molecule.hpp"

namespace bp = boost::python;

// ==> boost python <==
// functions and classes

void export_moldata()
{
    // CMolecule class
    // Note: CMolecule has two constructors

    bp::class_< CMolecule, std::shared_ptr<CMolecule> >
        (
            "CMolecule",
            bp::init<
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<double>&,
                const std::vector<std::string>&,
                const std::vector<int32_t>&
                >()
        )
        .def(bp::init<>())
        .def("print_geometry", &CMolecule::printGeometry)
        .def("get_sub_molecule", &CMolecule::getSubMolecule)
    ;
}
