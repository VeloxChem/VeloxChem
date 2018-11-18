//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "CubeGenerator.hpp"
#include "ExportVisualization.hpp"

namespace bp = boost::python;

namespace bp_visualization { // bp_visualization namespace

// Exports classes/functions in src/visualization to python

void export_visualization()
{
    bp::def("get_psi_molecular_orbital",
            &cubes::getPsiMolecularOrbital);

    bp::def("get_psi_density",
            &cubes::getPsiDensity);
}

} // bp_visualization namespace
