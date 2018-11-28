//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "VisualizationDriver.hpp"
#include "ExportVisualization.hpp"

namespace bp = boost::python;

namespace bp_visualization { // bp_visualization namespace

// Exports classes/functions in src/visualization to python

void export_visualization()
{
    // CVisualizationDriver class
    // Note: Need member function pointers for proper overloading

    double (CVisualizationDriver::*compute_mol_orb)(
            const CMolecule&          molecule,
            const CMolecularBasis&    basis,
            const CMolecularOrbitals& mo,
            const int32_t             moidx,
            const std::string&        mospin,
            const double              xp,
            const double              yp,
            const double              zp) const
        = &CVisualizationDriver::compute;

    double (CVisualizationDriver::*compute_density)(
            const CMolecule&        molecule,
            const CMolecularBasis&  basis,
            const CAODensityMatrix& density,
            const int32_t           denidx,
            const std::string&      denspin,
            const double            xp,
            const double            yp,
            const double            zp) const
        = &CVisualizationDriver::compute;

    bp::class_< CVisualizationDriver, std::shared_ptr<CVisualizationDriver> >
        (
            "VisualizationDriver",
            bp::init<>()
        )
        .def("compute", compute_mol_orb)
        .def("compute", compute_density)
    ;
}

} // bp_visualization namespace
