//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include "ExportGeneral.hpp"
#include "ExportStreams.hpp"
#include "ExportMolData.hpp"
#include "ExportOrbData.hpp"
#include "ExportOneInts.hpp"
#include "ExportTwoInts.hpp"
#include "ExportMath.hpp"
#include "ExportGpu.hpp"
#include "ExportSolvers.hpp"
#include "ExportExciton.hpp"
#include "ExportVisualization.hpp"

BOOST_PYTHON_MODULE(veloxchemlib)
{
    bp_general::export_general();

    bp_streams::export_streams();

    bp_moldata::export_moldata();

    bp_orbdata::export_orbdata();

    bp_oneints::export_oneints();

    bp_twoints::export_twoints();

    bp_math::export_math();

    bp_gpu::export_gpu();

    bp_solvers::export_solvers();

    bp_exciton::export_exciton();

    bp_visualization::export_visualization();
}
