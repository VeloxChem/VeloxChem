//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>

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

PYBIND11_MODULE(veloxchemlib, m)
{
    bp_general::export_general(m);

    bp_streams::export_streams(m);

    bp_moldata::export_moldata(m);

    bp_orbdata::export_orbdata(m);

    bp_oneints::export_oneints(m);

    bp_twoints::export_twoints(m);

    bp_math::export_math(m);

    bp_gpu::export_gpu(m);

    bp_solvers::export_solvers(m);

    bp_exciton::export_exciton(m);

    bp_visualization::export_visualization(m);
}
