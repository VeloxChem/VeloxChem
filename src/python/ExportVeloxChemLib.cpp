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
    vlx_general::export_general(m);

    vlx_streams::export_streams(m);

    vlx_moldata::export_moldata(m);

    vlx_orbdata::export_orbdata(m);

    vlx_oneints::export_oneints(m);

    vlx_twoints::export_twoints(m);

    vlx_math::export_math(m);

    vlx_gpu::export_gpu(m);

    vlx_solvers::export_solvers(m);

    vlx_exciton::export_exciton(m);

    vlx_visualization::export_visualization(m);
}
