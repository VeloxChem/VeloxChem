//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/pybind11.h>

#include "ExportDFT.hpp"
#include "ExportGeneral.hpp"
#include "ExportGpu.hpp"
#include "ExportMath.hpp"
#include "ExportMolData.hpp"
#include "ExportOneInts.hpp"
#include "ExportOrbData.hpp"
#include "ExportResponse.hpp"
#include "ExportTwoInts.hpp"
#include "ExportVisualization.hpp"

PYBIND11_MODULE(veloxchemlib, m)
{
    vlx_general::export_general(m);

    vlx_moldata::export_moldata(m);

    vlx_orbdata::export_orbdata(m);

    vlx_oneints::export_oneints(m);

    vlx_twoints::export_twoints(m);

    vlx_math::export_math(m);

    vlx_gpu::export_gpu(m);

    vlx_visualization::export_visualization(m);

    vlx_response::export_response(m);

    vlx_dft::export_dft(m);
}
