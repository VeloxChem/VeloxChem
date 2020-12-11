//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <pybind11/pybind11.h>

#include "ExportDFT.hpp"
#include "ExportXTB.hpp"
#include "ExportGeneral.hpp"
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

    vlx_visualization::export_visualization(m);

    vlx_response::export_response(m);

    vlx_dft::export_dft(m);

    vlx_xtb::export_xtb(m);
}
