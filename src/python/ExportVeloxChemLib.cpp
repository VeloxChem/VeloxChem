#include <pybind11/pybind11.h>

#include "ExportMath.hpp"
#include "ExportMoldata.hpp"
#include "ExportOrbdata.hpp"
#include "ExportGeneral.hpp"
#include "ExportT2CIntegrals.hpp"

PYBIND11_MODULE(veloxchemlib, m)
{
    vlx_general::export_general(m);
    
    vlx_math::export_math(m);
    
    vlx_moldata::export_moldata(m);
    
    vlx_orbdata::export_orbdata(m);
    
    vlx_t2cintegrals::export_t2cintegrals(m);
}
