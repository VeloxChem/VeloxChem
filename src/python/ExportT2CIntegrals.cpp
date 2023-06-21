#include "ExportT2CIntegrals.hpp"

#include <pybind11/operators.h>

#include "OverlapDriver.hpp"

namespace vlx_t2cintegrals {  // vlx_t2cintegrals namespace

// Exports classes/functions in src/t2c_* to python

void
export_t2cintegrals(py::module& m)
{
   
    // COverlapDriver class

    PyClass<COverlapDriver>(m, "OverlapDriver")
        .def(py::init<>())
        .def("compute",
             [](const COverlapDriver& ovl_drv,
                const CMolecularBasis& basis,
                const CMolecule&       molecule) -> std::shared_ptr<CMatrix>
             {
                return std::make_shared<CMatrix>(ovl_drv.compute(basis, molecule)); 
             },
             "Computes overlap matrix for given molecule and basis.");
    
    // ...
}

}  // namespace vlx_orbdata
