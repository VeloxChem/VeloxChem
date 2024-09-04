#include "ExportT4CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "FockDriver.hpp"
#include "T4CUtils.hpp"

namespace vlx_t4cintegrals {  // vlx_t4cintegrals namespace

// Exports classes/functions in src/t4c_* to python

void
export_t4cintegrals(py::module& m)
{
    // FockDriver class

    PyClass<CFockDriver>(m, "FockDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockDriver&     fock_drv,
               const CMolecularBasis& basis,
               const CMolecule&       molecule,
               const CMatrix&         density,
               const std::string&     label,
               const double           exchange_factor,
               const double           omega) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(fock_drv.compute(basis, molecule, density, label, exchange_factor, omega));
            },
             "Computes single Fock matrix of requested type for given molecule and basis.");
}

}  // namespace vlx_t4cintegrals
