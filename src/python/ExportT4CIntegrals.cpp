#include "ExportT4CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "FockDriver.hpp"
#include "T4CUtils.hpp"
#include "T4CScreener.hpp"

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
    
    // CT4CScreener class
    PyClass<CT4CScreener>(m, "T4CScreener")
        .def(py::init<>())
        .def("partition", &CT4CScreener::partition, "Partition basis funtion pairs blocks for given molecule and basis.")
        .def("gto_pair_blocks", &CT4CScreener::gto_pair_blocks, "Gets vector of blocked basis function pairs blocks.");

}

}  // namespace vlx_t4cintegrals
