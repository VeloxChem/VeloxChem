#include "ExportT3CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "T3FlatBuffer.hpp"
#include "RIFockDriver.hpp"

namespace vlx_t3cintegrals {

// Exports classes/functions in src/t3c_* to python

void
export_t3cintegrals(py::module& m)
{
    // CThreeCenterElectronRepulsionDriver class
    PyClass<CThreeCenterElectronRepulsionDriver>(m, "ThreeCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionDriver& eri_drv, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularBasis& aux_basis) -> CT3FlatBuffer<double> {
                return eri_drv.compute(basis, aux_basis, molecule);
            },
            "Computes electron repulsion integrals for given molecule, basis and auxilary basis.");
    
    // CRIFockDriver class
    PyClass<CRIFockDriver>(m, "RIFockDriver")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def("prepare_buffers", &CRIFockDriver::prepare_buffers, "Computes three center electron repulsion integral buffers.")
        .def("compute", &CRIFockDriver::compute, "Computes Coulomb Fock matrix for given density.");
    
}

}  // namespace vlx_t2cintegrals
