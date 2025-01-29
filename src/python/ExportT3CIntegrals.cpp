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
            [](const CThreeCenterElectronRepulsionDriver& eri_drv, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularBasis& aux_basis) -> std::shared_ptr<CT3FlatBuffer<double>> {
                return std::make_shared<CT3FlatBuffer<double>>(eri_drv.compute(basis, aux_basis, molecule));
            },
            "Computes electron repulsion integrals for given molecule, basis and auxilary basis.");
    
    // CRIFockDriver class
    PyClass<CRIFockDriver>(m, "RIFockDriver")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def(py::init<const CSubMatrix&, const CSubMatrix&>())
        .def("prepare_buffers", &CRIFockDriver::prepare_buffers, "Computes three center electron repulsion integral buffers.")
        .def("comp_gamma_vector", &CRIFockDriver::comp_gamma_vector, "Computes Gamma vector required for J fitting.")
        .def("trafo_gamma_vector", &CRIFockDriver::trafo_gamma_vector, "Transforms Gamma vector with J metric.")
        .def("comp_j_vector", &CRIFockDriver::comp_j_vector, "Computes effective flatened J vector.")
        .def("compute", &CRIFockDriver::compute, "Computes Fock matrix for given density."); 
    
}

}  // namespace vlx_t2cintegrals
