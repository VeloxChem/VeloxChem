#include "ExportT3CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
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
        .def("prepare_buffers", &CRIFockDriver::prepare_buffers, "Computes three center electron repulsion integral buffers.")
        .def("compute", &CRIFockDriver::compute, "Computes Coulomb Fock matrix for given density.");
    
    // ThreeCenterElectronRepulsionGeom100Driver class
    PyClass<CThreeCenterElectronRepulsionGeomX00Driver<1>>(m, "ThreeCenterElectronRepulsionGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, CMolecularBasis& aux_basis, const int iatom)
             -> TPoint<double> { return geom_drv.compute(basis, aux_basis, molecule, iatom); },
            "Computes gradient contribution for given molecule, basis, auxilary basis and selected atom.");
}

}  // namespace vlx_t2cintegrals
