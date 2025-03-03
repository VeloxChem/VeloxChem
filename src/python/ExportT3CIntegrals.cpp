#include "ExportT3CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
#include "ThreeCenterElectronRepulsionGeom0X0Driver.hpp"
#include "T3FlatBuffer.hpp"
#include "RIFockDriver.hpp"
#include "RIFockGradDriver.hpp"

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
        .def("compute", &CRIFockDriver::compute, "Computes Coulomb Fock matrix for given density.")
        .def("compute_bq_vector", &CRIFockDriver::compute_bq_vector, "Computes transformed Gamma vector for given density.");
    
    // CRIFockGradDriver class
    PyClass<CRIFockGradDriver>(m, "RIFockGradDriver")
        .def(py::init<>())
        .def("compute", &CRIFockGradDriver::compute, "Computes Coulomb Fock contribution to atom gradient.");
    
    // ThreeCenterElectronRepulsionGeom100Driver class
    PyClass<CThreeCenterElectronRepulsionGeomX00Driver<1>>(m, "ThreeCenterElectronRepulsionGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis,  CMolecularBasis& aux_basis, const int iatom) -> std::shared_ptr<CT3FlatBuffer<double>> {
                return std::make_shared<CT3FlatBuffer<double>>(geom_drv.compute(basis, aux_basis, molecule, iatom));
             },
            "Computes gradient integrals for given molecule, basis, auxilary basis and selected atom.");
    
    
    // ThreeCenterElectronRepulsionGeom010Driver class
    PyClass<CThreeCenterElectronRepulsionGeom0X0Driver<1>>(m, "ThreeCenterElectronRepulsionGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis,  CMolecularBasis& aux_basis, const int iatom) -> std::shared_ptr<CT3RectFlatBuffer<double>> {
                return std::make_shared<CT3RectFlatBuffer<double>>(geom_drv.compute(basis, aux_basis, molecule, iatom));
             },
            "Computes gradient integrals for given molecule, basis, auxilary basis and selected atom.");
}

}  // namespace vlx_t2cintegrals
