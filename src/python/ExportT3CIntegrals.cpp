#include "ExportT3CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
#include "ThreeCenterElectronRepulsionGeom0X0Driver.hpp"
#include "T3FlatBuffer.hpp"
#include "RIFockDriver.hpp"
#include "RIFockGradDriver.hpp"
#include "T4CScreener.hpp"

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
            "Computes electron repulsion integrals for given molecule, basis and auxilary basis.")
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionDriver& eri_drv, const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularBasis& aux_basis, const std::vector<int>& atoms) -> CT3FlatBuffer<double> {
                return eri_drv.compute(basis, aux_basis, molecule, atoms);
            },
            "Computes electron repulsion integrals for given molecule, basis, auxilary basis, and list of atoms.");
    
    // CRIFockDriver class
    PyClass<CRIFockDriver>(m, "RIFockDriver")
        .def(py::init<>())
        .def(py::init<const CSubMatrix&>())
        .def("prepare_buffers", py::overload_cast<const CMolecule&,
                                                  const CMolecularBasis&,
                                                  const CMolecularBasis&>
             (&CRIFockDriver::prepare_buffers),
             "Computes three center electron repulsion integral buffers.")
        .def("prepare_buffers", py::overload_cast<const CMolecule&,
                                                  const CMolecularBasis&,
                                                  const CMolecularBasis&,
                                                  const std::vector<int>&>
             (&CRIFockDriver::prepare_buffers),
             "Computes three center electron repulsion integral buffers.")
        .def("compute",  py::overload_cast<const CMatrix&, const std::string&>
             (&CRIFockDriver::compute, py::const_),
             "Computes Coulomb Fock matrix for given density.")
        .def("compute",  py::overload_cast<const CMatrix&, const std::vector<double>&, const std::string&>
             (&CRIFockDriver::compute, py::const_),
             "Computes Coulomb Fock matrix for given density.")
        .def("local_compute", &CRIFockDriver::local_compute, "Computes local Coulomb Fock matrix for given density.")
        .def("compute_bq_vector", &CRIFockDriver::compute_bq_vector, "Computes transformed Gamma vector for given density.")
        .def("compute_local_bq_vector", &CRIFockDriver::compute_local_bq_vector, "Computes transformed local Gamma vector for given density.");
    
    // CRIFockGradDriver class
    PyClass<CRIFockGradDriver>(m, "RIFockGradDriver")
        .def(py::init<>())
        .def("compute",
             [] (const CRIFockGradDriver&   grad_drv,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const int                  iatom) -> TPoint<double>
             {
                return grad_drv.compute(basis, aux_basis, molecule, gamma, density, iatom);
             },
          "Computes Coulomb Fock contribution to atom's gradient.")
        .def("compute",
             [] (const CRIFockGradDriver&   grad_drv,
                 const CT4CScreener&        screener,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const int                  iatom,
                 const int                  ithreshold) -> TPoint<double>
             {
                return grad_drv.compute(screener, basis, aux_basis, molecule, gamma, density, iatom, ithreshold);
             },
          "Computes Coulomb Fock contribution to atom's gradient.")
        .def("direct_compute",
             [] (const CRIFockGradDriver&   grad_drv,
                 const CT4CScreener&        screener,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const int                  iatom,
                 const int                  ithreshold) -> TPoint<double>
             {
                return grad_drv.direct_compute(screener, basis, aux_basis, molecule, gamma, density, iatom, ithreshold);
             },
          "Computes Coulomb Fock contribution to atom's gradient.")
        .def("compute",
             [] (const CRIFockGradDriver&   grad_drv,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const std::vector<int>&    atoms) -> std::vector<TPoint<double>>
             {
                return grad_drv.compute(basis, aux_basis, molecule, gamma, density, atoms);
             },
          "Computes Coulomb Fock contribution to atoms gradient.")
        .def("compute",
             [] (const CRIFockGradDriver&   grad_drv,
                 const CT4CScreener&        screener,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const CMatrix&             density,
                 const std::vector<int>&    atoms,
                 const int                  ithreshold) -> std::vector<TPoint<double>>
             {
                return grad_drv.compute(screener, basis, aux_basis, molecule, gamma, density, atoms, ithreshold);
             },
          "Computes Coulomb Fock contribution to atoms gradient.");
    
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
