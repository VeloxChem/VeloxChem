#include "ExportT2CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ElectricDipoleMomentumDriver.hpp"
#include "ElectricDipoleMomentumGeomX00Driver.hpp"
#include "KineticEnergyDriver.hpp"
#include "KineticEnergyGeomX00Driver.hpp"
#include "KineticEnergyGeomX0YDriver.hpp"
#include "NuclearPotentialDriver.hpp"
#include "NuclearPotentialErfDriver.hpp"
#include "NuclearPotentialGeom0X0Driver.hpp"
#include "NuclearPotentialGeomX00Driver.hpp"
#include "NuclearPotentialGeomX0YDriver.hpp"
#include "NuclearPotentialGeomXY0Driver.hpp"
#include "NuclearPotentialErfGeom0X0Driver.hpp"
#include "NuclearPotentialErfGeomX00Driver.hpp"
#include "OverlapDriver.hpp"
#include "OverlapGeomX00Driver.hpp"
#include "OverlapGeomX0YDriver.hpp"

namespace vlx_t2cintegrals {

// Exports classes/functions in src/t2c_* to python

void
export_t2cintegrals(py::module& m)
{
    // COverlapDriver class
    PyClass<COverlapDriver>(m, "OverlapDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(ovl_drv.compute(basis, molecule));
            },
            "Computes overlap matrix for given molecule and basis.")
        .def(
            "compute",
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis)
                -> std::shared_ptr<CMatrix> { return std::make_shared<CMatrix>(ovl_drv.compute(bra_basis, ket_basis, molecule)); },
            "Computes overlap matrix for given molecule and pair of bases.");

    // CKineticEnergyDriver class
    PyClass<CKineticEnergyDriver>(m, "KineticEnergyDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyDriver& kin_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(kin_drv.compute(basis, molecule));
            },
            "Computes kinetic energy matrix for given molecule and basis.");

    // CNuclearPotentialDriver class
    PyClass<CNuclearPotentialDriver>(m, "NuclearPotentialDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialDriver&            npot_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                charges,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrix> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrix>(npot_drv.compute(charges, points, basis, molecule));
            },
            "Computes nuclear potential matrix for given molecule, basis and vector of external charges.")
        .def(
            "compute",
            [](const CNuclearPotentialDriver& npot_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(npot_drv.compute(basis, molecule));
            },
            "Computes nuclear potential matrix for given molecule and basis.");

    // CNuclearPotentialErfDriver class
    PyClass<CNuclearPotentialErfDriver>(m, "NuclearPotentialErfDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialErfDriver&         npot_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                charges,
               const std::vector<std::array<double, 3>>& coords,
               const std::vector<double>&                omegas) -> std::shared_ptr<CMatrix> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrix>(npot_drv.compute(charges, points, omegas, basis, molecule));
            },
            "Computes range separated nuclear potential matrix for given molecule, basis and vector of external "
            "charges.")
        .def(
            "compute",
            [](const CNuclearPotentialErfDriver&         npot_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                charges,
               const std::vector<std::array<double, 3>>& coords,
               const double                              omega) -> std::shared_ptr<CMatrix> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrix>(npot_drv.compute(charges, points, omega, basis, molecule));
            },
            "Computes range separated nuclear potential matrix for given molecule, basis and vector of external "
            "charges.");

    // CElectricDipoleMomentumDriver class
    PyClass<CElectricDipoleMomentumDriver>(m, "ElectricDipoleMomentumDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CElectricDipoleMomentumDriver& dip_drv,
               const CMolecule&                     molecule,
               const CMolecularBasis&               basis,
               const std::array<double, 3>&         origin) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(dip_drv.compute(basis, molecule, TPoint<double>(origin)));
            },
            "Computes the electric dipole momentum matrices for a given molecule and basis.");

    // TODO: Replace Electric dipole, and higher multipoles code with templated single variant

    // CNuclearPotentialGeom010Driver class
    PyClass<CNuclearPotentialGeom0X0Driver<1>>(m, "NuclearPotentialGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<1>&  geom_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                dipoles,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrices> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(dipoles, points, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.")
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes nuclear potential derivatives matrices for given molecule, basis and selected atom.");
    
    // CNuclearPotentialGeom010Driver class
    PyClass<CNuclearPotentialErfGeom0X0Driver<1>>(m, "NuclearPotentialErfGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialErfGeom0X0Driver<1>&  geom_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                dipoles,
               const std::vector<std::array<double, 3>>& coords,
               const std::vector<double>&                omegas) -> std::shared_ptr<CMatrices> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(dipoles, points, omegas, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.")
        .def(
            "compute",
            [](const CNuclearPotentialErfGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const double omega, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, omega, iatom)); },
            "Computes nuclear potential derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom020Driver class
    PyClass<CNuclearPotentialGeom0X0Driver<2>>(m, "NuclearPotentialGeom020Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<2>&  geom_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                quadrupoles,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrices> {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(quadrupoles, points, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "quadrupoles.")
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes nuclear potential derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom100Driver class
    PyClass<CNuclearPotentialGeomX00Driver<1>>(m, "NuclearPotentialGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis and selected atom.");
    
    // CNuclearPotentialErfGeom100Driver class
    PyClass<CNuclearPotentialErfGeomX00Driver<1>>(m, "NuclearPotentialErfGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialErfGeomX00Driver<1>& geom_drv, const std::vector<double> &omegas, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(omegas, basis, molecule, iatom)); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom200Driver class
    PyClass<CNuclearPotentialGeomX00Driver<2>>(m, "NuclearPotentialGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom101Driver class
    PyClass<CNuclearPotentialGeomX0YDriver<1, 1>>(m, "NuclearPotentialGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomX0YDriver<1, 1>& geom_drv,
               const CMolecule&                            molecule,
               const CMolecularBasis&                      basis,
               const int                                   iatom,
               const int                                   jatom) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom, jatom));
            },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom110Driver class
    PyClass<CNuclearPotentialGeomXY0Driver<1, 1>>(m, "NuclearPotentialGeom110Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomXY0Driver<1, 1>& geom_drv,
               const CMolecule&                            molecule,
               const CMolecularBasis&                      basis,
               const int                                   iatom,
               const int                                   jatom) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom, jatom));
            },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom100Driver class
    PyClass<COverlapGeomX00Driver<1>>(m, "OverlapGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes overlap first derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom200Driver class
    PyClass<COverlapGeomX00Driver<2>>(m, "OverlapGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes overlap second derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom101Driver class
    PyClass<COverlapGeomX0YDriver<1, 1>>(m, "OverlapGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX0YDriver<1, 1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom, const int jatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom, jatom)); },
            "Computes overlap second derivatives matrices for given molecule, basis and selected atom.");

    // CKineticEnergyGeom100Driver class
    PyClass<CKineticEnergyGeomX00Driver<1>>(m, "KineticEnergyGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes kinetic energy first derivatives matrices for given molecule, basis and selected atom.");

    // CKineticEnergyGeom200Driver class
    PyClass<CKineticEnergyGeomX00Driver<2>>(m, "KineticEnergyGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom)); },
            "Computes kinetic energy first derivatives matrices for given molecule, basis and selected atom.");

    // CKineticEnergyGeom101Driver class
    PyClass<CKineticEnergyGeomX0YDriver<1, 1>>(m, "KineticEnergyGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyGeomX0YDriver<1, 1>& geom_drv,
               const CMolecule&                         molecule,
               const CMolecularBasis&                   basis,
               const int                                iatom,
               const int                                jatom) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(geom_drv.compute(basis, molecule, iatom, jatom));
            },
            "Computes kinetic energy second derivatives matrices for given molecule, basis and selected atom.");

    // CElectricDipoleMomentumGeom100Driver class
    PyClass<CElectricDipoleMomentumGeomX00Driver<1>>(m, "ElectricDipoleMomentumGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CElectricDipoleMomentumGeomX00Driver<1>& dip_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const std::array<double, 3>&                   origin,
               const int                                      iatom) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(dip_drv.compute(basis, molecule, TPoint<double>(origin), iatom));
            },
            "Computes the electric dipole momentum derivatives matrices for a given molecule, basis and selected atom.");
}

}  // namespace vlx_t2cintegrals
