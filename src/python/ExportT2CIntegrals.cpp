#include "ExportT2CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ElectricDipoleMomentumDriver.hpp"
#include "KineticEnergyDriver.hpp"
#include "NuclearPotentialDriver.hpp"
#include "NuclearPotentialErfDriver.hpp"
#include "NuclearPotentialGeom0X0Driver.hpp"
#include "OverlapDriver.hpp"
#include "OverlapGeomX00Driver.hpp"

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
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& basis)
                -> std::shared_ptr<CMatrix> { return std::make_shared<CMatrix>(ovl_drv.compute(basis, molecule)); },
            "Computes overlap matrix for given molecule and basis.")
        .def(
            "compute",
            [](const COverlapDriver&  ovl_drv,
               const CMolecule&       molecule,
               const CMolecularBasis& bra_basis,
               const CMolecularBasis& ket_basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(ovl_drv.compute(bra_basis, ket_basis, molecule));
            },
            "Computes overlap matrix for given molecule and pair of bases.");

    // CKineticEnergyDriver class
    PyClass<CKineticEnergyDriver>(m, "KineticEnergyDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyDriver& kin_drv, const CMolecule& molecule, const CMolecularBasis& basis)
                -> std::shared_ptr<CMatrix> { return std::make_shared<CMatrix>(kin_drv.compute(basis, molecule)); },
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
                std::ranges::transform(
                    coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrix>(npot_drv.compute(charges, points, basis, molecule));
            },
            "Computes nuclear potential matrix for given molecule, basis and vector of external charges.")
        .def(
            "compute",
            [](const CNuclearPotentialDriver& npot_drv, const CMolecule& molecule, const CMolecularBasis& basis)
                -> std::shared_ptr<CMatrix> { return std::make_shared<CMatrix>(npot_drv.compute(basis, molecule)); },
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
                std::ranges::transform(
                    coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
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
                std::ranges::transform(
                    coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
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
            [](const CNuclearPotentialGeom0X0Driver<1>& geom_drv,
               const CMolecule&                         molecule,
               const CMolecularBasis&                   basis,
               const std::vector<double>&               dipoles,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrices> {
                   auto points = std::vector<TPoint<double>>();
                   points.reserve(coords.size());
                   std::ranges::transform(
                       coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(dipoles, points, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.");

    // CNuclearPotentialGeom020Driver class
    PyClass<CNuclearPotentialGeom0X0Driver<2>>(m, "NuclearPotentialGeom020Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<2>& geom_drv,
               const CMolecule&                         molecule,
               const CMolecularBasis&                   basis,
               const std::vector<double>&               quadrupoles,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrices> {
                   auto points = std::vector<TPoint<double>>();
                   points.reserve(coords.size());
                   std::ranges::transform(
                       coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(quadrupoles, points, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "quadrupoles.");
    
    // COverlapGeom100Driver class
    PyClass<COverlapGeomX00Driver<1>>(m, "OverlapGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<1>& geom_drv,
               const CMolecule&                         molecule,
               const CMolecularBasis&                   basis,
               const std::vector<double>&               dipoles,
               const std::vector<std::array<double, 3>>& coords) -> std::shared_ptr<CMatrices> {
                   auto points = std::vector<TPoint<double>>();
                   points.reserve(coords.size());
                   std::ranges::transform(
                       coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return std::make_shared<CMatrices>(geom_drv.compute(dipoles, points, basis, molecule));
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.");
}

}  // namespace vlx_t2cintegrals
