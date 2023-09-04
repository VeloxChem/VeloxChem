#include "ExportT2CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "DipoleDriver.hpp"
#include "KineticEnergyDriver.hpp"
#include "NuclearPotentialDriver.hpp"
#include "OctupoleDriver.hpp"
#include "OverlapDriver.hpp"
#include "Point.hpp"
#include "QuadrupoleDriver.hpp"

namespace vlx_t2cintegrals {  // vlx_t2cintegrals namespace

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
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(ovl_drv.compute(bra_basis, ket_basis, molecule));
            },
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
            [](const CNuclearPotentialDriver& npot_drv,
               const CMolecule&               molecule,
               const CMolecularBasis&         basis,
               const std::vector<double>&     charges,
               const std::vector<TPoint3D>&   points) -> std::shared_ptr<CMatrices> {
                return std::make_shared<CMatrices>(npot_drv.compute(basis, molecule, charges, points));
            },
            "Computes nuclear potential matrix for given molecule, basis and vector of external charges.");

    // CDipoleDriver class

    PyClass<CDipoleDriver>(m, "DipoleDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CDipoleDriver& dip_drv, const CMolecule& molecule, const CMolecularBasis& basis, const TPoint3D& point)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(dip_drv.compute(basis, molecule, point)); },
            "Computes dipole matrix for given molecule, basis and origin.");

    // CQuadrupoleDriver class

    PyClass<CQuadrupoleDriver>(m, "QuadrupoleDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CQuadrupoleDriver& quad_drv, const CMolecule& molecule, const CMolecularBasis& basis, const TPoint3D& point)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(quad_drv.compute(basis, molecule, point)); },
            "Computes quadrupole matrix for given molecule, basis and origin.");

    // COctupoleDriver class

    PyClass<COctupoleDriver>(m, "OctupoleDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COctupoleDriver& octu_drv, const CMolecule& molecule, const CMolecularBasis& basis, const TPoint3D& point)
                -> std::shared_ptr<CMatrices> { return std::make_shared<CMatrices>(octu_drv.compute(basis, molecule, point)); },
            "Computes octupole matrix for given molecule, basis and origin.");

    // ...
}

}  // namespace vlx_t2cintegrals
