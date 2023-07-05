#include "ExportT2CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "QuadrupoleDriver.hpp"
#include "DipoleDriver.hpp"
#include "KineticEnergyDriver.hpp"
#include "OverlapDriver.hpp"
#include "Point.hpp"

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
            "Computes overlap matrix for given molecule and basis.");

    // CKineticEnergyDriver class

    PyClass<CKineticEnergyDriver>(m, "KineticEnergyDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyDriver& kin_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> std::shared_ptr<CMatrix> {
                return std::make_shared<CMatrix>(kin_drv.compute(basis, molecule));
            },
            "Computes kinetic energy matrix for given molecule and basis.");

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
            "Computes dipole matrix for given molecule, basis and origin.");

    // ...
}

}  // namespace vlx_t2cintegrals
