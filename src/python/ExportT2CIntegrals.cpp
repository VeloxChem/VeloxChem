//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
#include "TwoCenterElectronRepulsionDriver.hpp"
#include "TwoCenterElectronRepulsionGeomX00Driver.hpp"
#include "LocalCorePotentialDriver.hpp"
#include "ProjectedCorePotentialDriver.hpp"
#include "CorePotentialDriver.hpp"
#include "ProjectedCorePotentialGeomX00Driver.hpp"
#include "ProjectedCorePotentialGeom0X0Driver.hpp"
#include "ProjectedCorePotentialGeomXY0Driver.hpp"
#include "ProjectedCorePotentialGeomX0YDriver.hpp"
#include "LocalCorePotentialGeomX00Driver.hpp"
#include "LocalCorePotentialGeomX0YDriver.hpp"
#include "LocalCorePotentialGeom0X0Driver.hpp"
#include "LocalCorePotentialGeomXY0Driver.hpp"
#include "CorePotentialGradientDriver.hpp"
#include "CorePotentialHessianDriver.hpp"

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
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> CMatrix {
                return ovl_drv.compute(basis, molecule);
            },
            "Computes overlap matrix for given molecule and basis.")
        .def(
            "compute",
            [](const COverlapDriver& ovl_drv, const CMolecule& molecule, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis)
                -> CMatrix { return ovl_drv.compute(bra_basis, ket_basis, molecule); },
            "Computes overlap matrix for given molecule and pair of bases.")
        .def(
            "compute",
            [](const COverlapDriver& ovl_drv, const CMolecule& bra_molecule,  const CMolecule& ket_molecule, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis)
                -> CMatrix { return ovl_drv.compute(bra_basis, ket_basis, bra_molecule, ket_molecule); },
            "Computes overlap matrix for given pair of molecules and pair of bases.");

    // CKineticEnergyDriver class
    PyClass<CKineticEnergyDriver>(m, "KineticEnergyDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyDriver& kin_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> CMatrix {
                return kin_drv.compute(basis, molecule);
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
               const std::vector<std::array<double, 3>>& coords) -> CMatrix {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return npot_drv.compute(charges, points, basis, molecule);
            },
            "Computes nuclear potential matrix for given molecule, basis and vector of external charges.")
        .def(
            "compute",
            [](const CNuclearPotentialDriver& npot_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> CMatrix {
                return npot_drv.compute(basis, molecule);
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
               const std::vector<double>&                omegas) -> CMatrix {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return npot_drv.compute(charges, points, omegas, basis, molecule);
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
               const double                              omega) -> CMatrix {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return npot_drv.compute(charges, points, omega, basis, molecule);
            },
            "Computes range separated nuclear potential matrix for given molecule, basis and vector of external "
            "charges.");

    // CElectricDipoleMomentumDriver class
    // TODO: rename CElectricDipoleMomentumDriver to CElectricDipoleMomentDriver
    PyClass<CElectricDipoleMomentumDriver>(m, "ElectricDipoleMomentDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CElectricDipoleMomentumDriver& dip_drv,
               const CMolecule&                     molecule,
               const CMolecularBasis&               basis,
               const std::array<double, 3>&         origin) -> CMatrices {
                return dip_drv.compute(basis, molecule, TPoint<double>(origin));
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
               const std::vector<std::array<double, 3>>& coords) -> CMatrices {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return geom_drv.compute(dipoles, points, basis, molecule);
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.")
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes nuclear potential derivatives matrices for given molecule, basis and selected atom.")
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom, const double alternative_charge)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom, alternative_charge); },
            "Computes nuclear potential derivatives matrices for given molecule, basis, selected atom and alternative nuclear charge (e.g. with ECP).");
    
    // CNuclearPotentialErfGeom010Driver class
    PyClass<CNuclearPotentialErfGeom0X0Driver<1>>(m, "NuclearPotentialErfGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialErfGeom0X0Driver<1>&  geom_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                dipoles,
               const std::vector<std::array<double, 3>>& coords,
               const std::vector<double>&                omegas) -> CMatrices {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return geom_drv.compute(dipoles, points, omegas, basis, molecule);
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "dipoles.")
        .def(
            "compute",
            [](const CNuclearPotentialErfGeom0X0Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const double omega, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, omega, iatom); },
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
               const std::vector<std::array<double, 3>>& coords) -> CMatrices {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                return geom_drv.compute(quadrupoles, points, basis, molecule);
            },
            "Computes nuclear potential derivatives matrices for given molecule, basis and vector of external "
            "quadrupoles.")
        .def(
            "compute",
            [](const CNuclearPotentialGeom0X0Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes nuclear potential derivatives matrices for given molecule, basis and selected atom.");

    // CNuclearPotentialGeom100Driver class
    PyClass<CNuclearPotentialGeomX00Driver<1>>(m, "NuclearPotentialGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis and selected atom.")
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<1>& geom_drv,
               const CMolecule& molecule,
               const CMolecularBasis& basis,
               const int iatom,
               const std::vector<std::array<double, 3>>& coords_array,
               const std::vector<double>& charges) -> CMatrices { return geom_drv.compute(basis, molecule, iatom, coords_array, charges); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis, selected atom, and point charges.");
    
    // CNuclearPotentialErfGeom100Driver class
    PyClass<CNuclearPotentialErfGeomX00Driver<1>>(m, "NuclearPotentialErfGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialErfGeomX00Driver<1>& geom_drv, const std::vector<double> &omegas, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(omegas, basis, molecule, iatom); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis and selected atom.")
        .def(
            "compute",
            [](const CNuclearPotentialErfGeomX00Driver<1>& geom_drv,
               const CMolecule& molecule,
               const CMolecularBasis& basis,
               const int iatom,
               const std::vector<std::array<double, 3>>& coords_array,
               const std::vector<double>& charges,
               const std::vector<double>& omegas) -> CMatrices { return geom_drv.compute(basis, molecule, iatom, coords_array, charges, omegas); },
            "Computes nuclear potential first derivatives matrices for given molecule, basis, selected atom, and point charges.");

    // CNuclearPotentialGeom200Driver class
    PyClass<CNuclearPotentialGeomX00Driver<2>>(m, "NuclearPotentialGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.")
        .def(
            "compute",
            [](const CNuclearPotentialGeomX00Driver<2>& geom_drv,
               const CMolecule& molecule,
               const CMolecularBasis& basis,
               const int iatom,
               const std::vector<std::array<double, 3>>& coordinates,
               const std::vector<double>& charges) -> CMatrices {
                return geom_drv.compute(basis, molecule, iatom, coordinates, charges);
            },
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
               const int                                   jatom) -> CMatrices {
                return geom_drv.compute(basis, molecule, iatom, jatom);
            },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.")
        .def(
            "compute",
            [](const CNuclearPotentialGeomX0YDriver<1, 1>& geom_drv,
               const CMolecule&                            molecule,
               const CMolecularBasis&                      basis,
               const int                                   iatom,
               const int                                   jatom,
               const std::vector<std::array<double, 3>>&   coordinates,
               const std::vector<double>&                  charges) -> CMatrices {
                return geom_drv.compute(basis, molecule, iatom, jatom, coordinates, charges);
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
               const int                                   jatom) -> CMatrices {
                return geom_drv.compute(basis, molecule, iatom, jatom);
            },
            "Computes nuclear potential second derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom100Driver class
    PyClass<COverlapGeomX00Driver<1>>(m, "OverlapGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes overlap first derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom200Driver class
    PyClass<COverlapGeomX00Driver<2>>(m, "OverlapGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes overlap second derivatives matrices for given molecule, basis and selected atom.");

    // COverlapGeom101Driver class
    PyClass<COverlapGeomX0YDriver<1, 1>>(m, "OverlapGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const COverlapGeomX0YDriver<1, 1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom, const int jatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom, jatom); },
            "Computes overlap second derivatives matrices for given molecule, basis and selected atom.");

    // CKineticEnergyGeom100Driver class
    PyClass<CKineticEnergyGeomX00Driver<1>>(m, "KineticEnergyGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes kinetic energy first derivatives matrices for given molecule, basis and selected atom.");

    // CKineticEnergyGeom200Driver class
    PyClass<CKineticEnergyGeomX00Driver<2>>(m, "KineticEnergyGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CKineticEnergyGeomX00Driver<2>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
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
               const int                                jatom) -> CMatrices {
                return geom_drv.compute(basis, molecule, iatom, jatom);
            },
            "Computes kinetic energy second derivatives matrices for given molecule, basis and selected atom.");

    // CElectricDipoleMomentumGeom100Driver class
    // TODO: rename CElectricDipoleMomentumGeom100Driver to CElectricDipoleMomentGeom100Driver
    PyClass<CElectricDipoleMomentumGeomX00Driver<1>>(m, "ElectricDipoleMomentGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CElectricDipoleMomentumGeomX00Driver<1>& dip_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const std::array<double, 3>&                   origin,
               const int                                      iatom) -> CMatrices {
                return dip_drv.compute(basis, molecule, TPoint<double>(origin), iatom);
            },
            "Computes the electric dipole momentum derivatives matrices for a given molecule, basis and selected atom.");
    
    // CTwoCenterElectronRepulsionDriver class
    PyClass<CTwoCenterElectronRepulsionDriver>(m, "TwoCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CTwoCenterElectronRepulsionDriver& eri_drv, const CMolecule& molecule, const CMolecularBasis& basis) -> CMatrix {
                return eri_drv.compute(basis, molecule);
            },
            "Computes electron repulsion matrix for given molecule and basis.");
    
    // COverlapGeom100Driver class
    PyClass<CTwoCenterElectronRepulsionGeomX00Driver<1>>(m, "TwoCenterElectronRepulsionGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CTwoCenterElectronRepulsionGeomX00Driver<1>& geom_drv, const CMolecule& molecule, const CMolecularBasis& basis, const int iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, iatom); },
            "Computes overlap first derivatives matrices for given molecule, basis and selected atom.");
    
    // CLocalCorePotentialDriver class
    PyClass<CLocalCorePotentialDriver>(m, "LocalECPDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialDriver& ecp_drv,
               const CMolecule&                 molecule,
               const CMolecularBasis&           basis,
               const CBaseCorePotential&        ecp_potential,
               const int                        iatom) -> CMatrix {
               return ecp_drv.compute(basis, molecule, ecp_potential, iatom);
            },
             "Computes local ECP matrix for given molecule, basis, base core potential.");
    
    // CProjectedCorePotentialDriver class
    PyClass<CProjectedCorePotentialDriver>(m, "ProjectedECPDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialDriver& ecp_drv,
               const CMolecule&                     molecule,
               const CMolecularBasis&               basis,
               const CBaseCorePotential&            ecp_potential,
               const int                            ecp_momentum,
               const int                            iatom) -> CMatrix {
               return ecp_drv.compute(basis, molecule, ecp_potential, ecp_momentum, iatom);
            },
             "Computes projected ECP matrix for given molecule, basis, base core potential.");
    
    // CCorePotentialDriver class
    PyClass<CCorePotentialDriver>(m, "ECPDriver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CCorePotentialDriver& ecp_drv,
               const CMolecule&            molecule,
               const CMolecularBasis&      basis,
               const CAtomCorePotential&   ecp_potential) -> CMatrix {
               return ecp_drv.compute(basis, molecule, ecp_potential);
            },
             "Computes ECP matrix for given molecule, basis, atom core potential.")
        .def(
            "compute",
            [](const CCorePotentialDriver& ecp_drv,
               const CMolecule&            molecule,
               const CMolecularBasis&      basis,
               const std::vector<int>&     atoms) -> CMatrix {
               return ecp_drv.compute(basis, molecule, atoms);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.");
    
    // CProjectedCorePotentialGeom100Driver class
    PyClass<CProjectedCorePotentialGeomX00Driver<1>>(m, "ProjectedECPGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeomX00Driver<1>& geom_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const CBaseCorePotential&                      ecp_potential,
               const int                                      ecp_momentum,
               const int                                      iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atom.");
    
    // CProjectedCorePotentialGeom200Driver class
    PyClass<CProjectedCorePotentialGeomX00Driver<2>>(m, "ProjectedECPGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeomX00Driver<2>& geom_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const CBaseCorePotential&                      ecp_potential,
               const int                                      ecp_momentum,
               const int                                      iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atom.");
    
    // CLocalCorePotentialGeom100Driver class
    PyClass<CLocalCorePotentialGeomX00Driver<1>>(m, "LocalECPGeom100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeomX00Driver<1>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential,
               const int                                  iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atom.");
    
    // CLocalCorePotentialGeom200Driver class
    PyClass<CLocalCorePotentialGeomX00Driver<2>>(m, "LocalECPGeom200Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeomX00Driver<2>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential,
               const int                                  iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atom.");
    
    // CLocalCorePotentialGeom101Driver class
    PyClass<CLocalCorePotentialGeomX0YDriver<1, 1>>(m, "LocalECPGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeomX0YDriver<1, 1>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential,
               const int                                  iatom,
               const int                                  jatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, iatom, jatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atoms.");
    
    // CLocalCorePotentialGeom010Driver class
    PyClass<CLocalCorePotentialGeom0X0Driver<1>>(m, "LocalECPGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeom0X0Driver<1>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential); },
            "Computes projected core potential first derivatives matrices for given molecule, basis.");
    
    // CProjectedCorePotentialGeom010Driver class
    PyClass<CProjectedCorePotentialGeom0X0Driver<1>>(m, "ProjectedECPGeom010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeom0X0Driver<1>& geom_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const CBaseCorePotential&                      ecp_potential,
               const int                                      ecp_momentum)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum); },
            "Computes projected core potential first derivatives matrices for given molecule, basis.");
    
    // CProjectedCorePotentialGeom020Driver class
    PyClass<CProjectedCorePotentialGeom0X0Driver<2>>(m, "ProjectedECPGeom020Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeom0X0Driver<2>& geom_drv,
               const CMolecule&                               molecule,
               const CMolecularBasis&                         basis,
               const CBaseCorePotential&                      ecp_potential,
               const int                                      ecp_momentum)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum); },
            "Computes projected core potential first derivatives matrices for given molecule, basis.");
    
    // CLocalCorePotentialGeom020Driver class
    PyClass<CLocalCorePotentialGeom0X0Driver<2>>(m, "LocalECPGeom020Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeom0X0Driver<2>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential); },
            "Computes projected core potential first derivatives matrices for given molecule, basis.");
    
    // CLocalCorePotentialGeom110Driver class
    PyClass<CLocalCorePotentialGeomXY0Driver<1, 1>>(m, "LocalECPGeom110Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CLocalCorePotentialGeomXY0Driver<1, 1>& geom_drv,
               const CMolecule&                              molecule,
               const CMolecularBasis&                        basis,
               const CBaseCorePotential&                     ecp_potential,
               const int                                     iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atoms.");
    
    // CProjectedCorePotentialGeom110Driver class
    PyClass<CProjectedCorePotentialGeomXY0Driver<1, 1>>(m, "ProjectedECPGeom110Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeomXY0Driver<1, 1>& geom_drv,
               const CMolecule&                                  molecule,
               const CMolecularBasis&                            basis,
               const CBaseCorePotential&                         ecp_potential,
               const int                                         ecp_momentum,
               const int                                         iatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum, iatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atoms.");
    
    // CProjectedCorePotentialGeom101Driver class
    PyClass<CProjectedCorePotentialGeomX0YDriver<1, 1>>(m, "ProjectedECPGeom101Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CProjectedCorePotentialGeomX0YDriver<1, 1>& geom_drv,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const CBaseCorePotential&                  ecp_potential,
               const int                                  ecp_momentum,
               const int                                  iatom,
               const int                                  jatom)
                -> CMatrices { return geom_drv.compute(basis, molecule, ecp_potential, ecp_momentum, iatom, jatom); },
            "Computes projected core potential first derivatives matrices for given molecule, basis and selected atoms.");
    
    // CCorePotentialGradientDriver class
    PyClass<CCorePotentialGradientDriver>(m, "ECPGradientDriver")
        .def(py::init<>())
        .def(
            "compute_bra_grad",
            [](const CCorePotentialGradientDriver& ecp_drv,
               const CMolecule&            molecule,
               const CMolecularBasis&      basis,
               const CAtomCorePotential&   ecp_potential,
               const int                   iatom) -> CMatrices {
               return ecp_drv.compute_bra_grad(basis, molecule, ecp_potential, iatom);
            },
             "Computes ECP matrix for given molecule, basis, atom core potential.")
        .def(
            "compute_bra_grad",
            [](const CCorePotentialGradientDriver& ecp_drv,
               const CMolecule&                    molecule,
               const CMolecularBasis&              basis,
               const std::vector<int>&             atoms,
               const int                           iatom) -> CMatrices {
               return ecp_drv.compute_bra_grad(basis, molecule, atoms, iatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.")
        .def(
            "compute_pot_grad",
            [](const CCorePotentialGradientDriver& ecp_drv,
               const CMolecule&            molecule,
               const CMolecularBasis&      basis,
               const CAtomCorePotential&   ecp_potential) -> CMatrices {
               return ecp_drv.compute_pot_grad(basis, molecule, ecp_potential);
            },
             "Computes ECP matrix for given molecule, basis, atom core potential.")
        .def(
            "compute_pot_grad",
            [](const CCorePotentialGradientDriver& ecp_drv,
               const CMolecule&                    molecule,
               const CMolecularBasis&              basis,
               const int                           iatom) -> CMatrices {
               return ecp_drv.compute_pot_grad(basis, molecule, iatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.");
    
    // CCorePotentialHessianDriver class
    PyClass<CCorePotentialHessianDriver>(m, "ECPHessianDriver")
        .def(py::init<>())
        .def(
            "compute_geom_200",
            [](const CCorePotentialHessianDriver& ecp_drv,
               const CMolecule&                   molecule,
               const CMolecularBasis&             basis,
               const std::vector<int>&            atoms,
               const int                          iatom) -> CMatrices {
               return ecp_drv.compute_geom_200(basis, molecule, atoms, iatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.")
        .def(
            "compute_geom_101",
            [](const CCorePotentialHessianDriver& ecp_drv,
               const CMolecule&                   molecule,
               const CMolecularBasis&             basis,
               const std::vector<int>&            atoms,
               const int                          iatom,
               const int                          jatom) -> CMatrices {
               return ecp_drv.compute_geom_101(basis, molecule, atoms, iatom, jatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.")
        .def(
            "compute_geom_110",
            [](const CCorePotentialHessianDriver& ecp_drv,
               const CMolecule&                   molecule,
               const CMolecularBasis&             basis,
               const int                          iatom,
               const int                          jatom) -> CMatrices {
               return ecp_drv.compute_geom_110(basis, molecule, iatom, jatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.")
        .def(
            "compute_geom_020",
            [](const CCorePotentialHessianDriver& ecp_drv,
               const CMolecule&                   molecule,
               const CMolecularBasis&             basis,
               const int                          iatom) -> CMatrices {
               return ecp_drv.compute_geom_020(basis, molecule, iatom);
            },
             "Computes ECP matrix for given molecule, basis, vector of atoms.");
        
}

}  // namespace vlx_t2cintegrals
