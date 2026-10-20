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

#include "ExportT3CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#ifdef VLX_USE_OLD_RI
#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
#include "ThreeCenterElectronRepulsionGeom0X0Driver.hpp"
#endif
#include "ThreeCenterOverlapDriver.hpp"
#include "T3FlatBuffer.hpp"
#include "T4CScreener.hpp"

namespace vlx_t3cintegrals {

// Exports classes/functions in src/t3c_* to python

void
export_t3cintegrals(py::module& m)
{
#ifdef VLX_USE_OLD_RI
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
            "Computes electron repulsion integrals for given molecule, basis, auxilary basis, and list of atoms.")
        .def(
            "compute",
            [](const CThreeCenterElectronRepulsionDriver& eri_drv, const CMolecule& molecule, const CT4CScreener& screener, const CMolecularBasis& aux_basis, const std::vector<int>& atoms, const int ithreshold) -> CT3FlatBuffer<double> {
                return eri_drv.compute(screener, aux_basis, molecule, atoms, ithreshold);
            },
            "Computes electron repulsion integrals for given molecule, screener, auxilary basis, and list of atoms.");
#endif
    
    
    // CThreeCenterOverlapDriver class
    PyClass<CThreeCenterOverlapDriver>(m, "ThreeCenterOverlapDriver")
        .def(py::init<>())
        .def(
            "compute",
             [](const CThreeCenterOverlapDriver&         t3ovl_drv,
               const CMolecule&                          molecule,
               const CMolecularBasis&                    basis,
               const std::vector<double>&                exponents,
               const std::vector<double>&                factors,
               const std::vector<std::array<double, 3>>& coords) -> CMatrix {
                auto points = std::vector<TPoint<double>>();
                points.reserve(coords.size());
                std::ranges::transform(coords, std::back_inserter(points), [](auto rxyz) { return TPoint<double>(rxyz); });
                   return t3ovl_drv.compute(exponents, factors, points, basis, molecule);
            },
            "Computes overlap matrix for given molecule, basis and vector of external scaled Gaussians.");

#ifdef VLX_USE_OLD_RI
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
#endif
}

}  // namespace vlx_t2cintegrals
