//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "ExportResponse.hpp"

#include <mpi.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ExcitationVector.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ScreeningContainer.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_response {  // vlx_response namespace

// Exports classes/functions in src/response to python

void
export_response(py::module& m)
{
    // CExcitationVector class

    PyClass<CExcitationVector>(m, "ExcitationVector")
        .def(py::init<>())
        .def(py::init<const szblock,
                      const std::vector<int32_t>&,
                      const std::vector<int32_t>&,
                      const std::vector<double>&,
                      const std::vector<double>&>())
        .def(py::init<const szblock, const int32_t, const int32_t, const int32_t, const int32_t, const bool>())
        .def(py::init<const szblock, const int32_t, const int32_t, const int32_t, const int32_t>())
        .def(py::init<const std::vector<double>&, const std::vector<CExcitationVector>&>())
        .def(py::init<const CExcitationVector&>())
        .def("__str__", &CExcitationVector::getString)
        .def("number_excitations", &CExcitationVector::getNumberOfExcitations, "Gets number of one particle excitations in excitations vector.")
        .def("get_zmatrix", &CExcitationVector::getMatrixZ, "Transforms Z vector to matrix (occ, virt) format.")
        .def("get_ymatrix", &CExcitationVector::getMatrixY, "Transforms Y vector to matrix (virt, occ) format.")
        .def("get_zdensity", &CExcitationVector::getDensityZ, "Transforms Z vector to AO density matrix.", "molecularOrbitals"_a)
        .def("get_ydensity", &CExcitationVector::getDensityY, "Transforms Y vector to AO density matrix.", "molecularOrbitals"_a)
        .def("set_zcoefficient", &CExcitationVector::setCoefficientZ, "Sets specific element of Z coefficients vector.", "zValue"_a, "iCoefficient"_a)
        .def("set_ycoefficient", &CExcitationVector::setCoefficientY, "Sets specific element of Y coefficients vector.", "yValue"_a, "iCoefficient"_a)
        .def(
            "set_yzcoefficients",
            [](CExcitationVector& self, const std::vector<double>& z_coef, const std::vector<double>& y_coef) -> void {
                CMemBlock<double> zCoefficients(z_coef);
                CMemBlock<double> yCoefficients(y_coef);
                self.setCoefficientsZY(zCoefficients, yCoefficients);
            },
            "Sets Z and Y vectors.",
            "z_coef"_a,
            "y_coef"_a)
        .def("bra_unique_indexes",
             &CExcitationVector::getBraUniqueIndexes,
             "Gets vector of unique molecular orbitals indexes associated with creation operator in one particle excitations vector.")
        .def("ket_unique_indexes",
             &CExcitationVector::getKetUniqueIndexes,
             "Gets vector of unique molecular orbitals indexes associated with anihilation operator in one particle excitations vector.")
        .def(
            "bra_indexes",
            [](const CExcitationVector& self) -> py::array_t<int32_t> {
                return vlx_general::pointer_to_numpy(self.getBraIndexes(), self.getNumberOfExcitations());
            },
            "Gets indexes of molecular orbitals associated with anihilation operators.")
        .def(
            "ket_indexes",
            [](const CExcitationVector& self) -> py::array_t<int32_t> {
                return vlx_general::pointer_to_numpy(self.getKetIndexes(), self.getNumberOfExcitations());
            },
            "Gets indexes of molecular orbitals associated with creation operators.")
        .def("small_energy_identifiers",
             &CExcitationVector::getSmallEnergyIdentifiers,
             "Determines indexes of single particle excitation operators associated with smallest approximate excitation energies i.e. e_a - e_i.",
             "molecularOrbitals"_a,
             "nExcitations"_a)
        .def(
            "zvector_to_numpy",
            [](const CExcitationVector& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoefficientsZ(), self.getNumberOfExcitations());
            },
            "Converts Z vector to numpy array.")
        .def(
            "yvector_to_numpy",
            [](const CExcitationVector& self) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.getCoefficientsY(), self.getNumberOfExcitations());
            },
            "Converts Y vector to numpy array.")
        .def(
            "diagonal_to_numpy",
            [](const CExcitationVector& self, const CMolecularOrbitals& molecularOrbitals) -> py::array_t<double> {
                auto diagmat = self.getApproximateDiagonal(molecularOrbitals);
                return vlx_general::pointer_to_numpy(diagmat.data(), diagmat.size());
            },
            "Converts approximate diagonal of A matrix to numpy array.",
            "molecularOrbitals"_a);
}

}  // namespace vlx_response
