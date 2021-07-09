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

// Helper function for converting Z vector to numpy array

static py::array_t<double>
CExcitationVector_zvector_to_numpy(const CExcitationVector& self)
{
    return vlx_general::pointer_to_numpy(self.getCoefficientsZ(), self.getNumberOfExcitations(), 1);
}

// Helper function for converting Y vector to numpy array

static py::array_t<double>
CExcitationVector_yvector_to_numpy(const CExcitationVector& self)
{
    return vlx_general::pointer_to_numpy(self.getCoefficientsY(), self.getNumberOfExcitations(), 1);
}

// Helper function for setting Z and Y vectors

static void
CExcitationVector_set_yzcoefficients(CExcitationVector& self, const std::vector<double>& z_coef, const std::vector<double>& y_coef)
{
    CMemBlock<double> zCoefficients(z_coef);

    CMemBlock<double> yCoefficients(y_coef);

    self.setCoefficientsZY(zCoefficients, yCoefficients);
}

// Helper function for converting approximate diagonal of A matrix to numpy array

static py::array_t<double>
CExcitationVector_diagonal_to_numpy(const CExcitationVector& self, const CMolecularOrbitals& molecularOrbitals)
{
    auto diagmat = self.getApproximateDiagonal(molecularOrbitals);

    return vlx_general::pointer_to_numpy(diagmat.data(), diagmat.size(), 1);
}

static py::array_t<int32_t>
CExcitationVector_bra_indexes_to_numpy(const CExcitationVector& self)
{
    return vlx_general::pointer_to_numpy(self.getBraIndexes(), self.getNumberOfExcitations());
}

static py::array_t<int32_t>
CExcitationVector_ket_indexes_to_numpy(const CExcitationVector& self)
{
    return vlx_general::pointer_to_numpy(self.getKetIndexes(), self.getNumberOfExcitations());
}

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
        .def("set_zcoefficient", &CExcitationVector::setCoefficientZ, "Sets specific element of Z coefficients vector.", "zValue"_a, "iCoefficient"_a)
        .def("set_ycoefficient", &CExcitationVector::setCoefficientY, "Sets specific element of Y coefficients vector.", "yValue"_a, "iCoefficient"_a)
        .def("set_yzcoefficients", &CExcitationVector_set_yzcoefficients, "Sets Z and Y vectors.", "z_coef"_a, "coef"_a)
        .def("number_excitations", &CExcitationVector::getNumberOfExcitations, "Gets number of one particle excitations in excitations vector.")
        .def("bra_unique_indexes",
             &CExcitationVector::getBraUniqueIndexes,
             "Gets vector of unique molecular orbitals indexes associated with creation operator in one particle excitations vector.")
        .def("ket_unique_indexes",
             &CExcitationVector::getKetUniqueIndexes,
             "Gets vector of unique molecular orbitals indexes associated with anihilation operator in one particle excitations vector.")
        .def("bra_indexes", &CExcitationVector_bra_indexes_to_numpy, "Gets indexes of molecular orbitals associated with anihilation operators.")
        .def("ket_indexes", &CExcitationVector_ket_indexes_to_numpy, "Gets indexes of molecular orbitals associated with creation operators.")
        .def("get_zmatrix", &CExcitationVector::getMatrixZ, "Transforms Z vector to matrix (occ, virt) format.")
        .def("get_ymatrix", &CExcitationVector::getMatrixY, "Transforms Y vector to matrix (virt, occ) format.")
        .def("get_zdensity", &CExcitationVector::getDensityZ, "Transforms Z vector to AO density matrix.", "molecularOrbitals"_a)
        .def("get_ydensity", &CExcitationVector::getDensityY, "Transforms Y vector to AO density matrix.", "molecularOrbitals"_a)
        .def("small_energy_identifiers",
             &CExcitationVector::getSmallEnergyIdentifiers,
             "Determines indexes of single particle excitation operators associated with smallest approximate excitation energies i.e. e_a - e_i.",
             "molecularOrbitals"_a,
             "nExcitations"_a)
        .def("zvector_to_numpy", &CExcitationVector_zvector_to_numpy, "Converts Z vector to numpy array.")
        .def("yvector_to_numpy", &CExcitationVector_yvector_to_numpy, "Converts Y vector to numpy array.")
        .def("diagonal_to_numpy",
             &CExcitationVector_diagonal_to_numpy,
             "Converts approximate diagonal of A matrix to numpy array.",
             "molecularOrbitals"_a);
}

}  // namespace vlx_response
