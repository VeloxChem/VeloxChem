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



#include "ExportSimdIntegrals.hpp"

#include <pybind11/pybind11.h>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "SimdKineticEnergyDriver.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"
#include "SimdOverlapDriver.hpp"
#include "SparseMatrix.hpp"

namespace vlx_simdintegrals {  // vlx_simdintegrals namespace

auto
export_simdintegrals(py::module &m) -> void
{
    // CSimdOverlapDriver class

    PyClass<CSimdOverlapDriver>(m, "SimdOverlapDriver")
        .def(py::init<>())
        .def(py::init<const double>(), "Creates an overlap driver with given screening threshold.", py::arg("threshold"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &>(&CSimdOverlapDriver::compute, py::const_),
             "Computes sparse overlap matrix for given molecule and basis.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &, const CMolecularBasis &>(&CSimdOverlapDriver::compute, py::const_),
             "Computes sparse overlap matrix for given molecule and pair of bases.",
             py::arg("molecule"),
             py::arg("bra_basis"),
             py::arg("ket_basis"))
        .def("get_threshold", &CSimdOverlapDriver::get_threshold, "Gets screening threshold of the integrals.");

    // CSimdKineticEnergyDriver class

    PyClass<CSimdKineticEnergyDriver>(m, "SimdKineticEnergyDriver")
        .def(py::init<>())
        .def(py::init<const double>(), "Creates a kinetic energy driver with given screening threshold.", py::arg("threshold"))
        .def("compute",
             &CSimdKineticEnergyDriver::compute,
             "Computes sparse kinetic energy matrix for given molecule and basis.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("get_threshold", &CSimdKineticEnergyDriver::get_threshold, "Gets screening threshold of the integrals.");

    // CSimdTwoCenterElectronRepulsionDriver class

    PyClass<CSimdTwoCenterElectronRepulsionDriver>(m, "SimdTwoCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def("compute",
             &CSimdTwoCenterElectronRepulsionDriver::compute,
             "Computes two-center electron repulsion matrix for given molecule and molecular basis.",
             py::arg("molecule"),
             py::arg("basis"));
}

}  // namespace vlx_simdintegrals
