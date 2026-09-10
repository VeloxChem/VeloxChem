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
#include <pybind11/stl.h>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"
#include "SimdKineticEnergyDriver.hpp"
#include "SimdOverlapDriver.hpp"
#include "SimdRIFockDriver.hpp"
#include "SimdRIJKFockDriver.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdTwoCenterElectronRepulsionDriver.hpp"
#include "SparseMatrix.hpp"

namespace vlx_simdintegrals {  // vlx_simdintegrals namespace

auto
export_simdintegrals(py::module &m) -> void
{
    // CSimdOverlapDriver class

    PyClass<CSimdOverlapDriver>(m, "SimdOverlapDriver")
        .def(py::init<>())
        .def(py::init<const double, const size_t>(),
             "Creates an overlap driver with given screening threshold and target block size.",
             py::arg("threshold"),
             py::arg("block_size") = 0)
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &>(&CSimdOverlapDriver::compute_matrix, py::const_),
             "Computes sparse overlap matrix for given molecule and basis.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &, const CMolecularBasis &>(
                 &CSimdOverlapDriver::compute_matrix, py::const_),
             "Computes sparse overlap matrix for given molecule and pair of bases.",
             py::arg("molecule"),
             py::arg("bra_basis"),
             py::arg("ket_basis"))
        .def("get_threshold", &CSimdOverlapDriver::get_threshold, "Gets screening threshold of the integrals.")
        .def("get_block_size", &CSimdOverlapDriver::get_block_size, "Gets target number of atom pairs of a block.");

    // CSimdKineticEnergyDriver class

    PyClass<CSimdKineticEnergyDriver>(m, "SimdKineticEnergyDriver")
        .def(py::init<>())
        .def(py::init<const double, const size_t>(),
             "Creates a kinetic energy driver with given screening threshold and target block size.",
             py::arg("threshold"),
             py::arg("block_size") = 0)
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &>(&CSimdKineticEnergyDriver::compute_matrix, py::const_),
             "Computes sparse kinetic energy matrix for given molecule and basis.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("get_threshold", &CSimdKineticEnergyDriver::get_threshold, "Gets screening threshold of the integrals.")
        .def("get_block_size", &CSimdKineticEnergyDriver::get_block_size, "Gets target number of atom pairs of a block.");

    // CSimdTwoCenterElectronRepulsionDriver class

    PyClass<CSimdTwoCenterElectronRepulsionDriver>(m, "SimdTwoCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def(py::init<const size_t>(),
             "Creates a two-center electron repulsion driver with given target block size.",
             py::arg("block_size"))
        .def("compute",
             &CSimdTwoCenterElectronRepulsionDriver::compute,
             "Computes packed two-center electron repulsion matrix for given molecule and basis.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("get_block_size",
             &CSimdTwoCenterElectronRepulsionDriver::get_block_size,
             "Gets target number of atom pairs of a block.");

    // CSimdThreeCenterElectronRepulsionDriver class

    PyClass<CSimdThreeCenterElectronRepulsionDriver>(m, "SimdThreeCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def("compute",
             static_cast<CSparseTensor (CSimdThreeCenterElectronRepulsionDriver::*)(
                 const CMolecule &, const CMolecularBasis &, const CMolecularBasis &, const double) const>(
                 &CSimdThreeCenterElectronRepulsionDriver::compute),
             "Computes sparse tensor of three-center electron repulsion integrals.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"))
        .def("compute",
             static_cast<CSparseTensor (CSimdThreeCenterElectronRepulsionDriver::*)(
                 const CMolecule &, const CMolecularBasis &, const CMolecularBasis &, const double,
                 const std::vector<int> &) const>(&CSimdThreeCenterElectronRepulsionDriver::compute),
             "Computes sparse tensor of three-center electron repulsion integrals for given atoms on c side.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"),
             py::arg("atoms"));

    // CSimdRIFockDriver class

    PyClass<CSimdRIFockDriver>(m, "SimdRIFockDriver")
        .def(py::init<>())
        .def("compute_bq_vectors",
             &CSimdRIFockDriver::compute_bq_vectors,
             "Computes sparse tensor of the B vectors of the resolution of the identity.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("inverse_metric"),
             py::arg("threshold"),
             py::arg("aux_atoms") = std::vector<int>{})
        .def("compute_y_vector",
             &CSimdRIFockDriver::compute_y_vector,
             "Contracts the B vectors with a density matrix.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("density"))
        .def("compute_fock_matrix",
             static_cast<CPackedMatrix (CSimdRIFockDriver::*)(
                 const CSparseTensor &, const CMolecularBasis &, const CMolecularBasis &,
                 const CPackedMatrix &) const>(&CSimdRIFockDriver::compute_fock_matrix),
             "Computes the Coulomb matrix of the resolution of the identity for a density.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("density"))
        .def("compute_fock_matrix",
             static_cast<CPackedMatrix (CSimdRIFockDriver::*)(
                 const CSparseTensor &, const CMolecularBasis &, const CMolecularBasis &,
                 const std::vector<double> &) const>(&CSimdRIFockDriver::compute_fock_matrix),
             "Computes the Coulomb matrix of the resolution of the identity for a Y vector.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("y_vector"))
        .def("compute_w_vectors",
             static_cast<std::vector<CPackedMatrix> (CSimdRIFockDriver::*)(
                 const CSparseTensor &, const CMolecularBasis &, const CMolecularBasis &,
                 const CPackedMatrix &, const size_t, const size_t) const>(
                 &CSimdRIFockDriver::compute_w_vectors),
             "Transforms one index of the B vectors into the molecular orbitals.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("coefficients"),
             py::arg("qfirst"),
             py::arg("qlast"))
        .def("set_dense_threshold", &CSimdRIFockDriver::set_dense_threshold,
             "Sets the density at which the transformation expands the B vectors.",
             py::arg("threshold"))
        .def("get_dense_threshold", &CSimdRIFockDriver::get_dense_threshold,
             "Gets the density at which the transformation expands the B vectors.")
        .def("compute_exchange_matrix",
             &CSimdRIFockDriver::compute_exchange_matrix,
             "Adds the exchange contribution of a range of the auxiliary basis to a matrix.",
             py::arg("w_vectors"),
             py::arg("matrix"),
             py::arg("factor") = 1.0);

    // CSimdRIJKFockDriver class

    PyClass<CSimdRIJKFockDriver>(m, "SimdRIJKFockDriver")
        .def(py::init<>())
        .def("required_memory",
             &CSimdRIJKFockDriver::required_memory,
             "Gets the memory of the B vectors in bytes.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"))
        .def("prepare",
             &CSimdRIJKFockDriver::prepare,
             "Forms the inverted factor of the metric and the B vectors.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"),
             py::arg("memory_budget"),
             py::arg("metric_threshold") = 1.0e-12,
             py::arg("use_inverse_square_root") = false)
        .def("compute",
             &CSimdRIJKFockDriver::compute,
             "Computes the Fock matrix, twice the Coulomb less the scaled exchange.",
             py::arg("density"),
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"))
        .def("is_prepared", &CSimdRIJKFockDriver::is_prepared, "Checks that the driver has been prepared.")
        .def("get_bq_vectors", &CSimdRIJKFockDriver::get_bq_vectors,
             py::return_value_policy::reference_internal, "Gets the B vectors the driver holds.")
        .def("get_metric", &CSimdRIJKFockDriver::get_metric,
             py::return_value_policy::reference_internal, "Gets the inverted factor of the metric.");
}

}  // namespace vlx_simdintegrals
