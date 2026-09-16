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

#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "PackedMatrix.hpp"
#include "SparseTensor.hpp"
#include "TripleSparsityPattern.hpp"
#include "SimdKineticEnergyDriver.hpp"
#include "SimdNuclearPotentialDriver.hpp"
#include "SimdOverlapDriver.hpp"
#include "SimdRIFockDriver.hpp"
#include "SimdRIJKFockDriver.hpp"
#include "SimdRIJKGradientDriver.hpp"
#include "SimdRIJKResponseDriver.hpp"
#include "SimdTwoCenterElectronRepulsionGradientDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"
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

    // CSimdNuclearPotentialDriver class

    PyClass<CSimdNuclearPotentialDriver>(m, "SimdNuclearPotentialDriver")
        .def(py::init<>())
        .def(py::init<const double, const size_t>(),
             "Creates a nuclear potential driver with given screening threshold and target block size.",
             py::arg("threshold"),
             py::arg("block_size") = 0)
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &, const std::vector<double> &,
                               const std::vector<double> &>(&CSimdNuclearPotentialDriver::compute_matrix, py::const_),
             "Computes sparse nuclear potential matrix for given molecule, basis and set of point charges. "
             "The positions of the charges are a flat array of three coordinates each, in bohr.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("charges"),
             py::arg("points"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &>(
                 &CSimdNuclearPotentialDriver::compute_matrix, py::const_),
             "Computes sparse nuclear potential matrix for given molecule and basis, the point charges being "
             "the nuclei of the molecule.",
             py::arg("molecule"),
             py::arg("basis"))
        .def("get_threshold", &CSimdNuclearPotentialDriver::get_threshold, "Gets screening threshold of the integrals.")
        .def("get_block_size", &CSimdNuclearPotentialDriver::get_block_size, "Gets target number of atom pairs of a block.");

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

    // CTripleSparsityPattern, which the gradient of the three-center integrals is
    // asked for rather than forming: it must be the one the calculation already
    // holds, so that the derivative indexes the same atom pairs as the B vectors.

    PyClass<CTripleSparsityPattern>(m, "TripleSparsityPattern")
        .def("number_of_blocks", &CTripleSparsityPattern::number_of_blocks,
             "Gets the number of blocks of the pattern.")
        .def("get_threshold", &CTripleSparsityPattern::get_threshold,
             "Gets the screening threshold the pattern was described with.");

    PyClass<CSimdThreeCenterElectronRepulsionDriver>(m, "SimdThreeCenterElectronRepulsionDriver")
        .def(py::init<>())
        .def("make_pattern",
             static_cast<CTripleSparsityPattern (CSimdThreeCenterElectronRepulsionDriver::*)(
                 const CMolecule &, const CMolecularBasis &, const CMolecularBasis &, const double) const>(
                 &CSimdThreeCenterElectronRepulsionDriver::make_pattern),
             "Creates the sparsity pattern the integrals are computed in.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"))
        .def("make_pattern",
             static_cast<CTripleSparsityPattern (CSimdThreeCenterElectronRepulsionDriver::*)(
                 const CMolecule &, const CMolecularBasis &, const CMolecularBasis &, const double,
                 const std::vector<int> &) const>(&CSimdThreeCenterElectronRepulsionDriver::make_pattern),
             "Creates it for the given atoms on the auxiliary side.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"),
             py::arg("atoms"))
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

    // the way the resolution of the identity driver forms its Fock matrices

    py::enum_<rimode>(m, "rimode")
        .value("automatic", rimode::automatic)
        .value("in_memory", rimode::in_memory)
        .value("direct", rimode::direct);

    // CSimdRIJKFockDriver class

    PyClass<CSimdRIJKFockDriver>(m, "SimdRIJKFockDriver")
        .def(py::init<>())
        .def("required_memory",
             &CSimdRIJKFockDriver::required_memory,
             "Gets the memory of the B vectors of the given auxiliary atoms in bytes, or of all of them.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"),
             py::arg("aux_atoms") = std::vector<int>{})
        .def("make_metric",
             &CSimdRIJKFockDriver::make_metric,
             "Forms the metric a way of building asks for, and answers the way it is for, "
             "which is the way asked for unless a fallback has changed it.",
             py::arg("molecule"),
             py::arg("aux_basis"),
             py::arg("metric_threshold") = 1.0e-12,
             py::arg("use_inverse_square_root") = false,
             py::arg("mode") = rimode::in_memory)
        .def("prepare",
             &CSimdRIJKFockDriver::prepare,
             "Forms the inverted factor of the metric and the B vectors.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"),
             py::arg("memory_budget"),
             py::arg("metric_threshold") = 1.0e-12,
             py::arg("use_inverse_square_root") = false,
             py::arg("mode") = rimode::automatic,
             py::arg("aux_atoms") = std::vector<int>{},
             py::arg("metric") = CPackedMatrix(),
             py::arg("min_parts") = 1)
        .def("compute",
             &CSimdRIJKFockDriver::compute,
             "Computes the Fock matrix, twice the Coulomb less the scaled exchange.",
             py::arg("density"),
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"))
        .def("compute_exchange",
             &CSimdRIJKFockDriver::compute_exchange,
             "Computes the exchange of a range of the orbitals, and the right hand side of the fitting it "
             "closes on the way. The ranks of a communicator take a range each and add what they answer.",
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"),
             py::arg("ofirst"),
             py::arg("olast"))
        .def("solve_fitting",
             &CSimdRIJKFockDriver::solve_fitting,
             "Solves the metric against the right hand side of the fitting, which has to be complete.",
             py::arg("gamma"))
        .def("compute_coulomb",
             &CSimdRIJKFockDriver::compute_coulomb,
             "Adds the Coulomb matrix of the given parts of the auxiliary basis to a matrix. The ranks of a "
             "communicator take some parts each.",
             py::arg("gamma"),
             py::arg("parts"),
             py::arg("matrix"))
        .def("aux_atom_weights",
             &CSimdRIJKFockDriver::aux_atom_weights,
             "Gets the memory of the B vectors each atom of the auxiliary basis carries, in bytes, "
             "which is what a communicator should divide by rather than the count of the atoms.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("threshold"))
        .def("number_of_aux_functions", &CSimdRIJKFockDriver::number_of_aux_functions,
             "Gets the number of auxiliary basis functions a build sweeps, which is what says whether a "
             "division over the ranks divided the work and not only the memory.")
        .def("number_of_parts", &CSimdRIJKFockDriver::number_of_parts,
             "Gets the number of parts the Coulomb pass of the direct way is divided over, which is what "
             "compute_coulomb indexes into.")
        .def("number_of_sweep_parts", &CSimdRIJKFockDriver::number_of_sweep_parts,
             "Gets the number of parts the exchange pass sweeps the auxiliary basis in. Every rank sweeps "
             "every one of them, so asking for more parts of the Coulomb pass must not add any here.")
        .def("is_prepared", &CSimdRIJKFockDriver::is_prepared, "Checks that the driver has been prepared.")
        .def("get_mode", &CSimdRIJKFockDriver::get_mode, "Gets the way the driver forms the Fock matrices.")
        .def("get_bq_vectors", &CSimdRIJKFockDriver::get_bq_vectors,
             py::return_value_policy::reference_internal, "Gets the B vectors the driver holds.")
        .def("get_metric", &CSimdRIJKFockDriver::get_metric,
             py::return_value_policy::reference_internal, "Gets the inverted factor of the metric.");

    // CSimdThreeCenterElectronRepulsionGradientDriver class

    PyClass<CSimdThreeCenterElectronRepulsionGradientDriver>(m, "SimdThreeCenterElectronRepulsionGradientDriver")
        .def(py::init<>())
        .def(py::init<const size_t>(),
             "Creates a driver with given target block size.",
             py::arg("block_size"))
        .def("compute",
             &CSimdThreeCenterElectronRepulsionGradientDriver::compute,
             "Computes the derivative of the three-center electron repulsion integrals with respect to "
             "the two atoms on bra side, over a pattern the caller already holds. The tensor carries six "
             "components an element; the derivative of the auxiliary center is the negative of the two.",
             py::arg("pattern"),
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"))
        .def("get_block_size", &CSimdThreeCenterElectronRepulsionGradientDriver::get_block_size,
             "Gets target number of atom pairs of a block.");

    // CSimdTwoCenterElectronRepulsionGradientDriver class

    PyClass<CSimdTwoCenterElectronRepulsionGradientDriver>(m, "SimdTwoCenterElectronRepulsionGradientDriver")
        .def(py::init<>())
        .def(py::init<const size_t>(),
             "Creates a gradient driver with given target block size.",
             py::arg("block_size"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &, const CPackedMatrix &,
                               const std::vector<int> &>(
                 &CSimdTwoCenterElectronRepulsionGradientDriver::compute, py::const_),
             "Computes the gradient of the two-center electron repulsion integrals of an auxiliary basis "
             "contracted with Omega, for the given atoms. No sign is applied: the term of an RI-JK "
             "gradient is the negative of what is returned.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("omega"),
             py::arg("atoms"))
        .def("compute",
             py::overload_cast<const CMolecule &, const CMolecularBasis &, const CPackedMatrix &>(
                 &CSimdTwoCenterElectronRepulsionGradientDriver::compute, py::const_),
             "Computes it for every atom of the molecule.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("omega"))
        .def("get_block_size", &CSimdTwoCenterElectronRepulsionGradientDriver::get_block_size,
             "Gets target number of atom pairs of a block.");

    // TFittedDensities, what the first phase of the gradient forms

    py::class_<TFittedDensities>(m, "FittedDensities")
        .def_readonly("coefficients", &TFittedDensities::coefficients,
                      "The fitting coefficients, one per auxiliary basis function.")
        .def_readonly("orbital_densities", &TFittedDensities::orbital_densities,
                      "The fitted densities of the occupied orbitals, one matrix per auxiliary function.")
        .def_readonly("omega", &TFittedDensities::omega,
                      "The two-index fitted density the derivative of the metric is contracted against.");

    // CSimdRIJKGradientDriver class

    PyClass<CSimdRIJKGradientDriver>(m, "SimdRIJKGradientDriver")
        .def(py::init<>())
        .def(py::init<const double, const size_t, const size_t>(),
             "Creates a gradient driver with given screening threshold, target block size and memory "
             "budget. The budget bounds the batch of auxiliary functions the transformation into the "
             "occupied orbitals holds at a time, and nothing else.",
             py::arg("threshold"),
             py::arg("block_size")    = 0,
             py::arg("memory_budget") = size_t{4} * 1024 * 1024 * 1024)
        .def("get_memory_budget", &CSimdRIJKGradientDriver::get_memory_budget,
             "Gets the memory the driver may hold, in bytes.")
        .def("compute",
             py::overload_cast<const CMolecule &,
                               const CMolecularBasis &,
                               const CMolecularBasis &,
                               const CSparseTensor &,
                               const CPackedMatrix &,
                               const CPackedMatrix &,
                               const CPackedMatrix &,
                               const double,
                               const std::vector<int> &,
                               const std::vector<int> &>(&CSimdRIJKGradientDriver::compute, py::const_),
             "Computes the Coulomb and exchange contributions to the gradient of the given atoms, from the "
             "B vectors and the metric a Fock driver holds. The rows of the atoms not asked for are zero. "
             "aux_atoms names the share of the auxiliary basis the B vectors span, which under MPI makes "
             "the gradient a partial one for the ranks to reduce.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("bq_vectors"),
             py::arg("metric"),
             py::arg("density"),
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"),
             py::arg("atoms"),
             py::arg("aux_atoms") = std::vector<int>{})
        .def("compute",
             py::overload_cast<const CMolecule &,
                               const CMolecularBasis &,
                               const CMolecularBasis &,
                               const CSparseTensor &,
                               const CPackedMatrix &,
                               const CPackedMatrix &,
                               const CPackedMatrix &,
                               const double>(&CSimdRIJKGradientDriver::compute, py::const_),
             "Computes the Coulomb and exchange contributions to the gradient of every atom of the molecule.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("bq_vectors"),
             py::arg("metric"),
             py::arg("density"),
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"))
        .def("fitted_densities",
             &CSimdRIJKGradientDriver::fitted_densities,
             "Forms the fitted densities the derivative integrals are contracted against: the fitting "
             "coefficients, the fitted densities of the occupied orbitals, and the two-index fitted "
             "density. Reads no integrals, being a transformation of the B vectors alone.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("metric"),
             py::arg("density"),
             py::arg("coefficients"),
             py::arg("exchange_scaling_factor"))
        .def("get_threshold", &CSimdRIJKGradientDriver::get_threshold,
             "Gets screening threshold of the integrals.")
        .def("get_block_size", &CSimdRIJKGradientDriver::get_block_size,
             "Gets target number of atom pairs of a block.");

    // CSimdRIJKResponseDriver class

    PyClass<CSimdRIJKResponseDriver>(m, "SimdRIJKResponseDriver")
        .def(py::init<>())
        .def(py::init<const double, const size_t, const size_t>(),
             "Creates a response driver with given threshold, block size and memory budget.",
             py::arg("threshold"),
             py::arg("block_size")    = 0,
             py::arg("memory_budget") = size_t{4} * 1024 * 1024 * 1024)
        .def("compute_exchange",
             &CSimdRIJKResponseDriver::compute_exchange,
             "Computes the exchange matrices of densities given as their factors, the density of one "
             "of them being the left factor times the transpose of its right factor. The matrices are "
             "general and carry no scaling by the fraction of exact exchange.",
             py::arg("bq_vectors"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("left"),
             py::arg("rights"))
        .def("get_threshold", &CSimdRIJKResponseDriver::get_threshold,
             "Gets screening threshold of the integrals.")
        .def("get_block_size", &CSimdRIJKResponseDriver::get_block_size,
             "Gets target number of atom pairs of a block.")
        .def("get_memory_budget", &CSimdRIJKResponseDriver::get_memory_budget,
             "Gets the memory the driver may hold, in bytes.");
}

}  // namespace vlx_simdintegrals
