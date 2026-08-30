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


#include "ExportSparse.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "AtomBasisDiagonalSparsity.hpp"
#include "AtomBasisPairSparsity.hpp"
#include "AtomBasisTripleSparsity.hpp"
#include "ExportGeneral.hpp"
#include "MolecularBasis.hpp"
#include "ScreeningFunc.hpp"
#include "SparseMatrix.hpp"
#include "SparseTensor.hpp"

namespace vlx_sparse {  // vlx_sparse namespace

/// @brief Reconstructs the dense matrix of a sparse matrix as a numpy array.
/// @param matrix The sparse matrix to reconstruct the dense matrix of.
/// @param bra_basis The molecular basis on bra side.
/// @param ket_basis The molecular basis on ket side.
/// @param max_memory The maximum memory of the dense matrix in gigabytes.
/// @return The dense matrix.
static auto
dense_to_numpy(const CSparseMatrix   &matrix,
               const CMolecularBasis &bra_basis,
               const CMolecularBasis &ket_basis,
               const double           max_memory) -> py::array_t<double>
{
    const auto nrows = bra_basis.dimensions_of_basis();

    const auto ncols = ket_basis.dimensions_of_basis();

    // NOTE: the dense matrix is reconstructed without regard to the sparsity, so
    // its memory is checked against a limit before it is allocated. The limit is
    // enforced with an exception rather than with a critical error, as the latter
    // would terminate the Python interpreter.

    const auto memory = static_cast<double>(nrows) * static_cast<double>(ncols) * static_cast<double>(sizeof(double));

    if (const auto limit = max_memory * 1024.0 * 1024.0 * 1024.0; memory > limit)
    {
        throw std::runtime_error("SparseMatrix.to_numpy: Dense matrix of " + std::to_string(memory / (1024.0 * 1024.0 * 1024.0)) +
                                 " GB exceeds limit of " + std::to_string(max_memory) + " GB");
    }

    auto array = py::array_t<double>(std::vector<py::ssize_t>{static_cast<py::ssize_t>(nrows), static_cast<py::ssize_t>(ncols)});

    matrix.to_dense(bra_basis, ket_basis, array.mutable_data());

    return array;
}

auto
export_sparse(py::module &m) -> void
{
    // screener enumeration

    py::enum_<screener>(m, "screener")
        .value("overlap", screener::overlap)
        .value("kinetic_energy", screener::kinetic_energy)
        .value("nuclear_potential", screener::nuclear_potential)
        .value("electron_repulsion", screener::electron_repulsion);

    // diagstor enumeration

    py::enum_<diagstor>(m, "diagstor").value("scalar", diagstor::scalar).value("full", diagstor::full);

    // valstat enumeration

    py::enum_<valstat>(m, "valstat").value("empty", valstat::empty).value("allocated", valstat::allocated);

    // CAtomBasisPairSparsity class

    PyClass<CAtomBasisPairSparsity>(m, "AtomBasisPairSparsity")
        .def("bra_index", &CAtomBasisPairSparsity::bra_index, "Gets index of atom basis on bra side.")
        .def("ket_index", &CAtomBasisPairSparsity::ket_index, "Gets index of atom basis on ket side.")
        .def("bra_atoms", &CAtomBasisPairSparsity::bra_atoms, "Gets atomic indices of atom pairs on bra side.")
        .def("ket_atoms", &CAtomBasisPairSparsity::ket_atoms, "Gets atomic indices of atom pairs on ket side.")
        .def("number_of_pairs",
             py::overload_cast<>(&CAtomBasisPairSparsity::number_of_pairs, py::const_),
             "Gets number of atom pairs in sparsity pattern.")
        .def("number_of_pairs",
             py::overload_cast<const int, const size_t, const int, const size_t>(&CAtomBasisPairSparsity::number_of_pairs, py::const_),
             "Gets number of surviving atom pairs for combination of basis functions.",
             py::arg("bra_angular_momentum"),
             py::arg("bra_index"),
             py::arg("ket_angular_momentum"),
             py::arg("ket_index"))
        .def("number_of_elements",
             py::overload_cast<>(&CAtomBasisPairSparsity::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of all combinations of basis functions.")
        .def("number_of_elements",
             py::overload_cast<const int, const size_t, const int, const size_t>(&CAtomBasisPairSparsity::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of combination of basis functions.",
             py::arg("bra_angular_momentum"),
             py::arg("bra_index"),
             py::arg("ket_angular_momentum"),
             py::arg("ket_index"))
        .def("element_offset",
             &CAtomBasisPairSparsity::element_offset,
             "Gets offset of the values of combination of basis functions in the values block.",
             py::arg("bra_angular_momentum"),
             py::arg("bra_index"),
             py::arg("ket_angular_momentum"),
             py::arg("ket_index"))
        .def("number_of_bra_basis_functions",
             &CAtomBasisPairSparsity::number_of_bra_basis_functions,
             "Gets number of basis functions on bra side.")
        .def("number_of_ket_basis_functions",
             &CAtomBasisPairSparsity::number_of_ket_basis_functions,
             "Gets number of basis functions on ket side.");

    // CAtomBasisDiagonalSparsity class

    PyClass<CAtomBasisDiagonalSparsity>(m, "AtomBasisDiagonalSparsity")
        .def("bra_index", &CAtomBasisDiagonalSparsity::bra_index, "Gets index of atom basis on bra side.")
        .def("ket_index", &CAtomBasisDiagonalSparsity::ket_index, "Gets index of atom basis on ket side.")
        .def("atoms", &CAtomBasisDiagonalSparsity::atoms, "Gets atomic indices of atoms in diagonal atom pairs.")
        .def("number_of_atoms", &CAtomBasisDiagonalSparsity::number_of_atoms, "Gets number of atoms in diagonal atom pairs.")
        .def("get_storage", &CAtomBasisDiagonalSparsity::get_storage, "Gets storage layout of the diagonal atom pair blocks.")
        .def("is_triangular",
             &CAtomBasisDiagonalSparsity::is_triangular,
             "Checks if only the triangle of the combinations of basis functions is stored.")
        .def("number_of_elements",
             py::overload_cast<>(&CAtomBasisDiagonalSparsity::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of all combinations of basis functions.")
        .def("number_of_elements",
             py::overload_cast<const int, const size_t, const int, const size_t>(&CAtomBasisDiagonalSparsity::number_of_elements,
                                                                                 py::const_),
             "Gets number of values required to store the integrals of combination of basis functions.",
             py::arg("bra_angular_momentum"),
             py::arg("bra_index"),
             py::arg("ket_angular_momentum"),
             py::arg("ket_index"))
        .def("element_offset",
             &CAtomBasisDiagonalSparsity::element_offset,
             "Gets offset of the values of combination of basis functions in the values block.",
             py::arg("bra_angular_momentum"),
             py::arg("bra_index"),
             py::arg("ket_angular_momentum"),
             py::arg("ket_index"))
        .def("number_of_bra_basis_functions",
             &CAtomBasisDiagonalSparsity::number_of_bra_basis_functions,
             "Gets number of basis functions on bra side.")
        .def("number_of_ket_basis_functions",
             &CAtomBasisDiagonalSparsity::number_of_ket_basis_functions,
             "Gets number of basis functions on ket side.");

    // CAtomBasisTripleSparsity class

    PyClass<CAtomBasisTripleSparsity>(m, "AtomBasisTripleSparsity")
        .def("a_index", &CAtomBasisTripleSparsity::a_index, "Gets index of atom basis on a side.")
        .def("b_index", &CAtomBasisTripleSparsity::b_index, "Gets index of atom basis on b side.")
        .def("c_index", &CAtomBasisTripleSparsity::c_index, "Gets index of atom basis on c side.")
        .def("a_atoms", &CAtomBasisTripleSparsity::a_atoms, "Gets atomic indices of atom pairs on a side.")
        .def("b_atoms", &CAtomBasisTripleSparsity::b_atoms, "Gets atomic indices of atom pairs on b side.")
        .def("c_atoms", &CAtomBasisTripleSparsity::c_atoms, "Gets atomic indices of atoms on c side.")
        .def("number_of_pairs",
             py::overload_cast<>(&CAtomBasisTripleSparsity::number_of_pairs, py::const_),
             "Gets number of atom pairs in sparsity pattern.")
        .def("number_of_pairs",
             py::overload_cast<const int, const size_t, const int, const size_t, const int, const size_t>(
                 &CAtomBasisTripleSparsity::number_of_pairs, py::const_),
             "Gets number of surviving atom pairs for combination of basis functions.",
             py::arg("a_angular_momentum"),
             py::arg("a_index"),
             py::arg("b_angular_momentum"),
             py::arg("b_index"),
             py::arg("c_angular_momentum"),
             py::arg("c_index"))
        .def("number_of_elements",
             py::overload_cast<>(&CAtomBasisTripleSparsity::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of all combinations of basis functions.")
        .def("number_of_elements",
             py::overload_cast<const int, const size_t, const int, const size_t, const int, const size_t>(
                 &CAtomBasisTripleSparsity::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of combination of basis functions.",
             py::arg("a_angular_momentum"),
             py::arg("a_index"),
             py::arg("b_angular_momentum"),
             py::arg("b_index"),
             py::arg("c_angular_momentum"),
             py::arg("c_index"))
        .def("element_offset",
             &CAtomBasisTripleSparsity::element_offset,
             "Gets offset of the values of combination of basis functions in the values block.",
             py::arg("a_angular_momentum"),
             py::arg("a_index"),
             py::arg("b_angular_momentum"),
             py::arg("b_index"),
             py::arg("c_angular_momentum"),
             py::arg("c_index"))
        .def("number_of_c_atoms", &CAtomBasisTripleSparsity::number_of_c_atoms, "Gets number of atoms on c side.")
        .def("number_of_a_basis_functions", &CAtomBasisTripleSparsity::number_of_a_basis_functions, "Gets number of basis functions on a side.")
        .def("number_of_b_basis_functions", &CAtomBasisTripleSparsity::number_of_b_basis_functions, "Gets number of basis functions on b side.")
        .def("number_of_c_basis_functions", &CAtomBasisTripleSparsity::number_of_c_basis_functions, "Gets number of basis functions on c side.");

    // CSparseMatrix class

    PyClass<CSparseMatrix>(m, "SparseMatrix")
        .def(py::init<>())
        .def(py::init<const CMolecule &, const CMolecularBasis &, const screener, const double, const mat_t, const diagstor>(),
             "Creates a sparse matrix for a molecular basis.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("bound"),
             py::arg("threshold"),
             py::arg("matrix_type"),
             py::arg("storage"))
        .def(py::init<const CMolecule &,
                      const CMolecularBasis &,
                      const CMolecularBasis &,
                      const screener,
                      const double,
                      const mat_t,
                      const diagstor>(),
             "Creates a sparse matrix for a pair of molecular bases.",
             py::arg("molecule"),
             py::arg("bra_basis"),
             py::arg("ket_basis"),
             py::arg("bound"),
             py::arg("threshold"),
             py::arg("matrix_type"),
             py::arg("storage"))
        .def("get_type", &CSparseMatrix::get_type, "Gets type of matrix.")
        .def("set_type", &CSparseMatrix::set_type, "Sets type of matrix.", py::arg("matrix_type"))
        .def("get_values_state", &CSparseMatrix::get_values_state, "Gets state of the values blocks of matrix.")
        .def("allocate", &CSparseMatrix::allocate, "Allocates the values blocks of matrix.")
        .def("deallocate", &CSparseMatrix::deallocate, "Deallocates the values blocks of matrix.")
        .def("zero", &CSparseMatrix::zero, "Sets the values of all values blocks of matrix to zero.")
        .def("scale", &CSparseMatrix::scale, "Scales the values of all values blocks of matrix.", py::arg("factor"))
        .def("number_of_pair_blocks", &CSparseMatrix::number_of_pair_blocks, "Gets number of off-diagonal blocks in matrix.")
        .def("number_of_diagonal_blocks", &CSparseMatrix::number_of_diagonal_blocks, "Gets number of diagonal blocks in matrix.")
        .def("number_of_value_blocks", &CSparseMatrix::number_of_value_blocks, "Gets number of values blocks in matrix.")
        .def("pair_block", &CSparseMatrix::pair_block, "Gets sparsity pattern of the off-diagonal block.", py::arg("index"))
        .def("diagonal_block", &CSparseMatrix::diagonal_block, "Gets sparsity pattern of the diagonal block.", py::arg("index"))
        .def("number_of_elements",
             py::overload_cast<>(&CSparseMatrix::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of all values blocks.")
        .def("number_of_elements",
             py::overload_cast<const size_t>(&CSparseMatrix::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of values block.",
             py::arg("index"))
        .def("memory_size", &CSparseMatrix::memory_size, "Gets memory required to store the values blocks of matrix in bytes.")
        .def(
            "pair_block_to_numpy",
            [](const CSparseMatrix &self, const size_t index) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.pair_values(index),
                                                     {static_cast<int>(self.pair_block(index).number_of_elements())});
            },
            "Gets copy of the values of the off-diagonal block. The copy may be large.",
            py::arg("index"))
        .def(
            "diagonal_block_to_numpy",
            [](const CSparseMatrix &self, const size_t index) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.diagonal_values(index),
                                                     {static_cast<int>(self.diagonal_block(index).number_of_elements())});
            },
            "Gets copy of the values of the diagonal block. The copy may be large.",
            py::arg("index"))
        .def(
            "to_numpy",
            [](const CSparseMatrix &self, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis, const double max_memory) -> py::array_t<double> {
                return vlx_sparse::dense_to_numpy(self, bra_basis, ket_basis, max_memory);
            },
            "Gets the dense matrix in atomic orbital basis. The matrix may be large.",
            py::arg("bra_basis"),
            py::arg("ket_basis"),
            py::arg("max_memory") = 8.0)
        .def(
            "to_numpy",
            [](const CSparseMatrix &self, const CMolecularBasis &basis, const double max_memory) -> py::array_t<double> {
                return vlx_sparse::dense_to_numpy(self, basis, basis, max_memory);
            },
            "Gets the dense matrix in atomic orbital basis. The matrix may be large.",
            py::arg("basis"),
            py::arg("max_memory") = 8.0);

    // CSparseTensor class

    PyClass<CSparseTensor>(m, "SparseTensor")
        .def(py::init<>())
        .def(py::init<const CMolecule &,
                      const CMolecularBasis &,
                      const CMolecularBasis &,
                      const screener,
                      const double,
                      const mat_t,
                      const size_t,
                      const size_t>(),
             "Creates a sparse tensor for a batch of atoms on c side.",
             py::arg("molecule"),
             py::arg("basis"),
             py::arg("aux_basis"),
             py::arg("bound"),
             py::arg("threshold"),
             py::arg("matrix_type"),
             py::arg("natoms_per_batch"),
             py::arg("batch_index"))
        .def_static("number_of_batches",
                    &CSparseTensor::number_of_batches,
                    "Gets number of batches of atoms on c side.",
                    py::arg("aux_basis"),
                    py::arg("natoms_per_batch"))
        .def("get_type", &CSparseTensor::get_type, "Gets type of tensor.")
        .def("set_type", &CSparseTensor::set_type, "Sets type of tensor.", py::arg("matrix_type"))
        .def("get_values_state", &CSparseTensor::get_values_state, "Gets state of the values blocks of tensor.")
        .def("allocate", &CSparseTensor::allocate, "Allocates the values blocks of tensor.")
        .def("deallocate", &CSparseTensor::deallocate, "Deallocates the values blocks of tensor.")
        .def("zero", &CSparseTensor::zero, "Sets the values of all values blocks of tensor to zero.")
        .def("scale", &CSparseTensor::scale, "Scales the values of all values blocks of tensor.", py::arg("factor"))
        .def("number_of_blocks", &CSparseTensor::number_of_blocks, "Gets number of blocks in tensor.")
        .def("block", &CSparseTensor::block, "Gets sparsity pattern of the block.", py::arg("index"))
        .def("number_of_elements",
             py::overload_cast<>(&CSparseTensor::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of all blocks.")
        .def("number_of_elements",
             py::overload_cast<const size_t>(&CSparseTensor::number_of_elements, py::const_),
             "Gets number of values required to store the integrals of block.",
             py::arg("index"))
        .def("memory_size", &CSparseTensor::memory_size, "Gets memory required to store the values blocks of tensor in bytes.")
        .def(
            "block_to_numpy",
            [](const CSparseTensor &self, const size_t index) -> py::array_t<double> {
                return vlx_general::pointer_to_numpy(self.values(index), {static_cast<int>(self.block(index).number_of_elements())});
            },
            "Gets copy of the values of the block. The copy may be large.",
            py::arg("index"));
}

}  // namespace vlx_sparse
