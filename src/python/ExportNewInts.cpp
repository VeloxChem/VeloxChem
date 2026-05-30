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

#include "ExportNewInts.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <vector>

#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"
#include "Molecule.hpp"
#include "OverlapDriver.hpp"
#include "SparseMatrix.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_newints {  // vlx_newints namespace

auto
export_newints(py::module &m) -> void
{
    // dedicated submodule mirroring the C++ newints namespace
    auto sub = m.def_submodule("newints", "New integrals computation routines.");

    // newints::SymmetryType enum
    py::enum_<newints::SymmetryType>(sub, "SymmetryType")
        .value("symmetric", newints::SymmetryType::symmetric)
        .value("antisymmetric", newints::SymmetryType::antisymmetric)
        .value("general", newints::SymmetryType::general);

    // newints::Block struct
    PyClass<newints::Block>(sub, "Block")
        .def(py::init<>())
        .def(py::init([](const std::size_t nrows, const std::size_t ncols, const std::vector<double> &data) {
                 return newints::Block{nrows, ncols, data};
             }),
             "nrows"_a,
             "ncols"_a,
             "data"_a)
        .def_readwrite("nrows", &newints::Block::nrows, "Number of rows (2 l_a + 1).")
        .def_readwrite("ncols", &newints::Block::ncols, "Number of columns (2 l_b + 1).")
        .def_readwrite("data", &newints::Block::data, "Row-major block values.")
        .def("__eq__", [](const newints::Block &self, const newints::Block &other) { return self == other; });

    // newints::SparseMatrix class
    PyClass<newints::SparseMatrix>(sub, "SparseMatrix")
        .def(py::init<>())
        .def(py::init<const newints::SymmetryType>(), "symmetry"_a)
        .def("set_symmetry", &newints::SparseMatrix::set_symmetry, "Sets symmetry of matrix.", "symmetry"_a)
        .def("symmetry", &newints::SparseMatrix::symmetry, "Gets symmetry of matrix.")
        .def("add",
             py::overload_cast<const newints::SparseMatrix::Key &, const newints::Block &>(&newints::SparseMatrix::add),
             "Adds block at given (bra, ket) key.",
             "key"_a,
             "block"_a)
        .def("add",
             py::overload_cast<const int, const int, const newints::Block &>(&newints::SparseMatrix::add),
             "Adds block at (i, j) contracted GTO indices.",
             "i"_a,
             "j"_a,
             "block"_a)
        .def("contains", &newints::SparseMatrix::contains, "Checks if a block exists at given key.", "key"_a)
        .def("block",
             py::overload_cast<const newints::SparseMatrix::Key &>(&newints::SparseMatrix::block),
             "Gets the block at given key, or None if absent.",
             "key"_a,
             py::return_value_policy::reference_internal)
        .def("zero", &newints::SparseMatrix::zero, "Sets all block values to zero.")
        .def("number_of_blocks", &newints::SparseMatrix::number_of_blocks, "Gets number of stored blocks.")
        .def("keys", &newints::SparseMatrix::keys, "Gets keys of all stored blocks in ascending order.")
        .def("__eq__", [](const newints::SparseMatrix &self, const newints::SparseMatrix &other) { return self == other; });

    // newints::MolecularBasisOutline class
    PyClass<newints::MolecularBasisOutline>(sub, "MolecularBasisOutline")
        .def(py::init<>())
        .def(py::init<const CMolecularBasis &>(), "basis"_a)
        .def("number_of_atoms", &newints::MolecularBasisOutline::number_of_atoms, "Gets the number of atoms.")
        .def("number_of_basis_functions",
             &newints::MolecularBasisOutline::number_of_basis_functions,
             "Gets the number of contracted basis functions (shells).")
        .def("number_of_atomic_orbitals",
             &newints::MolecularBasisOutline::number_of_atomic_orbitals,
             "Gets the total angular-expanded dimension (number of atomic orbitals).")
        .def("basis_function_indices",
             &newints::MolecularBasisOutline::basis_function_indices,
             "Gets the contracted-GTO (shell) indices on a given atom.",
             "atom"_a)
        .def("indices", &newints::MolecularBasisOutline::indices, "Gets the contracted-GTO index of each shell (atom-major).")
        .def("angular_momenta", &newints::MolecularBasisOutline::angular_momenta, "Gets the angular momentum of each shell (atom-major).")
        .def("atom_indices", &newints::MolecularBasisOutline::atom_indices, "Gets the owning atom of each shell (atom-major).")
        .def("__eq__",
             [](const newints::MolecularBasisOutline &self, const newints::MolecularBasisOutline &other) { return self == other; });

    // newints::OverlapDriver class
    PyClass<newints::OverlapDriver>(sub, "OverlapDriver")
        .def(py::init<>())
        .def("compute",
             &newints::OverlapDriver::compute,
             "Computes two-center overlap matrix for given molecule and basis.",
             "molecule"_a,
             "basis"_a,
             "threshold"_a);
}

}  // namespace vlx_newints
