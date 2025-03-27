//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#include "ExportOrbdata.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "AtomBasis.hpp"
#include "AODensityMatrix.hpp"
#include "AOIndices.hpp"
#include "BasisFunction.hpp"
#include "BlockedGtoPairBlock.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "GtoPairBlock.hpp"
#include "GtoPairBlockFunc.hpp"
#include "MolecularBasis.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_orbdata {  // vlx_orbdata namespace

// Helper function for CAODensityMatrix constructor

static std::shared_ptr<CAODensityMatrix>
CAODensityMatrix_from_numpy_list(const std::vector<py::array_t<double>> &arrays, const denmat den_type)
{
    std::vector<CDenseMatrix> dmat;

    for (size_t i = 0; i < arrays.size(); i++)
    {
        dmat.push_back(*vlx_math::CDenseMatrix_from_numpy(arrays[i]));
    }

    return std::make_shared<CAODensityMatrix>(dmat, den_type);
}

void
export_orbdata(py::module &m)
{
    // exposing functions from GtoFunc.hpp
    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &>(&gtofunc::make_gto_blocks),
          "Creates vector of basis functions blocks for given basis and "
          "molecule.");
    m.def("make_gto_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &, const std::vector<int> &>(&gtofunc::make_gto_blocks),
          "Creates vector of basis functions blocks for selected atoms in given "
          "basis and molecule.");

    m.def("make_gto_pair_blocks",
          py::overload_cast<const CMolecularBasis &, const CMolecule &>(&gtofunc::make_gto_pair_blocks),
          "Creates vector of basis function pairs blocks for given basis and molecule.");
    m.def("make_gto_pair_blocks",
          py::overload_cast<const std::vector<CGtoBlock>&>(&gtofunc::make_gto_pair_blocks),
          "Creates vector of basis function pairs blocks for given vector of basis function blocks.");
    m.def("make_gto_pair_blocks",
          py::overload_cast<const std::vector<CGtoBlock> &, const std::vector<CGtoBlock> &>(&gtofunc::make_gto_pair_blocks),
          "Creates vector of basis function pairs blocks for given pair of vectors of basis function blocks.");

    // CBasisFunction class
    PyClass<CBasisFunction>(m, "BasisFunction")
        .def(py::init<>())
        .def(py::init<const CBasisFunction &>())
        .def(py::init<const std::vector<double> &, const std::vector<double> &, const int>())
        .def(py::pickle(
            [](const CBasisFunction &bf) { return py::make_tuple(bf.get_exponents(), bf.get_normalization_factors(), bf.get_angular_momentum()); },
            [](py::tuple t) { return CBasisFunction(t[0].cast<std::vector<double>>(), t[1].cast<std::vector<double>>(), t[2].cast<int>()); }))
        .def("set_exponents", &CBasisFunction::set_exponents, "Sets exponents of basis function.")
        .def("set_normalization_factors", &CBasisFunction::set_normalization_factors, "Gets name of chemical element.")
        .def("set_angular_momentum", &CBasisFunction::set_angular_momentum, "Sets angular momentum of basis function.")
        .def("add",
             &CBasisFunction::add,
             "Add primitive i.e. exponent and normalization factor to basis "
             "function.")
        .def("normalize", &CBasisFunction::normalize, "Normalizes primitive GTOs in basis function.")
        .def("get_exponents", &CBasisFunction::get_exponents, "Gets vector of exponents in basis function.")
        .def("get_normalization_factors", &CBasisFunction::get_normalization_factors, "Gets vector of normalization factors in basis function.")
        .def("get_angular_momentum", &CBasisFunction::get_angular_momentum, "Gets angular momentum of basis function.")
        .def("number_of_primitives", &CBasisFunction::number_of_primitive_functions, "Gets number of primitives in basis function.")
        .def("__eq__", [](const CBasisFunction &self, const CBasisFunction &other) { return self == other; })
        .def("__copy__", [](const CBasisFunction &self) { return CBasisFunction(self); })
        .def("__deepcopy__", [](const CBasisFunction &self, py::dict) { return CBasisFunction(self); });

    // CAtomBasis class
    PyClass<CAtomBasis>(m, "AtomBasis")
        .def(py::init<>())
        .def(py::init<const CAtomBasis &>())
        .def(py::init<const std::vector<CBasisFunction> &, const std::string &, const std::string &, const int>())
        .def(py::pickle(
            [](const CAtomBasis &abas) {
                return py::make_tuple(abas.basis_functions(), abas.get_name(), abas.get_ecp_label(), abas.get_identifier());
            },
            [](py::tuple t) {
                return CAtomBasis(t[0].cast<std::vector<CBasisFunction>>(), t[1].cast<std::string>(), t[2].cast<std::string>(), t[3].cast<int>());
            }))
        .def("set_identifier", &CAtomBasis::set_identifier, "Sets identifier of atom basis.")
        .def("set_name", &CAtomBasis::set_name, "Sets name of atom basis.")
        .def("set_ecp_label", &CAtomBasis::set_ecp_label, "Sets effective core potential label of atom basis.")
        .def("add", &CAtomBasis::add, "Adds basis function to atom basis.")
        .def("reduce_to_valence_basis", &CAtomBasis::reduce_to_valence_basis, "Reduces atom basis to it's valence only form.")
        .def("get_basis_functions", py::overload_cast<>(&CAtomBasis::basis_functions, py::const_), "Gets GTOs.")
        .def("get_basis_functions",
             py::overload_cast<const int>(&CAtomBasis::basis_functions, py::const_),
             "Gets GTOs with specific angular momentum.")
        .def("get_basis_functions",
             py::overload_cast<const int, const size_t>(&CAtomBasis::basis_functions, py::const_),
             "Gets GTOs with specific angular momentum and number of primitives.")
        .def("get_identifier", &CAtomBasis::get_identifier, "Gets identifier of atom basis.")
        .def("get_name", &CAtomBasis::get_name, "Gets name of atom basis.")
        .def("get_ecp_label", &CAtomBasis::get_ecp_label, "Gets effective core potential label of atom basis.")
        .def("need_ecp", &CAtomBasis::need_ecp, "Checks if atom basis requires effective core potential.")
        .def("max_angular_momentum", &CAtomBasis::max_angular_momentum, "Gets maximum angular momentum in atom basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int>(&CAtomBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum.")
        .def("number_of_basis_functions",
             py::overload_cast<const int, const size_t>(&CAtomBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of "
             "primitives.")
        .def("number_of_primitive_basis_functions",
             &CAtomBasis::number_of_primitive_functions,
             "Gets number of primitive GTOs with specific angular momentum.")
        .def("contraction_depths", &CAtomBasis::contraction_depths, "Gets contraction depths of GTOs with specific angular momentum.")
        .def("contraction_str", &CAtomBasis::contraction_string, "Gets contraction string of atom basis.")
        .def("primitives_str", &CAtomBasis::primitives_string, "Gets primitive GTOs string of atom basis.")
        .def("__eq__", [](const CAtomBasis &self, const CAtomBasis &other) { return self == other; })
        .def("__copy__", [](const CAtomBasis &self) { return CAtomBasis(self); })
        .def("__deepcopy__", [](const CAtomBasis &self, py::dict) { return CAtomBasis(self); });

    // CMolecularBasis class
    PyClass<CMolecularBasis>(m, "MolecularBasis")
        .def(py::init<>())
        .def(py::init<const CMolecularBasis &>())
        .def(py::init<const std::vector<CAtomBasis> &, const std::vector<int> &>())
        .def(py::pickle([](const CMolecularBasis &mbas) { return py::make_tuple(mbas.basis_sets(), mbas.basis_sets_indices()); },
                        [](py::tuple t) { return CMolecularBasis(t[0].cast<std::vector<CAtomBasis>>(), t[1].cast<std::vector<int>>()); }))
        .def("add", &CMolecularBasis::add, "Adds atomic basis to molecular basis.")
        .def("slice", &CMolecularBasis::slice, "Slices fraction of molecular basis for specific atoms.")
        .def("reduce_to_valence_basis", &CMolecularBasis::reduce_to_valence_basis, "Reduces molecular basis to it's valence only form.")
        .def("basis_sets", &CMolecularBasis::basis_sets, "Gets unique atomic basis sets in molecular basis")
        .def("basis_sets_indices", &CMolecularBasis::basis_sets_indices, "Gets vector of basis sets indices.")
        .def("get_label", &CMolecularBasis::get_label, "Gets name of molecular basis.")
        .def("get_ao_basis_map", &CMolecularBasis::get_ao_basis_map, "Creates string representation map of basis functions.", "molecule"_a)
        .def("max_angular_momentum",
             py::overload_cast<>(&CMolecularBasis::max_angular_momentum, py::const_),
             "Gets maximum angular momentum of molecular basis.")
        .def("max_angular_momentum",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::max_angular_momentum, py::const_),
             "Gets maximum angular momentum of molecular basis for list of "
             "specific atoms.")
        .def("basis_functions", py::overload_cast<>(&CMolecularBasis::basis_functions, py::const_), "Gets vector of GTOs from molecular basis.")
        .def("basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular "
             "basis.")
        .def("basis_functions",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of "
             "primitive GTOs from "
             "molecular basis.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs from molecular basis for list of specific "
             "atoms.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum from molecular "
             "basis for list of specific atoms.")
        .def("basis_functions",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::basis_functions, py::const_),
             "Gets vector of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis for list of specific atoms.")
        .def("atomic_indices",
             py::overload_cast<>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const int>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum and number of primitive GTOs from molecular basis.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs from molecular basis for "
             "list of specific "
             "atoms.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum from molecular basis for list of specific atoms.")
        .def("atomic_indices",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::atomic_indices, py::const_),
             "Gets vector of atomic indices for GTOs with specific angular "
             "momentum and number of primitive GTOs from molecular basis for "
             "list of "
             "specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular "
             "basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum from molecular "
             "basis for list of specific atoms.")
        .def("number_of_basis_functions",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::number_of_basis_functions, py::const_),
             "Gets number of GTOs with specific angular momentum and number of "
             "primitive GTOs from molecular basis for list of specific atoms.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const int>(&CMolecularBasis::number_of_primitive_functions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from "
             "molecular basis.")
        .def("number_of_primitive_basis_functions",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::number_of_primitive_functions, py::const_),
             "Gets number of primitive GTOs with specific angular momentum from "
             "molecular basis for list of specific atoms.")
        .def("contraction_depths",
             py::overload_cast<const int>(&CMolecularBasis::contraction_depths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum "
             "from molecular basis.")
        .def("contraction_depths",
             py::overload_cast<const std::vector<int> &, const int>(&CMolecularBasis::contraction_depths, py::const_),
             "Gets contraction depths of GTOs with specific angular momentum "
             "from molecular basis for list of specific atoms.")
        .def("get_dimension_of_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets full dimensions of basis of molecular basis.")
        .def("get_dimensions_of_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets full dimensions of basis of molecular basis.")
        .def("get_dimension_of_basis",
             py::overload_cast<const int>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets partial dimensions of basis of molecular basis up to specific "
             "angular momentum.")
        .def("get_dimensions_of_basis",
             py::overload_cast<const int>(&CMolecularBasis::dimensions_of_basis, py::const_),
             "Gets partial dimensions of basis of molecular basis up to specific "
             "angular momentum.")
        .def("get_dimension_of_primitive_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_primitive_basis, py::const_),
             "Gets full dimensions of primitive basis of molecular basis.")
        .def("get_dimensions_of_primitive_basis",
             py::overload_cast<>(&CMolecularBasis::dimensions_of_primitive_basis, py::const_),
             "Gets full dimensions of primitive basis of molecular basis.")
        .def("get_index_map",
             py::overload_cast<const int, const size_t>(&CMolecularBasis::index_map, py::const_),
             "Gets compressed global index map for basis sets of specific "
             "angular "
             "momentum and number of primitive GTOs from molecular basis.")
        .def("get_index_map",
             py::overload_cast<const std::vector<int> &, const int, const size_t>(&CMolecularBasis::index_map, py::const_),
             "Gets compressed global index map for basis sets of specific "
             "angular "
             "momentum and number of primitive GTOs from molecular basis for "
             "list of specific atoms.")
        .def("get_main_basis_label", &CMolecularBasis::main_basis_label, "Gets main basis set label in molecular basis.")
        .def("__eq__", [](const CMolecularBasis &self, const CMolecularBasis &other) { return self == other; })
        .def("__copy__", [](const CMolecularBasis &self) { return CMolecularBasis(self); })
        .def("__deepcopy__", [](const CMolecularBasis &self, py::dict) { return CMolecularBasis(self); });

    // CGtoBlock class
    PyClass<CGtoBlock>(m, "GtoBlock")
        .def(py::init<>())
        .def(py::init<const CGtoBlock &>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const int, const int>())
        .def(py::init<const CMolecularBasis &, const CMolecule &, const std::vector<int> &, const int, const int>())
        .def("reduce", &CGtoBlock::reduce, "Reduces basis functions block with given exclusion mask.")
        .def("coordinates", &CGtoBlock::coordinates, "Gets vector of basis function Cartesian coordinates.")
        .def("exponents", &CGtoBlock::exponents, "Gets vector of basis function exponents.")
        .def("normalization_factors", &CGtoBlock::normalization_factors, "Gets vector of basis function normalization factors.")
        .def("orbital_indices", &CGtoBlock::orbital_indices, "Gets vector of orbital indices of GTOs.")
        .def("atomic_indices", &CGtoBlock::atomic_indices, "Gets vector of atomic indices of GTOs.")
        .def("angular_momentum", &CGtoBlock::angular_momentum, "Gets angular momentum of GTOs block.")
        .def("number_of_primitives", &CGtoBlock::number_of_primitives, "Gets number of primitive GTOs in basis function.")
        .def("number_of_basis_functions", &CGtoBlock::number_of_basis_functions, "Gets number of GTOs in basis function block.")
        .def("__eq__", [](const CGtoBlock &self, const CGtoBlock &other) { return self == other; })
        .def("__copy__", [](const CGtoBlock &self) { return CGtoBlock(self); })
        .def("__deepcopy__", [](const CGtoBlock &self, py::dict) { return CGtoBlock(self); });

    // CGtoPairBlock class
    PyClass<CGtoPairBlock>(m, "GtoPairBlock")
        .def(py::init<>())
        .def(py::init<const CGtoPairBlock &>())
        .def(py::init<const CGtoBlock &>())
        .def(py::init<const CGtoBlock &, const CGtoBlock &>())
        .def(py::init<const std::vector<TPoint<double>> &,
                      const std::vector<TPoint<double>> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<size_t> &,
                      const std::vector<size_t> &,
                      const std::vector<int> &,
                      const std::vector<int> &,
                      const std::pair<int, int> &,
                      const int>())
        .def(py::pickle(
            [](const CGtoPairBlock &gp_block) {
                return py::make_tuple(gp_block.bra_coordinates(),
                                      gp_block.ket_coordinates(),
                                      gp_block.bra_exponents(),
                                      gp_block.ket_exponents(),
                                      gp_block.normalization_factors(),
                                      gp_block.overlap_factors(),
                                      gp_block.bra_orbital_indices(),
                                      gp_block.ket_orbital_indices(),
                                      gp_block.bra_atomic_indices(),
                                      gp_block.ket_atomic_indices(),
                                      gp_block.angular_momentums(),
                                      gp_block.number_of_primitive_pairs());
            },
            [](py::tuple t) {
                auto gp_block = CGtoPairBlock(t[0].cast<const std::vector<TPoint<double>>>(),
                                              t[1].cast<const std::vector<TPoint<double>>>(),
                                              t[2].cast<const std::vector<double>>(),
                                              t[3].cast<const std::vector<double>>(),
                                              t[4].cast<const std::vector<double>>(),
                                              t[5].cast<const std::vector<double>>(),
                                              t[6].cast<const std::vector<size_t>>(),
                                              t[7].cast<const std::vector<size_t>>(),
                                              t[8].cast<const std::vector<int>>(),
                                              t[9].cast<const std::vector<int>>(),
                                              t[10].cast<const std::pair<int, int>>(),
                                              t[11].cast<const int>());
                return gp_block;
            }))
        .def("bra_coordinates", &CGtoPairBlock::bra_coordinates, "Gets vector of Cartesian coordinates of center A.")
        .def("ket_coordinates", &CGtoPairBlock::ket_coordinates, "Gets vector of Cartesian coordinates of center B.")
        .def("bra_exponents", &CGtoPairBlock::bra_exponents, "Gets vector of basis function exponents of center A.")
        .def("ket_exponents", &CGtoPairBlock::ket_exponents, "Gets vector of basis function exponents of center B.")
        .def("normalization_factors", &CGtoPairBlock::normalization_factors, "Gets vector of basis functions normalization factors.")
        .def("overlap_factors", &CGtoPairBlock::overlap_factors, "Gets vector of basis functions overlap factors.")
        .def("bra_orbital_indices", &CGtoPairBlock::bra_orbital_indices, "Gets vector of orbital indices of center A.")
        .def("ket_orbital_indices", &CGtoPairBlock::ket_orbital_indices, "Gets vector of orbital indices of center B.")
        .def("bra_atomic_indices", &CGtoPairBlock::bra_atomic_indices, "Gets vector of atomic indices of center A.")
        .def("ket_atomic_indices", &CGtoPairBlock::ket_atomic_indices, "Gets vector of atomic indices of center B.")
        .def("angular_momentums", &CGtoPairBlock::angular_momentums, "Gets angular momentums of GTOs pair.")
        .def("number_of_primitive_pairs", &CGtoPairBlock::number_of_primitive_pairs, "Gets number of primitive GTO pairs in GTO pair.")
        .def("number_of_contracted_pairs",
             &CGtoPairBlock::number_of_contracted_pairs,
             "Gets number of contracted GTO pairs in basis function pairs block.")
        .def("__eq__", [](const CGtoPairBlock &self, const CGtoPairBlock &other) { return self == other; })
        .def("__neq__", [](const CGtoPairBlock &self, const CGtoPairBlock &other) { return self != other; })
        .def("__copy__", [](const CGtoPairBlock &self) { return CGtoPairBlock(self); })
        .def("__deepcopy__", [](const CGtoPairBlock &self, py::dict) { return CGtoPairBlock(self); });

    // CBlockedGtoPairBlock class
    PyClass<CBlockedGtoPairBlock>(m, "BlockedGtoPairBlock")
        .def(py::init<>())
        .def(py::init<const CBlockedGtoPairBlock &>())
        .def(py::init<const std::array<CGtoPairBlock, 16> &>())
        .def(py::init<const std::vector<CGtoPairBlock> &, const std::vector<int> &>())
        .def(py::init<const CGtoPairBlock &, const std::vector<double> &>())
        .def(py::pickle([](const CBlockedGtoPairBlock &bgto_block) { return py::make_tuple(bgto_block.gto_pair_blocks()); },
                        [](py::tuple t) { return CBlockedGtoPairBlock(t[0].cast<const std::array<CGtoPairBlock, 16>>()); }))
        .def("gto_pair_block", &CBlockedGtoPairBlock::gto_pair_block, "Gets specific basis function pairs block.")
        .def("is_empty_gto_pair_block", &CBlockedGtoPairBlock::is_empty_gto_pair_block, "Checks if specific basis function pairs block is empty.")
        .def("gto_pair_blocks", &CBlockedGtoPairBlock::gto_pair_blocks, "Gets array of basis function pairs blocks.")
        .def("__eq__", [](const CBlockedGtoPairBlock &self, const CBlockedGtoPairBlock &other) { return self == other; })
        .def("__neq__", [](const CBlockedGtoPairBlock &self, const CBlockedGtoPairBlock &other) { return self != other; })
        .def("__copy__", [](const CBlockedGtoPairBlock &self) { return CBlockedGtoPairBlock(self); })
        .def("__deepcopy__", [](const CBlockedGtoPairBlock &self, py::dict) { return CBlockedGtoPairBlock(self); });

    // CAODensityMatrix class

    // clang-format off
    py::enum_<denmat>(m, "denmat")
        .value("rest", denmat::rest)
        .value("unrest", denmat::unrest);
    // clang-format on

    PyClass<CAODensityMatrix>(m, "AODensityMatrix")
        .def(py::init<>())
        .def(py::init<const CAODensityMatrix &>())
        .def(py::init(&CAODensityMatrix_from_numpy_list))
        .def(
            "alpha_to_numpy",
            [](const CAODensityMatrix &self, const int iDensityMatrix) -> py::array_t<double> {
                auto numRows = self.getNumberOfRows(iDensityMatrix);
                auto numCols = self.getNumberOfColumns(iDensityMatrix);
                return vlx_general::pointer_to_numpy(self.alphaDensity(iDensityMatrix), {numRows, numCols});
            },
            "Converts alpha density matrix to numpy array.",
            "i_dens"_a)
        .def(
            "beta_to_numpy",
            [](const CAODensityMatrix &self, const int iDensityMatrix) -> py::array_t<double> {
                auto numRows = self.getNumberOfRows(iDensityMatrix);
                auto numCols = self.getNumberOfColumns(iDensityMatrix);
                return vlx_general::pointer_to_numpy(self.betaDensity(iDensityMatrix), {numRows, numCols});
            },
            "Converts beta density matrix to numpy array.",
            "i_dens"_a)
        .def("number_of_density_matrices", &CAODensityMatrix::getNumberOfDensityMatrices, "Gets number of density matrices.")
        .def("get_density_type", &CAODensityMatrix::getDensityType, "Gets type of density matrix.")
        .def(py::self == py::self);

    // exposing functions

    m.def("get_dimer_ao_indices",
          &aoindices::getDimerAOIndices,
          "Gets AO indices of the two molecules in a molecular dimer.",
          "mol_1"_a,
          "mol_2"_a,
          "basis_1"_a,
          "basis_2"_a);
}

}  // namespace vlx_orbdata
