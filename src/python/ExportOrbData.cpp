//
//                           VELOXCHEM 1.0-RC
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

#include "ExportOrbData.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi.h>

#include <fstream>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOIndices.hpp"
#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "GtoTransform.hpp"
#include "MolecularBasis.hpp"
#include "MolecularOrbitals.hpp"
#include "MolecularOrbitalsType.hpp"
#include "SADGuessDriver.hpp"
#include "StringFormat.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_orbdata {  // vlx_orbdata namespace

// Helper function for broadcasting CMolecularBasis object

static void
CMolecularBasis_broadcast(CMolecularBasis& self, int32_t rank, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, comm);
}

// Helper function for converting CAODensityMatrix to numpy array

static py::array_t<double>
CAODensityMatrix_alpha_density_to_numpy(const CAODensityMatrix& self, const int32_t iDensityMatrix)
{
    return vlx_general::pointer_to_numpy(
        self.alphaDensity(iDensityMatrix), self.getNumberOfRows(iDensityMatrix), self.getNumberOfColumns(iDensityMatrix));
}

static py::array_t<double>
CAODensityMatrix_beta_density_to_numpy(const CAODensityMatrix& self, const int32_t iDensityMatrix)
{
    return vlx_general::pointer_to_numpy(
        self.betaDensity(iDensityMatrix), self.getNumberOfRows(iDensityMatrix), self.getNumberOfColumns(iDensityMatrix));
}

// Helper function for CAODensityMatrix constructor

static std::shared_ptr<CAODensityMatrix>
CAODensityMatrix_from_numpy_list(const std::vector<py::array_t<double>>& arrays, const denmat den_type)
{
    std::vector<CDenseMatrix> dmat;

    for (size_t i = 0; i < arrays.size(); i++)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(arrays[i]);

        dmat.push_back(*mp);
    }

    return std::make_shared<CAODensityMatrix>(dmat, den_type);
}

// Helper function for broadcasting CAODensityMatrix object

static void
CAODensityMatrix_broadcast(CAODensityMatrix& self, int32_t rank, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, comm);
}

// Helper function for converting CMolecularOrbitals to numpy array

static py::array_t<double>
CMolecularOrbitals_alpha_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.alphaOrbitals(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CMolecularOrbitals_beta_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.betaOrbitals(), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CMolecularOrbitals_alpha_energies_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.alphaEnergies(), self.getNumberOfColumns());
}

static py::array_t<double>
CMolecularOrbitals_beta_energies_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.betaEnergies(), self.getNumberOfColumns());
}

// Helper function for CMolecularOrbitals constructor

static std::shared_ptr<CMolecularOrbitals>
CMolecularOrbitals_from_numpy_list(const std::vector<py::array_t<double>>& mol_orbs,
                                   const std::vector<py::array_t<double>>& eig_vals,
                                   const molorb                            orbs_type)
{
    std::vector<CDenseMatrix> cmos;

    for (size_t i = 0; i < mol_orbs.size(); i++)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(mol_orbs[i]);

        cmos.push_back(*mp);
    }

    std::vector<CMemBlock<double>> ceigs;

    for (size_t i = 0; i < eig_vals.size(); i++)
    {
        const py::array_t<double>& arr = eig_vals[i];

        std::string errdim("MolecularOrbitals: Expecting 1D numpy arrays for eigenvalues");

        errors::assertMsgCritical(arr.ndim() == 1, errdim);

        if (arr.data() == nullptr || arr.size() == 0)
        {
            return std::make_shared<CMolecularOrbitals>();
        }

        std::vector<double> vec(arr.data(), arr.data() + arr.size());

        ceigs.push_back(CMemBlock<double>(vec));
    }

    return std::make_shared<CMolecularOrbitals>(cmos, ceigs, orbs_type);
}

// Helper function for broadcasting CMolecularOrbitals object

static void
CMolecularOrbitals_broadcast(CMolecularOrbitals& self, int32_t rank, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, comm);
}

// Exports classes/functions in src/orbdata to python

void
export_orbdata(py::module& m)
{
    // CBasisFunction class

    PyClass<CBasisFunction>(m, "BasisFunction")
        .def(py::init<>())
        .def(py::init<const std::vector<double>&, const std::vector<double>&, const int32_t>())
        .def("normalize", &CBasisFunction::normalize)
        .def_property("angular_momentum", &CBasisFunction::getAngularMomentum, &CBasisFunction::setAngularMomentum)
        .def_property("exponents", &CBasisFunction::getExponents, &CBasisFunction::setExponents)
        .def_property("normalization_factors", &CBasisFunction::getNormalizationFactors, &CBasisFunction::setNormalizationFactors)
        .def_property_readonly("n_primitives", &CBasisFunction::getNumberOfPrimitiveFunctions)
        .def("__repr__", &CBasisFunction::repr);

    // CAtomBasis class

    PyClass<CAtomBasis>(m, "AtomBasis")
        .def(py::init<>())
        .def("add_basis_function", &CAtomBasis::addBasisFunction)
        .def("set_elemental_id", &CAtomBasis::setIdElemental)
        .def("__repr__", &CAtomBasis::repr)
        .def(
            "__iter__",
            [](const CAtomBasis& obj) { return py::make_iterator(obj.begin(), obj.end()); },
            py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

    // CMolecularBasis class

    PyClass<CMolecularBasis>(m, "MolecularBasis")
        .def(py::init<>())
        .def("get_string", vlx_general::overload_cast_<const std::string&, const CMolecule&>()(&CMolecularBasis::printBasis, py::const_), "title"_a, "molecule"_a)
        .def("get_string", vlx_general::overload_cast_<const CMolecule&>()(&CMolecularBasis::printBasis, py::const_), "molecule"_a)
        .def("set_label", &CMolecularBasis::setLabel)
        .def("get_label", &CMolecularBasis::getLabel)
        .def("get_ao_basis_map", &CMolecularBasis::getAOBasisMap)
        .def("broadcast", &CMolecularBasis_broadcast)
        .def("get_valence_basis", &CMolecularBasis::reduceToValenceBasis)
        .def("add_atom_basis", &CMolecularBasis::addAtomBasis)
        .def("get_dimensions_of_basis", &CMolecularBasis::getDimensionsOfBasis)
        .def("get_dimensions_of_primitive_basis", &CMolecularBasis::getDimensionsOfPrimitiveBasis)
        .def(py::self == py::self)
        .def("__repr__", &CMolecularBasis::repr)
        .def("n_basis_functions",
             vlx_general::overload_cast_<const CMolecule&, int32_t>()(&CMolecularBasis::getNumberOfBasisFunctions, py::const_),
             "molecule"_a,
             "ang_mom"_a)
        .def("n_primitive_basis_functions",
             vlx_general::overload_cast_<const CMolecule&, int32_t>()(&CMolecularBasis::getNumberOfPrimitiveBasisFunctions, py::const_),
             "molecule"_a,
             "ang_mom"_a)
        .def(
            "__iter__",
            [](const CMolecularBasis& obj) { return py::make_iterator(obj.begin(), obj.end()); },
            py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

    // denmat enum class
    // clang-format off
    py::enum_<denmat>(m, "denmat")
        .value("rest", denmat::rest)
        .value("unrest", denmat::unrest)
        .value("osrest", denmat::osrest)
        .value("rmoij", denmat::rmoij)
        .value("umoij", denmat::umoij)
        .value("rgen", denmat::rgen);
    // clang-format on

    // CAODensityMatrix class

    PyClass<CAODensityMatrix>(m, "AODensityMatrix")
        .def(py::init<>())
        .def(py::init<const CAODensityMatrix&>())
        .def(py::init(&CAODensityMatrix_from_numpy_list))
        .def("__str__", &CAODensityMatrix::getString)
        .def("alpha_to_numpy", &CAODensityMatrix_alpha_density_to_numpy)
        .def("beta_to_numpy", &CAODensityMatrix_beta_density_to_numpy)
        .def("number_of_density_matrices", &CAODensityMatrix::getNumberOfDensityMatrices)
        .def("get_density_type", &CAODensityMatrix::getDensityType)
        .def("sub", &CAODensityMatrix::sub)
        .def("append", &CAODensityMatrix::append)
        .def("broadcast", &CAODensityMatrix_broadcast)
        .def(py::self == py::self);

    // molorb enum class

    py::enum_<molorb>(m, "molorb").value("rest", molorb::rest).value("unrest", molorb::unrest);

    // CMolecularOrbitals class

    PyClass<CMolecularOrbitals>(m, "MolecularOrbitals")
        .def(py::init<>())
        .def(py::init<const CMolecularOrbitals&>())
        .def(py::init(&CMolecularOrbitals_from_numpy_list))
        .def("__str__", &CMolecularOrbitals::getString)
        .def("alpha_to_numpy", &CMolecularOrbitals_alpha_orbitals_to_numpy)
        .def("beta_to_numpy", &CMolecularOrbitals_beta_orbitals_to_numpy)
        .def("ea_to_numpy", &CMolecularOrbitals_alpha_energies_to_numpy)
        .def("eb_to_numpy", &CMolecularOrbitals_beta_energies_to_numpy)
        .def("get_orbitals_type", &CMolecularOrbitals::getOrbitalsType)
        .def("get_ao_density",
             vlx_general::overload_cast_<const int32_t>()(&CMolecularOrbitals::getAODensity, py::const_),
             "Compute spin-restricted AO density matrix for given number of electrons",
             "n_electrons"_a)
        .def("get_ao_density",
             vlx_general::overload_cast_<const int32_t, const int32_t>()(&CMolecularOrbitals::getAODensity, py::const_),
             "Compute spin-restricted AO density matrix for given number of alpha and beta electrons",
             "n_alpha_electrons"_a,
             "n_beta_electrons"_a)
        .def("get_pair_density",
             vlx_general::overload_cast_<const std::vector<int32_t>&, const std::vector<int32_t>&>()(
                 &CMolecularOrbitals::getRestrictedPairDensity, py::const_),
             "Computes set of restricted pair C_i C_j^T density matrices in AO basis",
             "i"_a,
             "j"_a)
        .def("get_pair_density",
             vlx_general::overload_cast_<const int32_t, const int32_t>()(
                 &CMolecularOrbitals::getRestrictedPairDensity, py::const_),
             "Computes restricted pair C_i C_j^T density matrix in AO basis",
             "i"_a,
             "j"_a)
        .def("insert", &CMolecularOrbitals::insert)
        .def("number_mos", &CMolecularOrbitals::getNumberOfColumns)
        .def("number_aos", &CMolecularOrbitals::getNumberOfRows)
        .def("alpha_orbitals",
             vlx_general::overload_cast_<const int32_t, const int32_t>()(&CMolecularOrbitals::alphaOrbitals, py::const_),
             "Get alpha orbitals within specific range",
             "i"_a,
             "n_orbitals"_a)
        .def("beta_orbitals",
             vlx_general::overload_cast_<const int32_t, const int32_t>()(&CMolecularOrbitals::betaOrbitals, py::const_),
             "Get beta orbitals within specific range",
             "i"_a,
             "n_orbitals"_a)
        .def("broadcast", &CMolecularOrbitals_broadcast)
        .def(py::self == py::self);

    // CSADGuessDriver class

    PyClass<CSADGuessDriver>(m, "SADGuessDriver")
        .def(py::init(&vlx_general::create<CSADGuessDriver>), "comm"_a = py::none())
        .def("compute", &CSADGuessDriver::compute);

    // exposing functions

    m.def("get_dimer_ao_indices", &aoindices::getDimerAOIndices);

    m.def("ao_matrix_to_veloxchem", &gtotra::to_veloxchem);

    m.def("ao_matrix_to_dalton", &gtotra::to_dalton);
}

}  // namespace vlx_orbdata
