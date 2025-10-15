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

#include "ExportGeneral.hpp"

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AtomicRadii.hpp"
#include "BatchFunc.hpp"
#include "Codata.hpp"
#include "CpcmUtils.hpp"
#include "ErrorHandler.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "SphericalMomentum.hpp"
#include "StringFormat.hpp"
#include "TensorComponents.hpp"
#include "TensorLabels.hpp"
#include "T3FlatBuffer.hpp"
#include "T3RectFlatBuffer.hpp"
#include "MathFunc.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_general {

// Gets shape and strides from dimension

static auto
get_shape_and_strides(const std::vector<int>& dimension) -> std::tuple<std::vector<py::ssize_t>, std::vector<py::ssize_t>>
{
    std::vector<py::ssize_t> shape, strides;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        shape.push_back(static_cast<py::ssize_t>(dimension[i]));

        size_t strd = 1;

        for (size_t j = i + 1; j < dimension.size(); j++)
        {
            strd *= dimension[j];
        }

        strides.push_back(static_cast<py::ssize_t>(strd * sizeof(double)));
    }

    return {shape, strides};
}

// Creates numpy ndarray from pointer and dimension
// Not a static function; used in other files

auto
pointer_to_numpy(const double* ptr, const std::vector<int>& dimension) -> py::array_t<double>
{
    if (dimension.size() == 0)
    {
        return py::array_t<double>();
    }
    else
    {
        if (ptr == nullptr)
        {
            // double check that total number of elements is zero
            int num_elems = 1;
            for (const auto& n : dimension)
            {
                num_elems *= n;
            }
            errors::assertMsgCritical(num_elems == 0, std::string("pointer_to_numpy: Invalid dimension for nullptr"));
        }

        const auto [shape, strides] = get_shape_and_strides(dimension);

        return py::array_t<double>(shape, strides, ptr);
    }
}

auto
export_general(py::module &m) -> void
{
    /// mpi master
    m.def(
        "mpi_master",
        []() -> int { return 0; },
        "Gets the rank of MPI master process.");

    /// mpi size limit
    /// here we choose a number that is smaller than the actual limit
    m.def(
        "mpi_size_limit",
        []() -> int { return static_cast<int>(1 << 30) / 5 * 9; },
        "Gets the size limit in MPI communication (below 2^31-1).");

    /// exposing functions from StringFormat.hpp
    m.def("upper_case", &format::upper_case, "Creates upper cased copy of given string.");
    m.def("lower_case", &format::lower_case, "Creates lower cased copy of given string.");

    // exposing functions from Codata.hpp
    m.def("bohr_in_angstrom", &units::bohr_in_angstrom, "Gets Bohr value in Angstroms.");
    m.def("bohr_in_angstroms", &units::bohr_in_angstrom, "Gets Bohr value in Angstroms.");
    m.def("hartree_in_ev", &units::hartree_in_ev, "Gets Hartree value in electronvolts.");
    m.def("hartree_in_kcalpermol", &units::getHartreeValueInKiloCaloriePerMole, "Gets Hartree value in kcal/mol.");
    m.def("hartree_in_kjpermol", &units::getHartreeValueInKiloJoulePerMole, "Gets Hartree value in kJ/mol.");
    m.def("hartree_in_inverse_nm", &units::getHartreeValueInInverseNanometer, "Gets Hartree value in inverse nanometer.");
    m.def("hartree_in_wavenumber", &units::getHartreeValueInWavenumbers, "Gets Hartree value in reciprocal cm.");
    m.def("hartree_in_wavenumbers", &units::getHartreeValueInWavenumbers, "Gets Hartree value in reciprocal cm.");
    m.def("electron_mass_in_amu", &units::getElectronMassInAtomicMassUnit, "Gets electron mass in amu.");
    m.def("amu_in_electron_mass", &units::getAtomicMassUnitInElectronMasses, "Gets atomic mass unit in electron masses.");
    m.def("amu_in_electron_masses", &units::getAtomicMassUnitInElectronMasses, "Gets atomic mass unit in electron masses.");
    m.def("amu_in_kg", &units::getAtomicMassUnitInKg, "Gets atomic mass unit in kg.");
    m.def("speed_of_light_in_vacuum_in_SI", &units::getSpeedOfLightInVacuumInSI, "Gets speed of light in vacuum in SI.");
    m.def("avogadro_constant", &units::getAvogadroConstant, "Gets Avogadro constant.");
    m.def("boltzmann_in_evperkelvin", &units::getBoltzmannConstantInElectronVoltsPerKelvin, "Gets Boltzmann constant in eV/K.");
    m.def("boltzmann_in_hartreeperkelvin", &units::getBoltzmannConstantInHartreePerKelvin, "Gets Boltzmann constant in Hartree/K.");

    m.def("dipole_in_debye", &units::getDipoleInDebye, "Gets convertion factor for dipole moment (a.u. -> Debye).");
    m.def("rotatory_strength_in_cgs", &units::getRotatoryStrengthInCGS, "Gets convertion factor for rotatory strength (a.u. -> 10^-40 cgs).");
    m.def("extinction_coefficient_from_beta",
          &units::getExtinctionCoefficientFromBeta,
          "Gets factor needed for the calculation of the extinction coefficent from the electric-dipole magnetic-dipole polarizability beta.");
    m.def("fine_structure_constant", &units::getFineStructureConstant, "Gets fine-structure constant.");

    // exposing functions from AtomicRadii.hpp
    m.def("get_vdw_radii_data_in_bohr", &atomicradii::buildVdwRadii, "Gets VDW radii data in atomic unit.");

    // exposing functions from TensorLabels.hpp
    m.def("tensor_cartesian_labels", &tensor::cartesian_labels, "Gets all Cartesian component labels of tensor.");
    m.def("tensor_spherical_labels", &tensor::spherical_labels, "Gets all spherical component labels of tensor.");
    m.def("tensor_cartesian_index", &tensor::cartesian_index, "Gets index of Cartesian tensor component.");
    m.def("tensor_label", &tensor::label, "Gets label of tensor.");
    m.def("tensor_order", &tensor::order, "Gets order of tensor.");

    // exposing functions from TensorComponents.hpp
    m.def(
        "number_of_cartesian_components",
        [](const int order) -> int { return tensor::number_of_cartesian_components(std::array<int, 1>{order}); },
        "Gets number of Cartesian components in tensor.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 2> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 3> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_cartesian_components",
        [](const std::array<int, 4> orders) -> int { return tensor::number_of_cartesian_components(orders); },
        "Gets number of Cartesian components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const int order) -> int { return tensor::number_of_spherical_components(std::array<int, 1>{order}); },
        "Gets number of spherical components in tensor.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 2> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 3> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");
    m.def(
        "number_of_spherical_components",
        [](const std::array<int, 4> orders) -> int { return tensor::number_of_spherical_components(orders); },
        "Gets number of spherical components in array of tensors.");

    // exposing functions from BatchFunc.hpp
    m.def("number_of_batches", &batch::number_of_batches<size_t>, "Gets number of batches.");
    m.def("batch_range",
          py::overload_cast<const size_t, const size_t, const size_t>(&batch::batch_range<size_t>),
          "Gets [first, last) range for requested batch.");
    m.def("batch_range",
          py::overload_cast<const size_t, const size_t, const size_t, const size_t>(&batch::batch_range<size_t>),
          "Gets [first, last) range for requested batch.");

    // exposing functions from OpenMPFunc.hpp
    m.def("set_number_of_threads", &omp::set_number_of_threads, "Sets number of OMP threads to requested value.");
    m.def("get_number_of_threads", &omp::get_number_of_threads, "Gets number of OMP threads available.");
    m.def("make_work_tasks",
          py::overload_cast<const std::vector<CGtoBlock> &>(&omp::make_work_tasks),
          "Gets work tasks for given vector of basis function blocks.");
    m.def("make_work_tasks",
          py::overload_cast<const std::vector<CGtoBlock> &, const std::vector<CGtoBlock> &>(&omp::make_work_tasks),
          "Gets work tasks for given two vectors of basis function blocks.");
    m.def("make_diag_work_tasks", &omp::make_diag_work_group, "Gets work tasks for diagonal integrals.");
    m.def("make_work_group",
          py::overload_cast<const std::vector<CBlockedGtoPairBlock>&, const int, const int>(&omp::make_work_group),
          "Gets work group for ERIs.");
    m.def("make_work_group",
          py::overload_cast<const std::vector<CBlockedGtoPairBlock>&, const int, const int, const int>(&omp::make_work_group),
          "Gets work group for ERIs.");
    m.def("make_bra_ket_work_group", &omp::make_bra_ket_work_group, "Gets work group for ERIs.");
    m.def("partition_atoms", &omp::partition_atoms, "Get atomic indices of partitioned atoms list.");
    m.def("partition_flat_buffer", &omp::partition_flat_buffer, "Creates vector of indices required for partitioning of flat buffer.");
    
    // exposing functions from SphericalMomentum.hpp
    m.def("spherical_momentum_s_factors", spher_mom::transformation_factors<0>, "Gets transformation factors for S type spherical momentum.");
    m.def("spherical_momentum_p_factors", spher_mom::transformation_factors<1>, "Gets transformation factors for P type spherical momentum.");
    m.def("spherical_momentum_d_factors", spher_mom::transformation_factors<2>, "Gets transformation factors for D type spherical momentum.");
    m.def("spherical_momentum_f_factors", spher_mom::transformation_factors<3>, "Gets transformation factors for F type spherical momentum.");
    m.def("spherical_momentum_g_factors", spher_mom::transformation_factors<4>, "Gets transformation factors for G type spherical momentum.");

    // exposing functions from CpcmUtils.hpp
    m.def("cpcm_local_matrix_A_diagonals",
          [](const py::array_t<double>& grid_data,
             const py::array_t<double>& sw_func,
             const int                  grid_index_start,
             const int                  grid_index_end) -> py::array_t<double> {
              std::string errstyle("cpcm_form_matrix_A: Expecting contiguous numpy arrays");
              auto        c_style_1 = py::detail::check_flags(grid_data.ptr(), py::array::c_style);
              auto        c_style_2 = py::detail::check_flags(sw_func.ptr(), py::array::c_style);
              errors::assertMsgCritical((c_style_1 && c_style_2), errstyle);
              std::string errsize("cpcm_form_matrix_A: Inconsistent sizes");
              errors::assertMsgCritical(grid_data.shape(0) == sw_func.shape(0), errsize);
              const auto npoints = static_cast<int>(grid_data.shape(0));
              const auto ncols = static_cast<int>(grid_data.shape(1));
              std::string errindex("cpcm_form_matrix_A: Invalid indices");
              errors::assertMsgCritical((0 <= grid_index_start) && (grid_index_start < npoints), errindex);
              errors::assertMsgCritical((0 < grid_index_end) && (grid_index_end <= npoints), errindex);
              auto Adiag = cpcm::local_matrix_A_diagonals(grid_data.data(), grid_index_start, grid_index_end, ncols, sw_func.data());
              return vlx_general::pointer_to_numpy(Adiag.data(), {static_cast<int>(Adiag.size())});
          },
          "Form local diagonals of C-PCM matrix A.",
          "grid_data"_a,
          "sw_func"_a,
          "grid_index_start"_a,
          "grid_index_end"_a);

    m.def("cpcm_local_matrix_A_dot_vector",
          [](const py::array_t<double>& grid_data,
             const py::array_t<double>& sw_func,
             const py::array_t<double>& vector,
             const int                  grid_index_start,
             const int                  grid_index_end) -> py::array_t<double> {
              std::string errstyle("cpcm_local_matrix_A_dot_vector: Expecting contiguous numpy arrays");
              auto        c_style_1 = py::detail::check_flags(grid_data.ptr(), py::array::c_style);
              auto        c_style_2 = py::detail::check_flags(sw_func.ptr(), py::array::c_style);
              auto        c_style_3 = py::detail::check_flags(vector.ptr(), py::array::c_style);
              errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3), errstyle);
              std::string errsize("cpcm_local_matrix_A_dot_vector: Inconsistent sizes");
              errors::assertMsgCritical(grid_data.shape(0) == sw_func.shape(0), errsize);
              errors::assertMsgCritical(grid_data.shape(0) == vector.shape(0), errsize);
              const auto npoints = static_cast<int>(grid_data.shape(0));
              const auto ncols = static_cast<int>(grid_data.shape(1));
              std::string errindex("cpcm_local_matrix_A_dot_vector: Invalid indices");
              errors::assertMsgCritical((0 <= grid_index_start) && (grid_index_start < npoints), errindex);
              errors::assertMsgCritical((0 < grid_index_end) && (grid_index_end <= npoints), errindex);
              auto prod = cpcm::local_matrix_A_dot_vector(npoints, grid_data.data(), grid_index_start, grid_index_end, ncols, sw_func.data(), vector.data());
              return vlx_general::pointer_to_numpy(prod.data(), {static_cast<int>(prod.size())});
          },
          "Form local product of C-PCM matrix A and a vector.",
          "grid_data"_a,
          "sw_func"_a,
          "vector"_a,
          "grid_index_start"_a,
          "grid_index_end"_a);

    m.def("cpcm_comp_grad_Aij",
          [](const py::array_t<double>& grid_coords,
             const py::array_t<double>& zeta,
             const py::array_t<int>&    atom_indices,
             const py::array_t<double>& q,
             const int                  grid_index_start,
             const int                  grid_index_end,
             const int                  natoms) -> py::array_t<double> {
              std::string errstyle("cpcm_comp_grad_Aij: Expecting contiguous numpy arrays");
              auto        c_style_1 = py::detail::check_flags(grid_coords.ptr(), py::array::c_style);
              auto        c_style_2 = py::detail::check_flags(zeta.ptr(), py::array::c_style);
              auto        c_style_3 = py::detail::check_flags(atom_indices.ptr(), py::array::c_style);
              auto        c_style_4 = py::detail::check_flags(q.ptr(), py::array::c_style);
              errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
              std::string errsize("cpcm_comp_grad_Aij: Inconsistent sizes");
              errors::assertMsgCritical(grid_coords.shape(0) == zeta.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(0) == atom_indices.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(0) == q.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(1) == 3, errsize);
              const auto npoints = static_cast<int>(grid_coords.shape(0));
              std::string errindex("cpcm_comp_grad_Aij: Invalid indices");
              errors::assertMsgCritical((0 <= grid_index_start) && (grid_index_start < npoints), errindex);
              errors::assertMsgCritical((0 < grid_index_end) && (grid_index_end <= npoints), errindex);
              auto grad_Amat = cpcm::comp_grad_Aij(grid_coords.data(), zeta.data(), atom_indices.data(), q.data(),
                                                   grid_index_start, grid_index_end, npoints, natoms);
              return vlx_general::pointer_to_numpy(grad_Amat.data(), {natoms, 3});
          },
          "Compute C-PCM gradient for Aij.",
          "grid_coords"_a,
          "zeta"_a,
          "atom_indices"_a,
          "q"_a,
          "grid_index_start"_a,
          "grid_index_end"_a,
          "natoms"_a);

    m.def("cpcm_comp_grad_Aii",
          [](const py::array_t<double>& grid_coords,
             const py::array_t<double>& zeta,
             const py::array_t<double>& sw_f,
             const py::array_t<int>&    atom_indices,
             const py::array_t<double>& q,
             const int                  grid_index_start,
             const int                  grid_index_end,
             const py::array_t<double>& atom_coords,
             const py::array_t<double>& atom_radii) -> py::array_t<double> {
              std::string errstyle("cpcm_comp_grad_Aii: Expecting contiguous numpy arrays");
              auto        c_style_1 = py::detail::check_flags(grid_coords.ptr(), py::array::c_style);
              auto        c_style_2 = py::detail::check_flags(zeta.ptr(), py::array::c_style);
              auto        c_style_3 = py::detail::check_flags(sw_f.ptr(), py::array::c_style);
              auto        c_style_4 = py::detail::check_flags(atom_indices.ptr(), py::array::c_style);
              auto        c_style_5 = py::detail::check_flags(q.ptr(), py::array::c_style);
              auto        c_style_6 = py::detail::check_flags(atom_coords.ptr(), py::array::c_style);
              auto        c_style_7 = py::detail::check_flags(atom_radii.ptr(), py::array::c_style);
              errors::assertMsgCritical((c_style_1 && c_style_2 && c_style_3 && c_style_4), errstyle);
              errors::assertMsgCritical((c_style_5 && c_style_6 && c_style_7), errstyle);
              std::string errsize("cpcm_comp_grad_Aii: Inconsistent sizes");
              errors::assertMsgCritical(grid_coords.shape(0) == zeta.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(0) == sw_f.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(0) == atom_indices.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(0) == q.shape(0), errsize);
              errors::assertMsgCritical(grid_coords.shape(1) == 3, errsize);
              errors::assertMsgCritical(atom_coords.shape(0) == atom_radii.shape(0), errsize);
              errors::assertMsgCritical(atom_coords.shape(1) == 3, errsize);
              const auto npoints = static_cast<int>(grid_coords.shape(0));
              const auto natoms = static_cast<int>(atom_coords.shape(0));
              std::string errindex("cpcm_comp_grad_Aii: Invalid indices");
              errors::assertMsgCritical((0 <= grid_index_start) && (grid_index_start < npoints), errindex);
              errors::assertMsgCritical((0 < grid_index_end) && (grid_index_end <= npoints), errindex);
              auto grad_Amat = cpcm::comp_grad_Aii(grid_coords.data(), zeta.data(), sw_f.data(), atom_indices.data(), q.data(),
                                                   atom_coords.data(), atom_radii.data(), grid_index_start, grid_index_end, npoints, natoms);
              return vlx_general::pointer_to_numpy(grad_Amat.data(), {natoms, 3});
          },
          "Compute C-PCM gradient for Aii.",
          "grid_coords"_a,
          "zeta"_a,
          "sw_f"_a,
          "atom_indices"_a,
          "q"_a,
          "grid_index_start"_a,
          "grid_index_end"_a,
          "atom_coords"_a,
          "atom_radii"_a);

    // TPoint class
    PyClass<TPoint<double>>(m, "Point")
        .def(py::init<>())
        .def(py::init<const std::array<double, 3> &>())
        .def(py::init<const TPoint<double> &>())
        .def(py::pickle([](const TPoint<double> &pnt) { return py::make_tuple(pnt.coordinates()); },
                        [](py::tuple t) { return TPoint<double>(t[0].cast<std::array<double, 3>>()); }))
        .def("coordinates", &TPoint<double>::coordinates, "Getter for Cartesian coordinates.")
        .def("scale", &TPoint<double>::scale, "Scales Cartesian coordinates by factor.")
        .def("length_square", &TPoint<double>::length_square, "Computes square of length for vector given by point.")
        .def("length", &TPoint<double>::length, "Computes length for vector given by point.")
        .def("distance_square", &TPoint<double>::distance_square, "Computes square of distance between two points.")
        .def("distance", &TPoint<double>::distance, "Computes distance between two points.")
        .def("__eq__", [](const TPoint<double> &self, const TPoint<double> &other) { return self == other; })
        .def("__ne__", [](const TPoint<double> &self, const TPoint<double> &other) { return self != other; })
        .def("__copy__", [](const TPoint<double> &self) { return TPoint<double>(self); })
        .def("__deepcopy__", [](const TPoint<double> &self, py::dict) { return TPoint<double>(self); });
    
    // CT3FlatBuffer class
    PyClass<CT3FlatBuffer<double>>(m, "T3FlatBuffer")
        .def(py::init<>())
        .def(py::init<const std::vector<size_t>&, const size_t>())
        .def(py::init<const CT3FlatBuffer<double> &>())
        .def("indices", &CT3FlatBuffer<double>::indices, "Gets indices vector along x axis of tensor.")
        .def("mask_indices", &CT3FlatBuffer<double>::mask_indices, "Gets mask of indices along x axis of tensor.")
        .def("width", &CT3FlatBuffer<double>::width, "Gets width along x axis of tensor.")
        .def("aux_width", &CT3FlatBuffer<double>::aux_width, "Gets width along y,x axes of tensor.")
        .def(
            "values",
             [](const CT3FlatBuffer<double> &self, const size_t index) -> py::array_t<double> {
                 const auto ndim = self.width();
                 const auto nelems = static_cast<py::ssize_t>(ndim * (ndim + 1) / 2);
                 const auto tdim  = static_cast<py::ssize_t>(sizeof(double));
                 return py::array_t<double>(
                                            std::vector<py::ssize_t>({nelems,}), std::vector<py::ssize_t>({tdim, }), self.data(index));
            },
            "Gets slice of tensor values along y,z axes.")
        .def("value",
             [](const CT3FlatBuffer<double> &self, const size_t index, const size_t i, const size_t j) -> double {
                if (i <= j)
                {
                    return self.data(index)[mathfunc::uplo_rm_index(i, j, self.width())];
                }
                else
                {
                    return self.data(index)[mathfunc::uplo_rm_index(j, i, self.width())];
                }
             },
             "Gets tensor element value.")
        .def("__copy__", [](CT3FlatBuffer<double> &self) { return CT3FlatBuffer<double>(self); })
        .def("__deepcopy__", [](const CT3FlatBuffer<double> &self, py::dict) { return CT3FlatBuffer<double>(self); });
    
    // CT3RectFlatBuffer class
    PyClass<CT3RectFlatBuffer<double>>(m, "T3RectFlatBuffer")
        .def(py::init<>())
        .def(py::init<const std::vector<size_t>&, const std::map<size_t, size_t>&, const size_t>())
        .def(py::init<const CT3RectFlatBuffer<double> &>())
        .def("indices", &CT3RectFlatBuffer<double>::indices, "Gets indices vector along x axis of tensor.")
        .def("mask_indices", &CT3RectFlatBuffer<double>::mask_indices, "Gets masked indices along y axis of tensor.")
        .def("width", &CT3RectFlatBuffer<double>::width, "Gets width along z axis of tensor.")
        .def(
            "values",
             [](const CT3RectFlatBuffer<double> &self, const size_t index) -> py::array_t<double> {
                 const auto nrows = (self.mask_indices().size() == 0) ? self.width() : self.mask_indices().size();
                 const auto ncols = self.width();
                 const auto nelems = static_cast<py::ssize_t>(nrows * ncols);
                 const auto tdim  = static_cast<py::ssize_t>(sizeof(double));
                 return py::array_t<double>(
                                            std::vector<py::ssize_t>({nelems,}), std::vector<py::ssize_t>({tdim, }), self.data(index));
            },
            "Gets slice of tensor values along y,z axes.")
        .def("value",
             [](const CT3RectFlatBuffer<double> &self, const size_t index, const size_t i, const size_t j) -> double {
                const auto mask = self.mask_indices();
                return self.data(index)[mask.at(i) * self.width() + j];
             },
             "Gets tensor element value.")
        .def("__copy__", [](CT3RectFlatBuffer<double> &self) { return CT3RectFlatBuffer<double>(self); })
        .def("__deepcopy__", [](const CT3RectFlatBuffer<double> &self, py::dict) { return CT3RectFlatBuffer<double>(self); });
}

}  // namespace vlx_general
