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

#include "ExportFika.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "ExportGeneral.hpp"
#include "fika_dipole_potential_driver.hpp"
#include "fika_nuclear_attraction_driver.hpp"
#include "fika_quadrupole_potential_driver.hpp"
#include "fika_summation.hpp"
#include "fika_veloxchem.hpp"
#include "fika_veloxchem_order.hpp"

using namespace py::literals;

namespace vlx_fika {

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// Rows of an (n, width) array, or of a flat array of n * width values; `what` names it in errors.
auto rows(const Array& array, std::size_t width, const std::string& what) -> std::size_t {
    const auto size = static_cast<std::size_t>(array.size());
    const bool shaped = array.ndim() == 2 && static_cast<std::size_t>(array.shape(1)) == width;
    const bool flat = array.ndim() == 1 && size % width == 0;
    if (!shaped && !flat)
    {
        throw std::invalid_argument("Fika: " + what + " must have shape (n, " + std::to_string(width) + ")");
    }
    return size / width;
}

auto points(const Array& coordinates, std::size_t count) -> std::vector<fika::Point3D<double>> {
    if (rows(coordinates, 3, "coordinates") != count)
    {
        throw std::invalid_argument("Fika: " + std::to_string(count) + " sources but " +
                                    std::to_string(rows(coordinates, 3, "coordinates")) + " coordinates");
    }
    const double* data = coordinates.data();
    std::vector<fika::Point3D<double>> result(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        result[i] = {data[3 * i], data[3 * i + 1], data[3 * i + 2]};
    }
    return result;
}

template <int Rank>
auto tensors(const Array& values, const std::string& what) -> std::vector<fika::SymmetricTensor<Rank>> {
    constexpr std::size_t width = fika::SymmetricTensor<Rank>::size;
    const std::size_t count = rows(values, width, what);
    const double* data = values.data();
    std::vector<fika::SymmetricTensor<Rank>> result(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t k = 0; k < width; ++k)
        {
            result[i].components[k] = data[width * i + k];
        }
    }
    return result;
}

/// Full n x n numpy matrix in VeloxChem's AO order.
auto to_numpy(const fika::BlockSparseMatrix& matrix, const fika::MolecularBasis& basis) -> py::array_t<double> {
    const auto full = fika::fika_to_veloxchem(matrix.to_dense_matrix(), basis).to_full();
    const auto n = static_cast<py::ssize_t>(basis.function_count());
    return vlx_general::pointer_to_numpy(full.data(), {n, n});
}

}  // namespace

auto
export_fika(py::module& m) -> void
{
    py::enum_<fika::ChargeSummation>(m, "FikaChargeSummation")
        .value("automatic", fika::ChargeSummation::automatic)
        .value("direct", fika::ChargeSummation::direct)
        .value("multipole", fika::ChargeSummation::multipole);

    // Integral drivers: potential of the point sources, sum_C <a| phi_C(r) |b> (no electron charge
    // sign), as numpy matrices in VeloxChem's AO order.

    py::class_<fika::NuclearAttractionDriver>(m, "FikaNuclearAttractionDriver")
        .def(py::init<std::size_t, fika::ChargeSummation>(),
             py::arg("block_size") = 0,
             py::arg("summation")  = fika::ChargeSummation::automatic)
        .def(
            "compute",
            [](const fika::NuclearAttractionDriver& self,
               const CMolecule&                     molecule,
               const CMolecularBasis&               basis,
               const Array&                         charges,
               const Array&                         coordinates,
               const double                         threshold) {
                const auto fmol   = fika::from_veloxchem(molecule);
                const auto fbasis = fika::from_veloxchem(basis, molecule);
                const auto count  = rows(charges, 1, "charges");
                const std::vector<double> q(charges.data(), charges.data() + count);
                const auto where = points(coordinates, count);
                return to_numpy(self.compute(fmol, fbasis, q, where, threshold), fbasis);
            },
            "Computes sum_C q_C <a| 1/|r - C| |b> for point charges q_C at coordinates C (bohr).",
            "molecule"_a,
            "basis"_a,
            "charges"_a,
            "coordinates"_a,
            "threshold"_a = 1.0e-12);

    py::class_<fika::DipolePotentialDriver>(m, "FikaDipolePotentialDriver")
        .def(py::init<std::size_t, fika::ChargeSummation>(),
             py::arg("block_size") = 0,
             py::arg("summation")  = fika::ChargeSummation::automatic)
        .def(
            "compute",
            [](const fika::DipolePotentialDriver& self,
               const CMolecule&                   molecule,
               const CMolecularBasis&             basis,
               const Array&                       dipoles,
               const Array&                       coordinates,
               const double                       threshold) {
                const auto fmol   = fika::from_veloxchem(molecule);
                const auto fbasis = fika::from_veloxchem(basis, molecule);
                const auto mu     = tensors<1>(dipoles, "dipoles");
                const auto where  = points(coordinates, mu.size());
                return to_numpy(self.compute(fmol, fbasis, mu, where, threshold), fbasis);
            },
            "Computes sum_D <a| mu_D . (r - D) / |r - D|^3 |b> for point dipoles (n, 3) at coordinates D (bohr).",
            "molecule"_a,
            "basis"_a,
            "dipoles"_a,
            "coordinates"_a,
            "threshold"_a = 1.0e-12);

    py::class_<fika::QuadrupolePotentialDriver>(m, "FikaQuadrupolePotentialDriver")
        .def(py::init<std::size_t, fika::ChargeSummation>(),
             py::arg("block_size") = 0,
             py::arg("summation")  = fika::ChargeSummation::automatic)
        .def(
            "compute",
            [](const fika::QuadrupolePotentialDriver& self,
               const CMolecule&                       molecule,
               const CMolecularBasis&                 basis,
               const Array&                           quadrupoles,
               const Array&                           coordinates,
               const double                           threshold) {
                const auto fmol   = fika::from_veloxchem(molecule);
                const auto fbasis = fika::from_veloxchem(basis, molecule);
                const auto theta  = tensors<2>(quadrupoles, "quadrupoles");
                const auto where  = points(coordinates, theta.size());
                return to_numpy(self.compute(fmol, fbasis, theta, where, threshold), fbasis);
            },
            "Computes the potential of point quadrupoles (n, 6: primitive Cartesian moments xx, xy, xz, yy, yz, zz) at "
            "coordinates E (bohr): sum_E <a| 1/2 sum_ij Q_ij (3 x_i x_j - r^2 delta_ij) / r^5 |b>, x = r - E.",
            "molecule"_a,
            "basis"_a,
            "quadrupoles"_a,
            "coordinates"_a,
            "threshold"_a = 1.0e-12);
}

}  // namespace vlx_fika
