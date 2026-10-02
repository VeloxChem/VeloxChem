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

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ExportGeneral.hpp"
#include "fika_classical_system.hpp"
#include "fika_dipole_potential_driver.hpp"
#include "fika_embedding.hpp"
#include "fika_force_field.hpp"
#include "fika_label.hpp"
#include "fika_mm_induced_dipoles.hpp"
#include "fika_nuclear_attraction_driver.hpp"
#include "fika_point_sources.hpp"
#include "fika_polarizable_sites.hpp"
#include "fika_quadrupole_potential_driver.hpp"
#include "fika_residue.hpp"
#include "fika_summation.hpp"
#include "fika_veloxchem.hpp"
#include "fika_veloxchem_order.hpp"

using namespace py::literals;

namespace vlx_fika {

namespace {

using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// Rows of a C-contiguous array whose last dimension is `width` ((n, width), or (m, k, width) as n = m k
/// rows), or of a flat array of n * width values; `what` names it in errors.
auto rows(const Array& array, std::size_t width, const std::string& what) -> std::size_t {
    const auto size = static_cast<std::size_t>(array.size());
    const bool shaped = array.ndim() >= 2 && static_cast<std::size_t>(array.shape(array.ndim() - 1)) == width;
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

/// (n, 3) numpy array of points.
auto to_numpy(const std::vector<fika::Point3D<double>>& values) -> py::array_t<double> {
    py::array_t<double> result({static_cast<py::ssize_t>(values.size()), py::ssize_t{3}});
    auto out = result.mutable_unchecked<2>();
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        const auto row = static_cast<py::ssize_t>(i);
        out(row, 0) = values[i].x;
        out(row, 1) = values[i].y;
        out(row, 2) = values[i].z;
    }
    return result;
}

/// Full n x n numpy matrix (any storage), in the order it is stored.
auto to_numpy(const fika::DenseMatrix& matrix) -> py::array_t<double> {
    const auto full = matrix.to_full();
    return vlx_general::pointer_to_numpy(
        full.data(), {static_cast<py::ssize_t>(matrix.rows()), static_cast<py::ssize_t>(matrix.columns())});
}

/// Molecule of a residue: atomic numbers and (n, 3) coordinates in bohr.
auto residue_molecule(const std::vector<int>& elements, const Array& coordinates) -> fika::Molecule<double> {
    const auto where = points(coordinates, elements.size());
    fika::Molecule<double> molecule;
    molecule.reserve(elements.size());
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        molecule.add_atom(fika::Element(elements[i]), where[i]);
    }
    return molecule;
}

/// Python-side drivers: options plus compute().
struct InducedDipolesDriver {
    fika::InducedDipoleOptions options;
};

struct QmmmEmbeddingDriver {
    fika::QmmmEmbeddingOptions options;
};

struct InducedDipoleFockDriver {
    double                threshold = 1.0e-12;
    fika::ChargeSummation summation = fika::ChargeSummation::automatic;
};

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

    // Classical (MM) system: force fields, residues and the system holding both.

    py::class_<fika::TholeDamping>(m, "FikaTholeDamping")
        .def(py::init([](const double a) { return fika::TholeDamping{a}; }), "a"_a = 2.1304)
        .def_readwrite("a", &fika::TholeDamping::a);

    py::class_<fika::ForceField>(m, "FikaForceField")
        .def(py::init([](const std::string&                          label,
                         const std::string&                          residue,
                         const std::vector<double>&                  charges,
                         const std::optional<std::map<int, double>>& polarizabilities) {
                 fika::ForceFieldParameters parameters;
                 parameters.label   = label;
                 parameters.residue = residue;
                 parameters.charges = charges;
                 if (polarizabilities)
                 {
                     std::vector<fika::AtomParameter<double>> isotropic;
                     for (const auto& [atom, alpha] : *polarizabilities)  // ordered by atom
                     {
                         if (atom < 0) throw std::invalid_argument("FikaForceField: negative atom index");
                         isotropic.push_back({static_cast<std::size_t>(atom), alpha});
                     }
                     parameters.polarizabilities = std::move(isotropic);
                 }
                 return fika::ForceField(std::move(parameters));
             }),
             "Force field of one residue type: a charge on every atom (in atom order) and isotropic "
             "polarizabilities {atom index (0-based): alpha} (a.u.) of a polarizable force field.",
             "label"_a,
             "residue"_a,
             "charges"_a,
             "polarizabilities"_a = py::none())
        .def_property_readonly("label", &fika::ForceField::label)
        .def_property_readonly("residue", &fika::ForceField::residue)
        .def_property_readonly("atom_count", &fika::ForceField::atom_count)
        .def_property_readonly("polarizable", &fika::ForceField::polarizable)
        .def_property_readonly("charges",
                               [](const fika::ForceField& self) {
                                   const auto q = self.charges();
                                   return std::vector<double>(q.begin(), q.end());
                               })
        .def_property_readonly("polarizabilities", [](const fika::ForceField& self) {
            std::map<std::size_t, double> result;
            for (const auto& [atom, alpha] : self.isotropic_polarizabilities()) result[atom] = alpha;
            return result;
        });

    py::class_<fika::Residue>(m, "FikaResidue")
        .def(py::init([](const std::string&      name,
                         const std::size_t       index,
                         const std::string&      force_field,
                         const std::vector<int>& elements,
                         const Array&            coordinates,
                         const bool              polarizable) {
                 return fika::Residue(name, index, force_field, residue_molecule(elements, coordinates), polarizable);
             }),
             "Residue `name` number `index` (0-based within its name) with force field label `force_field`: "
             "atomic numbers and coordinates (n, 3, bohr) in force-field atom order.",
             "name"_a,
             "index"_a,
             "force_field"_a,
             "elements"_a,
             "coordinates"_a,
             "polarizable"_a = false)
        .def_property_readonly("name", &fika::Residue::name)
        .def_property_readonly("index", &fika::Residue::index)
        .def_property_readonly("force_field", &fika::Residue::force_field)
        .def_property_readonly("polarizable", &fika::Residue::polarizable)
        .def_property_readonly("elements",
                               [](const fika::Residue& self) {
                                   std::vector<int> z;
                                   for (const auto& e : self.molecule().elements()) z.push_back(e.atomic_number());
                                   return z;
                               })
        .def_property_readonly("coordinates", [](const fika::Residue& self) {
            const auto c = self.molecule().coordinates();
            return to_numpy(std::vector<fika::Point3D<double>>(c.begin(), c.end()));
        });

    py::class_<fika::ClassicalSystem>(m, "FikaClassicalSystem")
        .def(py::init<>())
        .def("add_force_field", &fika::ClassicalSystem::add_force_field, "Adds (or replaces) a force field.", "force_field"_a)
        .def("add_residue", &fika::ClassicalSystem::add_residue, "Adds a residue to the region of its kind.", "residue"_a)
        .def(
            "add_residues",
            [](fika::ClassicalSystem&  self,
               const std::string&      name,
               const std::string&      force_field,
               const std::vector<int>& elements,
               const Array&            coordinates,
               const bool              polarizable) {
                const std::size_t atoms = elements.size();
                if (atoms == 0) throw std::invalid_argument("FikaClassicalSystem.add_residues: no atoms");
                const std::size_t total = rows(coordinates, 3, "coordinates");
                if (total % atoms != 0)
                {
                    throw std::invalid_argument("FikaClassicalSystem.add_residues: coordinates of " +
                                                std::to_string(total) + " atoms for residues of " +
                                                std::to_string(atoms));
                }
                // Indices continue those of the residues with this name already present.
                const auto label = fika::detail::normalized_label(name, "FikaClassicalSystem.add_residues: residue name");
                std::size_t index = 0;
                for (const auto* region : {&self.polarizable_region(), &self.nonpolarizable_region()})
                {
                    for (const auto& residue : region->residues()) index += residue.name() == label ? 1 : 0;
                }
                const double* data = coordinates.data();
                for (std::size_t r = 0; r < total / atoms; ++r)
                {
                    fika::Molecule<double> molecule;
                    molecule.reserve(atoms);
                    for (std::size_t a = 0; a < atoms; ++a)
                    {
                        const double* p = data + 3 * (r * atoms + a);
                        molecule.add_atom(fika::Element(elements[a]), {p[0], p[1], p[2]});
                    }
                    self.add_residue(fika::Residue(name, index++, force_field, std::move(molecule), polarizable));
                }
            },
            "Adds residues of one kind: atomic numbers of one residue and coordinates (residues, atoms, 3) or "
            "(residues * atoms, 3) in bohr.",
            "name"_a,
            "force_field"_a,
            "elements"_a,
            "coordinates"_a,
            "polarizable"_a = false)
        .def(
            "number_of_residues",
            [](const fika::ClassicalSystem& self, const bool polarizable) {
                return polarizable ? self.polarizable_region().size() : self.nonpolarizable_region().size();
            },
            "Residues in the polarizable or the nonpolarizable region.",
            "polarizable"_a)
        .def("force_field", [](const fika::ClassicalSystem& self, const std::string& label, const std::string& residue) {
            const auto* field = self.force_field(label, residue);
            return field == nullptr ? std::optional<fika::ForceField>() : std::optional<fika::ForceField>(*field);
        }, "label"_a, "residue"_a)
        .def("empty", &fika::ClassicalSystem::empty);

    m.def(
        "fika_polarizable_sites",
        [](const fika::ClassicalSystem& system) {
            const auto sites = fika::polarizable_sites(system);
            std::vector<double> alpha(sites.polarizabilities.size());
            for (std::size_t i = 0; i < alpha.size(); ++i) alpha[i] = sites.polarizabilities[i].components[0];
            return py::make_tuple(to_numpy(sites.positions), alpha, sites.owners);
        },
        "Polarizable sites: positions (n, 3, bohr), isotropic polarizabilities and owning residues (polarizable "
        "region index).",
        "system"_a);

    m.def(
        "fika_classical_charges",
        [](const fika::ClassicalSystem& system) {
            const auto sources = fika::classical_charges(system);
            return py::make_tuple(sources.charges, to_numpy(sources.coordinates), sources.residue_offsets);
        },
        "Permanent charges of both regions: charges, coordinates (n, 3, bohr) and residue offsets.",
        "system"_a);

    // Induced dipoles and the QM/MM embedding.

    py::class_<fika::InducedDipoleOptions>(m, "FikaInducedDipoleOptions")
        .def(py::init<>())
        .def_readwrite("tolerance", &fika::InducedDipoleOptions::tolerance)
        .def_readwrite("max_iterations", &fika::InducedDipoleOptions::max_iterations)
        .def_property(
            "initial_guess",
            [](const fika::InducedDipoleOptions& self) { return to_numpy(self.initial_guess); },
            [](fika::InducedDipoleOptions& self, const std::optional<Array>& guess) {
                self.initial_guess = guess ? points(*guess, rows(*guess, 3, "initial guess")) : std::vector<fika::Point3D<double>>{};
            })
        .def_readwrite("scale_initial_guess", &fika::InducedDipoleOptions::scale_initial_guess)
        .def_readwrite("summation", &fika::InducedDipoleOptions::summation)
        .def_readwrite("field_accuracy", &fika::InducedDipoleOptions::field_accuracy)
        .def_readwrite("permanent_field", &fika::InducedDipoleOptions::permanent_field);

    py::class_<fika::InducedDipoles>(m, "FikaInducedDipoles")
        .def_property_readonly("dipoles", [](const fika::InducedDipoles& self) { return to_numpy(self.dipoles); })
        .def_property_readonly("field", [](const fika::InducedDipoles& self) { return to_numpy(self.field); })
        .def_readonly("iterations", &fika::InducedDipoles::iterations)
        .def_readonly("residual", &fika::InducedDipoles::residual)
        .def_readonly("rms_residual", &fika::InducedDipoles::rms_residual)
        .def_readonly("converged", &fika::InducedDipoles::converged)
        .def_readonly("guess_used", &fika::InducedDipoles::guess_used)
        .def_readonly("guess_scale", &fika::InducedDipoles::guess_scale)
        .def_readonly("field_summation", &fika::InducedDipoles::field_summation)
        .def_readonly("coupling_summation", &fika::InducedDipoles::coupling_summation)
        .def_readonly("field_order", &fika::InducedDipoles::field_order)
        .def_readonly("coupling_order", &fika::InducedDipoles::coupling_order);

    py::class_<InducedDipolesDriver>(m, "FikaInducedDipolesDriver")
        .def(py::init([](const fika::InducedDipoleOptions& options) { return InducedDipolesDriver{options}; }),
             "options"_a = fika::InducedDipoleOptions{})
        .def_readwrite("options", &InducedDipolesDriver::options)
        .def(
            "compute",
            [](const InducedDipolesDriver& self, const fika::ClassicalSystem& system, const std::optional<fika::TholeDamping>& damping) {
                py::gil_scoped_release release;
                return fika::induced_dipoles(system, damping, self.options);
            },
            "Induced dipoles of the polarizable region in the field of the permanent charges (damping None: undamped).",
            "system"_a,
            "damping"_a = fika::TholeDamping{});

    py::enum_<fika::QmmmSources>(m, "FikaQmmmSources")
        .value("all", fika::QmmmSources::all)
        .value("electrons_only", fika::QmmmSources::electrons_only);

    py::class_<fika::QmmmEmbeddingOptions>(m, "FikaQmmmEmbeddingOptions")
        .def(py::init<>())
        .def_readwrite("induced", &fika::QmmmEmbeddingOptions::induced)
        .def_readwrite("sources", &fika::QmmmEmbeddingOptions::sources)
        .def_readwrite("fock_threshold", &fika::QmmmEmbeddingOptions::fock_threshold)
        .def_readwrite("fock_summation", &fika::QmmmEmbeddingOptions::fock_summation);

    py::class_<fika::QmmmEmbedding>(m, "FikaQmmmEmbedding")
        .def_readonly("induced", &fika::QmmmEmbedding::induced)
        .def_readonly("electron_field_summation", &fika::QmmmEmbedding::electron_field_summation)
        .def_readonly("electron_field_order", &fika::QmmmEmbedding::electron_field_order)
        .def_readonly("electron_permanent_energy", &fika::QmmmEmbedding::electron_permanent_energy)
        .def_readonly("nuclear_permanent_energy", &fika::QmmmEmbedding::nuclear_permanent_energy)
        .def_readonly("polarization_energy", &fika::QmmmEmbedding::polarization_energy)
        .def_property_readonly("permanent_fock", [](const fika::QmmmEmbedding& self) { return to_numpy(self.permanent_fock); })
        .def_property_readonly("induced_fock", [](const fika::QmmmEmbedding& self) { return to_numpy(self.induced_fock); })
        .def("energy", &fika::QmmmEmbedding::energy, "Sum of the energy components (Hartree).")
        .def("fock", [](const fika::QmmmEmbedding& self) { return to_numpy(self.fock()); }, "Sum of the Fock contributions.");

    py::class_<QmmmEmbeddingDriver>(m, "FikaQmmmEmbeddingDriver")
        .def(py::init([](const fika::QmmmEmbeddingOptions& options) { return QmmmEmbeddingDriver{options}; }),
             "options"_a = fika::QmmmEmbeddingOptions{})
        .def_readwrite("options", &QmmmEmbeddingDriver::options)
        .def(
            "compute",
            [](const QmmmEmbeddingDriver&                 self,
               const CMolecule&                           molecule,
               const CMolecularBasis&                     basis,
               const Array&                               density,
               const fika::ClassicalSystem&               system,
               const std::optional<fika::TholeDamping>&   damping) {
                const auto fmol   = fika::from_veloxchem(molecule);
                const auto fbasis = fika::from_veloxchem(basis, molecule);
                const auto n      = fbasis.function_count();
                if (density.ndim() != 2 || static_cast<std::size_t>(density.shape(0)) != n ||
                    static_cast<std::size_t>(density.shape(1)) != n)
                {
                    throw std::invalid_argument("FikaQmmmEmbeddingDriver: the density must be " + std::to_string(n) + " x " +
                                                std::to_string(n));
                }
                fika::DenseMatrix matrix(n, n, fika::MatrixSymmetry::general);
                std::copy(density.data(), density.data() + n * n, matrix.values().begin());
                py::gil_scoped_release release;
                return fika::qmmm_embedding(fmol, fbasis, matrix, system, damping, self.options);
            },
            "QM/MM embedding of the density (n_ao x n_ao, VeloxChem AO order, alpha + beta; any symmetry, its symmetric "
            "part acts): induced dipoles, energy components and Fock contributions in VeloxChem AO order.",
            "molecule"_a,
            "basis"_a,
            "density"_a,
            "system"_a,
            "damping"_a = fika::TholeDamping{});

    py::class_<InducedDipoleFockDriver>(m, "FikaInducedDipoleFockDriver")
        .def(py::init([](const double threshold, const fika::ChargeSummation summation) {
                 return InducedDipoleFockDriver{threshold, summation};
             }),
             "threshold"_a = 1.0e-12,
             "summation"_a = fika::ChargeSummation::automatic)
        .def_readwrite("threshold", &InducedDipoleFockDriver::threshold)
        .def_readwrite("summation", &InducedDipoleFockDriver::summation)
        .def(
            "compute",
            [](const InducedDipoleFockDriver& self, const CMolecule& molecule, const CMolecularBasis& basis,
               const Array& positions, const Array& dipoles) {
                const auto fmol   = fika::from_veloxchem(molecule);
                const auto fbasis = fika::from_veloxchem(basis, molecule);
                const auto count  = rows(dipoles, 3, "dipoles");
                const auto mu     = points(dipoles, count);
                const auto where  = points(positions, count);
                fika::DenseMatrix fock(0, fika::MatrixSymmetry::symmetric);
                {
                    py::gil_scoped_release release;
                    fock = fika::induced_dipole_fock(fmol, fbasis, where, mu, self.threshold, self.summation);
                }
                return to_numpy(fock);
            },
            "Fock contribution -sum_s mu_s . <a| (r - s) / |r - s|^3 |b> of dipoles (n, 3) at positions (n, 3, bohr), "
            "VeloxChem AO order.",
            "molecule"_a,
            "basis"_a,
            "positions"_a,
            "dipoles"_a);
}

}  // namespace vlx_fika
