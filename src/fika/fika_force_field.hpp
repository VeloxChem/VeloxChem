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

#ifndef fika_force_field_hpp
#define fika_force_field_hpp

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <variant>
#include <vector>

#include "fika_symmetric_tensor.hpp"
#include "fika_molecule.hpp"

namespace fika {

/// A value assigned to one atom, by its index in the residue's atom order.
template <typename T>
struct AtomParameter {
  std::size_t atom = 0;
  T value{};
};

/// Anisotropic polarizability: symmetric 3 x 3 tensor (xx, xy, xz, yy, yz, zz).
using Polarizability = SymmetricTensor<2>;

/// Input of a ForceField (atomic units; reference geometry in bohr). Every atom carries a charge;
/// dipoles, quadrupoles and polarizabilities may skip atoms (listed by increasing atom index).
struct ForceFieldParameters {
  std::string label;            // force-field label, e.g. "loprop-mm3"
  std::string residue;          // standard label of the residue it describes, e.g. "hoh"
  std::vector<double> charges;  // one per atom: defines the atom count
  std::vector<AtomParameter<Dipole>> dipoles;
  std::vector<AtomParameter<Quadrupole>> quadrupoles;  // primitive Cartesian moments
  // Not polarizable, isotropic or anisotropic polarizabilities (one kind per force field).
  std::variant<std::monostate, std::vector<AtomParameter<double>>,
               std::vector<AtomParameter<Polarizability>>>
      polarizabilities;
  // Geometry the dipoles, quadrupoles and anisotropic polarizabilities are given in; required
  // when any of them is present.
  std::optional<Molecule<double>> reference_geometry;
};

/// Parameters of one residue type in a force field with a strict atom order: permanent charges
/// on every atom, permanent dipoles and quadrupoles and (if polarizable) isotropic or anisotropic
/// polarizabilities on some atoms, and the reference geometry of the orientation-dependent
/// quantities. Both labels are stored lowercased.
class ForceField {
 public:
  /// Throws std::invalid_argument if a label is empty or contains whitespace, there are no
  /// charges, an atom index is out of range or indices of a quantity are not strictly
  /// increasing, a value is not finite, an isotropic polarizability is negative, a
  /// polarizability block is empty, or the reference geometry is missing when needed or has a
  /// different number of atoms.
  explicit ForceField(ForceFieldParameters parameters);

  auto label() const noexcept -> const std::string& { return parameters_.label; }
  auto residue() const noexcept -> const std::string& { return parameters_.residue; }
  auto atom_count() const noexcept -> std::size_t { return parameters_.charges.size(); }

  auto charges() const noexcept -> std::span<const double> { return parameters_.charges; }
  auto dipoles() const noexcept -> std::span<const AtomParameter<Dipole>> {
    return parameters_.dipoles;
  }
  auto quadrupoles() const noexcept -> std::span<const AtomParameter<Quadrupole>> {
    return parameters_.quadrupoles;
  }

  auto polarizable() const noexcept -> bool {
    return !std::holds_alternative<std::monostate>(parameters_.polarizabilities);
  }
  /// Isotropic polarizabilities (empty unless the force field has them).
  auto isotropic_polarizabilities() const noexcept -> std::span<const AtomParameter<double>>;
  /// Anisotropic polarizabilities (empty unless the force field has them).
  auto anisotropic_polarizabilities() const noexcept
      -> std::span<const AtomParameter<Polarizability>>;

  auto reference_geometry() const noexcept -> const std::optional<Molecule<double>>& {
    return parameters_.reference_geometry;
  }

 private:
  ForceFieldParameters parameters_;
};

}  // namespace fika

#endif  // fika_force_field_hpp
