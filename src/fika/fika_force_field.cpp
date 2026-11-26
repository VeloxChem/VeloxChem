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

#include "fika_force_field.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "fika_label.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::ForceField: " + reason);
}

auto finite(double value) -> bool {
  return std::isfinite(value);
}

template <int Rank>
auto finite(const SymmetricTensor<Rank>& tensor) -> bool {
  return std::ranges::all_of(tensor.components, [](double v) { return std::isfinite(v); });
}

/// Checks the atom indices (in range, strictly increasing) and values of one quantity.
template <typename T>
void check(const std::vector<AtomParameter<T>>& entries, std::size_t atoms, const char* what) {
  for (std::size_t i = 0; i < entries.size(); ++i) {
    const std::size_t atom = entries[i].atom;
    if (atom >= atoms) {
      fail(std::string(what) + " on atom " + std::to_string(atom) + " of " + std::to_string(atoms));
    }
    if (i > 0 && atom <= entries[i - 1].atom) {
      fail(std::string(what) + " atom indices must be strictly increasing");
    }
    if (!finite(entries[i].value)) {
      fail(std::string(what) + " on atom " + std::to_string(atom) + " is not finite");
    }
  }
}

}  // namespace

ForceField::ForceField(ForceFieldParameters parameters) : parameters_(std::move(parameters)) {
  parameters_.label =
      detail::normalized_label(std::move(parameters_.label), "fika::ForceField: force-field label");
  parameters_.residue =
      detail::normalized_label(std::move(parameters_.residue), "fika::ForceField: residue label");
  const std::size_t atoms = parameters_.charges.size();
  if (atoms == 0) {
    fail("no charges (every atom needs one)");
  }
  if (!std::ranges::all_of(parameters_.charges, [](double q) { return std::isfinite(q); })) {
    fail("charges must be finite");
  }
  check(parameters_.dipoles, atoms, "dipole");
  check(parameters_.quadrupoles, atoms, "quadrupole");
  bool oriented = !parameters_.dipoles.empty() || !parameters_.quadrupoles.empty();
  if (const auto* isotropic =
          std::get_if<std::vector<AtomParameter<double>>>(&parameters_.polarizabilities)) {
    if (isotropic->empty()) {
      fail("empty isotropic polarizability block");
    }
    check(*isotropic, atoms, "isotropic polarizability");
    if (std::ranges::any_of(*isotropic, [](const auto& entry) { return entry.value < 0.0; })) {
      fail("isotropic polarizabilities must be non-negative");
    }
  }
  if (const auto* anisotropic =
          std::get_if<std::vector<AtomParameter<Polarizability>>>(&parameters_.polarizabilities)) {
    if (anisotropic->empty()) {
      fail("empty anisotropic polarizability block");
    }
    check(*anisotropic, atoms, "anisotropic polarizability");
    oriented = true;
  }
  const auto& geometry = parameters_.reference_geometry;
  if (oriented && !geometry) {
    fail("dipoles, quadrupoles and anisotropic polarizabilities need a reference geometry");
  }
  if (geometry && geometry->size() != atoms) {
    fail("reference geometry has " + std::to_string(geometry->size()) + " atoms, charges " +
         std::to_string(atoms));
  }
  if (geometry && !std::ranges::all_of(geometry->coordinates(), [](const auto& p) {
        return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
      })) {
    fail("reference geometry coordinates must be finite");
  }
}

auto ForceField::isotropic_polarizabilities() const noexcept
    -> std::span<const AtomParameter<double>> {
  const auto* isotropic =
      std::get_if<std::vector<AtomParameter<double>>>(&parameters_.polarizabilities);
  return isotropic != nullptr ? std::span<const AtomParameter<double>>(*isotropic)
                              : std::span<const AtomParameter<double>>();
}

auto ForceField::anisotropic_polarizabilities() const noexcept
    -> std::span<const AtomParameter<Polarizability>> {
  const auto* anisotropic =
      std::get_if<std::vector<AtomParameter<Polarizability>>>(&parameters_.polarizabilities);
  return anisotropic != nullptr ? std::span<const AtomParameter<Polarizability>>(*anisotropic)
                                : std::span<const AtomParameter<Polarizability>>();
}

}  // namespace fika
