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

#include "fika_veloxchem.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BasisFunction.hpp"
#include "fika_basis_shell.hpp"
#include "fika_element.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::from_veloxchem: " + reason);
}

auto element_of(int identifier, const std::string& what) -> Element {
  if (identifier < 1 || identifier > detail::max_atomic_number) {
    fail(what + " has identifier " + std::to_string(identifier) +
         " (ghost atoms are not supported)");
  }
  return Element(identifier);
}

}  // namespace

auto from_veloxchem(const CMolecule& molecule) -> Molecule<double> {
  const auto identifiers = molecule.identifiers();
  const auto& coordinates = molecule.coordinates();  // bohr
  Molecule<double> result;
  result.reserve(identifiers.size());
  for (std::size_t i = 0; i < identifiers.size(); ++i) {
    const auto xyz = coordinates[i].coordinates();
    result.add_atom(element_of(identifiers[i], "atom " + std::to_string(i)),
                    {xyz[0], xyz[1], xyz[2]});
  }
  return result;
}

auto from_veloxchem(const CAtomBasis& basis) -> AtomBasis {
  if (basis.has_ecp()) {
    fail("atom basis " + basis.get_name() + " has an effective core potential");
  }
  std::vector<BasisShell> shells;
  for (const CBasisFunction& function : basis.functions()) {
    const int l = function.get_angular_momentum();
    const auto& exponents = function.exponents();
    const auto& norms = function.normalization_factors();
    std::size_t nonzero = 0, last = 0;
    for (std::size_t i = 0; i < norms.size(); ++i) {
      if (norms[i] != 0.0) {
        ++nonzero;
        last = i;
      }
    }
    if (nonzero == 1) {
      shells.emplace_back(UncontractedShell(l, exponents[last], norms[last], EffectiveCoefficients{}));
    } else {
      shells.emplace_back(SegmentedShell(l, exponents, norms, EffectiveCoefficients{}));
    }
  }
  return AtomBasis(element_of(basis.get_identifier(), "atom basis " + basis.get_name()),
                   basis.get_name(), std::move(shells));
}

auto from_veloxchem(const CMolecularBasis& basis, const CMolecule& molecule) -> MolecularBasis {
  const auto indices = basis.basis_sets_indices();
  const auto identifiers = molecule.identifiers();
  if (indices.size() != identifiers.size()) {
    fail("basis for " + std::to_string(indices.size()) + " atoms, molecule of " +
         std::to_string(identifiers.size()));
  }
  const auto& sets = basis.basis_sets();
  std::vector<AtomBasis> unique;
  unique.reserve(sets.size());
  for (const CAtomBasis& set : sets) {
    unique.push_back(from_veloxchem(set));
  }
  std::vector<std::uint16_t> atom_indices(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const int index = indices[i];
    if (index < 0 || static_cast<std::size_t>(index) >= sets.size()) {
      fail("atom " + std::to_string(i) + " has basis index " + std::to_string(index));
    }
    if (sets[static_cast<std::size_t>(index)].get_identifier() != identifiers[i]) {
      fail("atom " + std::to_string(i) + " (identifier " + std::to_string(identifiers[i]) +
           ") has the basis of identifier " +
           std::to_string(sets[static_cast<std::size_t>(index)].get_identifier()));
    }
    atom_indices[i] = static_cast<std::uint16_t>(index);
  }
  return MolecularBasis(std::move(unique), std::move(atom_indices));
}

}  // namespace fika
