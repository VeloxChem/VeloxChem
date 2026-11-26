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

#ifndef fika_classical_system_hpp
#define fika_classical_system_hpp

#include <map>
#include <span>
#include <string>
#include <utility>

#include "fika_classical_region.hpp"
#include "fika_force_field.hpp"
#include "fika_residue.hpp"

namespace fika {

/// The classical (MM) environment: a polarizable and a non-polarizable region, either or both of
/// which may be empty. The regions are read-only from outside, so each keeps its kind; residues
/// are added through add_residue and edited through the residue spans.
///
/// The system also holds the force fields of its residues, keyed by force-field label and residue
/// name (both lowercased): every residue's force field must be added through add_force_field.
class ClassicalSystem {
 public:
  /// Both regions empty.
  ClassicalSystem() = default;

  /// Throws std::invalid_argument if `polarizable` is not a polarizable region or
  /// `nonpolarizable` is.
  ClassicalSystem(ClassicalRegion polarizable, ClassicalRegion nonpolarizable);

  auto polarizable_region() const noexcept -> const ClassicalRegion& { return polarizable_; }
  auto nonpolarizable_region() const noexcept -> const ClassicalRegion& { return nonpolarizable_; }

  /// Appends a residue to the region of its kind (its polarizable flag).
  auto add_residue(Residue residue) -> void;

  /// Editable residues of each region; residues cannot be added or removed through these views.
  auto polarizable_residues() noexcept -> std::span<Residue> { return polarizable_.residues(); }
  auto nonpolarizable_residues() noexcept -> std::span<Residue> {
    return nonpolarizable_.residues();
  }

  /// Adds `force_field` for the residues named force_field.residue() with force-field label
  /// force_field.label(), replacing one added before for the same pair.
  auto add_force_field(ForceField force_field) -> void;

  /// Force field `label` of residue `residue` (lowercase labels, as Residue and ForceField store
  /// them); nullptr if none was added.
  auto force_field(const std::string& label, const std::string& residue) const noexcept
      -> const ForceField*;

  /// Whether both regions are empty.
  auto empty() const noexcept -> bool { return polarizable_.empty() && nonpolarizable_.empty(); }

 private:
  ClassicalRegion polarizable_{true};
  ClassicalRegion nonpolarizable_{false};
  std::map<std::pair<std::string, std::string>, ForceField> force_fields_;  // (label, residue)
};

}  // namespace fika

#endif  // fika_classical_system_hpp
