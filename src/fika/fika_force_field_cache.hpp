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

#ifndef fika_force_field_cache_hpp
#define fika_force_field_cache_hpp

// Internal: force fields of residue types, from the classical system.

#include <stdexcept>
#include <string>
#include <utility>

#include "fika_classical_system.hpp"
#include "fika_force_field.hpp"
#include "fika_residue.hpp"

namespace fika::detail {

class ForceFieldCache {
 public:
  /// `caller` prefixes error messages, e.g. "fika::classical_charges".
  ForceFieldCache(const ClassicalSystem& system, std::string caller)
      : system_(&system), caller_(std::move(caller)) {}

  /// Force field of the residue's type; throws std::runtime_error if the system holds none for it
  /// and std::invalid_argument if the residue's atom count differs from it.
  auto of(const Residue& residue) const -> const ForceField& {
    const ForceField* field = system_->force_field(residue.force_field(), residue.name());
    if (field == nullptr) {
      throw std::runtime_error(caller_ + ": no force field " + residue.force_field() +
                               " for residue " + residue.name());
    }
    if (residue.molecule().size() != field->atom_count()) {
      throw std::invalid_argument(
          caller_ + ": residue " + residue.name() + " " + std::to_string(residue.index()) +
          " has " + std::to_string(residue.molecule().size()) + " atoms, force field " +
          field->label() + " " + std::to_string(field->atom_count()));
    }
    return *field;
  }

 private:
  const ClassicalSystem* system_;
  std::string caller_;
};

}  // namespace fika::detail

#endif  // fika_force_field_cache_hpp
