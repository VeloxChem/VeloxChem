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

#include "fika_classical_system.hpp"

#include <stdexcept>
#include <utility>

namespace fika {

ClassicalSystem::ClassicalSystem(ClassicalRegion polarizable, ClassicalRegion nonpolarizable)
    : polarizable_(std::move(polarizable)), nonpolarizable_(std::move(nonpolarizable)) {
  if (!polarizable_.polarizable()) {
    throw std::invalid_argument("fika::ClassicalSystem: the polarizable region is not polarizable");
  }
  if (nonpolarizable_.polarizable()) {
    throw std::invalid_argument("fika::ClassicalSystem: the non-polarizable region is polarizable");
  }
}

auto ClassicalSystem::add_residue(Residue residue) -> void {
  (residue.polarizable() ? polarizable_ : nonpolarizable_).add_residue(std::move(residue));
}

auto ClassicalSystem::add_force_field(ForceField force_field) -> void {
  auto key = std::pair{force_field.label(), force_field.residue()};
  force_fields_.insert_or_assign(std::move(key), std::move(force_field));
}

auto ClassicalSystem::force_field(const std::string& label, const std::string& residue) const noexcept
    -> const ForceField* {
  const auto found = force_fields_.find(std::pair{label, residue});
  return found == force_fields_.end() ? nullptr : &found->second;
}

}  // namespace fika
