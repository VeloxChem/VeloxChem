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

#include "fika_classical_region.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace fika {

ClassicalRegion::ClassicalRegion(bool polarizable, std::vector<Residue> residues)
    : polarizable_(polarizable), residues_(std::move(residues)) {
  for (const Residue& residue : residues_) {
    check(residue);
  }
}

auto ClassicalRegion::add_residue(Residue residue) -> void {
  check(residue);
  residues_.push_back(std::move(residue));
}

void ClassicalRegion::check(const Residue& residue) const {
  if (residue.polarizable() != polarizable_) {
    throw std::invalid_argument(
        "fika::ClassicalRegion: residue " + residue.name() + " " + std::to_string(residue.index()) +
        (residue.polarizable() ? " is polarizable" : " is not polarizable") + ", the region " +
        (polarizable_ ? "is" : "is not"));
  }
}

}  // namespace fika
