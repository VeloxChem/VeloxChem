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

#ifndef fika_residue_hpp
#define fika_residue_hpp

#include <cstddef>
#include <string>

#include "fika_molecule.hpp"

namespace fika {

/// A residue of the MM region: its standard residue label (e.g. "hoh"), its 0-based index among
/// the residues with that label, the label of the force field assigned to it, its atoms
/// (elements and coordinates in bohr), and whether that force field is polarizable. Both labels
/// are stored lowercased.
class Residue {
 public:
  /// Lowercases both labels (ASCII). Throws std::invalid_argument if a label is empty or contains
  /// whitespace.
  Residue(std::string name, std::size_t index, std::string force_field, Molecule<double> molecule,
          bool polarizable = false);

  /// Standard residue label, lowercased.
  auto name() const noexcept -> const std::string& { return name_; }

  /// 0-based index among the residues with this label.
  auto index() const noexcept -> std::size_t { return index_; }

  /// Label of the force field assigned to the residue, lowercased.
  auto force_field() const noexcept -> const std::string& { return force_field_; }

  /// Whether the residue's force field is polarizable.
  auto polarizable() const noexcept -> bool { return polarizable_; }

  auto molecule() const noexcept -> const Molecule<double>& { return molecule_; }
  auto molecule() noexcept -> Molecule<double>& { return molecule_; }

 private:
  std::string name_;
  std::size_t index_;
  std::string force_field_;
  Molecule<double> molecule_;
  bool polarizable_;
};

}  // namespace fika

#endif  // fika_residue_hpp
