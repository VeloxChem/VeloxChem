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

#ifndef fika_molecule_hpp
#define fika_molecule_hpp

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_element.hpp"

namespace fika {

/// Molecule stored as parallel arrays of elements and atomic positions (bohr).
///
/// Both arrays always have the same length; atom i is (elements()[i], coordinates()[i]).
template <real_scalar T>
class Molecule {
 public:
  Molecule() = default;

  /// Takes ownership of both arrays; throws std::invalid_argument if their sizes differ.
  Molecule(std::vector<Element> elements, std::vector<Point3D<T>> coordinates)
      : elements_(std::move(elements)), coordinates_(std::move(coordinates)) {
    if (elements_.size() != coordinates_.size()) {
      throw std::invalid_argument("fika::Molecule: " + std::to_string(elements_.size()) +
                                  " elements but " + std::to_string(coordinates_.size()) +
                                  " coordinates");
    }
  }

  /// Appends an atom; position in bohr. Leaves the molecule unchanged if allocation fails.
  auto add_atom(Element element, Point3D<T> position) -> void {
    elements_.push_back(element);
    try {
      coordinates_.push_back(position);
    } catch (...) {
      elements_.pop_back();
      throw;
    }
  }

  auto reserve(std::size_t count) -> void {
    elements_.reserve(count);
    coordinates_.reserve(count);
  }

  auto size() const noexcept -> std::size_t { return elements_.size(); }

  auto empty() const noexcept -> bool { return elements_.empty(); }

  auto elements() const noexcept -> std::span<const Element> { return elements_; }

  /// Atomic positions in bohr.
  auto coordinates() const noexcept -> std::span<const Point3D<T>> { return coordinates_; }

  /// Centre of mass in bohr, weighting atoms by Element::mass() (most abundant isotope). Throws
  /// std::logic_error for an empty molecule.
  auto centre_of_mass() const -> Point3D<T> {
    if (empty()) {
      throw std::logic_error("fika::Molecule: centre of mass of an empty molecule");
    }
    T total{0};
    Point3D<T> weighted{T{0}, T{0}, T{0}};
    for (std::size_t atom = 0; atom < size(); ++atom) {
      const auto mass = static_cast<T>(elements_[atom].mass());
      total += mass;
      weighted.x += mass * coordinates_[atom].x;
      weighted.y += mass * coordinates_[atom].y;
      weighted.z += mass * coordinates_[atom].z;
    }
    return {weighted.x / total, weighted.y / total, weighted.z / total};
  }

  /// Editable atomic positions in bohr; the number of atoms cannot change through this view.
  auto coordinates() noexcept -> std::span<Point3D<T>> { return coordinates_; }

 private:
  std::vector<Element> elements_;
  std::vector<Point3D<T>> coordinates_;
};

}  // namespace fika

#endif  // fika_molecule_hpp
