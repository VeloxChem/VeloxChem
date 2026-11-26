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

#ifndef fika_element_hpp
#define fika_element_hpp

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

#include "fika_element_data.hpp"

namespace fika {

namespace detail {

inline constexpr int max_atomic_number = static_cast<int>(element_data.size());

[[noreturn]] inline void throw_invalid_atomic_number(int atomic_number) {
  throw std::out_of_range("fika::Element: atomic number " + std::to_string(atomic_number) +
                          " outside 1.." + std::to_string(max_atomic_number));
}

[[noreturn]] inline void throw_unknown_label(std::string_view label) {
  throw std::invalid_argument("fika::Element: unknown element label '" + std::string(label) + "'");
}

/// Maps 'a'..'z' and 'A'..'Z' to 0..25, anything else to 26.
constexpr auto letter_index(char c) noexcept -> std::size_t {
  if (c >= 'a' && c <= 'z') {
    return static_cast<std::size_t>(c - 'a');
  }
  if (c >= 'A' && c <= 'Z') {
    return static_cast<std::size_t>(c - 'A');
  }
  return 26;
}

/// Keys 0..26*27-1 encode one- and two-letter labels; the last key marks invalid labels.
inline constexpr std::size_t invalid_label_key = 26 * 27;

/// Case-insensitive key of a label, invalid_label_key if it cannot be an element symbol.
constexpr auto label_key(std::string_view label) noexcept -> std::size_t {
  if (label.empty() || label.size() > 2) {
    return invalid_label_key;
  }
  const std::size_t first = letter_index(label[0]);
  const std::size_t second = label.size() == 2 ? letter_index(label[1]) : 0;
  if (first > 25 || (label.size() == 2 && second > 25)) {
    return invalid_label_key;
  }
  return first * 27 + (label.size() == 2 ? second + 1 : 0);
}

/// Atomic number for each label key, 0 where no element exists (including invalid_label_key).
inline constexpr auto atomic_number_by_label_key = [] {
  std::array<std::uint8_t, invalid_label_key + 1> table{};
  for (std::size_t i = 0; i < element_data.size(); ++i) {
    table[label_key(element_data[i].label)] = static_cast<std::uint8_t>(i + 1);
  }
  return table;
}();

}  // namespace detail

/// Chemical element identified by its atomic number (1..118).
///
/// Stored in a single byte and trivially copyable, so containers of elements are compact.
/// Properties are read from a compile-time table.
class Element {
 public:
  constexpr explicit Element(int atomic_number) : atomic_number_(checked(atomic_number)) {}

  /// Element from its IUPAC symbol, ignoring case ("he", "HE" and "He" all give helium).
  static constexpr auto from_label(std::string_view label) -> Element {
    const int atomic_number = detail::atomic_number_by_label_key[detail::label_key(label)];
    if (atomic_number == 0) {
      detail::throw_unknown_label(label);
    }
    return Element(atomic_number);
  }

  constexpr auto atomic_number() const noexcept -> int { return atomic_number_; }

  /// IUPAC symbol, e.g. "He".
  constexpr auto label() const noexcept -> std::string_view { return data().label; }

  /// Nuclear charge in atomic units (equal to the atomic number).
  constexpr auto charge() const noexcept -> double { return static_cast<double>(atomic_number_); }

  /// Atomic mass in Da of the most abundant isotope, or of the longest-lived isotope for
  /// elements without naturally occurring isotopes (AME2020 / NUBASE2020).
  constexpr auto mass() const noexcept -> double { return data().mass; }

  friend constexpr auto operator<=>(Element, Element) = default;

 private:
  static constexpr auto checked(int atomic_number) -> std::uint8_t {
    if (atomic_number < 1 || atomic_number > detail::max_atomic_number) {
      detail::throw_invalid_atomic_number(atomic_number);
    }
    return static_cast<std::uint8_t>(atomic_number);
  }

  constexpr auto data() const noexcept -> const detail::ElementData& {
    return detail::element_data[atomic_number_ - 1];
  }

  std::uint8_t atomic_number_;
};

}  // namespace fika

#endif  // fika_element_hpp
