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

#ifndef fika_text_parsing_hpp
#define fika_text_parsing_hpp

// Internal helpers for line-oriented text formats (XYZ, basis-set files).

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace fika::detail {

/// Splits text into lines, dropping a trailing '\r', and counts line numbers from 1.
class LineReader {
 public:
  explicit LineReader(std::string_view text) : text_(text) {}

  /// Next line, or false at the end of the text.
  auto next(std::string_view& line) -> bool {
    if (position_ > text_.size()) {
      return false;
    }
    const std::size_t end = std::min(text_.find('\n', position_), text_.size());
    line = text_.substr(position_, end - position_);
    if (!line.empty() && line.back() == '\r') {
      line.remove_suffix(1);
    }
    position_ = end + 1;
    ++line_number_;
    return true;
  }

  auto line_number() const noexcept -> std::size_t { return line_number_; }

 private:
  std::string_view text_;
  std::size_t position_ = 0;
  std::size_t line_number_ = 0;
};

constexpr auto is_blank(char c) noexcept -> bool {
  return c == ' ' || c == '\t';
}

inline auto is_blank_line(std::string_view line) -> bool {
  return std::ranges::all_of(line, is_blank);
}

/// Calls field(std::string_view) for each whitespace-separated field; stops early and returns
/// false as soon as field returns false.
template <typename F>
auto for_each_field(std::string_view line, F&& field) -> bool {
  std::size_t position = 0;
  while (true) {
    while (position < line.size() && is_blank(line[position])) {
      ++position;
    }
    if (position == line.size()) {
      return true;
    }
    const std::size_t start = position;
    while (position < line.size() && !is_blank(line[position])) {
      ++position;
    }
    if (!field(line.substr(start, position - start))) {
      return false;
    }
  }
}

/// Splits a line into at most N whitespace-separated fields; returns the number found, or N + 1
/// if the line holds more than N fields.
template <std::size_t N>
auto split_fields(std::string_view line, std::array<std::string_view, N>& fields) -> std::size_t {
  std::size_t count = 0;
  const bool fits = for_each_field(line, [&](std::string_view field) {
    if (count == N) {
      return false;
    }
    fields[count++] = field;
    return true;
  });
  return fits ? count : N + 1;
}

/// Splits a line into all its whitespace-separated fields, reusing the vector's storage.
inline void split_fields(std::string_view line, std::vector<std::string_view>& fields) {
  fields.clear();
  for_each_field(line, [&](std::string_view field) {
    fields.push_back(field);
    return true;
  });
}

/// ASCII lower case (independent of the locale; other bytes unchanged).
inline auto lower_case(std::string_view text) -> std::string {
  std::string result(text);
  std::ranges::transform(result, result.begin(), [](char c) {
    return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
  });
  return result;
}

/// Non-negative integer, or nullopt if the field is not one.
inline auto parse_count(std::string_view field) -> std::optional<std::size_t> {
  std::size_t value = 0;
  const auto [end, error] = std::from_chars(field.data(), field.data() + field.size(), value);
  if (error != std::errc{} || end != field.data() + field.size()) {
    return std::nullopt;
  }
  return value;
}

/// Decimal exponent of the leading digit of a number from_chars matched (sign, digits, optional
/// point and fraction, optional exponent), e.g. 2 for "123.4", -3 for "0.0012e0"; the exponent
/// part saturates instead of overflowing.
inline auto decimal_magnitude(std::string_view number) -> long long {
  std::size_t i = number.empty() || number.front() != '-' ? 0 : 1;
  long long magnitude = -1;  // digits before the point, less one
  bool leading = true;       // still in leading zeros
  for (; i < number.size() && number[i] >= '0' && number[i] <= '9'; ++i) {
    if (!leading || number[i] != '0') {
      leading = false;
      ++magnitude;
    }
  }
  if (i < number.size() && number[i] == '.') {
    for (++i; i < number.size() && number[i] >= '0' && number[i] <= '9'; ++i) {
      if (leading && number[i] == '0') {
        --magnitude;  // 0.00x: leading fractional zeros
      } else {
        leading = false;
      }
    }
  }
  long long exponent = 0;
  if (i < number.size() && (number[i] == 'e' || number[i] == 'E')) {
    ++i;
    const bool negative = i < number.size() && number[i] == '-';
    if (i < number.size() && (number[i] == '-' || number[i] == '+')) {
      ++i;
    }
    for (; i < number.size() && number[i] >= '0' && number[i] <= '9'; ++i) {
      exponent = std::min(exponent * 10 + (number[i] - '0'), 1'000'000'000LL);
    }
    exponent = negative ? -exponent : exponent;
  }
  return magnitude + exponent;
}

/// Finite floating-point number (a leading '+' is accepted), or nullopt if the field is not one.
/// Values below the smallest subnormal read as zero of their sign (as strtod gives); values
/// beyond the largest double are rejected.
inline auto parse_real(std::string_view field) -> std::optional<double> {
  // from_chars rejects a leading '+', which some programs write.
  if (field.size() > 1 && field.front() == '+' && field[1] != '-') {
    field.remove_prefix(1);
  }
  double value = 0.0;
  const auto [end, error] = std::from_chars(field.data(), field.data() + field.size(), value);
  if (error == std::errc::result_out_of_range && end == field.data() + field.size() &&
      decimal_magnitude(field) < 0) {
    return field.front() == '-' ? -0.0 : 0.0;  // underflow
  }
  if (error != std::errc{} || end != field.data() + field.size() || !std::isfinite(value)) {
    return std::nullopt;
  }
  return value;
}

/// Whole file as a string, or nullopt if it cannot be opened.
inline auto read_text_file(const std::filesystem::path& path) -> std::optional<std::string> {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    return std::nullopt;
  }
  return std::string{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

}  // namespace fika::detail

#endif  // fika_text_parsing_hpp
