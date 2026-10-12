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
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef fika_switching_hpp
#define fika_switching_hpp

#include <cmath>
#include <stdexcept>
#include <string>

namespace fika {

/// Inner switching function of the short-range coupling of pairs of different residues, e.g. of
/// a permanent charge with an induced dipole or of two induced dipoles:
///   S(r) = 0 for r <= r_on, 1 for r >= r_off, and
///   S(r) = 1/2 [1 - cos(pi (r - r_on) / (r_off - r_on))] in between,
/// so S and its first derivative are continuous and S = 0 removes the pair while S = 1 leaves
/// the bare Coulomb coupling. Distances in bohr; the radii exclude short-range inter-residue
/// couplings by distance (typical values: 1.0 and 1.8 Angstrom).
struct InnerSwitch {
  double r_on = 0.0;
  double r_off = 0.0;
};

/// Throws std::invalid_argument (prefixed with `caller`) unless 0 <= r_on < r_off and both are
/// finite.
inline auto check_inner_switch(const InnerSwitch& switching, const char* caller) -> void {
  if (!std::isfinite(switching.r_on) || !std::isfinite(switching.r_off) ||
      switching.r_on < 0.0 || !(switching.r_off > switching.r_on)) {
    throw std::invalid_argument(std::string(caller) +
                                ": switching radii need 0 <= r_on < r_off");
  }
}

/// Value S(r); `switching` must satisfy check_inner_switch.
inline auto inner_switch(double r, const InnerSwitch& switching) noexcept -> double {
  if (r <= switching.r_on) {
    return 0.0;
  }
  if (r >= switching.r_off) {
    return 1.0;
  }
  constexpr double pi = 3.14159265358979323846;
  return 0.5 * (1.0 - std::cos(pi * (r - switching.r_on) / (switching.r_off - switching.r_on)));
}

}  // namespace fika

#endif  // fika_switching_hpp
