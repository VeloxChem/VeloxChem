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

#include "fika_quadrupole_potential_driver.hpp"

#include "fika_point_source_operator.hpp"

namespace fika {

QuadrupolePotentialDriver::QuadrupolePotentialDriver(std::size_t block_size,
                                                     ChargeSummation summation)
    : block_size_(block_size), summation_(summation) {}

auto QuadrupolePotentialDriver::compute(const Molecule<double>& molecule,
                                        const MolecularBasis& basis,
                                        std::span<const Quadrupole> quadrupoles,
                                        std::span<const Point3D<double>> quadrupole_coordinates,
                                        double threshold) const -> BlockSparseMatrix {
  return detail::point_source_matrix(
      molecule, basis,
      detail::PointSources{.quadrupoles = quadrupoles,
                           .quadrupole_coordinates = quadrupole_coordinates},
      detail::PointSourceKind::quadrupoles, threshold,
      {.caller = "fika::QuadrupolePotentialDriver",
       .block_size = block_size_,
       .summation = summation_,
       .multipole_count = multipole_quadrupole_count});
}

}  // namespace fika
