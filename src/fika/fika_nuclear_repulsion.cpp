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

#include "fika_nuclear_repulsion.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_compensated_sum.hpp"
#include "fika_point3d.hpp"
#include "fika_element.hpp"

namespace fika {

namespace {

/// Below this many atom pairs threading costs more than it saves (measured on Apple M-series).
constexpr std::size_t parallel_threshold_pairs = 100'000;

/// Atoms of the larger molecule per intermolecular work block: aims at ~256 blocks for load
/// balance, bounded to keep per-block overhead low and blocks cache sized. Depends only on the
/// molecule size, so the summation order is independent of the thread count.
auto intermolecular_block_size(std::size_t large_size) -> std::size_t {
  return std::clamp<std::size_t>((large_size + 255) / 256, 64, 4096);
}

template <real_scalar T>
auto distance(const Point3D<T>& a, const Point3D<T>& b) -> double {
  const double dx = static_cast<double>(a.x) - static_cast<double>(b.x);
  const double dy = static_cast<double>(a.y) - static_cast<double>(b.y);
  const double dz = static_cast<double>(a.z) - static_cast<double>(b.z);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// Coincident nuclei give Z / 0 = inf in the sum; checked once instead of per pair.
auto checked(double energy) -> double {
  if (!std::isfinite(energy)) {
    throw std::domain_error("fika::nuclear_repulsion_energy: coincident nuclei");
  }
  return energy;
}

/// Z_i * sum over j < i of Z_j / r_ij.
template <real_scalar T>
auto row_energy(const Molecule<T>& molecule, std::size_t i) -> double {
  const auto elements = molecule.elements();
  const auto coordinates = molecule.coordinates();
  double row = 0.0;
  for (std::size_t j = 0; j < i; ++j) {
    row += elements[j].charge() / distance(coordinates[i], coordinates[j]);
  }
  return elements[i].charge() * row;
}

/// Interaction of all atoms of `small` with atoms [begin, end) of `large`.
template <real_scalar T>
auto block_energy(const Molecule<T>& small, const Molecule<T>& large, std::size_t begin,
                  std::size_t end) -> double {
  const auto elements_small = small.elements();
  const auto coordinates_small = small.coordinates();
  const auto elements_large = large.elements();
  const auto coordinates_large = large.coordinates();
  double energy = 0.0;
  for (std::size_t i = 0; i < small.size(); ++i) {
    double row = 0.0;
    for (std::size_t j = begin; j < end; ++j) {
      row += elements_large[j].charge() / distance(coordinates_small[i], coordinates_large[j]);
    }
    energy += elements_small[i].charge() * row;
  }
  return energy;
}

}  // namespace

// Row energies are summed in row order after the parallel loop, so the result does not depend
// on the number of threads and equals the serial sum.
template <real_scalar T>
auto nuclear_repulsion_energy(const Molecule<T>& molecule) -> double {
  molecule.check_finite_coordinates("fika::nuclear_repulsion_energy");
  const std::size_t n = molecule.size();
  if (n < 2) {
    return 0.0;
  }

  double energy = 0.0;
  if (n * (n - 1) / 2 < parallel_threshold_pairs) {
    for (std::size_t i = 1; i < n; ++i) {
      energy += row_energy(molecule, i);
    }
    return checked(energy);
  }

  std::vector<double> rows(n);
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t i = 1; i < n; ++i) {
    rows[i] = row_energy(molecule, i);
  }
  for (std::size_t i = 1; i < n; ++i) {
    energy += rows[i];
  }
  return checked(energy);
}

// The larger molecule is split into fixed-size blocks whose energies are summed in block order,
// so the result does not depend on the number of threads.
template <real_scalar T>
auto nuclear_repulsion_energy(const Molecule<T>& a, const Molecule<T>& b) -> double {
  a.check_finite_coordinates("fika::nuclear_repulsion_energy");
  b.check_finite_coordinates("fika::nuclear_repulsion_energy");
  const bool a_is_small = a.size() <= b.size();
  const Molecule<T>& small = a_is_small ? a : b;
  const Molecule<T>& large = a_is_small ? b : a;

  const std::size_t block_size = intermolecular_block_size(large.size());
  const std::size_t block_count = (large.size() + block_size - 1) / block_size;
  if (block_count <= 1) {
    return checked(block_energy(small, large, 0, large.size()));
  }

  std::vector<double> blocks(block_count);
  const bool parallel = small.size() * large.size() >= parallel_threshold_pairs;
#pragma omp parallel for schedule(dynamic, 1) if (parallel)
  for (std::size_t block = 0; block < block_count; ++block) {
    const std::size_t begin = block * block_size;
    const std::size_t end = std::min(begin + block_size, large.size());
    blocks[block] = block_energy(small, large, begin, end);
  }

  double energy = 0.0;
  for (const double block : blocks) {
    energy += block;
  }
  return checked(energy);
}

namespace {

/// Charges per block of nuclear_point_charge_energy (fixed, so the summation order is too).
constexpr std::size_t point_charge_block_size = 4096;

/// Compensated sum of q_C sum_A Z_A / |R_A - C| over charges [begin, end); the inner sum has
/// only positive terms and is summed plainly.
auto point_charge_block_energy(std::span<const double> x, std::span<const double> y,
                               std::span<const double> z, std::span<const double> nuclear_charges,
                               std::span<const double> charges,
                               std::span<const Point3D<double>> coordinates, std::size_t begin,
                               std::size_t end) -> double {
  detail::NeumaierSum energy;
  const std::size_t atoms = nuclear_charges.size();
  for (std::size_t c = begin; c < end; ++c) {
    if (charges[c] == 0.0) {
      continue;  // a zero charge on a nucleus would give 0 * inf
    }
    const Point3D<double>& position = coordinates[c];
    double potential = 0.0;  // sum_A Z_A / |R_A - C|
    for (std::size_t a = 0; a < atoms; ++a) {
      const double dx = x[a] - position.x;
      const double dy = y[a] - position.y;
      const double dz = z[a] - position.z;
      potential += nuclear_charges[a] / std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    energy.add(charges[c] * potential);
  }
  return energy.value();
}

}  // namespace

auto nuclear_point_charge_energy(const Molecule<double>& molecule, std::span<const double> charges,
                                 std::span<const Point3D<double>> coordinates) -> double {
  if (charges.size() != coordinates.size()) {
    throw std::invalid_argument(
        "fika::nuclear_point_charge_energy: " + std::to_string(charges.size()) + " charges but " +
        std::to_string(coordinates.size()) + " coordinates");
  }
  const auto finite = [](const Point3D<double>& p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
  };
  if (!std::ranges::all_of(charges, [](double q) { return std::isfinite(q); }) ||
      !std::ranges::all_of(coordinates, finite)) {
    throw std::invalid_argument(
        "fika::nuclear_point_charge_energy: charges and coordinates must be finite");
  }
  molecule.check_finite_coordinates("fika::nuclear_point_charge_energy");
  // Nuclei as contiguous arrays for the inner loop.
  const std::size_t atoms = molecule.size();
  std::vector<double> x(atoms), y(atoms), z(atoms), nuclear_charges(atoms);
  for (std::size_t a = 0; a < atoms; ++a) {
    x[a] = molecule.coordinates()[a].x;
    y[a] = molecule.coordinates()[a].y;
    z[a] = molecule.coordinates()[a].z;
    nuclear_charges[a] = molecule.elements()[a].charge();
  }
  const std::size_t block_count =
      (charges.size() + point_charge_block_size - 1) / point_charge_block_size;
  std::vector<double> blocks(block_count);
  const bool parallel = atoms * charges.size() >= parallel_threshold_pairs;
#pragma omp parallel for schedule(dynamic, 1) if (parallel)
  for (std::size_t block = 0; block < block_count; ++block) {
    const std::size_t begin = block * point_charge_block_size;
    const std::size_t end = std::min(begin + point_charge_block_size, charges.size());
    blocks[block] =
        point_charge_block_energy(x, y, z, nuclear_charges, charges, coordinates, begin, end);
  }
  detail::NeumaierSum energy;
  for (const double block : blocks) {
    energy.add(block);
  }
  const double result = energy.value();
  if (!std::isfinite(result)) {
    throw std::domain_error("fika::nuclear_point_charge_energy: a charge sits on a nucleus");
  }
  return result;
}

template auto nuclear_repulsion_energy(const Molecule<float>&) -> double;
template auto nuclear_repulsion_energy(const Molecule<double>&) -> double;
template auto nuclear_repulsion_energy(const Molecule<float>&, const Molecule<float>&) -> double;
template auto nuclear_repulsion_energy(const Molecule<double>&, const Molecule<double>&) -> double;

}  // namespace fika
