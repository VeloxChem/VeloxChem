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

#include "fika_far_field_expansion.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "fika_overlap_screening.hpp"
#include "fika_multipole_expansion.hpp"
#include "fika_octree_keys.hpp"

namespace fika::detail {

namespace {

constexpr int max_target_depth = 8;

auto distance(const Point3D<double>& a, const Point3D<double>& b) -> double {
  return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                   (a.z - b.z) * (a.z - b.z));
}

// Whether the segment a..b meets the box [low, high].
auto segment_meets_box(const Point3D<double>& a, const Point3D<double>& b,
                       const Point3D<double>& low, const Point3D<double>& high) -> bool {
  double t_low = 0.0;
  double t_high = 1.0;
  const double start[3] = {a.x, a.y, a.z};
  const double step[3] = {b.x - a.x, b.y - a.y, b.z - a.z};
  const double box_low[3] = {low.x, low.y, low.z};
  const double box_high[3] = {high.x, high.y, high.z};
  for (int axis = 0; axis < 3; ++axis) {
    if (step[axis] == 0.0) {
      if (start[axis] < box_low[axis] || start[axis] > box_high[axis]) {
        return false;
      }
      continue;
    }
    double t0 = (box_low[axis] - start[axis]) / step[axis];
    double t1 = (box_high[axis] - start[axis]) / step[axis];
    if (t0 > t1) {
      std::swap(t0, t1);
    }
    t_low = std::max(t_low, t0);
    t_high = std::min(t_high, t1);
    if (t_low > t_high) {
      return false;
    }
  }
  return true;
}

// Smallest order p (rank <= p <= max_order) at which the error bound per unit charge of a
// translation (radii sum `radii`, centre distance d) is at most tolerance[Lambda] for every
// Lambda <= rank, or -1 if the translation is not admissible.
auto required_order(double radii, double d, double penetration_radius,
                    std::span<const double> tolerances, int max_order, double tolerance_scale = 1.0)
    -> int {
  if (d - radii < penetration_radius || radii >= d) {
    return -1;
  }
  const int rank = static_cast<int>(tolerances.size()) - 1;
  const double theta = radii / d;
  const double scale = 1.0 / ((1.0 - theta) * d);
  if (std::pow(theta, max_order + 1) * scale > tolerance_scale * tolerances[0]) {
    return -1;  // Lambda = 0 fails even at max_order
  }
  int order = rank;
  for (int big_l = 0; big_l <= rank; ++big_l) {
    // bound(p) = C(p + 1, Lambda) theta^(p - Lambda + 1) scale^(Lambda + 1) from p = rank.
    double binomial = 1.0;
    for (int i = 1; i <= big_l; ++i) {
      binomial = binomial * (rank + 2 - i) / i;
    }
    double bound = binomial * std::pow(theta, rank - big_l + 1) * std::pow(scale, big_l + 1);
    int p = rank;
    while (bound > tolerance_scale * tolerances[static_cast<std::size_t>(big_l)]) {
      if (p == max_order) {
        return -1;
      }
      bound *= theta * (p + 2) / (p + 2 - big_l);
      ++p;
    }
    order = std::max(order, p);
  }
  return order;
}

// Smallest order at which point sources of kind `kind` (1 dipoles, 2 quadrupoles) within
// `source_radius` of O are admissible for a cell of radius `target_radius` at centre distance d
// (tolerances per unit dipole or quadrupole), or -1: the smallest over the bound's
// delta = gap / 2^k of the charge requirement for radius source_radius + delta with the
// tolerances scaled by delta / 3 (dipole_field_tensor_error_bound) or
// delta^2 3 sqrt 3 / 20 (quadrupole_field_tensor_error_bound).
auto derivative_required_order(std::size_t kind, double source_radius, double target_radius,
                               double d, double penetration_radius,
                               std::span<const double> tolerances, int max_order) -> int {
  const double gap = d - source_radius - target_radius;
  if (gap < penetration_radius || gap <= 0.0) {
    return -1;
  }
  int best = -1;
  for (int k = 1; k <= 10; ++k) {
    const double delta = gap / static_cast<double>(1 << k);
    const double scale = kind == 1 ? delta / 3.0 : delta * delta * 3.0 * std::sqrt(3.0) / 20.0;
    const int order =
        required_order(source_radius + delta + target_radius, d, 0.0, tolerances, max_order, scale);
    if (order >= 0 && (best < 0 || order < best)) {
      best = order;
    }
  }
  return best;
}

// Distances from which a single site (r_S = 0) is admissible for a cell of radius `radius`, per
// order: entry p - rank is the smallest distance d (found by bisection, rounded up) with
// order_at(d) in [0, p], or infinity if none; order_at(d) (the site's required order, -1 if not
// admissible) decreases with d.
template <typename OrderAt>
auto admissible_distances(double radius, double penetration_radius, int rank, int max_order,
                          OrderAt&& order_at) -> std::vector<double> {
  std::vector<double> distances;
  for (int p = rank; p <= max_order; ++p) {
    const auto admissible = [&](double d) {
      const int order = order_at(d);
      return order >= 0 && order <= p;
    };
    double low = radius + penetration_radius;
    double high = std::max(2.0 * low, 1.0);
    while (!admissible(high) && high < 1e12) {
      low = high;
      high *= 2.0;
    }
    if (!admissible(high)) {
      distances.push_back(std::numeric_limits<double>::infinity());
      continue;
    }
    for (int iteration = 0; iteration < 100 && low < high; ++iteration) {
      const double middle = 0.5 * (low + high);
      if (middle <= low || middle >= high) {
        break;
      }
      (admissible(middle) ? high : low) = middle;
    }
    distances.push_back(high);
  }
  return distances;
}

// CSR from (key, value) pairs in insertion order.
void to_csr(const std::vector<std::pair<std::size_t, std::uint32_t>>& pairs, std::size_t keys,
            std::vector<std::size_t>& offsets, std::vector<std::uint32_t>& values) {
  offsets.assign(keys + 1, 0);
  for (const auto& [key, value] : pairs) {
    ++offsets[key + 1];
  }
  for (std::size_t k = 0; k < keys; ++k) {
    offsets[k + 1] += offsets[k];
  }
  values.resize(pairs.size());
  std::vector<std::size_t> next(offsets.begin(), offsets.end() - 1);
  for (const auto& [key, value] : pairs) {
    values[next[key]++] = value;
  }
}

/// Norm sums sum |q|, sum |mu| and sum ||Theta||_2 of the sources.
auto norm_sums(const PointSources& sources) -> std::array<double, 3> {
  std::array<double, 3> sums{};
  for (const double q : sources.charges) {
    sums[0] += std::abs(q);
  }
  for (const Dipole& mu : sources.dipoles) {
    sums[1] += std::hypot(mu(0), mu(1), mu(2));
  }
  for (const Quadrupole& q : sources.quadrupoles) {
    sums[2] += traceless_quadrupole_norm(q);
  }
  return sums;
}

}  // namespace

FarFieldExpansion::FarFieldExpansion(std::span<const double> charges,
                                     std::span<const Point3D<double>> coordinates,
                                     std::span<const Point3D<double>> atoms,
                                     std::span<const AtomPair> pairs,
                                     const FarFieldOptions& options)
    : FarFieldExpansion(PointSources{.charges = charges, .charge_coordinates = coordinates}, atoms,
                        pairs, options) {}

FarFieldExpansion::FarFieldExpansion(const PointSources& sources,
                                     std::span<const Point3D<double>> atoms,
                                     std::span<const AtomPair> pairs,
                                     const FarFieldOptions& options) {
  if (sources.charges.size() != sources.charge_coordinates.size()) {
    throw std::invalid_argument("FarFieldExpansion: charge and coordinate counts differ");
  }
  if (sources.dipoles.size() != sources.dipole_coordinates.size()) {
    throw std::invalid_argument("FarFieldExpansion: dipole and coordinate counts differ");
  }
  if (sources.quadrupoles.size() != sources.quadrupole_coordinates.size()) {
    throw std::invalid_argument("FarFieldExpansion: quadrupole and coordinate counts differ");
  }
  if (sources.charges.size() + sources.dipoles.size() + sources.quadrupoles.size() >
      std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("FarFieldExpansion: too many sources");
  }
  if (options.rank < 0 || options.tolerances.size() != static_cast<std::size_t>(options.rank) + 1 ||
      std::ranges::any_of(
          options.tolerances, [](double t) { return !(t > 0.0) || !std::isfinite(t); })) {
    throw std::invalid_argument("FarFieldExpansion: need a positive tolerance per rank");
  }
  const int highest_order = !sources.quadrupoles.empty() ? max_expansion_order - 2
                            : !sources.dipoles.empty()   ? max_expansion_order - 1
                                                         : max_expansion_order;
  if (!(options.penetration_radius >= 0.0) || !std::isfinite(options.penetration_radius) ||
      !(options.leaf_edge > 0.0) || !std::isfinite(options.leaf_edge) ||
      options.leaf_capacity == 0 || options.max_order < options.rank ||
      options.max_order > highest_order) {
    throw std::invalid_argument("FarFieldExpansion: invalid options");
  }
  for (const auto& pair : pairs) {
    if (pair.bra >= atoms.size() || pair.ket >= atoms.size()) {
      throw std::invalid_argument("FarFieldExpansion: atom pair out of range");
    }
  }
  if (options.source_norms && std::ranges::any_of(*options.source_norms, [](double v) {
        return !(v >= 0.0) || !std::isfinite(v);
      })) {
    throw std::invalid_argument("FarFieldExpansion: source norms must be finite and non-negative");
  }
  norms_ = options.source_norms.value_or(norm_sums(sources));
  counts_ = {sources.charges.size(), sources.dipoles.size(), sources.quadrupoles.size()};
  build_sources(sources, options.leaf_capacity);
  build_targets(atoms, pairs, options.leaf_edge);
  traverse(options, norms_);
  upward_pass();
  downward_pass();
}

void FarFieldExpansion::update(const PointSources& sources) {
  const std::array<std::size_t, kind_count> counts{sources.charges.size(), sources.dipoles.size(),
                                                   sources.quadrupoles.size()};
  const std::array<std::span<const Point3D<double>>, kind_count> coordinates{
      sources.charge_coordinates, sources.dipole_coordinates, sources.quadrupole_coordinates};
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    if (counts[kind] != counts_[kind] || coordinates[kind].size() != counts_[kind]) {
      throw std::invalid_argument("FarFieldExpansion::update: source counts differ");
    }
  }
  for (std::size_t i = 0; i < sorted_coordinates_.size(); ++i) {
    const std::size_t kind = mixed_ ? static_cast<std::size_t>(sorted_kind_[i]) : charge_kind;
    if (coordinates[kind][sorted_index_[i]] != sorted_coordinates_[i]) {
      throw std::invalid_argument("FarFieldExpansion::update: source positions differ");
    }
  }
  const auto sums = norm_sums(sources);
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    if (norms_[kind] == 0.0 && sums[kind] > 0.0) {
      throw std::invalid_argument(
          "FarFieldExpansion::update: a source kind had zero norm when the lists were built");
    }
  }
  set_values(sources);
  upward_pass();
  downward_pass();
}

void FarFieldExpansion::set_values(const PointSources& sources) {
  dipoles_.assign(sources.dipoles.begin(), sources.dipoles.end());
  quadrupoles_.assign(sources.quadrupoles.begin(), sources.quadrupoles.end());
  for (std::size_t i = 0; i < sorted_charges_.size(); ++i) {
    const bool charge = !mixed_ || static_cast<std::size_t>(sorted_kind_[i]) == charge_kind;
    sorted_charges_[i] = charge ? sources.charges[sorted_index_[i]] : 0.0;
  }
}

void FarFieldExpansion::build_sources(const PointSources& sources, std::size_t capacity) {
  // Sites: the charges, then the dipoles, then the quadrupoles.
  const std::size_t charge_count = sources.charges.size();
  const std::size_t dipole_end = charge_count + sources.dipoles.size();
  const std::size_t n = dipole_end + sources.quadrupoles.size();
  mixed_ = n > charge_count;
  if (n == 0) {
    return;
  }
  dipoles_.assign(sources.dipoles.begin(), sources.dipoles.end());
  quadrupoles_.assign(sources.quadrupoles.begin(), sources.quadrupoles.end());
  std::vector<Point3D<double>> all;
  std::span<const Point3D<double>> coordinates = sources.charge_coordinates;
  if (mixed_) {
    all.reserve(n);
    all.insert(all.end(), sources.charge_coordinates.begin(), sources.charge_coordinates.end());
    all.insert(all.end(), sources.dipole_coordinates.begin(), sources.dipole_coordinates.end());
    all.insert(all.end(), sources.quadrupole_coordinates.begin(),
               sources.quadrupole_coordinates.end());
    coordinates = all;
  }
  const auto cube = bounding_cube(coordinates);  // no structured binding: OpenMP captures these
  const Point3D<double> centre = cube.first;
  const double half = cube.second;
  source_half_edge_ = half;
  const double scale = static_cast<double>(1 << source_depth) / (2.0 * half);
  const auto cell = [&](double value, double low) {
    const double index = std::floor((value - low) * scale);
    return static_cast<std::uint64_t>(
        std::clamp(index, 0.0, static_cast<double>((1 << source_depth) - 1)));
  };
  std::vector<std::pair<std::uint64_t, std::uint32_t>> keys(n);
#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < n; ++c) {
    keys[c] = {
        morton(cell(coordinates[c].x, centre.x - half), cell(coordinates[c].y, centre.y - half),
               cell(coordinates[c].z, centre.z - half)),
        static_cast<std::uint32_t>(c)};
  }
  radix_sort(keys);
  sorted_charges_.resize(n);
  sorted_coordinates_.resize(n);
  sorted_index_.resize(n);
  if (mixed_) {
    sorted_kind_.resize(n);
  }
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t site = keys[i].second;
    const std::size_t kind = site < charge_count ? charge_kind
                             : site < dipole_end ? dipole_kind
                                                 : quadrupole_kind;
    const std::size_t first = kind == charge_kind   ? 0
                              : kind == dipole_kind ? charge_count
                                                    : dipole_end;
    sorted_index_[i] = static_cast<std::uint32_t>(site - first);
    sorted_charges_[i] = kind == charge_kind ? sources.charges[site] : 0.0;
    sorted_coordinates_[i] = coordinates[site];
    if (mixed_) {
      sorted_kind_[i] = static_cast<char>(kind);
    }
  }

  // Cells in creation order: a cell's children are appended together, after it.
  Cell root;
  root.centre = centre;
  root.count = n;
  sources_.push_back(root);
  std::vector<double> half_edges{half};
  for (std::size_t index = 0; index < sources_.size(); ++index) {
    const Cell cell_copy = sources_[index];
    const double child_half = 0.5 * half_edges[index];
    if (cell_copy.count <= capacity || cell_copy.level == source_depth) {
      continue;
    }
    const int shift = 3 * (source_depth - cell_copy.level - 1);
    sources_[index].first_child = sources_.size();
    std::size_t begin = cell_copy.first;
    const std::size_t end = cell_copy.first + cell_copy.count;
    while (begin < end) {
      const std::uint64_t octant = (keys[begin].first >> shift) & 7;
      std::size_t stop = begin;
      while (stop < end && ((keys[stop].first >> shift) & 7) == octant) {
        ++stop;
      }
      Cell child;
      child.centre = {cell_copy.centre.x + ((octant & 4) != 0 ? child_half : -child_half),
                      cell_copy.centre.y + ((octant & 2) != 0 ? child_half : -child_half),
                      cell_copy.centre.z + ((octant & 1) != 0 ? child_half : -child_half)};
      child.first = begin;
      child.count = stop - begin;
      child.parent = index;
      child.level = cell_copy.level + 1;
      sources_.push_back(child);
      half_edges.push_back(child_half);
      ++sources_[index].children;
      begin = stop;
    }
  }
#pragma omp parallel for schedule(dynamic, 16)
  for (std::size_t s = 0; s < sources_.size(); ++s) {
    Cell& cell_entry = sources_[s];
    double radius = 0.0;
    for (std::size_t i = cell_entry.first; i < cell_entry.first + cell_entry.count; ++i) {
      radius = std::max(radius, distance(sorted_coordinates_[i], cell_entry.centre));
    }
    cell_entry.radius = radius;
    cell_entry.kinds[charge_kind] = !mixed_;
    if (mixed_) {
      for (std::size_t i = cell_entry.first; i < cell_entry.first + cell_entry.count; ++i) {
        cell_entry.kinds[static_cast<std::size_t>(sorted_kind_[i])] = true;
      }
    }
  }
  statistics_.source_cells = sources_.size();
}

void FarFieldExpansion::build_targets(std::span<const Point3D<double>> atoms,
                                      std::span<const AtomPair> pairs, double leaf_edge) {
  if (atoms.empty() || pairs.empty()) {
    return;
  }
  const auto [centre, half] = bounding_cube(atoms);
  int depth = 0;
  while (depth < max_target_depth && 2.0 * half / static_cast<double>(1 << depth) > leaf_edge) {
    ++depth;
  }
  grid_size_ = std::size_t{1} << depth;
  grid_edge_ = 2.0 * half / static_cast<double>(grid_size_);
  grid_origin_ = {centre.x - half, centre.y - half, centre.z - half};
  const auto grid = static_cast<std::int64_t>(grid_size_);

  // Active leaves: grid cells (boxes padded against rounding) crossed by an atom-pair segment.
  std::vector<char> active(grid_size_ * grid_size_ * grid_size_, 0);
  const double pad = 1e-9 * grid_edge_;
  const auto index_of = [&](double value, double origin) {
    return std::clamp(static_cast<std::int64_t>(std::floor((value - origin) / grid_edge_)),
                      std::int64_t{0}, grid - 1);
  };
  for (const auto& pair : pairs) {
    const Point3D<double>& a = atoms[pair.bra];
    const Point3D<double>& b = atoms[pair.ket];
    const std::int64_t low[3] = {index_of(std::min(a.x, b.x) - pad, grid_origin_.x),
                                 index_of(std::min(a.y, b.y) - pad, grid_origin_.y),
                                 index_of(std::min(a.z, b.z) - pad, grid_origin_.z)};
    const std::int64_t high[3] = {index_of(std::max(a.x, b.x) + pad, grid_origin_.x),
                                  index_of(std::max(a.y, b.y) + pad, grid_origin_.y),
                                  index_of(std::max(a.z, b.z) + pad, grid_origin_.z)};
    for (std::int64_t ix = low[0]; ix <= high[0]; ++ix) {
      for (std::int64_t iy = low[1]; iy <= high[1]; ++iy) {
        for (std::int64_t iz = low[2]; iz <= high[2]; ++iz) {
          const std::size_t linear = static_cast<std::size_t>((ix * grid + iy) * grid + iz);
          if (active[linear] != 0) {
            continue;
          }
          const Point3D<double> box_low{
              grid_origin_.x + static_cast<double>(ix) * grid_edge_ - pad,
              grid_origin_.y + static_cast<double>(iy) * grid_edge_ - pad,
              grid_origin_.z + static_cast<double>(iz) * grid_edge_ - pad};
          const Point3D<double> box_high{box_low.x + grid_edge_ + 2 * pad,
                                         box_low.y + grid_edge_ + 2 * pad,
                                         box_low.z + grid_edge_ + 2 * pad};
          if (segment_meets_box(a, b, box_low, box_high)) {
            active[linear] = 1;
          }
        }
      }
    }
  }

  // Octree cells by level, root first; keys are Morton codes of the integer cell coordinates.
  std::vector<std::vector<std::uint64_t>> levels(static_cast<std::size_t>(depth) + 1);
  for (std::size_t ix = 0; ix < grid_size_; ++ix) {
    for (std::size_t iy = 0; iy < grid_size_; ++iy) {
      for (std::size_t iz = 0; iz < grid_size_; ++iz) {
        if (active[(ix * grid_size_ + iy) * grid_size_ + iz] != 0) {
          levels.back().push_back(morton(ix, iy, iz));
        }
      }
    }
  }
  std::sort(levels.back().begin(), levels.back().end());
  for (int level = depth - 1; level >= 0; --level) {
    auto& keys = levels[static_cast<std::size_t>(level)];
    for (const std::uint64_t key : levels[static_cast<std::size_t>(level) + 1]) {
      if (keys.empty() || keys.back() != key >> 3) {
        keys.push_back(key >> 3);
      }
    }
  }
  std::vector<std::size_t> level_start(levels.size() + 1, 0);
  for (std::size_t level = 0; level < levels.size(); ++level) {
    level_start[level + 1] = level_start[level] + levels[level].size();
  }
  targets_.resize(level_start.back());
  for (std::size_t level = 0; level < levels.size(); ++level) {
    const double edge = 2.0 * half / static_cast<double>(std::size_t{1} << level);
    std::size_t child = 0;
    for (std::size_t i = 0; i < levels[level].size(); ++i) {
      const std::uint64_t key = levels[level][i];
      Cell& cell = targets_[level_start[level] + i];
      cell.level = static_cast<int>(level);
      cell.centre = {grid_origin_.x + (static_cast<double>(compact_bits(key >> 2)) + 0.5) * edge,
                     grid_origin_.y + (static_cast<double>(compact_bits(key >> 1)) + 0.5) * edge,
                     grid_origin_.z + (static_cast<double>(compact_bits(key)) + 0.5) * edge};
      cell.radius = 0.5 * std::sqrt(3.0) * edge;
      if (level + 1 < levels.size()) {
        const auto& next = levels[level + 1];
        cell.first_child = level_start[level + 1] + child;
        while (child < next.size() && next[child] >> 3 == key) {
          targets_[level_start[level + 1] + child].parent = level_start[level] + i;
          ++child;
          ++cell.children;
        }
      }
    }
  }
  grid_leaf_.assign(active.size(), -1);
  for (std::size_t i = 0; i < levels.back().size(); ++i) {
    const std::uint64_t key = levels.back()[i];
    const std::size_t linear =
        (compact_bits(key >> 2) * grid_size_ + compact_bits(key >> 1)) * grid_size_ +
        compact_bits(key);
    grid_leaf_[linear] = static_cast<std::int32_t>(i);
    leaf_nodes_.push_back(level_start.back() - levels.back().size() + i);
  }
  statistics_.target_cells = targets_.size();
  statistics_.target_leaves = leaf_nodes_.size();
}

void FarFieldExpansion::traverse(const FarFieldOptions& options,
                                 const std::array<double, kind_count>& sums) {
  // Tolerances per unit source of each kind; the budget is split equally among the kinds
  // present.
  const auto present =
      static_cast<double>(std::ranges::count_if(sums, [](double v) { return v > 0.0; }));
  const double share = present > 0.0 ? 1.0 / present : 1.0;
  std::array<std::vector<double>, kind_count> tolerances;
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    tolerances[kind].resize(options.tolerances.size());
    for (std::size_t i = 0; i < options.tolerances.size(); ++i) {
      tolerances[kind][i] = sums[kind] > 0.0 ? share * options.tolerances[i] / sums[kind]
                                             : std::numeric_limits<double>::max();
    }
  }
  const double penetration_radius = options.penetration_radius;
  const int max_order = options.max_order;
  order_ = options.rank;
  std::vector<std::size_t> leaf_of(targets_.size(), 0);
  for (std::size_t leaf = 0; leaf < leaf_nodes_.size(); ++leaf) {
    leaf_of[leaf_nodes_[leaf]] = leaf;
  }
  if (targets_.empty() || sources_.empty()) {
    to_csr({}, targets_.size(), m2l_offsets_, m2l_sources_);
    to_csr({}, leaf_nodes_.size(), p2l_offsets_, p2l_charges_);
    for (std::size_t kind = 0; kind < kind_count; ++kind) {
      to_csr({}, leaf_nodes_.size(), near_offsets_[kind], near_[kind]);
    }
    return;
  }
  // Required order of sites of `kind` within `source_radius` of a centre at distance d from a
  // target cell of radius `target_radius`, or -1.
  const auto kind_order = [&](std::size_t kind, double source_radius, double target_radius,
                              double d) {
    return kind == charge_kind
               ? required_order(target_radius + source_radius, d, penetration_radius,
                                tolerances[kind], max_order)
               : derivative_required_order(kind, source_radius, target_radius, d,
                                           penetration_radius, tolerances[kind], max_order);
  };
  // Required order of an M2L from a source cell (every kind it holds), or -1.
  const auto cell_order = [&](const Cell& target, const Cell& source) {
    const double d = distance(target.centre, source.centre);
    int order = options.rank;
    for (std::size_t kind = 0; kind < kind_count; ++kind) {
      if (source.kinds[kind]) {
        const int kind_required = kind_order(kind, source.radius, target.radius, d);
        if (kind_required < 0) {
          return -1;
        }
        order = std::max(order, kind_required);
      }
    }
    return order;
  };
  // Single sites are tested against the leaf radius by distance alone: admissible from
  // distances.back(), at order rank + (first index whose distance they reach).
  const double leaf_radius = targets_[leaf_nodes_.front()].radius;
  std::array<std::vector<double>, kind_count> distances;
  std::array<double, kind_count> nearest{};
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    if (kind == charge_kind || sums[kind] > 0.0) {
      distances[kind] =
          admissible_distances(leaf_radius, penetration_radius, options.rank, max_order,
                               [&](double d) { return kind_order(kind, 0.0, leaf_radius, d); });
      nearest[kind] = distances[kind].back();
    }
  }

  struct Lists {
    std::vector<std::pair<std::size_t, std::uint32_t>> m2l, p2l;
    std::array<std::vector<std::pair<std::size_t, std::uint32_t>>, kind_count> near;
    int order = 0;
    std::array<double, kind_count> closest{std::numeric_limits<double>::infinity(),
                                           std::numeric_limits<double>::infinity(),
                                           std::numeric_limits<double>::infinity()};
  };
  // Handles the pair (t, s): M2L if admissible, sites of two leaves one by one, otherwise
  // `split` gets the pairs of the larger cell's children.
  const auto visit = [&](std::size_t t, std::size_t s, Lists& lists, auto&& split) {
    const Cell& target = targets_[t];
    const Cell& source = sources_[s];
    const int order = cell_order(target, source);
    if (order >= 0) {
      lists.m2l.emplace_back(t, static_cast<std::uint32_t>(s));
      lists.order = std::max(lists.order, order);
      return;
    }
    const bool target_leaf = target.children == 0;
    const bool source_leaf = source.children == 0;
    if (target_leaf && source_leaf) {
      const std::size_t leaf = leaf_of[t];
      for (std::size_t c = source.first; c < source.first + source.count; ++c) {
        const double d = distance(target.centre, sorted_coordinates_[c]);
        const std::size_t kind = mixed_ ? static_cast<std::size_t>(sorted_kind_[c]) : charge_kind;
        if (d >= nearest[kind]) {
          lists.p2l.emplace_back(leaf, static_cast<std::uint32_t>(c));
          lists.closest[kind] = std::min(lists.closest[kind], d);
        } else {
          lists.near[kind].emplace_back(leaf, sorted_index_[c]);
        }
      }
      return;
    }
    if (!source_leaf && (target_leaf || source.radius > target.radius)) {
      for (std::size_t child = source.first_child; child < source.first_child + source.children;
           ++child) {
        split(t, child);
      }
    } else {
      for (std::size_t child = target.first_child; child < target.first_child + target.children;
           ++child) {
        split(child, s);
      }
    }
  };

  // Breadth-first expansion (serial) into independent tasks, then depth-first traversal of the
  // tasks in parallel; the lists are merged in task order, independent of the thread count.
  constexpr std::size_t task_count = 4096;
  Lists first;
  std::vector<std::pair<std::size_t, std::size_t>> frontier{{0, 0}}, next;
  while (!frontier.empty() && frontier.size() < task_count) {
    next.clear();
    bool split_any = false;
    for (const auto& [t, s] : frontier) {
      const bool leaves = targets_[t].children == 0 && sources_[s].children == 0;
      if (leaves) {
        next.emplace_back(t, s);  // leaf pairs become tasks
        continue;
      }
      visit(t, s, first, [&](std::size_t ct, std::size_t cs) { next.emplace_back(ct, cs); });
      split_any = true;
    }
    std::swap(frontier, next);
    if (!split_any) {
      break;
    }
  }
  std::vector<Lists> tasks(frontier.size());
#pragma omp parallel for schedule(dynamic, 8)
  for (std::size_t i = 0; i < frontier.size(); ++i) {
    std::vector<std::pair<std::size_t, std::size_t>> stack{frontier[i]};
    while (!stack.empty()) {
      const auto [t, s] = stack.back();
      stack.pop_back();
      const std::size_t mark = stack.size();
      visit(t, s, tasks[i], [&](std::size_t ct, std::size_t cs) { stack.emplace_back(ct, cs); });
      std::reverse(stack.begin() + static_cast<std::ptrdiff_t>(mark), stack.end());
    }
  }
  std::vector<std::pair<std::size_t, std::uint32_t>> m2l = std::move(first.m2l), p2l;
  std::array<std::vector<std::pair<std::size_t, std::uint32_t>>, kind_count> near;
  order_ = std::max(order_, first.order);
  std::array<double, kind_count> closest = first.closest;
  for (Lists& task : tasks) {
    m2l.insert(m2l.end(), task.m2l.begin(), task.m2l.end());
    p2l.insert(p2l.end(), task.p2l.begin(), task.p2l.end());
    for (std::size_t kind = 0; kind < kind_count; ++kind) {
      near[kind].insert(near[kind].end(), task.near[kind].begin(), task.near[kind].end());
      closest[kind] = std::min(closest[kind], task.closest[kind]);
    }
    order_ = std::max(order_, task.order);
  }
  // The closest P2L site of each kind sets the order that kind needs.
  const auto order_for = [&](double closest, const std::vector<double>& distances) {
    int p = options.rank;
    while (closest < distances[static_cast<std::size_t>(p - options.rank)]) {
      ++p;
    }
    return p;
  };
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    if (closest[kind] < std::numeric_limits<double>::infinity()) {
      order_ = std::max(order_, order_for(closest[kind], distances[kind]));
    }
  }
  to_csr(m2l, targets_.size(), m2l_offsets_, m2l_sources_);
  to_csr(p2l, leaf_nodes_.size(), p2l_offsets_, p2l_charges_);
  for (std::size_t kind = 0; kind < kind_count; ++kind) {
    to_csr(near[kind], leaf_nodes_.size(), near_offsets_[kind], near_[kind]);
  }
  statistics_.multipole_to_local = m2l.size();
  statistics_.charge_to_local = p2l.size();
  statistics_.near_charges = near[charge_kind].size();
  statistics_.near_dipoles = near[dipole_kind].size();
  statistics_.near_quadrupoles = near[quadrupole_kind].size();
  for (std::size_t leaf = 0; leaf < leaf_nodes_.size(); ++leaf) {
    std::size_t count = 0;
    for (std::size_t kind = 0; kind < kind_count; ++kind) {
      count += near_offsets_[kind][leaf + 1] - near_offsets_[kind][leaf];
    }
    statistics_.largest_near = std::max(statistics_.largest_near, count);
  }
}

void FarFieldExpansion::upward_pass() {
  const std::size_t size = expansion_size(order_);
  multipoles_.assign(sources_.size() * size, 0.0);
  if (sources_.empty()) {
    return;
  }
  int depth = 0;
  for (const auto& cell : sources_) {
    depth = std::max(depth, cell.level);
  }
  // Only multipoles that some M2L reads, directly or through an ancestor, are computed (cells
  // near the target region never are: their charges go to P2L or near lists).
  std::vector<char> needed(sources_.size(), 0);
  for (const std::uint32_t s : m2l_sources_) {
    needed[s] = 1;
  }
  for (std::size_t s = 1; s < sources_.size(); ++s) {  // parents before children
    needed[s] = needed[s] != 0 || needed[sources_[s].parent] != 0 ? 1 : 0;
  }
  std::vector<std::vector<std::size_t>> levels(static_cast<std::size_t>(depth) + 1);
  for (std::size_t s = 0; s < sources_.size(); ++s) {
    if (needed[s] != 0) {
      levels[static_cast<std::size_t>(sources_[s].level)].push_back(s);
    }
  }
  // A child's centre is its parent's +- the child half edge in each coordinate, so the M2M shifts
  // are 8 per level (octant bits x, y, z as in the Morton code).
  std::vector<MultipoleShift> shifts(8 * (static_cast<std::size_t>(depth) + 1));
  for (int level = 1; level <= depth; ++level) {
    const double h = source_half_edge_ / static_cast<double>(std::size_t{1} << level);
    for (std::size_t octant = 0; octant < 8; ++octant) {
      const Point3D<double> offset{(octant & 4) != 0 ? h : -h, (octant & 2) != 0 ? h : -h,
                                   (octant & 1) != 0 ? h : -h};
      shifts[8 * static_cast<std::size_t>(level) + octant] =
          make_multipole_shift(offset, {0.0, 0.0, 0.0}, order_);
    }
  }
  for (int level = depth; level >= 0; --level) {
    const auto& cells = levels[static_cast<std::size_t>(level)];
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
      Sites sites;
#pragma omp for schedule(dynamic, 4)
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const Cell& cell = sources_[cells[i]];
        const std::span<std::complex<double>> multipole(multipoles_.data() + cells[i] * size, size);
        if (cell.children == 0) {
          if (!mixed_) {
            add_charges_to_multipole(std::span(sorted_charges_).subspan(cell.first, cell.count),
                                     std::span(sorted_coordinates_).subspan(cell.first, cell.count),
                                     cell.centre, order_, multipole, workspace);
            continue;
          }
          gather(cell.first, cell.first + cell.count, [](std::size_t i) { return i; }, sites);
          add_charges_to_multipole(sites.charges, sites.charge_coordinates, cell.centre, order_,
                                   multipole, workspace);
          add_dipoles_to_multipole(sites.dipoles, sites.dipole_coordinates, cell.centre, order_,
                                   multipole, workspace);
          add_quadrupoles_to_multipole(sites.quadrupoles, sites.quadrupole_coordinates, cell.centre,
                                       order_, multipole, workspace);
          continue;
        }
        for (std::size_t child = cell.first_child; child < cell.first_child + cell.children;
             ++child) {
          const Point3D<double>& position = sources_[child].centre;
          const std::size_t octant = (position.x > cell.centre.x ? 4 : 0) |
                                     (position.y > cell.centre.y ? 2 : 0) |
                                     (position.z > cell.centre.z ? 1 : 0);
          translate_multipole(std::span(multipoles_).subspan(child * size, size),
                              shifts[8 * static_cast<std::size_t>(sources_[child].level) + octant],
                              multipole, workspace);
        }
      }
    }
  }
}

void FarFieldExpansion::downward_pass() {
  const std::size_t size = expansion_size(order_);
  locals_.assign(targets_.size() * size, 0.0);
  if (targets_.empty()) {
    return;
  }
  std::vector<std::size_t> leaf_of(targets_.size(), 0);
  for (std::size_t leaf = 0; leaf < leaf_nodes_.size(); ++leaf) {
    leaf_of[leaf_nodes_[leaf]] = leaf;
  }
  // Cells are stored by level, root first. Per level, the M2L and P2L work is cut into chunks of
  // fixed size (independent of the thread count) that accumulate into their own buffers; each
  // cell then adds L2L from its parent and its chunks in order, so results do not depend on the
  // number of threads.
  constexpr std::size_t m2l_chunk = 32;
  constexpr std::size_t p2l_chunk = 1024;
  struct Chunk {
    std::size_t cell;
    bool charges;  // P2L (else M2L)
    std::size_t first;
    std::size_t last;
  };
  std::vector<Chunk> chunks;
  std::vector<std::size_t> chunk_offsets;  // per cell of the level
  std::vector<std::complex<double>> buffers;
  std::size_t begin = 0;
  while (begin < targets_.size()) {
    std::size_t end = begin;
    while (end < targets_.size() && targets_[end].level == targets_[begin].level) {
      ++end;
    }
    chunks.clear();
    chunk_offsets.assign(1, 0);
    for (std::size_t t = begin; t < end; ++t) {
      for (std::size_t i = m2l_offsets_[t]; i < m2l_offsets_[t + 1]; i += m2l_chunk) {
        chunks.push_back({t, false, i, std::min(i + m2l_chunk, m2l_offsets_[t + 1])});
      }
      if (targets_[t].children == 0) {
        const std::size_t leaf = leaf_of[t];
        for (std::size_t i = p2l_offsets_[leaf]; i < p2l_offsets_[leaf + 1]; i += p2l_chunk) {
          chunks.push_back({t, true, i, std::min(i + p2l_chunk, p2l_offsets_[leaf + 1])});
        }
      }
      chunk_offsets.push_back(chunks.size());
    }
    buffers.assign(chunks.size() * size, 0.0);
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
      Sites sites;
#pragma omp for schedule(dynamic, 1)
      for (std::size_t k = 0; k < chunks.size(); ++k) {
        const Chunk& chunk = chunks[k];
        const Cell& cell = targets_[chunk.cell];
        const std::span<std::complex<double>> buffer(buffers.data() + k * size, size);
        if (chunk.charges) {
          gather(chunk.first, chunk.last, [&](std::size_t i) { return p2l_charges_[i]; }, sites);
          add_charges_to_local(sites.charges, sites.charge_coordinates, cell.centre, order_, buffer,
                               workspace);
          if (mixed_) {
            add_dipoles_to_local(sites.dipoles, sites.dipole_coordinates, cell.centre, order_,
                                 buffer, workspace);
            add_quadrupoles_to_local(sites.quadrupoles, sites.quadrupole_coordinates, cell.centre,
                                     order_, buffer, workspace);
          }
          continue;
        }
        for (std::size_t i = chunk.first; i < chunk.last; ++i) {
          const std::size_t s = m2l_sources_[i];
          multipole_to_local(std::span(multipoles_).subspan(s * size, size), sources_[s].centre,
                             cell.centre, order_, buffer, workspace);
        }
      }
#pragma omp for schedule(dynamic, 1)
      for (std::size_t t = begin; t < end; ++t) {
        const Cell& cell = targets_[t];
        const std::span<std::complex<double>> local(locals_.data() + t * size, size);
        if (cell.level > 0) {
          translate_local(std::span(locals_).subspan(cell.parent * size, size),
                          targets_[cell.parent].centre, cell.centre, order_, local, workspace);
        }
        for (std::size_t k = chunk_offsets[t - begin]; k < chunk_offsets[t - begin + 1]; ++k) {
          for (std::size_t i = 0; i < size; ++i) {
            local[i] += buffers[k * size + i];
          }
        }
      }
    }
    begin = end;
  }
}

auto FarFieldExpansion::contains(const Point3D<double>& point) const -> bool {
  if (grid_size_ == 0) {
    return false;
  }
  const double coordinates[3] = {(point.x - grid_origin_.x) / grid_edge_,
                                 (point.y - grid_origin_.y) / grid_edge_,
                                 (point.z - grid_origin_.z) / grid_edge_};
  const auto limit = static_cast<double>(grid_size_);
  for (const double value : coordinates) {
    if (!(value >= 0.0 && value < limit)) {
      return false;
    }
  }
  const auto index = [&](double value) { return static_cast<std::size_t>(value); };
  return grid_leaf_[(index(coordinates[0]) * grid_size_ + index(coordinates[1])) * grid_size_ +
                    index(coordinates[2])] >= 0;
}

auto FarFieldExpansion::leaf(const Point3D<double>& point) const -> std::size_t {
  assert(grid_size_ > 0);
  const auto index = [&](double value, double origin) {
    return static_cast<std::size_t>(std::clamp(std::floor((value - origin) / grid_edge_), 0.0,
                                               static_cast<double>(grid_size_ - 1)));
  };
  const std::int32_t leaf =
      grid_leaf_[(index(point.x, grid_origin_.x) * grid_size_ + index(point.y, grid_origin_.y)) *
                     grid_size_ +
                 index(point.z, grid_origin_.z)];
  assert(leaf >= 0);
  return static_cast<std::size_t>(leaf);
}

auto FarFieldExpansion::centre(std::size_t leaf) const -> const Point3D<double>& {
  return targets_[leaf_nodes_[leaf]].centre;
}

auto FarFieldExpansion::local(std::size_t leaf) const -> std::span<const std::complex<double>> {
  const std::size_t size = expansion_size(order_);
  return std::span(locals_).subspan(leaf_nodes_[leaf] * size, size);
}

}  // namespace fika::detail
