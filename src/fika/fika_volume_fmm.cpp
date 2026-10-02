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

#include "fika_volume_fmm.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "fika_multipole_expansion.hpp"
#include "fika_octree_keys.hpp"

namespace fika::detail {

namespace {

using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;

/// Target cells per M2L chunk (fixed, so the evaluation order is independent of the threads).
constexpr std::size_t chunk_size = 256;

constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

// Calibrated error model (see volume_fmm.hpp).
constexpr double error_ratio = 0.55;
constexpr double charge_error_constant = 1.8e-4;
constexpr double dipole_error_constant = 7e-3;
constexpr double error_margin = 1.0;
constexpr int smallest_automatic_order = 14;  // the calibrated range starts here
constexpr std::size_t order_sample_size = 64;

auto shift_of(int level) -> int {
  return 3 * (source_depth - level);
}

auto seconds(Clock::time_point from) -> double {
  return std::chrono::duration<double>(Clock::now() - from).count();
}

/// Integer cell coordinates of a Morton key.
struct Coordinates {
  std::int64_t x, y, z;
};

auto coordinates_of(std::uint64_t key) -> Coordinates {
  return {static_cast<std::int64_t>(compact_bits(key >> 2)),
          static_cast<std::int64_t>(compact_bits(key >> 1)),
          static_cast<std::int64_t>(compact_bits(key))};
}

}  // namespace

VolumeFmm::VolumeFmm(std::span<const Point3D<double>> charge_positions,
                     std::span<const Point3D<double>> dipole_positions,
                     std::span<const Point3D<double>> targets, const VolumeFmmOptions& options)
    : separation_(options.separation), sample_size_(options.sample_size) {
  if (options.leaf_capacity == 0) {
    throw std::invalid_argument("fika::VolumeFmm: the leaf capacity must be positive");
  }
  if (options.separation < 1 || options.separation > 2) {
    throw std::invalid_argument("fika::VolumeFmm: separation must be 1 or 2");
  }
  if (options.order < 0 || options.order > max_expansion_order) {
    throw std::invalid_argument("fika::VolumeFmm: order outside 0.." +
                                std::to_string(max_expansion_order));
  }
  if (options.order == 0 && options.separation != 2) {
    // The calibrated model holds for separation 2 only (at 1 it underestimates ~1000x).
    throw std::invalid_argument("fika::VolumeFmm: the automatic order needs separation 2");
  }
  if (options.order == 0 &&
      (!(options.absolute_accuracy > 0.0) || !std::isfinite(options.absolute_accuracy) ||
       (!charge_positions.empty() && !(options.charge_scale > 0.0)) ||
       (!dipole_positions.empty() && !(options.dipole_scale > 0.0)))) {
    throw std::invalid_argument(
        "fika::VolumeFmm: the automatic order needs a positive accuracy and source scales");
  }
  const std::size_t limit = std::numeric_limits<std::uint32_t>::max();
  if (charge_positions.size() > limit || dipole_positions.size() > limit ||
      targets.size() > limit) {
    throw std::invalid_argument("fika::VolumeFmm: too many points");
  }
  build(charge_positions, dipole_positions, targets, options.leaf_capacity);

  int order = options.order;
  if (order == 0) {
    // Largest model estimate over a sample of targets (evenly spaced in Morton order).
    const std::size_t n_targets = targets_.index.size();
    const std::size_t n_sample = std::min(order_sample_size, n_targets);
    std::vector<double> estimates(n_sample, 0.0);
    // Sources within separation leaf edges are (mostly) summed directly.
    const double near = static_cast<double>(separation_) * statistics_.leaf_edge;
    const double near2 = near * near;
#pragma omp parallel for schedule(dynamic, 1)
    for (std::size_t k = 0; k < n_sample; ++k) {
      const Point3D<double>& x = targets_.positions[k * n_targets / n_sample];
      double charge_sum = 0.0, dipole_sum = 0.0;
      for (std::size_t j = 0; j < charges_.x.size(); ++j) {
        const double rx = x.x - charges_.x[j], ry = x.y - charges_.y[j], rz = x.z - charges_.z[j];
        const double r2 = rx * rx + ry * ry + rz * rz;
        charge_sum += r2 > near2 ? 1.0 / r2 : 0.0;
      }
      for (std::size_t j = 0; j < dipoles_.x.size(); ++j) {
        const double rx = x.x - dipoles_.x[j], ry = x.y - dipoles_.y[j], rz = x.z - dipoles_.z[j];
        const double r2 = rx * rx + ry * ry + rz * rz;
        dipole_sum += r2 > near2 ? 2.0 / (r2 * std::sqrt(r2)) : 0.0;
      }
      estimates[k] = charge_error_constant * options.charge_scale * charge_sum +
                     dipole_error_constant * options.dipole_scale * dipole_sum;
    }
    const double scale =
        estimates.empty() ? 0.0 : *std::max_element(estimates.begin(), estimates.end());
    order = smallest_automatic_order;
    double predicted = scale * std::pow(error_ratio, order);
    while (order < largest_automatic_fmm_order &&
           error_margin * predicted > options.absolute_accuracy) {
      ++order;
      predicted *= error_ratio;
    }
    statistics_.predicted_error = predicted;
  }
  statistics_.order = order;
  // The operators (tens to hundreds of MB at high orders) only when a translation uses them.
  if (statistics_.multipole_to_local > 0) {
    m2l_.emplace(order, separation_);
    statistics_.operator_bytes = m2l_->operator_bytes();
  }
}

void VolumeFmm::build(std::span<const Point3D<double>> charges,
                      std::span<const Point3D<double>> dipoles,
                      std::span<const Point3D<double>> targets, std::size_t capacity) {
  std::vector<Point3D<double>> all(charges.begin(), charges.end());
  all.insert(all.end(), dipoles.begin(), dipoles.end());
  all.insert(all.end(), targets.begin(), targets.end());
  if (all.empty()) {
    return;
  }
  const auto [centre, half] = bounding_cube(all);
  origin_ = {centre.x - half, centre.y - half, centre.z - half};
  edge_ = 2.0 * half;
  const double scale = static_cast<double>(1 << source_depth) / edge_;
  const auto cell = [&](double value, double low) {
    return static_cast<std::uint64_t>(std::clamp(std::floor((value - low) * scale), 0.0,
                                                 static_cast<double>((1 << source_depth) - 1)));
  };
  const auto sort_points = [&](std::span<const Point3D<double>> points, Sorted& sorted) {
    std::vector<std::pair<std::uint64_t, std::uint32_t>> keys(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      keys[i] = {morton(cell(points[i].x, origin_.x), cell(points[i].y, origin_.y),
                        cell(points[i].z, origin_.z)),
                 static_cast<std::uint32_t>(i)};
    }
    if (!keys.empty()) {
      radix_sort(keys);
    }
    for (const auto& [key, index] : keys) {
      sorted.keys.push_back(key);
      sorted.index.push_back(index);
      sorted.positions.push_back(points[index]);
      sorted.x.push_back(points[index].x);
      sorted.y.push_back(points[index].y);
      sorted.z.push_back(points[index].z);
    }
  };
  sort_points(charges, charges_);
  sort_points(dipoles, dipoles_);
  sort_points(targets, targets_);
  std::vector<std::uint64_t> merged, keys;
  std::merge(charges_.keys.begin(), charges_.keys.end(), dipoles_.keys.begin(), dipoles_.keys.end(),
             std::back_inserter(merged));
  std::merge(merged.begin(), merged.end(), targets_.keys.begin(), targets_.keys.end(),
             std::back_inserter(keys));

  // Leaf level: the shallowest (>= 2) with at most `capacity` points per non-empty leaf.
  const auto distinct = [&](int level) {
    std::size_t count = 0;
    std::uint64_t previous = std::numeric_limits<std::uint64_t>::max();
    for (const std::uint64_t key : keys) {
      const std::uint64_t prefix = key >> shift_of(level);
      count += prefix != previous ? 1 : 0;
      previous = prefix;
    }
    return count;
  };
  int depth = 2;
  while (depth < source_depth && keys.size() > capacity * distinct(depth)) {
    ++depth;
  }

  // Non-empty cells of every level.
  const auto range = [&](const Sorted& sorted, std::uint64_t prefix, int level, std::size_t* out) {
    const int shift = shift_of(level);
    out[0] = static_cast<std::size_t>(
        std::lower_bound(sorted.keys.begin(), sorted.keys.end(), prefix << shift) -
        sorted.keys.begin());
    out[1] = static_cast<std::size_t>(
        std::lower_bound(sorted.keys.begin(), sorted.keys.end(), (prefix + 1) << shift) -
        sorted.keys.begin());
  };
  levels_.assign(static_cast<std::size_t>(depth) + 1, {});
  for (int level = 0; level <= depth; ++level) {
    auto& cells = levels_[static_cast<std::size_t>(level)];
    const double cell_edge = edge_ / static_cast<double>(std::uint64_t{1} << level);
    std::uint64_t previous = std::numeric_limits<std::uint64_t>::max();
    for (const std::uint64_t key : keys) {
      const std::uint64_t prefix = key >> shift_of(level);
      if (prefix == previous) {
        continue;
      }
      previous = prefix;
      Cell entry;
      entry.key = prefix;
      const auto c = coordinates_of(prefix);
      entry.centre = {origin_.x + (static_cast<double>(c.x) + 0.5) * cell_edge,
                      origin_.y + (static_cast<double>(c.y) + 0.5) * cell_edge,
                      origin_.z + (static_cast<double>(c.z) + 0.5) * cell_edge};
      range(charges_, prefix, level, entry.charges);
      range(dipoles_, prefix, level, entry.dipoles);
      range(targets_, prefix, level, entry.targets);
      cells.push_back(entry);
    }
  }
  const auto find = [&](int level, std::uint64_t key) {
    const auto& cells = levels_[static_cast<std::size_t>(level)];
    const auto it = std::lower_bound(cells.begin(), cells.end(), key,
                                     [](const Cell& c, std::uint64_t k) { return c.key < k; });
    return it != cells.end() && it->key == key ? static_cast<std::size_t>(it - cells.begin())
                                               : none;
  };
  for (int level = 0; level <= depth; ++level) {
    for (Cell& entry : levels_[static_cast<std::size_t>(level)]) {
      if (level > 0) {
        entry.parent = find(level - 1, entry.key >> 3);
      }
      if (level < depth) {
        const auto& below = levels_[static_cast<std::size_t>(level) + 1];
        const auto lower = [](const Cell& c, std::uint64_t k) { return c.key < k; };
        entry.children[0] = static_cast<std::size_t>(
            std::lower_bound(below.begin(), below.end(), entry.key << 3, lower) - below.begin());
        entry.children[1] = static_cast<std::size_t>(
            std::lower_bound(below.begin(), below.end(), (entry.key + 1) << 3, lower) -
            below.begin());
      }
    }
  }
  // Cell at integer coordinates (none outside the grid or if empty).
  const auto cell_at = [&](int level, std::int64_t x, std::int64_t y, std::int64_t z) {
    const auto grid = static_cast<std::int64_t>(std::uint64_t{1} << level);
    if (x < 0 || y < 0 || z < 0 || x >= grid || y >= grid || z >= grid) {
      return none;
    }
    return find(level, morton(static_cast<std::uint64_t>(x), static_cast<std::uint64_t>(y),
                              static_cast<std::uint64_t>(z)));
  };

  // Offset index of v = target - source (components in -reach..reach).
  const int reach = 2 * separation_ + 1;
  const auto width = static_cast<std::size_t>(2 * reach + 1);
  std::vector<std::uint32_t> offset_index(width * width * width,
                                          std::numeric_limits<std::uint32_t>::max());
  const auto offset_key = [&](std::int64_t x, std::int64_t y, std::int64_t z) {
    return static_cast<std::size_t>(((x + reach) * static_cast<std::int64_t>(width) + y + reach) *
                                        static_cast<std::int64_t>(width) +
                                    z + reach);
  };
  const auto offsets = UniformM2L::interaction_offsets(separation_);
  for (std::size_t o = 0; o < offsets.size(); ++o) {
    offset_index[offset_key(offsets[o][0], offsets[o][1], offsets[o][2])] =
        static_cast<std::uint32_t>(o);
  }

  // M2L pairs per level, in chunks of target cells grouped by offset.
  const std::int64_t s = separation_;
  chunks_.assign(levels_.size(), {});
  for (int level = 2; level <= depth; ++level) {
    const auto& cells = levels_[static_cast<std::size_t>(level)];
    const auto& above = levels_[static_cast<std::size_t>(level) - 1];
    std::vector<std::uint32_t> target_cells;
    for (std::size_t t = 0; t < cells.size(); ++t) {
      if (cells[t].has_targets()) {
        target_cells.push_back(static_cast<std::uint32_t>(t));
      }
    }
    using Pairs = std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>;
    std::vector<Pairs> chunks((target_cells.size() + chunk_size - 1) / chunk_size,
                              Pairs(offsets.size()));
#pragma omp parallel for schedule(dynamic, 1)
    for (std::size_t k = 0; k < chunks.size(); ++k) {
      Pairs& by_offset = chunks[k];
      const std::size_t last = std::min(target_cells.size(), (k + 1) * chunk_size);
      for (std::size_t i = k * chunk_size; i < last; ++i) {
        const std::uint32_t t = target_cells[i];
        const auto c = coordinates_of(cells[t].key);
        for (std::int64_t dx = -s; dx <= s; ++dx) {
          for (std::int64_t dy = -s; dy <= s; ++dy) {
            for (std::int64_t dz = -s; dz <= s; ++dz) {
              const std::size_t n =
                  cell_at(level - 1, (c.x >> 1) + dx, (c.y >> 1) + dy, (c.z >> 1) + dz);
              if (n == none) {
                continue;
              }
              for (std::size_t child = above[n].children[0]; child < above[n].children[1];
                   ++child) {
                if (!cells[child].has_sources()) {
                  continue;
                }
                const auto d = coordinates_of(cells[child].key);
                const std::int64_t vx = c.x - d.x, vy = c.y - d.y, vz = c.z - d.z;
                if (std::max({std::abs(vx), std::abs(vy), std::abs(vz)}) <= s) {
                  continue;  // a neighbour: handled below this level or directly
                }
                by_offset[offset_index[offset_key(vx, vy, vz)]].emplace_back(
                    t, static_cast<std::uint32_t>(child));
              }
            }
          }
        }
      }
    }
    auto& level_chunks = chunks_[static_cast<std::size_t>(level)];
    level_chunks.resize(chunks.size());
    for (std::size_t k = 0; k < chunks.size(); ++k) {
      Chunk& chunk = level_chunks[k];
      chunk.group_begin.push_back(0);
      for (std::size_t o = 0; o < offsets.size(); ++o) {
        if (chunks[k][o].empty()) {
          continue;
        }
        chunk.group_offset.push_back(static_cast<std::uint32_t>(o));
        for (const auto& [t, source] : chunks[k][o]) {
          chunk.targets.push_back(t);
          chunk.sources.push_back(source);
        }
        chunk.group_begin.push_back(chunk.targets.size());
      }
      statistics_.multipole_to_local += chunk.targets.size();
    }
  }

  // P2P source leaves of every target leaf.
  const auto& leaves = levels_[static_cast<std::size_t>(depth)];
  near_begin_.assign(leaves.size() + 1, 0);
  std::vector<std::vector<std::uint32_t>> near(leaves.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t t = 0; t < leaves.size(); ++t) {
    if (!leaves[t].has_targets()) {
      continue;
    }
    const auto c = coordinates_of(leaves[t].key);
    for (std::int64_t dx = -s; dx <= s; ++dx) {
      for (std::int64_t dy = -s; dy <= s; ++dy) {
        for (std::int64_t dz = -s; dz <= s; ++dz) {
          const std::size_t n = cell_at(depth, c.x + dx, c.y + dy, c.z + dz);
          if (n != none && leaves[n].has_sources()) {
            near[t].push_back(static_cast<std::uint32_t>(n));
          }
        }
      }
    }
  }
  for (std::size_t t = 0; t < leaves.size(); ++t) {
    near_begin_[t + 1] = near_begin_[t] + near[t].size();
    near_leaves_.insert(near_leaves_.end(), near[t].begin(), near[t].end());
    statistics_.target_leaves += leaves[t].has_targets() ? 1 : 0;
    const std::size_t targets_here = leaves[t].targets[1] - leaves[t].targets[0];
    for (const std::uint32_t n : near[t]) {
      statistics_.near_interactions += targets_here * (leaves[n].charges[1] - leaves[n].charges[0] +
                                                       leaves[n].dipoles[1] - leaves[n].dipoles[0]);
    }
  }
  statistics_.near_leaf_pairs = near_leaves_.size();
  statistics_.depth = depth;
  statistics_.leaf_edge = edge_ / static_cast<double>(std::uint64_t{1} << depth);
  for (const auto& cells : levels_) {
    statistics_.cells += cells.size();
  }
}

void VolumeFmm::field(std::span<const double> charges, std::span<const Dipole> dipoles,
                      std::span<Point3D<double>> out, VolumeFmmReport* report) const {
  if (charges.size() != charges_.index.size() || dipoles.size() != dipoles_.index.size() ||
      out.size() != targets_.index.size()) {
    throw std::invalid_argument(
        "fika::VolumeFmm::field: " + std::to_string(charges.size()) + " charges, " +
        std::to_string(dipoles.size()) + " dipoles and " + std::to_string(out.size()) +
        " outputs for a tree of " + std::to_string(charges_.index.size()) + ", " +
        std::to_string(dipoles_.index.size()) + " and " + std::to_string(targets_.index.size()));
  }
  if (levels_.empty()) {
    return;
  }
  VolumeFmmReport local_report;
  VolumeFmmReport& time = report != nullptr ? *report : local_report;
  auto start = Clock::now();
  const int depth = static_cast<int>(levels_.size()) - 1;
  const int order = statistics_.order;
  const std::size_t size = expansion_size(order);
  std::vector<double> q(charges_.index.size());
  for (std::size_t i = 0; i < q.size(); ++i) {
    q[i] = charges[charges_.index[i]];
  }
  std::vector<Dipole> mu(dipoles_.index.size());
  for (std::size_t i = 0; i < mu.size(); ++i) {
    mu[i] = dipoles[dipoles_.index[i]];
  }
  const auto cell_edge = [&](int level) {
    return edge_ / static_cast<double>(std::uint64_t{1} << level);
  };

  // Upward: P2M in the leaves, M2M to level 2.
  std::vector<std::vector<Complex>> multipoles(levels_.size());
  for (int level = 2; level <= depth; ++level) {
    multipoles[static_cast<std::size_t>(level)].assign(
        levels_[static_cast<std::size_t>(level)].size() * size, Complex{});
  }
  {
    const auto& leaves = levels_[static_cast<std::size_t>(depth)];
    auto& target = multipoles[static_cast<std::size_t>(depth)];
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
#pragma omp for schedule(dynamic, 16)
      for (std::size_t c = 0; c < leaves.size(); ++c) {
        const Cell& leaf = leaves[c];
        const std::span<Complex> multipole(target.data() + c * size, size);
        const std::size_t nq = leaf.charges[1] - leaf.charges[0];
        if (nq > 0) {
          add_charges_to_multipole(std::span(q).subspan(leaf.charges[0], nq),
                                   std::span(charges_.positions).subspan(leaf.charges[0], nq),
                                   leaf.centre, order, multipole, workspace);
        }
        const std::size_t nd = leaf.dipoles[1] - leaf.dipoles[0];
        if (nd > 0) {
          add_dipoles_to_multipole(std::span(mu).subspan(leaf.dipoles[0], nd),
                                   std::span(dipoles_.positions).subspan(leaf.dipoles[0], nd),
                                   leaf.centre, order, multipole, workspace);
        }
      }
    }
  }
  for (int level = depth - 1; level >= 2; --level) {
    const auto& cells = levels_[static_cast<std::size_t>(level)];
    const auto& below = levels_[static_cast<std::size_t>(level) + 1];
    const double h = 0.5 * cell_edge(level + 1);
    std::vector<MultipoleShift> shifts(8);
    for (std::size_t octant = 0; octant < 8; ++octant) {
      shifts[octant] = make_multipole_shift(
          {(octant & 4) != 0 ? h : -h, (octant & 2) != 0 ? h : -h, (octant & 1) != 0 ? h : -h},
          {0.0, 0.0, 0.0}, order);
    }
    auto& target = multipoles[static_cast<std::size_t>(level)];
    const auto& source = multipoles[static_cast<std::size_t>(level) + 1];
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
#pragma omp for schedule(dynamic, 16)
      for (std::size_t c = 0; c < cells.size(); ++c) {
        if (!cells[c].has_sources()) {
          continue;
        }
        for (std::size_t child = cells[c].children[0]; child < cells[c].children[1]; ++child) {
          if (!below[child].has_sources()) {
            continue;
          }
          const std::uint64_t key = below[child].key;
          const std::size_t octant = static_cast<std::size_t>(key & 7);
          translate_multipole(std::span<const Complex>(source.data() + child * size, size),
                              shifts[octant], std::span<Complex>(target.data() + c * size, size),
                              workspace);
        }
      }
    }
  }
  time.upward = seconds(start);

  // M2L: chunks of target cells in parallel, one batched product per offset.
  start = Clock::now();
  std::vector<std::vector<Complex>> locals(levels_.size());
  for (int level = 2; level <= depth; ++level) {
    locals[static_cast<std::size_t>(level)].assign(
        levels_[static_cast<std::size_t>(level)].size() * size, Complex{});
  }
  const auto offsets = m2l_ ? m2l_->offsets() : std::span<const CellOffset>{};
  for (int level = 2; m2l_ && level <= depth; ++level) {
    const auto& chunks = chunks_[static_cast<std::size_t>(level)];
    const auto& source = multipoles[static_cast<std::size_t>(level)];
    auto& target = locals[static_cast<std::size_t>(level)];
    const double edge = cell_edge(level);
#pragma omp parallel
    {
      UniformM2LWorkspace workspace;
      std::vector<Complex> gathered, translated;
#pragma omp for schedule(dynamic, 1)
      for (std::size_t k = 0; k < chunks.size(); ++k) {
        const Chunk& chunk = chunks[k];
        for (std::size_t g = 0; g < chunk.group_offset.size(); ++g) {
          const std::size_t begin = chunk.group_begin[g], end = chunk.group_begin[g + 1];
          const std::size_t n = end - begin;
          gathered.resize(n * size);
          translated.assign(n * size, Complex{});
          for (std::size_t i = 0; i < n; ++i) {
            std::copy_n(source.data() + chunk.sources[begin + i] * size, size,
                        gathered.data() + i * size);
          }
          m2l_->apply(offsets[chunk.group_offset[g]], edge, gathered, translated, n, workspace);
          for (std::size_t i = 0; i < n; ++i) {
            Complex* local = target.data() + chunk.targets[begin + i] * size;
            const Complex* add = translated.data() + i * size;
            for (std::size_t j = 0; j < size; ++j) {
              local[j] += add[j];
            }
          }
        }
      }
    }
  }
  time.multipole_to_local = seconds(start);

  // Downward: L2L from level 2 to the leaves.
  start = Clock::now();
  for (int level = 3; level <= depth; ++level) {
    const auto& cells = levels_[static_cast<std::size_t>(level)];
    const auto& above = levels_[static_cast<std::size_t>(level) - 1];
    const auto& parent_locals = locals[static_cast<std::size_t>(level) - 1];
    auto& target = locals[static_cast<std::size_t>(level)];
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
#pragma omp for schedule(dynamic, 16)
      for (std::size_t c = 0; c < cells.size(); ++c) {
        if (!cells[c].has_targets()) {
          continue;
        }
        const std::size_t parent = cells[c].parent;
        translate_local(std::span<const Complex>(parent_locals.data() + parent * size, size),
                        above[parent].centre, cells[c].centre, order,
                        std::span<Complex>(target.data() + c * size, size), workspace);
      }
    }
  }
  time.downward = seconds(start);

  // Leaves: L2P and P2P.
  start = Clock::now();
  const auto& leaves = levels_[static_cast<std::size_t>(depth)];
  const auto& leaf_locals = locals[static_cast<std::size_t>(depth)];
  const double* cx = charges_.x.data();
  const double* cy = charges_.y.data();
  const double* cz = charges_.z.data();
  const double* dx = dipoles_.x.data();
  const double* dy = dipoles_.y.data();
  const double* dz = dipoles_.z.data();
#pragma omp parallel
  {
    ExpansionWorkspace workspace;
    std::vector<double> phi(4);
#pragma omp for schedule(dynamic, 8)
    for (std::size_t t = 0; t < leaves.size(); ++t) {
      const Cell& leaf = leaves[t];
      if (!leaf.has_targets()) {
        continue;
      }
      const std::span<const Complex> local(leaf_locals.data() + t * size, size);
      for (std::size_t i = leaf.targets[0]; i < leaf.targets[1]; ++i) {
        const Point3D<double>& x = targets_.positions[i];
        // Phi_1,M = sum q S_1M(C - x) / |C - x|^3 with S_1 = (y, z, x): minus the field.
        local_field_tensor(local, leaf.centre, order, x, 1, phi, workspace);
        double ex = -phi[3], ey = -phi[1], ez = -phi[2];
        for (std::size_t k = near_begin_[t]; k < near_begin_[t + 1]; ++k) {
          const Cell& source = leaves[near_leaves_[k]];
          double fx = 0.0, fy = 0.0, fz = 0.0;
          for (std::size_t j = source.charges[0]; j < source.charges[1]; ++j) {
            const double rx = x.x - cx[j], ry = x.y - cy[j], rz = x.z - cz[j];
            const double r2 = rx * rx + ry * ry + rz * rz;
            const double scale = r2 > 0.0 ? q[j] / (r2 * std::sqrt(r2)) : 0.0;
            fx += scale * rx;
            fy += scale * ry;
            fz += scale * rz;
          }
          for (std::size_t j = source.dipoles[0]; j < source.dipoles[1]; ++j) {
            const double rx = x.x - dx[j], ry = x.y - dy[j], rz = x.z - dz[j];
            const double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 == 0.0) {
              continue;
            }
            const double inverse3 = 1.0 / (r2 * std::sqrt(r2));
            const double projection =
                3.0 * (rx * mu[j](0) + ry * mu[j](1) + rz * mu[j](2)) * inverse3 / r2;
            fx += projection * rx - inverse3 * mu[j](0);
            fy += projection * ry - inverse3 * mu[j](1);
            fz += projection * rz - inverse3 * mu[j](2);
          }
          ex += fx;
          ey += fy;
          ez += fz;
        }
        out[targets_.index[i]] = {ex, ey, ez};
      }
    }
  }
  time.leaves = seconds(start);

  // Exact field at a fixed sample of targets (evenly spaced in Morton order).
  if (report == nullptr || sample_size_ == 0 || targets_.index.empty()) {
    return;
  }
  start = Clock::now();
  const std::size_t n_targets = targets_.index.size();
  const std::size_t n_sample = std::min(sample_size_, n_targets);
  std::vector<double> errors(n_sample), relative(n_sample);
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t k = 0; k < n_sample; ++k) {
    const std::size_t i = k * n_targets / n_sample;
    const Point3D<double>& x = targets_.positions[i];
    double ex = 0.0, ey = 0.0, ez = 0.0, absolute = 0.0;
    for (std::size_t j = 0; j < q.size(); ++j) {
      const double rx = x.x - cx[j], ry = x.y - cy[j], rz = x.z - cz[j];
      const double r2 = rx * rx + ry * ry + rz * rz;
      if (r2 == 0.0) {
        continue;
      }
      const double scale = q[j] / (r2 * std::sqrt(r2));
      ex += scale * rx;
      ey += scale * ry;
      ez += scale * rz;
      absolute += std::abs(q[j]) / r2;
    }
    for (std::size_t j = 0; j < mu.size(); ++j) {
      const double rx = x.x - dx[j], ry = x.y - dy[j], rz = x.z - dz[j];
      const double r2 = rx * rx + ry * ry + rz * rz;
      if (r2 == 0.0) {
        continue;
      }
      const double inverse3 = 1.0 / (r2 * std::sqrt(r2));
      const double projection =
          3.0 * (rx * mu[j](0) + ry * mu[j](1) + rz * mu[j](2)) * inverse3 / r2;
      ex += projection * rx - inverse3 * mu[j](0);
      ey += projection * ry - inverse3 * mu[j](1);
      ez += projection * rz - inverse3 * mu[j](2);
      absolute += 2.0 * std::hypot(mu[j](0), mu[j](1), mu[j](2)) * inverse3;
    }
    const Point3D<double>& e = out[targets_.index[i]];
    errors[k] = std::hypot(e.x - ex, e.y - ey, e.z - ez);
    relative[k] = absolute > 0.0 ? errors[k] / absolute : 0.0;
  }
  report->sampled_error = *std::max_element(errors.begin(), errors.end());
  report->sampled_relative_error = *std::max_element(relative.begin(), relative.end());
  report->sampled_targets = n_sample;
  report->sample = seconds(start);
}

}  // namespace fika::detail
