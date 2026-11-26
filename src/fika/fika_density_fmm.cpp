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

#include "fika_density_fmm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "fika_multipole_expansion.hpp"

namespace fika::detail {

namespace {

using Complex = std::complex<double>;

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::add_density_field_fmm: " + reason);
}

auto distance(const Point3D<double>& a, const Point3D<double>& b) -> double {
  return std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}

/// Adaptive octree over points: cells in depth-first preorder (parents before children), each
/// holding a contiguous range of `order` (point indices).
struct Octree {
  struct Cell {
    Point3D<double> centre;
    double half = 0.0;    // half edge of the cube
    double radius = 0.0;  // largest distance of a point from the centre
    std::size_t first = 0;
    std::size_t count = 0;
    std::int32_t parent = -1;
    int level = 0;
    std::vector<std::uint32_t> children;
    auto leaf() const noexcept -> bool { return children.empty(); }
  };
  std::vector<Cell> cells;
  std::vector<std::uint32_t> order;
  std::vector<std::vector<std::uint32_t>> levels;  // cells per level
};

auto build_octree(std::span<const Point3D<double>> points, std::size_t leaf_size) -> Octree {
  constexpr int max_level = 30;
  Octree tree;
  if (points.empty()) {
    return tree;
  }
  tree.order.resize(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    tree.order[i] = static_cast<std::uint32_t>(i);
  }
  Point3D<double> low = points[0], high = points[0];
  for (const auto& p : points) {
    low = {std::min(low.x, p.x), std::min(low.y, p.y), std::min(low.z, p.z)};
    high = {std::max(high.x, p.x), std::max(high.y, p.y), std::max(high.z, p.z)};
  }
  Octree::Cell root;
  root.centre = {0.5 * (low.x + high.x), 0.5 * (low.y + high.y), 0.5 * (low.z + high.z)};
  root.half = 0.5 * std::max({high.x - low.x, high.y - low.y, high.z - low.z, 1e-12}) * (1 + 1e-9);
  root.count = points.size();
  tree.cells.push_back(root);
  std::vector<std::uint32_t> stack{0};
  std::vector<std::uint32_t> scratch;
  while (!stack.empty()) {
    const std::uint32_t index = stack.back();
    stack.pop_back();
    Octree::Cell cell = tree.cells[index];
    const auto range = std::span(tree.order).subspan(cell.first, cell.count);
    bool coincident = true;  // all points of the cell at one position (one-centre sources)
    for (const std::uint32_t i : range) {
      cell.radius = std::max(cell.radius, distance(points[i], cell.centre));
      coincident = coincident && points[i] == points[range.front()];
    }
    // Coincident points stay together: splitting them would only add single-child levels.
    if (cell.count > leaf_size && cell.level < max_level && !coincident) {
      // Stable partition by octant.
      std::array<std::vector<std::uint32_t>, 8> octants;
      for (const std::uint32_t i : range) {
        const auto& p = points[i];
        const int octant = (p.x >= cell.centre.x ? 1 : 0) + (p.y >= cell.centre.y ? 2 : 0) +
                           (p.z >= cell.centre.z ? 4 : 0);
        octants[static_cast<std::size_t>(octant)].push_back(i);
      }
      std::size_t first = cell.first;
      for (int octant = 0; octant < 8; ++octant) {
        const auto& members = octants[static_cast<std::size_t>(octant)];
        if (members.empty()) {
          continue;
        }
        std::ranges::copy(members, tree.order.begin() + static_cast<std::ptrdiff_t>(first));
        Octree::Cell child;
        const double h = 0.5 * cell.half;
        child.centre = {cell.centre.x + ((octant & 1) != 0 ? h : -h),
                        cell.centre.y + ((octant & 2) != 0 ? h : -h),
                        cell.centre.z + ((octant & 4) != 0 ? h : -h)};
        child.half = h;
        child.first = first;
        child.count = members.size();
        child.parent = static_cast<std::int32_t>(index);
        child.level = cell.level + 1;
        first += members.size();
        cell.children.push_back(static_cast<std::uint32_t>(tree.cells.size()));
        tree.cells.push_back(child);
      }
      // Children are visited in octant order (pushed in reverse).
      for (auto it = cell.children.rbegin(); it != cell.children.rend(); ++it) {
        stack.push_back(*it);
      }
    }
    tree.cells[index] = std::move(cell);
  }
  // Renumber in depth-first preorder so that parents precede children.
  std::vector<std::uint32_t> preorder;
  preorder.reserve(tree.cells.size());
  stack.assign(1, 0);
  while (!stack.empty()) {
    const std::uint32_t index = stack.back();
    stack.pop_back();
    preorder.push_back(index);
    const auto& children = tree.cells[index].children;
    for (auto it = children.rbegin(); it != children.rend(); ++it) {
      stack.push_back(*it);
    }
  }
  std::vector<std::uint32_t> rank(tree.cells.size());
  for (std::size_t i = 0; i < preorder.size(); ++i) {
    rank[preorder[i]] = static_cast<std::uint32_t>(i);
  }
  std::vector<Octree::Cell> cells(tree.cells.size());
  for (std::size_t i = 0; i < preorder.size(); ++i) {
    Octree::Cell cell = std::move(tree.cells[preorder[i]]);
    for (auto& child : cell.children) {
      child = rank[child];
    }
    if (cell.parent >= 0) {
      cell.parent = static_cast<std::int32_t>(rank[static_cast<std::size_t>(cell.parent)]);
    }
    cells[i] = std::move(cell);
  }
  tree.cells = std::move(cells);
  for (std::size_t i = 0; i < tree.cells.size(); ++i) {
    const auto level = static_cast<std::size_t>(tree.cells[i].level);
    if (tree.levels.size() <= level) {
      tree.levels.resize(level + 1);
    }
    tree.levels[level].push_back(static_cast<std::uint32_t>(i));
  }
  return tree;
}

/// Interaction lists of the dual traversal.
struct Interactions {
  std::vector<std::vector<std::uint32_t>> expansions;  // per target cell: source cells (M2L)
  std::vector<std::vector<std::uint32_t>> direct;      // per target cell (leaves): source leaves
};

auto traverse(const Octree& targets, const Octree& sources, std::span<const double> reach,
              double theta) -> Interactions {
  Interactions lists;
  lists.expansions.resize(targets.cells.size());
  lists.direct.resize(targets.cells.size());
  std::vector<std::pair<std::uint32_t, std::uint32_t>> stack{{0, 0}};
  while (!stack.empty()) {
    const auto [t, s] = stack.back();
    stack.pop_back();
    const auto& target = targets.cells[t];
    const auto& source = sources.cells[s];
    const double d = distance(target.centre, source.centre);
    if (target.radius + source.radius <= theta * d && d - target.radius >= reach[s]) {
      lists.expansions[t].push_back(s);
    } else if (target.leaf() && source.leaf()) {
      lists.direct[t].push_back(s);
    } else if (!source.leaf() && (target.leaf() || source.radius >= target.radius)) {
      for (auto it = source.children.rbegin(); it != source.children.rend(); ++it) {
        stack.push_back({t, *it});
      }
    } else {
      for (auto it = target.children.rbegin(); it != target.children.rend(); ++it) {
        stack.push_back({*it, s});
      }
    }
  }
  return lists;
}

}  // namespace

auto add_density_field_fmm(const DensityMultipoles& sources, std::span<const Point3D<double>> sites,
                           std::span<Point3D<double>> field, const DensityFmmOptions& options)
    -> DensityFmmReport {
  if (field.size() != sites.size()) {
    fail("field has " + std::to_string(field.size()) + " entries for " +
         std::to_string(sites.size()) + " sites");
  }
  if (!(options.accuracy > 0.0) || !std::isfinite(options.accuracy)) {
    fail("the accuracy must be positive and finite");
  }
  if (!(options.theta > 0.0 && options.theta < 1.0)) {
    fail("theta must lie in (0, 1)");
  }
  if (options.source_leaf == 0 || options.target_leaf == 0 || options.max_order < 1 ||
      options.max_order > max_expansion_order || options.order < 0) {
    fail("leaf sizes must be positive and max_order within 1.." +
         std::to_string(max_expansion_order));
  }
  DensityFmmReport report;
  if (sources.size() == 0 || sites.empty()) {
    return report;
  }
  const Octree source_tree = build_octree(sources.centres, options.source_leaf);
  const Octree target_tree = build_octree(sites, options.target_leaf);
  report.source_cells = source_tree.cells.size();
  report.target_cells = target_tree.cells.size();

  // Per source cell: reach of the penetration spheres and moment sums per rank.
  int max_rank = 0;
  for (std::size_t s = 0; s < sources.size(); ++s) {
    max_rank = std::max(max_rank, sources.rank(s));
  }
  const auto ranks = static_cast<std::size_t>(max_rank) + 1;
  std::vector<double> reach(source_tree.cells.size(), 0.0);
  std::vector<double> moment_sums(source_tree.cells.size() * ranks, 0.0);
  for (std::size_t c = 0; c < source_tree.cells.size(); ++c) {
    const auto& cell = source_tree.cells[c];
    for (std::size_t i = cell.first; i < cell.first + cell.count; ++i) {
      const std::uint32_t s = source_tree.order[i];
      reach[c] = std::max(reach[c], distance(sources.centres[s], cell.centre) +
                                        std::sqrt(sources.penetration_squared[s]));
      const auto moments = sources.moments_of(s);
      for (int k = 0; k <= sources.rank(s); ++k) {
        double norm = 0.0;
        for (int m = -k; m <= k; ++m) {
          const double q = moments[static_cast<std::size_t>(k * k + m + k)];
          norm += q * q;
        }
        moment_sums[c * ranks + static_cast<std::size_t>(k)] += std::sqrt(norm);
      }
    }
  }
  const Interactions lists = traverse(target_tree, source_tree, reach, options.theta);
  for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
    report.expansion_pairs += lists.expansions[t].size();
    for (const std::uint32_t s : lists.direct[t]) {
      report.direct_pairs += target_tree.cells[t].count * source_tree.cells[s].count;
    }
  }

  // Smallest order whose bound, summed over the expansions reaching each target leaf, meets the
  // accuracy (the bound decreases with the order). The bound of every expansion pair is computed
  // in parallel, then summed down the tree in preorder (the same order as one serial pass).
  std::vector<std::size_t> expansion_offsets(target_tree.cells.size() + 1, 0);
  for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
    expansion_offsets[t + 1] = expansion_offsets[t] + lists.expansions[t].size();
  }
  std::vector<double> pair_bounds(expansion_offsets.back());
  const auto max_leaf_bound = [&](int order) {
#pragma omp parallel for schedule(dynamic, 16)
    for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
      const auto& target = target_tree.cells[t];
      std::size_t k = expansion_offsets[t];
      for (const std::uint32_t s : lists.expansions[t]) {
        const auto& source = source_tree.cells[s];
        pair_bounds[k++] = real_multipole_field_tensor_error_bound(
            order, 1, std::span(moment_sums).subspan(s * ranks, ranks), source.radius,
            target.radius, distance(target.centre, source.centre));
      }
    }
    std::vector<double> total(target_tree.cells.size(), 0.0);
    double worst = 0.0;
    for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
      const auto& target = target_tree.cells[t];
      double sum = target.parent >= 0 ? total[static_cast<std::size_t>(target.parent)] : 0.0;
      for (std::size_t k = expansion_offsets[t]; k < expansion_offsets[t + 1]; ++k) {
        sum += pair_bounds[k];
      }
      total[t] = sum;
      if (target.leaf()) {
        worst = std::isnan(sum) ? std::numeric_limits<double>::infinity() : std::max(worst, sum);
      }
    }
    return worst;
  };
  // NaN bounds give an infinite worst bound, which misses the accuracy.
  const auto meets = [&](int order) { return max_leaf_bound(order) <= options.accuracy; };
  int low = 1, high = options.max_order;
  if (options.order > 0) {
    low = high = std::min(options.order, options.max_order);
  } else if (!meets(high)) {
    throw std::runtime_error("fika::add_density_field_fmm: order " +
                             std::to_string(options.max_order) + " misses the accuracy");
  }
  while (low < high) {
    const int middle = (low + high) / 2;
    if (meets(middle)) {
      high = middle;
    } else {
      low = middle + 1;
    }
  }
  int order = low;

  // Sites gathered per target leaf, and the sample checked against the direct sum.
  std::vector<Point3D<double>> gathered(sites.size());
  for (std::size_t i = 0; i < sites.size(); ++i) {
    gathered[i] = sites[target_tree.order[i]];
  }
  const std::size_t samples = std::min(options.sample_size, sites.size());
  std::vector<Point3D<double>> sample_sites(samples);
  std::vector<std::size_t> sample_indices(samples);
  for (std::size_t k = 0; k < samples; ++k) {
    sample_indices[k] = k * sites.size() / samples;
    sample_sites[k] = sites[sample_indices[k]];
  }
  std::vector<Point3D<double>> exact(samples);
  add_density_field(sources, sample_sites, exact);

  std::vector<Point3D<double>> result(sites.size());
  while (true) {
    const std::size_t size = expansion_size(order);
    // Upward pass: P2M at the source leaves, then M2M in reverse preorder (children first).
    std::vector<Complex> multipoles(source_tree.cells.size() * size, 0.0);
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
#pragma omp for schedule(dynamic, 1)
      for (std::size_t c = 0; c < source_tree.cells.size(); ++c) {
        const auto& cell = source_tree.cells[c];
        if (!cell.leaf()) {
          continue;
        }
        const auto multipole = std::span(multipoles).subspan(c * size, size);
        for (std::size_t i = cell.first; i < cell.first + cell.count; ++i) {
          const std::uint32_t s = source_tree.order[i];
          add_real_multipole_to_multipole(sources.moments_of(s), sources.rank(s),
                                          sources.centres[s], cell.centre, order, multipole,
                                          workspace);
        }
      }
      for (std::size_t level = source_tree.levels.size(); level-- > 1;) {
        const auto& cells = source_tree.levels[level - 1];
#pragma omp for schedule(dynamic, 1)
        for (std::size_t k = 0; k < cells.size(); ++k) {
          const auto& parent = source_tree.cells[cells[k]];
          const auto target = std::span(multipoles).subspan(cells[k] * size, size);
          for (const std::uint32_t child : parent.children) {
            translate_multipole(std::span(multipoles).subspan(child * size, size),
                                source_tree.cells[child].centre, parent.centre, order, target,
                                workspace);
          }
        }
      }
    }
    // M2L into each target cell, then L2L down the tree level by level.
    std::vector<Complex> locals(target_tree.cells.size() * size, 0.0);
#pragma omp parallel
    {
      ExpansionWorkspace workspace;
#pragma omp for schedule(dynamic, 1)
      for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
        const auto local = std::span(locals).subspan(t * size, size);
        for (const std::uint32_t s : lists.expansions[t]) {
          multipole_to_local(std::span(multipoles).subspan(s * size, size),
                             source_tree.cells[s].centre, target_tree.cells[t].centre, order, local,
                             workspace);
        }
      }
      for (std::size_t level = 1; level < target_tree.levels.size(); ++level) {
        const auto& cells = target_tree.levels[level];
#pragma omp for schedule(dynamic, 1)
        for (std::size_t k = 0; k < cells.size(); ++k) {
          const auto& cell = target_tree.cells[cells[k]];
          const auto parent = static_cast<std::size_t>(cell.parent);
          translate_local(std::span(locals).subspan(parent * size, size),
                          target_tree.cells[parent].centre, cell.centre, order,
                          std::span(locals).subspan(cells[k] * size, size), workspace);
        }
      }
      // L2P and the exact leaf pairs.
      DensityFieldWorkspace direct;
      std::vector<double> sums;
      std::vector<double> phi(4);
#pragma omp for schedule(dynamic, 1)
      for (std::size_t t = 0; t < target_tree.cells.size(); ++t) {
        const auto& cell = target_tree.cells[t];
        if (!cell.leaf()) {
          continue;
        }
        const auto leaf_sites = std::span(gathered).subspan(cell.first, cell.count);
        sums.assign(3 * cell.count, 0.0);
        for (const std::uint32_t s : lists.direct[t]) {
          const auto& source = source_tree.cells[s];
          add_density_field_block(sources,
                                  std::span(source_tree.order).subspan(source.first, source.count),
                                  leaf_sites, sums, direct);
        }
        const auto local = std::span(locals).subspan(t * size, size);
        for (std::size_t i = 0; i < cell.count; ++i) {
          local_field_tensor(local, cell.centre, order, leaf_sites[i], 1, phi, workspace);
          result[target_tree.order[cell.first + i]] = {
              phi[3] + sums[3 * i + 2], phi[1] + sums[3 * i], phi[2] + sums[3 * i + 1]};
        }
      }
    }
    // A NaN deviation counts as an infinite error.
    double error = 0.0;
    for (std::size_t k = 0; k < samples; ++k) {
      const double deviation = distance(result[sample_indices[k]], exact[k]);
      error = std::isnan(deviation) ? std::numeric_limits<double>::infinity()
                                    : std::max(error, deviation);
    }
    report.sampled_error = error;
    if (error <= options.accuracy) {
      break;
    }
    if (order >= options.max_order) {
      throw std::runtime_error("fika::add_density_field_fmm: sampled error " +
                               std::to_string(error) + " above the accuracy at order " +
                               std::to_string(order));
    }
    order = std::min(order + 2, options.max_order);
    ++report.retries;
  }
  report.order = order;
  report.bound = max_leaf_bound(order);
  for (std::size_t i = 0; i < sites.size(); ++i) {
    field[i].x += result[i].x;
    field[i].y += result[i].y;
    field[i].z += result[i].z;
  }
  return report;
}

}  // namespace fika::detail
