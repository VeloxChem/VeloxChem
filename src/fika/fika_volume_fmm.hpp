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

#ifndef fika_volume_fmm_hpp
#define fika_volume_fmm_hpp

// Internal: fast multipole evaluation of the electric field of many point charges and dipoles at
// many target points spread through the same volume (the permanent field at, and the coupling
// of, polarizable MM sites).
//
// A uniform octree covers all points; its leaf level is the shallowest (>= 2) whose non-empty
// leaves hold at most leaf_capacity points on average, and only non-empty cells are stored. Two
// cells of one level interact through their expansions (M2L) once separated by more than
// `separation` cells, their parents not; target leaves sum the source leaves within
// `separation` cells directly (P2P). Multipole-to-local translations use precomputed operators
// as batched matrix products (UniformM2L); all expansions have one order, which sets the
// accuracy. The tree and lists depend only on the positions: they are built once, and the field
// can be evaluated for any source values.
//
// Accuracy (calibrated, not bounded): with separation 2 and the default leaf capacity, the
// largest field error over molecular systems (osimertinib water droplets of 5000-80000 waters,
// every or 3000 targets; atoms at least ~1.8 bohr apart) stays below rho^p (C_q A_q + C_d A_d),
// rho = 0.55, C_q = 1.8e-4, C_d = 7e-3, with A_q = max|q| sum 1 / r^2 and A_d = max|mu| sum
// 2 / r^3 over the sources farther than separation leaf edges (only those enter through
// expansions). Without an explicit order, the order is the smallest with this estimate (largest
// over a sample of targets) at most absolute_accuracy, within the calibrated range 14..26. Far
// denser point sets (e.g. random points 0.3 bohr apart) can exceed it by up to ~3x; the sampled
// error of each evaluation (which underestimates the largest one by up to ~10x) lets callers
// raise the order.

#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_uniform_m2l.hpp"

namespace fika::detail {

/// Highest order of the calibrated automatic order, and of the order retries of its users.
inline constexpr int largest_automatic_fmm_order = 26;

struct VolumeFmmOptions {
  double absolute_accuracy = 1e-9;  // field error target (a.u.) of the automatic order
  int order = 0;                    // expansion order p; 0: from absolute_accuracy
  int separation = 2;               // cells between interacting cells (1 or 2)
  std::size_t leaf_capacity = 512;  // average points (sources and targets) per non-empty leaf
  std::size_t sample_size = 64;     // targets whose exact field checks each evaluation (0: none)
  double charge_scale = 0.0;        // largest |q| expected (automatic order with charges)
  double dipole_scale = 0.0;        // largest |mu| expected (automatic order with dipoles)
};

struct VolumeFmmStatistics {
  int order = 0;                       // expansion order
  double predicted_error = 0.0;        // error estimate of the calibrated model (automatic order)
  std::size_t operator_bytes = 0;      // memory of the M2L operators
  int depth = 0;                       // leaf level
  double leaf_edge = 0.0;              // bohr
  std::size_t cells = 0;               // non-empty cells of all levels
  std::size_t target_leaves = 0;       // leaves holding targets
  std::size_t multipole_to_local = 0;  // M2L translations
  std::size_t near_leaf_pairs = 0;     // P2P pairs of leaves
  std::size_t near_interactions = 0;   // source-target point pairs summed directly
};

/// Wall times (s) of the stages of one field evaluation, and its error on the sample.
struct VolumeFmmReport {
  double upward = 0.0;              // P2M and M2M
  double multipole_to_local = 0.0;  // M2L
  double downward = 0.0;            // L2L
  double leaves = 0.0;              // L2P and P2P
  double sample = 0.0;              // exact field at the sampled targets
  // Largest |E_fmm - E_exact| over the sampled targets (a.u.), and largest ratio to the absolute
  // field A = sum_j |q_j| / r^2 + sum_k 2 |mu_k| / r^3 (0 without a sample).
  double sampled_error = 0.0;
  double sampled_relative_error = 0.0;
  std::size_t sampled_targets = 0;
};

class VolumeFmm {
 public:
  /// Tree and lists for charges at `charge_positions`, dipoles at `dipole_positions` and fields
  /// at `targets` (bohr; points of different sets may coincide). Throws std::invalid_argument
  /// for an order outside 0..max_expansion_order, a separation other than 1 or 2, a zero leaf
  /// capacity, the automatic order with separation 1, a nonpositive accuracy, a missing scale for
  /// the automatic order, or too many points.
  VolumeFmm(std::span<const Point3D<double>> charge_positions,
            std::span<const Point3D<double>> dipole_positions,
            std::span<const Point3D<double>> targets, const VolumeFmmOptions& options);

  /// Sets out[i] to the field at target i of the charges, sum_j q_j r / r^3, and the dipoles,
  /// sum_k (3 r (r . mu_k) - r^2 mu_k) / r^5 (r = x_i - source), without sources at the target
  /// point itself. Independent of the thread count. With a report, also the stage timings and the
  /// error on a fixed sample of targets (evenly spaced in Morton order), computed exactly.
  /// Throws std::invalid_argument if the sizes do not match.
  void field(std::span<const double> charges, std::span<const Dipole> dipoles,
             std::span<Point3D<double>> out, VolumeFmmReport* report = nullptr) const;

  auto order() const noexcept -> int { return statistics_.order; }
  auto statistics() const noexcept -> const VolumeFmmStatistics& { return statistics_; }

 private:
  /// Points of one set in Morton order.
  struct Sorted {
    std::vector<std::uint64_t> keys;
    std::vector<std::uint32_t> index;  // original index
    std::vector<Point3D<double>> positions;
    std::vector<double> x, y, z;
  };

  struct Cell {
    std::uint64_t key = 0;  // Morton key at its level
    Point3D<double> centre;
    std::size_t charges[2] = {0, 0};  // [begin, end) in the sorted charges
    std::size_t dipoles[2] = {0, 0};
    std::size_t targets[2] = {0, 0};
    std::size_t parent = 0;            // index in the level above
    std::size_t children[2] = {0, 0};  // [begin, end) in the level below

    auto has_sources() const noexcept -> bool {
      return charges[1] > charges[0] || dipoles[1] > dipoles[0];
    }
    auto has_targets() const noexcept -> bool { return targets[1] > targets[0]; }
  };

  /// M2L pairs of a chunk of target cells of one level, grouped by offset (each target at
  /// most once per group).
  struct Chunk {
    std::vector<std::uint32_t> group_offset;  // offset index of each group
    std::vector<std::size_t> group_begin;     // groups + 1
    std::vector<std::uint32_t> targets;       // cell indices (level), group by group
    std::vector<std::uint32_t> sources;
  };

  void build(std::span<const Point3D<double>> charges, std::span<const Point3D<double>> dipoles,
             std::span<const Point3D<double>> targets, std::size_t capacity);

  std::optional<UniformM2L> m2l_;
  int separation_ = 1;
  std::size_t sample_size_ = 0;
  Point3D<double> origin_{};  // lowest corner of the root cube
  double edge_ = 0.0;         // root edge
  Sorted charges_, dipoles_, targets_;
  std::vector<std::vector<Cell>> levels_;   // non-empty cells by level, by key
  std::vector<std::vector<Chunk>> chunks_;  // per level
  std::vector<std::size_t> near_begin_;     // per leaf: P2P source leaves (CSR)
  std::vector<std::uint32_t> near_leaves_;
  VolumeFmmStatistics statistics_;
};

}  // namespace fika::detail

#endif  // fika_volume_fmm_hpp
