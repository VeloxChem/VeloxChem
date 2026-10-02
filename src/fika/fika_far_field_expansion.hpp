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

#ifndef fika_far_field_expansion_hpp
#define fika_far_field_expansion_hpp

// Internal: fast multipole evaluation of the field tensor of many point sources (charges,
// dipoles and quadrupoles) at points of a small target region (the product centres of a QM basis
// among MM sources). A source octree over all sites carries multipole expansions; a target octree
// over the region carries local expansions. For a point x in an active target leaf, the field
// tensor
//   Phi_Lambda,M(x) = sum_C q_C S_Lambda,M(C - x) / |C - x|^(2 Lambda + 1)
//                     + sum_D mu_D . grad_D [S_Lambda,M(D - x) / |D - x|^(2 Lambda + 1)]
//                     + sum_E 1/2 Theta_E : grad_E grad_E [S_Lambda,M(E - x) / |E - x|^(2 Lambda +
//                     1)]
// is the leaf's local expansion evaluated at x (L2P) plus the direct sums over the leaf's near
// charges, dipoles and quadrupoles, with error at most tolerance[Lambda] (Euclidean norm over M)
// for every Lambda <= rank.
//
// Every translation into a target cell (M2L of a source cell, or P2L of a single site) is
// admissible only if (a) its gap d - r_S - r_T is at least the penetration radius, so every
// source closer than that to any point of the cell stays near, and (b) for the kinds of source
// it holds, the error bound per unit charge (field_tensor_error_bound) is at most
// tolerance_q[Lambda] / sum_C |q_C|, the bound per unit dipole
// (dipole_field_tensor_error_bound) at most tolerance_mu[Lambda] / sum_D |mu_D| and the bound per
// unit quadrupole (quadrupole_field_tensor_error_bound) at most
// tolerance_theta[Lambda] / sum_E ||Theta_E||_2, for every Lambda <= rank at some order
// <= max_order; the tolerance is split equally among the kinds present. Each point sees every
// source exactly once (through one ancestor's M2L, the leaf's P2L or a near list), so the errors
// add up to at most tolerance[Lambda]. All expansions use the largest order any translation needs.

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_atom_pair_group.hpp"

namespace fika::detail {

/// Point sources of a potential (any of the kinds may be empty); quadrupoles are primitive
/// Cartesian moments, of which only the traceless part acts.
struct PointSources {
  std::span<const double> charges = {};
  std::span<const Point3D<double>> charge_coordinates = {};
  std::span<const Dipole> dipoles = {};
  std::span<const Point3D<double>> dipole_coordinates = {};
  std::span<const Quadrupole> quadrupoles = {};
  std::span<const Point3D<double>> quadrupole_coordinates = {};
};

struct FarFieldOptions {
  int rank = 0;                    // highest field-tensor rank Lambda
  std::vector<double> tolerances;  // per Lambda <= rank, absolute
  // Norm sums (sum |q|, sum |mu|, sum ||Theta||_2) the tolerances refer to; by default those of
  // the sources passed to the constructor. The interaction lists hold the error below
  // tolerance[Lambda] for any source values whose norm sums are at most these; with other values
  // (update()) the bound scales with the largest ratio of norm sums.
  std::optional<std::array<double, 3>> source_norms;
  double penetration_radius = 0.0;
  double leaf_edge = 8.0;          // largest target leaf edge (bohr)
  std::size_t leaf_capacity = 64;  // charges per source leaf
  int max_order = 28;
};

struct FarFieldStatistics {
  std::size_t source_cells = 0;
  std::size_t target_cells = 0;
  std::size_t target_leaves = 0;
  std::size_t multipole_to_local = 0;
  std::size_t charge_to_local = 0;   // P2L of single sites (of every kind)
  std::size_t near_charges = 0;      // sum over leaves
  std::size_t near_dipoles = 0;      // sum over leaves
  std::size_t near_quadrupoles = 0;  // sum over leaves
  std::size_t largest_near = 0;      // sites of every kind of one leaf
};

class FarFieldExpansion {
 public:
  /// `atoms` and `pairs` define the target region: the leaves crossed by the segments between
  /// the atoms of each pair (a pair of equal atoms marks its atom's leaf). Throws
  /// std::invalid_argument for inconsistent options or charge counts.
  FarFieldExpansion(std::span<const double> charges, std::span<const Point3D<double>> coordinates,
                    std::span<const Point3D<double>> atoms, std::span<const AtomPair> pairs,
                    const FarFieldOptions& options);

  /// Sources of every kind (dipoles need max_order <= max_expansion_order - 1, quadrupoles
  /// max_order <= max_expansion_order - 2).
  FarFieldExpansion(const PointSources& sources, std::span<const Point3D<double>> atoms,
                    std::span<const AtomPair> pairs, const FarFieldOptions& options);

  /// Replaces the source values (charges, dipoles, quadrupoles) at unchanged positions and
  /// recomputes the expansions; trees, interaction lists and order are kept, so the error bound
  /// scales with the norm sums of the new values relative to those the lists were built for.
  /// Throws std::invalid_argument if the source counts or positions differ from the
  /// constructor's, or a kind whose norm sum was zero at construction gets nonzero values.
  void update(const PointSources& sources);

  /// Expansion order p of all multipole and local expansions.
  auto order() const noexcept -> int { return order_; }

  /// Leaf containing `point` (which must lie in an active leaf).
  auto leaf(const Point3D<double>& point) const -> std::size_t;

  /// Whether `point` lies in an active leaf.
  auto contains(const Point3D<double>& point) const -> bool;

  auto leaf_count() const noexcept -> std::size_t { return leaf_nodes_.size(); }

  auto centre(std::size_t leaf) const -> const Point3D<double>&;

  /// Local expansion of the leaf (expansion_size(order()) coefficients).
  auto local(std::size_t leaf) const -> std::span<const std::complex<double>>;

  /// Indices (into the charges passed in) of the leaf's near charges.
  auto near(std::size_t leaf) const -> std::span<const std::uint32_t> {
    return near_sites(charge_kind, leaf);
  }

  /// Indices (into the dipoles passed in) of the leaf's near dipoles.
  auto near_dipoles(std::size_t leaf) const -> std::span<const std::uint32_t> {
    return near_sites(dipole_kind, leaf);
  }

  /// Indices (into the quadrupoles passed in) of the leaf's near quadrupoles.
  auto near_quadrupoles(std::size_t leaf) const -> std::span<const std::uint32_t> {
    return near_sites(quadrupole_kind, leaf);
  }

  auto statistics() const noexcept -> const FarFieldStatistics& { return statistics_; }

 private:
  static constexpr std::size_t charge_kind = 0;
  static constexpr std::size_t dipole_kind = 1;
  static constexpr std::size_t quadrupole_kind = 2;
  static constexpr std::size_t kind_count = 3;

  auto near_sites(std::size_t kind, std::size_t leaf) const -> std::span<const std::uint32_t> {
    return std::span(near_[kind])
        .subspan(near_offsets_[kind][leaf],
                 near_offsets_[kind][leaf + 1] - near_offsets_[kind][leaf]);
  }

  struct Cell {
    Point3D<double> centre;
    double radius = 0.0;    // sources: largest charge distance; targets: half-diagonal
    std::size_t first = 0;  // sources: first sorted charge
    std::size_t count = 0;  // sources: charges
    std::size_t first_child = 0;
    std::size_t children = 0;  // 0 for leaves
    std::size_t parent = 0;
    int level = 0;
    std::array<bool, kind_count> kinds{};  // sources: holds sites of each kind
  };

  /// Sites gathered from sorted sites by kind (per-thread scratch).
  struct Sites {
    std::vector<double> charges;
    std::vector<Point3D<double>> charge_coordinates;
    std::vector<Dipole> dipoles;
    std::vector<Point3D<double>> dipole_coordinates;
    std::vector<Quadrupole> quadrupoles;
    std::vector<Point3D<double>> quadrupole_coordinates;
  };

  /// Gathers the sorted sites site_of(first..last - 1) into `sites`, split by kind.
  template <typename SiteOf>
  void gather(std::size_t first, std::size_t last, SiteOf&& site_of, Sites& sites) const {
    sites.charges.clear();
    sites.charge_coordinates.clear();
    sites.dipoles.clear();
    sites.dipole_coordinates.clear();
    sites.quadrupoles.clear();
    sites.quadrupole_coordinates.clear();
    for (std::size_t i = first; i < last; ++i) {
      const std::size_t site = site_of(i);
      switch (mixed_ ? static_cast<std::size_t>(sorted_kind_[site]) : charge_kind) {
        case dipole_kind:
          sites.dipoles.push_back(dipoles_[sorted_index_[site]]);
          sites.dipole_coordinates.push_back(sorted_coordinates_[site]);
          break;
        case quadrupole_kind:
          sites.quadrupoles.push_back(quadrupoles_[sorted_index_[site]]);
          sites.quadrupole_coordinates.push_back(sorted_coordinates_[site]);
          break;
        default:
          sites.charges.push_back(sorted_charges_[site]);
          sites.charge_coordinates.push_back(sorted_coordinates_[site]);
      }
    }
  }

  void build_sources(const PointSources& sources, std::size_t capacity);
  void build_targets(std::span<const Point3D<double>> atoms, std::span<const AtomPair> pairs,
                     double leaf_edge);
  void traverse(const FarFieldOptions& options, const std::array<double, kind_count>& sums);
  void set_values(const PointSources& sources);
  void upward_pass();
  void downward_pass();

  int order_ = 0;
  // Norm sums the interaction lists were built for, and the source count of each kind.
  std::array<double, kind_count> norms_{};
  std::array<std::size_t, kind_count> counts_{};
  // Sources of every kind in Morton order.
  bool mixed_ = false;                  // dipoles or quadrupoles present
  std::vector<double> sorted_charges_;  // 0 for other kinds
  std::vector<char> sorted_kind_;       // with mixed_
  std::vector<Point3D<double>> sorted_coordinates_;
  std::vector<std::uint32_t> sorted_index_;  // original index per sorted site (within its kind)
  std::vector<Dipole> dipoles_;              // copies of the dipoles passed in
  std::vector<Quadrupole> quadrupoles_;      // copies of the quadrupoles passed in
  double source_half_edge_ = 0.0;            // of the root cell
  std::vector<Cell> sources_;                // children after their parent
  std::vector<std::complex<double>> multipoles_;
  // Targets: uniform leaves of a cube grid, levels above them; cells by level (root first).
  Point3D<double> grid_origin_{};
  double grid_edge_ = 0.0;               // leaf edge
  std::size_t grid_size_ = 0;            // leaves per axis
  std::vector<std::int32_t> grid_leaf_;  // active leaf per grid cell, -1 if inactive
  std::vector<Cell> targets_;
  std::vector<std::size_t> leaf_nodes_;  // target cell of each leaf
  std::vector<std::complex<double>> locals_;
  // Interaction lists per target cell (CSR).
  std::vector<std::size_t> m2l_offsets_;
  std::vector<std::uint32_t> m2l_sources_;
  std::vector<std::size_t> p2l_offsets_;                           // per leaf
  std::vector<std::uint32_t> p2l_charges_;                         // sorted site indices
  std::array<std::vector<std::size_t>, kind_count> near_offsets_;  // per kind and leaf
  std::array<std::vector<std::uint32_t>, kind_count> near_;
  FarFieldStatistics statistics_;
};

}  // namespace fika::detail

#endif  // fika_far_field_expansion_hpp
