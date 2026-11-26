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

#ifndef fika_nuclear_attraction_kernels_hpp
#define fika_nuclear_attraction_kernels_hpp

// Internal: nuclear-attraction (point-charge potential) kernels, V_ab = sum_C q_C <a|1/|r - C||b>
// (nuclear-attraction notes). Same-atom blocks (A|V|A) use ordering A (Section 4.3): the charges
// are accumulated at the primitive level, then contracted and assembled once; charges far from
// the atom enter through its multipole field Omega_LM (Section 4.7).

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "fika_basis_shell.hpp"
#include "fika_molecular_basis.hpp"
#include "fika_point3d.hpp"
#include "fika_symmetric_tensor.hpp"
#include "fika_atom_pair_group.hpp"
#include "fika_two_centre_kernels.hpp"
#include "fika_solid_harmonics.hpp"
#include "fika_molecule.hpp"
#include "fika_far_field_expansion.hpp"

namespace fika::detail {

/// Smallest T from which the asymptotic Boys functions F_k(T) ~ Gamma(k + 1/2) / (2 T^(k+1/2)),
/// k <= n, have relative error <= 2^-53 (bound Q(n + 1/2, T) on the regularized upper incomplete
/// gamma function).
auto far_field_threshold(int n) -> double;

/// The point sources as seen from each atom. Charges with p_min(A) |C - A|^2 >= T_far(2 l_max(A)),
/// dipoles with p_min(A) |D - A|^2 >= T_far(2 l_max(A) + 1) and quadrupoles with
/// p_min(A) |E - A|^2 >= T_far(2 l_max(A) + 2) (one and two Boys orders more; p_min twice the
/// atom's smallest exponent, so the criterion holds for every shell pair on A) are far and enter
/// through the field tensor of the far sources at A, L <= 2 l_max(A), W = C - A:
///   Omega_LM(A) = sum_far q S_LM(W) / W^(2L+1) - (2L + 1) sum_far A+_LM(W, mu) / W^(2L+3)
///                 + (2L + 1)(2L + 3) / 2 sum_far B^(L+2)_LM(W, theta) / W^(2L+5)
/// (the dipole and quadrupole terms are mu . grad_W and 1/2 Theta : grad_W grad_W of the charge
/// term; dipole notes Eqs. 14 and 23, quadrupole notes Eqs. 14 and 18); the others are listed as
/// near. With a far-field expansion of the sources (whose near lists cover every source near an
/// atom), each atom visits only the near sources of its leaf and takes the others' field from the
/// leaf's local expansion.
class SourceFields {
 public:
  SourceFields(const Molecule<double>& molecule, const MolecularBasis& basis,
               const PointSources& sources, const FarFieldExpansion* far_field = nullptr);

  auto charges() const noexcept -> std::span<const double> { return sources_.charges; }
  auto charge_coordinates() const noexcept -> std::span<const Point3D<double>> {
    return sources_.charge_coordinates;
  }
  auto dipoles() const noexcept -> std::span<const Dipole> { return sources_.dipoles; }
  auto dipole_coordinates() const noexcept -> std::span<const Point3D<double>> {
    return sources_.dipole_coordinates;
  }
  auto quadrupoles() const noexcept -> std::span<const Quadrupole> { return sources_.quadrupoles; }
  auto quadrupole_coordinates() const noexcept -> std::span<const Point3D<double>> {
    return sources_.quadrupole_coordinates;
  }
  /// theta_mu (quadrupole_components) of each quadrupole.
  auto thetas() const noexcept -> std::span<const std::array<double, 5>> { return thetas_; }
  auto atom_position(std::size_t atom) const -> const Point3D<double>& { return atoms_[atom]; }

  /// Indices of the charges near `atom`.
  auto near(std::size_t atom) const -> std::span<const std::uint32_t> {
    return std::span(near_).subspan(near_offsets_[atom],
                                    near_offsets_[atom + 1] - near_offsets_[atom]);
  }

  /// Indices of the dipoles near `atom`.
  auto near_dipoles(std::size_t atom) const -> std::span<const std::uint32_t> {
    return std::span(near_dipoles_)
        .subspan(near_dipole_offsets_[atom],
                 near_dipole_offsets_[atom + 1] - near_dipole_offsets_[atom]);
  }

  /// Indices of the quadrupoles near `atom`.
  auto near_quadrupoles(std::size_t atom) const -> std::span<const std::uint32_t> {
    return std::span(near_quadrupoles_)
        .subspan(near_quadrupole_offsets_[atom],
                 near_quadrupole_offsets_[atom + 1] - near_quadrupole_offsets_[atom]);
  }

  /// Omega_LM of `atom`: row L^2 + M + L, L <= its order (2 l_max).
  auto omega(std::size_t atom) const -> std::span<const double> {
    return std::span(omega_).subspan(omega_offsets_[atom],
                                     omega_offsets_[atom + 1] - omega_offsets_[atom]);
  }

 private:
  PointSources sources_;
  std::vector<std::array<double, 5>> thetas_;
  std::vector<Point3D<double>> atoms_;
  std::vector<std::size_t> near_offsets_;  // atoms + 1
  std::vector<std::uint32_t> near_;
  std::vector<std::size_t> near_dipole_offsets_;  // atoms + 1
  std::vector<std::uint32_t> near_dipoles_;
  std::vector<std::size_t> near_quadrupole_offsets_;  // atoms + 1
  std::vector<std::uint32_t> near_quadrupoles_;
  std::vector<std::size_t> omega_offsets_;  // atoms + 1
  std::vector<double> omega_;
};

/// Geometry-independent data of a shell pair on one atom (bra shell l with K_A primitives and
/// N_A contracted functions, ket shell l' with K_B and N_B): p_ab (a-major), the far-field
/// moments Gamma(k + L + 3/2) / (2L + 1) p_ab^-(k + L + 3/2) per coupling L = |l - l'|..l + l'
/// (step 2, k = (l + l' - L)/2), and the dense coefficient matrices c^T (N_A x K_A), d (K_B x N_B).
struct SameAtomShellPair {
  int l = 0;
  int l_prime = 0;
  std::size_t bra_primitives = 0;
  std::size_t ket_primitives = 0;
  std::size_t bra_contractions = 0;
  std::size_t ket_contractions = 0;
  std::vector<double> p;
  std::vector<double> far_moments;       // couplings x K_A K_B
  std::vector<double> bra_coefficients;  // c^T: N_A x K_A
  std::vector<double> ket_coefficients;  // d: K_B x N_B
};

auto make_same_atom_shell_pair(const BasisShell& bra, const BasisShell& ket) -> SameAtomShellPair;

/// Far-field expansion of the point sources for the kernels of `basis`: rank 2 l_max, penetration
/// radius sqrt(T_far(2 l_max + n) / p_min) (n = 2 with quadrupoles, 1 with dipoles, else 0), and
/// tolerances eps_Lambda = threshold /
/// ((rank + 1) K_Lambda) with K_Lambda the largest far-field moment bound of a shell pair, so the
/// expansion changes no element of V by more than the threshold. The target region covers the
/// segments of the atom pairs within the largest screening cutoff (bounds of every kind present) of
/// their bases.
auto make_source_far_field(const Molecule<double>& molecule, const MolecularBasis& basis,
                           const PointSources& sources, double threshold) -> FarFieldExpansion;

/// Coulomb kernel of the one-centre block (nuclear-attraction notes):
///   G_kappa,Lambda(p, U^2) = kappa! sum_{i=0..kappa} U^(2i) / i! F_(Lambda+i)(p U^2) /
///   p^(kappa-i+1),
/// with `boys` holding F_0..F_(Lambda+kappa) at T = p U^2.
auto coulomb_kernel(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> double;

/// Radial kernels of the one-centre dipole block (dipole notes, Eq. 12):
///   H+ = -2 U^(2 kappa) F_(Lambda+kappa+1)(p U^2),
///   H- = (2 Lambda - 1) [G_kappa,Lambda - 2 U^(2 kappa + 2) F_(Lambda+kappa+1) / (2 Lambda + 1)]
/// (H- = 0 for Lambda = 0), with `boys` holding F_0..F_(Lambda+kappa+1).
struct DipoleKernels {
  double plus = 0.0;
  double minus = 0.0;
};

auto dipole_kernels(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> DipoleKernels;

/// One term of an angular dipole coupling of rank Lambda (dipole notes, Eq. 13):
/// A_Lambda,big_m += value d_(m) S_(Lambda +- 1),harmonic_m(U), with the Racah components of the
/// dipole d_(1) = d_x, d_(0) = d_z, d_(-1) = d_y.
struct DipoleCouplingTerm {
  int big_m;
  int m;
  int harmonic_m;
  double value;
};

/// Terms of A+ (C^{Lambda+1,M'}_{1m,Lambda M}, harmonics of degree Lambda + 1) and A-
/// (C^{Lambda M}_{1m,(Lambda-1)M'}, degree Lambda - 1; none for Lambda = 0) of rank Lambda <=
/// 2 max_angular_momentum (built once per rank, thread-safe).
struct DipoleCouplings {
  std::vector<DipoleCouplingTerm> plus;
  std::vector<DipoleCouplingTerm> minus;
};

auto dipole_couplings(int big_l) -> const DipoleCouplings&;

/// Racah component d_(m) (m = -1, 0, 1) of a dipole: y, z, x.
constexpr auto dipole_component(const Dipole& dipole, int m) noexcept -> double {
  return dipole.components[m == 1 ? 0U : (m == 0 ? 2U : 1U)];
}

/// Spherical components theta_mu (index mu + 2) of the traceless part of a primitive quadrupole
/// (quadrupole notes, Eq. 6): Theta : U U = sum_mu theta_mu S_2mu(U) = U . Theta . U.
auto quadrupole_components(const Quadrupole& quadrupole) -> std::array<double, 5>;

/// Radial kernels of the one-centre quadrupole block (quadrupole notes, Eqs. 11-12), with
/// x = U^2, G = G_kappa,Lambda(p, x), G' = -x^kappa F_(Lambda+kappa+1)(p x),
/// G'' = -kappa x^(kappa-1) F_(Lambda+kappa+1) + p x^kappa F_(Lambda+kappa+2):
///   K0 = G'',  K1 = x G'' + (Lambda + 3/2) G',
///   K2 = x^2 G'' + (2 Lambda + 1) x G' + (Lambda + 1/2)(Lambda - 1/2) G,
/// with `boys` holding F_0..F_(Lambda+kappa+2).
struct QuadrupoleKernels {
  double k0 = 0.0;
  double k1 = 0.0;
  double k2 = 0.0;
};

auto quadrupole_kernels(int kappa, int big_l, double p, double u2, std::span<const double> boys)
    -> QuadrupoleKernels;

/// One term of an angular quadrupole coupling of rank Lambda (quadrupole notes, Eq. 13):
/// B^(L')_Lambda,big_m += value theta_mu S_L',harmonic_m(U).
struct QuadrupoleCouplingTerm {
  int big_m;
  int mu;
  int harmonic_m;
  double value;
};

/// Terms of B^(Lambda+2) (with K0), B^(Lambda) (with K1) and B^(Lambda-2) (with K2) of rank
/// Lambda <= 2 max_angular_momentum, from C^{L'M'}_{2mu,Lambda M} (built once per rank,
/// thread-safe). The one-centre quadrupole block is
///   1/2 Theta : grad_U grad_U I_kappa,Lambda,M(p, U) = 4 pi (K0 B^(Lambda+2) + K1 B^(Lambda)
///                                                              + K2 B^(Lambda-2))_M.
struct QuadrupoleCouplings {
  std::vector<QuadrupoleCouplingTerm> raised;
  std::vector<QuadrupoleCouplingTerm> same;
  std::vector<QuadrupoleCouplingTerm> lowered;
};

auto quadrupole_couplings(int big_l) -> const QuadrupoleCouplings&;

/// Translation (addition-theorem) coefficient of S_lm(a + b) = sum T^lm_{lambda mu, nu}
/// S_lambda,mu(a) S_(l-lambda),nu(b): (l over lambda)!! C^{lm}_{lambda mu,(l-lambda) nu}.
struct TranslationTerm {
  int m;
  int lambda;
  int mu;
  int nu;
  double value;
};

/// All nonzero translation coefficients of order l (built once per l, thread-safe).
auto translation_terms(int l) -> std::span<const TranslationTerm>;

/// How the (A|V|B) kernel assembles V from the contracted fields: the flat list of Gamma terms,
/// or the factorized form V = 2 pi sum G^(lambda)(R) Z^(lambda lambda') H^(lambda')(R)^T
/// (geometric translation matrices and a Gaunt contraction); automatic picks whichever needs
/// fewer operations for the pair (l, l').
enum class AssemblyForm { automatic, flat, factorized };

/// Geometry-independent data of a shell pair on different atoms (bra shell l on A with K_A
/// primitives and N_A contracted functions, ket shell l' on B with K_B, N_B): exponents, the
/// dense coefficient matrices c^T (N_A x K_A) and d (K_B x N_B), the far-field thresholds
/// T_far(l + l') of charges, T_far(l + l' + 1) of dipoles and T_far(l + l' + 2) of quadrupoles,
/// and the primitive-pair screening data: a primitive pair (a, b) is skipped when
///   |c_a| |d_b| (L + 1) [(2 pi / p) sum|q| + (4 pi s_1 / sqrt(p)) sum|mu|
///                        + 4 pi s_2 sum||Theta||_2] max(1, R^L) exp(-mu R^2)
/// < threshold / (K_A K_B) (the charge, dipole and quadrupole shell-pair bounds at the primitive
/// level, so a block drops at most `threshold` in total), with |c_a|, |d_b| the largest over the
/// contracted functions.
struct TwoAtomShellPair {
  int l = 0;
  int l_prime = 0;
  std::vector<double> bra_exponents;
  std::vector<double> ket_exponents;
  std::size_t bra_contractions = 0;
  std::size_t ket_contractions = 0;
  std::vector<double> bra_coefficients;  // c^T: N_A x K_A
  std::vector<double> ket_coefficients;  // d: K_B x N_B
  double far_threshold = 0.0;
  double dipole_far_threshold = 0.0;
  double quadrupole_far_threshold = 0.0;
  std::vector<double> bra_largest;    // max_I |c_aI| per primitive
  std::vector<double> ket_largest;    // max_J |d_bJ| per primitive
  double bound_prefactor = 0.0;       // (L + 1) 2 pi sum|q|
  double dipole_prefactor = 0.0;      // (L + 1) 2 pi 2 s_1 sum|mu| (over sqrt(p))
  double quadrupole_prefactor = 0.0;  // (L + 1) 4 pi s_2 sum ||Theta||_2
  double primitive_threshold = 0.0;   // threshold / (K_A K_B); 0 disables the screening
  bool factorized = false;            // assembly form (resolved from AssemblyForm)
};

/// Total magnitudes of the point sources: sum |q_C|, sum |mu_D| and sum ||Theta_E||_2.
struct SourceSums {
  double charges = 0.0;
  double dipoles = 0.0;
  double quadrupoles = 0.0;
};

/// `sums` of the sources and the screening `threshold` of the matrix.
auto make_two_atom_shell_pair(const BasisShell& bra, const BasisShell& ket, const SourceSums& sums,
                              double threshold, AssemblyForm form = AssemblyForm::automatic)
    -> TwoAtomShellPair;

/// Bound on the far-field moments of the shell pair: sup over R of sum_ab of
///   W_Lambda(ab; R) = |c_a| |d_b| exp(-mu R^2) sum_ij C(l, i) C(l', j) (b R / p)^(l - i)
///                     (a R / p)^(l' - j) 2 pi Gamma((i + j + Lambda + 3) / 2)
///                     p^-((i+j+Lambda+3)/2)
/// >= integral |rho_ab| |r - P|^Lambda >= ||Q_Lambda(ab)|| (Euclidean norm over M of the moments
/// Q_Lambda,M = integral rho_ab S_Lambda,M(r - P) of the primitive product rho_ab about its
/// centre P), with |c_a|, |d_b| the largest over the contracted functions and the supremum
/// taken term by term (sup R^k exp(-mu R^2) = (k / (2 mu e))^(k/2)). An error delta Phi_Lambda
/// of the far-field tensor at every P changes each element of the shell-pair block by at most
/// sum_Lambda bound_Lambda ||delta Phi_Lambda||. Zero for Lambda > l + l'.
auto far_field_moment_bound(const BasisShell& bra, const BasisShell& ket, int rank) -> double;

/// Field rows of a shell pair (l, l') in Scheme II: each distinct (Lambda, kappa), Lambda <= l +
/// l', owns rows first_row .. first_row + 2 Lambda (one per M); the potential of any source is
/// linear in these rows (Y^(Lambda kappa)_M of the kernels).
struct FieldRowBlock {
  int big_l;
  int kappa;
  std::size_t first_row;
};

/// The field-row blocks of (l, l') in row order (built once per pair, thread-safe).
auto field_row_blocks(int l, int l_prime) -> std::span<const FieldRowBlock>;

/// Number of field rows of (l, l').
auto field_row_count(int l, int l_prime) -> std::size_t;

/// Weights of the field rows in the (A|V|B) block of primitive a (bra, on A) and primitive b (ket,
/// on B) of `pair`, without contraction coefficients: V_(m m') = sum_rows W_(m m'),row Y_row with
///   W = 2 pi exp(-mu R^2) sum (-t)^(l - lambda) (1 - t)^(l' - lambda') Gamma-term
/// (t = beta_b / p, product centre P = A - t R), for any point source. `harmonics` holds S up to
/// order l + l' of the single separation R = A - B (zero for one atom). Far sources give Y_row =
/// Gamma(kappa + Lambda + 3/2) / (2 Lambda + 1) p^-(kappa + Lambda + 3/2) times their field tensor
/// at P, so sum_kappa of W times these factors are the multipole moments
/// Q_LM(m, m') = integral chi_a,m chi_b,m' S_LM(r - P) d^3r of the product (rank <= l + l'); a near
/// dipole mu at P + U gives Y = H+_kappa,Lambda(p, U^2) A+_Lambda,M(U, mu) + H- A-. Layout:
/// component (m + l)(2l' + 1) + m' + l', then field row.
void primitive_pair_field_weights(const TwoAtomShellPair& pair, const SolidHarmonics& harmonics,
                                  std::size_t a, std::size_t b, std::span<double> weights);

/// Blocks of different atoms (A|V|B) of the potential of the point sources (charges, dipoles and
/// quadrupoles; Scheme II, dipole notes Eq. 22, quadrupole notes Eq. 23) for the n atom pairs
/// `pairs` (A != B) whose separations R = A - B are those of `harmonics`: layout as for the
/// same-atom blocks (all components). Uncontracted and segmented pairs fold c_a d_b into the
/// primitive prefactors; pairs with a general contraction contract c^T (prefactor * Y) d. Far
/// sources of a product centre enter through field tensors: charges S/U^(2L+1), dipoles
/// -(2L+1) A+/U^(2L+3), quadrupoles (2L+1)(2L+3)/2 B^(L+2)/U^(2L+5). With `far_field`, each
/// primitive pair sums only the near sources of the leaf holding its product centre and takes the
/// others from the leaf's local expansion.
void two_atom_potential_values(const TwoAtomShellPair& pair, const SourceFields& fields,
                               const SolidHarmonics& harmonics, std::span<const AtomPair> pairs,
                               std::size_t n, KernelWorkspace& workspace, std::span<double> values,
                               const FarFieldExpansion* far_field = nullptr);

/// Same-atom blocks (A|V|A) of the potential of the point sources (charges, dipoles and
/// quadrupoles; dipole notes Eq. 17, quadrupole notes Eq. 17), for the n atoms pairs[i] = (A, A):
/// N_A N_B (2l + 1)(2l' + 1) rows of n entries, row ((I N_B + J)(2l + 1) + m + l)(2l' + 1) + m' +
/// l' (all components).
void same_atom_potential_values(const SameAtomShellPair& pair, const SourceFields& fields,
                                std::span<const AtomPair> pairs, std::size_t n,
                                KernelWorkspace& workspace, std::span<double> values);

}  // namespace fika::detail

#endif  // fika_nuclear_attraction_kernels_hpp
