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

#include "ElectronRepulsionABDriver.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <utility>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "MathConst.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

#include "ElectronRepulsionABRecSS.hpp"
#include "ElectronRepulsionABRecPS.hpp"
#include "ElectronRepulsionABRecPP.hpp"
#include "ElectronRepulsionABRecDS.hpp"
#include "ElectronRepulsionABRecDP.hpp"
#include "ElectronRepulsionABRecDD.hpp"
#include "ElectronRepulsionABRecFS.hpp"
#include "ElectronRepulsionABRecFP.hpp"
#include "ElectronRepulsionABRecFD.hpp"
#include "ElectronRepulsionABRecFF.hpp"
#include "ElectronRepulsionABRecGS.hpp"
#include "ElectronRepulsionABRecGP.hpp"
#include "ElectronRepulsionABRecGD.hpp"
#include "ElectronRepulsionABRecGF.hpp"
#include "ElectronRepulsionABRecGG.hpp"
#include "ElectronRepulsionABRecHS.hpp"
#include "ElectronRepulsionABRecHP.hpp"
#include "ElectronRepulsionABRecHD.hpp"
#include "ElectronRepulsionABRecHF.hpp"
#include "ElectronRepulsionABRecHG.hpp"
#include "ElectronRepulsionABRecHH.hpp"
#include "ElectronRepulsionABRecIS.hpp"
#include "ElectronRepulsionABRecIP.hpp"
#include "ElectronRepulsionABRecID.hpp"
#include "ElectronRepulsionABRecIF.hpp"
#include "ElectronRepulsionABRecIG.hpp"
#include "ElectronRepulsionABRecIH.hpp"
#include "ElectronRepulsionABRecII.hpp"

namespace newints {

namespace {

/// @brief Integer power base^exp by repeated multiplication (cheaper than
/// std::pow for the small exponents that occur here, l = 0..6).
inline auto
ipow(const double base, const int exp) -> double
{
    auto result = 1.0;

    for (int k = 0; k < exp; k++) result *= base;

    return result;
}

/// @brief Screening test for an off-diagonal (two-center) Coulomb shell pair.
/// Returns true if the block must be computed.
///
/// Implements the per-pair pre-screener of per_pair_predicate.pdf (eq. 9): skip
/// the pair if
///   2 pi^{5/2} / (a_min b_min sqrt(a_min + b_min))
///     * min(1, sqrt(pi) / (2 sqrt(rho_min) R_AB)) * C(l_a, l_b) < tau,
/// where the first factor is the short-range (s|s) value, the min(...) is the
/// Boys-function envelope F0(rho R^2) <= min(1, 1/2 sqrt(pi/T)), and
///   C(l_a,l_b) = sqrt((2 l_a + 1)!! (2 l_b + 1)!!)
///                * (1 + a_min R^2)^{l_a/2} (1 + b_min R^2)^{l_b/2}.
/// a_min, b_min are the smallest exponents of each shell (stored in decreasing
/// order, so the smallest is the last one); rho_min = a_min b_min / (a_min + b_min).
/// This is a Level-1 scaling pre-screen using only exponents (no contraction
/// coefficients); it never underestimates the magnitude.
auto
coulomb_screener(const CBasisFunction &bra,
                 const CBasisFunction &ket,
                 const TPoint<double> &bra_center,
                 const TPoint<double> &ket_center,
                 const double          threshold) -> bool
{
    // smallest exponents (exponents are stored in decreasing order -> last element)
    const auto alpha = bra.exponents().back();

    const auto beta = ket.exponents().back();

    const auto pmin = alpha + beta;

    const auto rho = alpha * beta / pmin;

    // squared bra-ket separation R_AB^2
    const auto ra = bra_center.coordinates();

    const auto rb = ket_center.coordinates();

    const auto dx = ra[0] - rb[0];

    const auto dy = ra[1] - rb[1];

    const auto dz = ra[2] - rb[2];

    const auto r2 = dx * dx + dy * dy + dz * dz;

    const auto la = bra.get_angular_momentum();

    const auto lb = ket.get_angular_momentum();

    const auto pi = mathconst::pi_value();

    // short-range (s|s) value 2 pi^{5/2} / (a_min b_min sqrt(p_min))
    const auto two_pi52 = 2.0 * pi * pi * std::sqrt(pi);

    const auto short_range = two_pi52 / (alpha * beta * std::sqrt(pmin));

    // Boys-function envelope F0(rho R^2) <= min(1, 1/2 sqrt(pi / (rho R^2)))
    auto boys_env = 1.0;

    const auto t = rho * r2;

    if (t > 0.0)
    {
        const auto bound = 0.5 * std::sqrt(pi / t);

        if (bound < boys_env) boys_env = bound;
    }

    // combinatorial prefactor sqrt((2 l_a + 1)!! (2 l_b + 1)!!)
    auto dfact = 1.0;

    for (int k = 2 * la + 1; k > 0; k -= 2) dfact *= static_cast<double>(k);

    for (int k = 2 * lb + 1; k > 0; k -= 2) dfact *= static_cast<double>(k);

    const auto gamma = std::sqrt(dfact);

    // angular envelope (1 + a_min R^2)^{l_a/2} (1 + b_min R^2)^{l_b/2}. The
    // exponents are half-integer, so (as in the kinetic screener) compute them as
    // integer powers plus a single sqrt rather than via std::pow's slow path.
    const auto base_a = 1.0 + alpha * r2;

    const auto base_b = 1.0 + beta * r2;

    auto ang = ipow(base_a, la / 2) * ipow(base_b, lb / 2);

    if (la & 1) ang *= std::sqrt(base_a);

    if (lb & 1) ang *= std::sqrt(base_b);

    const auto estimate = short_range * boys_env * gamma * ang;

    return estimate >= threshold;
}

/// @brief Same-center two-center Coulomb value for two shells of equal angular
/// momentum l.
///
/// On a common center (R_AB = 0) the Coulomb block is diagonal in (l, m) and
/// independent of m, like the overlap and kinetic-energy blocks. For concentric
/// solid-harmonic Gaussians the momentum-space Coulomb integral gives, in
/// VeloxChem's unnormalized (Racah, int Y^2 dOmega = 4 pi / (2l+1)) convention,
/// the per-primitive-pair value
///   2 pi^{5/2} (2l-1)!! / ((2l+1) 2^l) * 1 / (alpha beta p^l sqrt(p)),  p = alpha + beta,
/// placed on every diagonal (m, m) entry by the caller. (Reduces to
/// 2 pi^{5/2} / (alpha beta sqrt(p)) for l = 0, matching the (s|s) seed.)
auto
coulomb_diagonal_value(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto l = bra.get_angular_momentum();  // == ket angular momentum (caller guarantees l_a == l_b)

    const auto &exps_a = bra.exponents();

    const auto &coefs_a = bra.normalization_factors();

    const auto &exps_b = ket.exponents();

    const auto &coefs_b = ket.normalization_factors();

    // (2l - 1)!! with (-1)!! = 1
    auto dfact = 1.0;

    for (int k = 2 * l - 1; k > 0; k -= 2) dfact *= static_cast<double>(k);

    const auto two_l = static_cast<double>(1 << l);

    const auto pi = mathconst::pi_value();

    const auto two_pi52 = 2.0 * pi * pi * std::sqrt(pi);

    // 2 pi^{5/2} (2l-1)!! / ((2l+1) 2^l)
    const auto prefac = two_pi52 * dfact / (static_cast<double>(2 * l + 1) * two_l);

    auto vab = 0.0;

    for (std::size_t i = 0; i < exps_a.size(); i++)
    {
        for (std::size_t j = 0; j < exps_b.size(); j++)
        {
            const auto alpha = exps_a[i];

            const auto beta = exps_b[j];

            const auto p = alpha + beta;

            // 1 / (alpha beta p^l sqrt(p))
            const auto denom = alpha * beta * ipow(p, l) * std::sqrt(p);

            vab += coefs_a[i] * coefs_b[j] * prefac / denom;
        }
    }

    return vab;
}

/// @brief Function-pointer type of a Tabula-generated spherical two-center Coulomb
/// kernel: it writes the (2 l_a + 1) x (2 l_b + 1) row-major block into the buffer.
using kernel_fn = void (*)(const CBasisFunction &, const CBasisFunction &, const TPoint<double> &, const TPoint<double> &, double *);

/// @brief Lower-triangular [l_a][l_b] dispatch table (l_a >= l_b) to the spherical
/// (l_a | l_b) Coulomb kernel, l = 0 (S) .. 6 (I). Upper-triangle entries are null:
/// an (l_a < l_b) block is obtained by transposing the (l_b, l_a) kernel, since the
/// two-center Coulomb matrix is symmetric.
constexpr kernel_fn t_kernels[7][7] = {
    {eri2cab::eri2c_s_s, nullptr,             nullptr,             nullptr,             nullptr,             nullptr,             nullptr},
    {eri2cab::eri2c_p_s, eri2cab::eri2c_p_p,  nullptr,             nullptr,             nullptr,             nullptr,             nullptr},
    {eri2cab::eri2c_d_s, eri2cab::eri2c_d_p,  eri2cab::eri2c_d_d,  nullptr,             nullptr,             nullptr,             nullptr},
    {eri2cab::eri2c_f_s, eri2cab::eri2c_f_p,  eri2cab::eri2c_f_d,  eri2cab::eri2c_f_f,  nullptr,             nullptr,             nullptr},
    {eri2cab::eri2c_g_s, eri2cab::eri2c_g_p,  eri2cab::eri2c_g_d,  eri2cab::eri2c_g_f,  eri2cab::eri2c_g_g,  nullptr,             nullptr},
    {eri2cab::eri2c_h_s, eri2cab::eri2c_h_p,  eri2cab::eri2c_h_d,  eri2cab::eri2c_h_f,  eri2cab::eri2c_h_g,  eri2cab::eri2c_h_h,  nullptr},
    {eri2cab::eri2c_i_s, eri2cab::eri2c_i_p,  eri2cab::eri2c_i_d,  eri2cab::eri2c_i_f,  eri2cab::eri2c_i_g,  eri2cab::eri2c_i_h,  eri2cab::eri2c_i_i},
};

/// @brief Two-center Coulomb block for two shells on different atoms, written
/// row-major into out (which must hold (2 l_a + 1) * (2 l_b + 1) doubles). For
/// l_a >= l_b dispatches to the matching eri2cab kernel directly; for l_a < l_b it
/// evaluates the (l_b, l_a) kernel (which exists) into a temporary and writes the
/// transpose, using the symmetry (a|b) = (b|a)^T. Pairs outside the S..I grid
/// (l > 6) leave out as an explicit zero block (the arena is not pre-zeroed).
auto
coulomb_kernel(const CBasisFunction &bra,
               const CBasisFunction &ket,
               const TPoint<double> &bra_center,
               const TPoint<double> &ket_center,
               double                *out) -> void
{
    const auto la = bra.get_angular_momentum();

    const auto lb = ket.get_angular_momentum();

    if (la < 0 || la > 6 || lb < 0 || lb > 6)
    {
        const auto n = (2 * la + 1) * (2 * lb + 1);

        for (int x = 0; x < n; x++) out[x] = 0.0;

        return;
    }

    if (la >= lb)
    {
        t_kernels[la][lb](bra, ket, bra_center, ket_center, out);
    }
    else
    {
        // evaluate the (l_b, l_a) kernel into a temporary, then transpose into out.
        // largest possible block is (2*6+1)^2 = 169 doubles.
        const auto nr = 2 * la + 1;  // rows of out

        const auto nc = 2 * lb + 1;  // cols of out

        double tmp[169];

        // tmp holds the (l_b | l_a) block: nc rows x nr cols, row-major
        t_kernels[lb][la](ket, bra, ket_center, bra_center, tmp);

        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
            {
                out[i * nc + j] = tmp[j * nr + i];
            }
        }
    }
}

}  // namespace

auto
ElectronRepulsionDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis, const double threshold) const -> SparseMatrix
{
    // the two-center Coulomb matrix is symmetric
    SparseMatrix matrix(SymmetryType::symmetric);

    const auto outline = MolecularBasisOutline(basis);

    const auto coords = molecule.coordinates("au");

    const auto basis_indices = basis.basis_sets_indices();

    const auto atom_bases = basis.basis_sets();

    const auto natoms = static_cast<int>(outline.number_of_atoms());

    // cache each atom's shells and contracted-GTO indices once, so the atom-pair
    // loop reads them by reference instead of reallocating these vectors on every
    // iteration (an atom appears in many pairs)
    std::vector<std::vector<CBasisFunction>> atom_shells(natoms);

    std::vector<std::vector<int>> atom_indices(natoms);

    for (int a = 0; a < natoms; a++)
    {
        atom_shells[a]  = atom_bases[basis_indices[a]].basis_functions();
        atom_indices[a] = outline.basis_function_indices(a);
    }

    // row offsets of the triangular atom-pair index (bra atom a <= ket atom b):
    // row_off[a] = number of pairs in rows before a. This O(N) table replaces the
    // O(N^2) materialized atom-pair list; the flat index t maps to (a, b) by a
    // binary search in the loop below.
    std::vector<std::size_t> row_off(static_cast<std::size_t>(natoms) + 1);

    row_off[0] = 0;

    for (int a = 0; a < natoms; a++) row_off[a + 1] = row_off[a] + static_cast<std::size_t>(natoms - a);

    // each thread writes block payloads directly into its own data arena (the
    // kernels take a raw double* output), recording a SparseMatrix::RawBlock per
    // block. No per-block allocation and no shared-container mutation in the
    // parallel region; the arenas are merged in bulk afterwards. The data arena
    // uses a non-zeroing allocator: off-diagonal kernels write every entry, and the
    // (few) diagonal blocks are zeroed explicitly, so resize() need not memset.
    struct ThreadArena
    {
        std::vector<double, DefaultInitAllocator<double>> data;
        std::vector<SparseMatrix::RawBlock> meta;
    };

    std::vector<ThreadArena> arenas(static_cast<std::size_t>(omp_get_max_threads()));

    // light warmup of each thread's arena to skip the first few reallocations
    for (auto &arena : arenas)
    {
        arena.data.reserve(1u << 14);

        arena.meta.reserve(1u << 10);
    }

    const auto npairs = static_cast<int>(row_off[natoms]);

    // dynamic scheduling balances the load: a heavy-heavy atom pair does far more
    // work than a light-light one
#pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < npairs; t++)
    {
        // map the flat triangular index t -> (a, b) with a <= b
        const auto row = std::upper_bound(row_off.begin(), row_off.end(), static_cast<std::size_t>(t));

        const int a = static_cast<int>(row - row_off.begin()) - 1;

        const int b = a + static_cast<int>(static_cast<std::size_t>(t) - row_off[a]);

        auto &arena = arenas[static_cast<std::size_t>(omp_get_thread_num())];

        const auto &bra_shells = atom_shells[a];

        const auto &bra_idx = atom_indices[a];

        if (a == b)
        {
            // same atom: only l == l' shell pairs contribute; triangular over shells.
            // The block is diagonal in m with value v_ab on every (m, m) entry.
            for (std::size_t p = 0; p < bra_shells.size(); p++)
            {
                for (std::size_t q = p; q < bra_shells.size(); q++)
                {
                    if (bra_shells[p].get_angular_momentum() != bra_shells[q].get_angular_momentum()) continue;

                    const auto vab = coulomb_diagonal_value(bra_shells[p], bra_shells[q]);

                    const auto n = static_cast<std::size_t>(2 * bra_shells[p].get_angular_momentum() + 1);

                    const auto off = arena.data.size();

                    if (p == q)
                    {
                        // diagonal block (i == j): packed lower-triangular, only the
                        // (m, m) entries non-zero -> zero the slot, then write them
                        const auto np = n * (n + 1) / 2;

                        arena.data.resize(off + np);

                        std::fill(arena.data.begin() + off, arena.data.begin() + off + np, 0.0);

                        for (std::size_t m = 0; m < n; m++) arena.data[off + m * (m + 1) / 2 + m] = vab;

                        arena.meta.push_back(SparseMatrix::RawBlock{bra_idx[p], bra_idx[q], n, n, off, Kind::lower_triangular});
                    }
                    else
                    {
                        // same-atom off-diagonal (i < j): full n x n block, diagonal in m
                        arena.data.resize(off + n * n);

                        std::fill(arena.data.begin() + off, arena.data.begin() + off + n * n, 0.0);

                        for (std::size_t m = 0; m < n; m++) arena.data[off + m * n + m] = vab;

                        arena.meta.push_back(SparseMatrix::RawBlock{bra_idx[p], bra_idx[q], n, n, off, Kind::full});
                    }
                }
            }
        }
        else
        {
            // different atoms: all (l, l') shell pairs, screened
            const auto &ket_shells = atom_shells[b];

            const auto &ket_idx = atom_indices[b];

            for (std::size_t p = 0; p < bra_shells.size(); p++)
            {
                for (std::size_t q = 0; q < ket_shells.size(); q++)
                {
                    if (!coulomb_screener(bra_shells[p], ket_shells[q], coords[a], coords[b], threshold)) continue;

                    const auto nr = static_cast<std::size_t>(2 * bra_shells[p].get_angular_momentum() + 1);

                    const auto nc = static_cast<std::size_t>(2 * ket_shells[q].get_angular_momentum() + 1);

                    const auto off = arena.data.size();

                    arena.data.resize(off + nr * nc);  // uninitialized; the kernel writes every entry

                    // pointer taken after resize, so a reallocation cannot dangle mid-kernel
                    coulomb_kernel(bra_shells[p], ket_shells[q], coords[a], coords[b], arena.data.data() + off);

                    arena.meta.push_back(SparseMatrix::RawBlock{bra_idx[p], ket_idx[q], nr, nc, off, Kind::full});
                }
            }
        }
    }

    // per-thread base offsets into the merged storage (serial, O(nthreads));
    // the heavy payload/descriptor copy is parallelized below
    const std::size_t nt = arenas.size();

    std::vector<std::size_t> data_base(nt + 1, 0);

    std::vector<std::size_t> meta_base(nt + 1, 0);

    for (std::size_t t = 0; t < nt; t++)
    {
        data_base[t + 1] = data_base[t] + arenas[t].data.size();

        meta_base[t + 1] = meta_base[t] + arenas[t].meta.size();
    }

    matrix.prepare(meta_base[nt], data_base[nt]);

    double *const dst = matrix.data();

    // parallel merge: each thread copies its arena into a disjoint slice of the
    // matrix data and fills its disjoint range of block descriptors (no serial
    // bottleneck beyond the prefix sums above)
#pragma omp parallel for schedule(static)
    for (int t = 0; t < static_cast<int>(nt); t++)
    {
        const auto &arena = arenas[static_cast<std::size_t>(t)];

        if (!arena.data.empty()) std::memcpy(dst + data_base[t], arena.data.data(), arena.data.size() * sizeof(double));

        auto m = meta_base[t];

        for (const auto &rb : arena.meta)
        {
            matrix.set_block(m++, SparseMatrix::RawBlock{rb.i, rb.j, rb.nrows, rb.ncols, data_base[t] + rb.offset, rb.kind});
        }
    }

    return matrix;
}

}  // namespace newints
