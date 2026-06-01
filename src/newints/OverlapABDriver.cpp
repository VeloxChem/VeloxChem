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

#include "OverlapABDriver.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "MathConst.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

#include "OverlapABRecSS.hpp"
#include "OverlapABRecSP.hpp"
#include "OverlapABRecSD.hpp"
#include "OverlapABRecSF.hpp"
#include "OverlapABRecSG.hpp"
#include "OverlapABRecSH.hpp"
#include "OverlapABRecSI.hpp"
#include "OverlapABRecPS.hpp"
#include "OverlapABRecPP.hpp"
#include "OverlapABRecPD.hpp"
#include "OverlapABRecPF.hpp"
#include "OverlapABRecPG.hpp"
#include "OverlapABRecPH.hpp"
#include "OverlapABRecPI.hpp"
#include "OverlapABRecDS.hpp"
#include "OverlapABRecDP.hpp"
#include "OverlapABRecDD.hpp"
#include "OverlapABRecDF.hpp"
#include "OverlapABRecDG.hpp"
#include "OverlapABRecDH.hpp"
#include "OverlapABRecDI.hpp"
#include "OverlapABRecFS.hpp"
#include "OverlapABRecFP.hpp"
#include "OverlapABRecFD.hpp"
#include "OverlapABRecFF.hpp"
#include "OverlapABRecFG.hpp"
#include "OverlapABRecFH.hpp"
#include "OverlapABRecFI.hpp"
#include "OverlapABRecGS.hpp"
#include "OverlapABRecGP.hpp"
#include "OverlapABRecGD.hpp"
#include "OverlapABRecGF.hpp"
#include "OverlapABRecGG.hpp"
#include "OverlapABRecGH.hpp"
#include "OverlapABRecGI.hpp"
#include "OverlapABRecHS.hpp"
#include "OverlapABRecHP.hpp"
#include "OverlapABRecHD.hpp"
#include "OverlapABRecHF.hpp"
#include "OverlapABRecHG.hpp"
#include "OverlapABRecHH.hpp"
#include "OverlapABRecHI.hpp"
#include "OverlapABRecIS.hpp"
#include "OverlapABRecIP.hpp"
#include "OverlapABRecID.hpp"
#include "OverlapABRecIF.hpp"
#include "OverlapABRecIG.hpp"
#include "OverlapABRecIH.hpp"
#include "OverlapABRecII.hpp"

namespace newints {

namespace {

/// @brief Integer power base^exp by repeated multiplication (cheaper than
/// std::pow for the small exponents that occur here, l = 0..6 and l_a + l_b).
inline auto
ipow(const double base, const int exp) -> double
{
    auto result = 1.0;

    for (int k = 0; k < exp; k++) result *= base;

    return result;
}

/// @brief STUB: screening test for an off-diagonal (two-center) shell pair.
/// Returns true if the shell-pair block must be computed. Its form will be
/// defined in a later step; for now nothing is screened out.
auto
overlap_screener(const CBasisFunction &bra,
                 const CBasisFunction &ket,
                 const TPoint<double> &bra_center,
                 const TPoint<double> &ket_center,
                 const double          threshold) -> bool
{
    // smallest exponents (exponents are stored in decreasing order -> last element)
    const auto alpha = bra.exponents().back();

    const auto beta = ket.exponents().back();

    const auto gamma = alpha + beta;

    // squared bra-ket separation R_AB^2
    const auto ra = bra_center.coordinates();

    const auto rb = ket_center.coordinates();

    const auto dx = ra[0] - rb[0];

    const auto dy = ra[1] - rb[1];

    const auto dz = ra[2] - rb[2];

    const auto r2 = dx * dx + dy * dy + dz * dz;

    const auto rab = std::sqrt(r2);

    const auto lsum = bra.get_angular_momentum() + ket.get_angular_momentum();

    // coarse overlap estimate using the smallest-exponent Gaussian pair
    const auto pg = mathconst::pi_value() / gamma;

    const auto estimate = pg * std::sqrt(pg) * std::max(1.0, ipow(rab, lsum)) * std::exp(-alpha * beta / gamma * r2);

    return estimate >= threshold;
}

/// @brief Same-center overlap value for two shells of equal angular momentum l.
///
/// On a common center the overlap is diagonal in (l, m) and independent of m
/// (eq. 8 of concentric_spherical_overlap.pdf). In VeloxChem's unnormalized
/// solid-harmonic convention (the (l,0) normalization absorbed into the
/// contraction coefficients) the per-primitive-pair value is
/// (pi/p)^{3/2} (2l-1)!! / (2^l p^l), p = alpha + beta, so the contracted value
/// is S_ab = sum_ij c_i c_j (pi/p)^{3/2} (2l-1)!! / (2^l p^l), placed on every
/// diagonal (m, m) entry by the caller.
auto
overlap_diagonal_value(const CBasisFunction &bra, const CBasisFunction &ket) -> double
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

    auto sab = 0.0;

    for (std::size_t i = 0; i < exps_a.size(); i++)
    {
        for (std::size_t j = 0; j < exps_b.size(); j++)
        {
            const auto p = exps_a[i] + exps_b[j];

            const auto pq = pi / p;  // (pi/p)^{3/2} = (pi/p) * sqrt(pi/p)

            sab += coefs_a[i] * coefs_b[j] * pq * std::sqrt(pq) * dfact / (two_l * ipow(p, l));
        }
    }

    return sab;
}

/// @brief Function-pointer type of a Tabula-generated spherical overlap kernel:
/// it writes the (2 l_a + 1) x (2 l_b + 1) row-major block into the buffer.
using kernel_fn = void (*)(const CBasisFunction &, const CBasisFunction &, const TPoint<double> &, const TPoint<double> &, double *);

/// @brief [l_a][l_b] dispatch table to the spherical (l_a | l_b) kernel, l = 0 (S) .. 6 (I).
constexpr kernel_fn s_kernels[7][7] = {
    {ovlab::overlap_s_s, ovlab::overlap_s_p, ovlab::overlap_s_d, ovlab::overlap_s_f, ovlab::overlap_s_g, ovlab::overlap_s_h, ovlab::overlap_s_i},
    {ovlab::overlap_p_s, ovlab::overlap_p_p, ovlab::overlap_p_d, ovlab::overlap_p_f, ovlab::overlap_p_g, ovlab::overlap_p_h, ovlab::overlap_p_i},
    {ovlab::overlap_d_s, ovlab::overlap_d_p, ovlab::overlap_d_d, ovlab::overlap_d_f, ovlab::overlap_d_g, ovlab::overlap_d_h, ovlab::overlap_d_i},
    {ovlab::overlap_f_s, ovlab::overlap_f_p, ovlab::overlap_f_d, ovlab::overlap_f_f, ovlab::overlap_f_g, ovlab::overlap_f_h, ovlab::overlap_f_i},
    {ovlab::overlap_g_s, ovlab::overlap_g_p, ovlab::overlap_g_d, ovlab::overlap_g_f, ovlab::overlap_g_g, ovlab::overlap_g_h, ovlab::overlap_g_i},
    {ovlab::overlap_h_s, ovlab::overlap_h_p, ovlab::overlap_h_d, ovlab::overlap_h_f, ovlab::overlap_h_g, ovlab::overlap_h_h, ovlab::overlap_h_i},
    {ovlab::overlap_i_s, ovlab::overlap_i_p, ovlab::overlap_i_d, ovlab::overlap_i_f, ovlab::overlap_i_g, ovlab::overlap_i_h, ovlab::overlap_i_i},
};

/// @brief Two-center overlap block for two shells on different atoms, written
/// row-major into out (which must hold (2 l_a + 1) * (2 l_b + 1) doubles).
/// Dispatches on (l_a, l_b) to the matching ovlab kernel; pairs outside the
/// S..I grid (l > 6) leave out as supplied by the caller (a zero block).
auto
overlap_kernel(const CBasisFunction &bra,
               const CBasisFunction &ket,
               const TPoint<double> &bra_center,
               const TPoint<double> &ket_center,
               double                *out) -> void
{
    const auto la = bra.get_angular_momentum();

    const auto lb = ket.get_angular_momentum();

    if (la < 0 || la > 6 || lb < 0 || lb > 6) return;

    s_kernels[la][lb](bra, ket, bra_center, ket_center, out);
}

}  // namespace

auto
OverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis, const double threshold) const -> SparseMatrix
{
    // the two-center overlap matrix is symmetric
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

    // flatten the triangular atom-pair loop (bra atom a <= ket atom b) into a work list
    std::vector<std::pair<int, int>> atom_pairs;

    atom_pairs.reserve(static_cast<std::size_t>(natoms) * (natoms + 1) / 2);

    for (int a = 0; a < natoms; a++)
    {
        for (int b = a; b < natoms; b++) atom_pairs.push_back({a, b});
    }

    // each thread writes block payloads directly into its own data arena (the
    // kernels take a raw double* output), recording a SparseMatrix::RawBlock per
    // block. No per-block allocation and no shared-container mutation in the
    // parallel region; the arenas are merged in bulk afterwards.
    struct ThreadArena
    {
        std::vector<double> data;
        std::vector<SparseMatrix::RawBlock> meta;
    };

    std::vector<ThreadArena> arenas(static_cast<std::size_t>(omp_get_max_threads()));

    // light warmup of each thread's arena to skip the first few reallocations
    for (auto &arena : arenas)
    {
        arena.data.reserve(1u << 14);

        arena.meta.reserve(1u << 10);
    }

    const auto npairs = static_cast<int>(atom_pairs.size());

    // dynamic scheduling balances the load: a heavy-heavy atom pair does far more
    // work than a light-light one
#pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < npairs; t++)
    {
        const auto a = atom_pairs[t].first;

        const auto b = atom_pairs[t].second;

        auto &arena = arenas[static_cast<std::size_t>(omp_get_thread_num())];

        const auto &bra_shells = atom_shells[a];

        const auto &bra_idx = atom_indices[a];

        if (a == b)
        {
            // same atom: only l == l' shell pairs contribute; triangular over shells.
            // The block is diagonal in m with value S_ab on every (m, m) entry.
            for (std::size_t p = 0; p < bra_shells.size(); p++)
            {
                for (std::size_t q = p; q < bra_shells.size(); q++)
                {
                    if (bra_shells[p].get_angular_momentum() != bra_shells[q].get_angular_momentum()) continue;

                    const auto sab = overlap_diagonal_value(bra_shells[p], bra_shells[q]);

                    const auto n = static_cast<std::size_t>(2 * bra_shells[p].get_angular_momentum() + 1);

                    const auto off = arena.data.size();

                    if (p == q)
                    {
                        // diagonal block (i == j): store packed lower-triangular
                        arena.data.resize(off + n * (n + 1) / 2);  // resize zero-fills

                        for (std::size_t m = 0; m < n; m++) arena.data[off + m * (m + 1) / 2 + m] = sab;

                        arena.meta.push_back(SparseMatrix::RawBlock{bra_idx[p], bra_idx[q], n, n, off, Kind::lower_triangular});
                    }
                    else
                    {
                        // same-atom off-diagonal (i < j): full n x n block, diagonal in m
                        arena.data.resize(off + n * n);  // resize zero-fills

                        for (std::size_t m = 0; m < n; m++) arena.data[off + m * n + m] = sab;

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
                    if (!overlap_screener(bra_shells[p], ket_shells[q], coords[a], coords[b], threshold)) continue;

                    const auto nr = static_cast<std::size_t>(2 * bra_shells[p].get_angular_momentum() + 1);

                    const auto nc = static_cast<std::size_t>(2 * ket_shells[q].get_angular_momentum() + 1);

                    const auto off = arena.data.size();

                    arena.data.resize(off + nr * nc);  // resize zero-fills (also the l > 6 fallback)

                    // pointer taken after resize, so a reallocation cannot dangle mid-kernel
                    overlap_kernel(bra_shells[p], ket_shells[q], coords[a], coords[b], arena.data.data() + off);

                    arena.meta.push_back(SparseMatrix::RawBlock{bra_idx[p], ket_idx[q], nr, nc, off, Kind::full});
                }
            }
        }
    }

    // serial bulk merge: append each thread's blocks (keys are disjoint across atom
    // pairs, so insertion order does not affect the result)
    std::size_t total = 0;

    std::size_t total_data = 0;

    for (const auto &arena : arenas)
    {
        total += arena.meta.size();

        total_data += arena.data.size();
    }

    matrix.reserve(total);

    matrix.reserve_data(total_data);

    for (const auto &arena : arenas) matrix.append_arena(arena.data, arena.meta);

    return matrix;
}

}  // namespace newints
