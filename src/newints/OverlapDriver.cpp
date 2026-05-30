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

#include "OverlapDriver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "MathConst.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

namespace newints {

namespace {

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
    const auto alpha = bra.get_exponents().back();

    const auto beta = ket.get_exponents().back();

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

    const auto estimate = pg * std::sqrt(pg) * std::max(1.0, std::pow(rab, lsum)) * std::exp(-alpha * beta / gamma * r2);

    return estimate >= threshold;
}

/// @brief Same-center overlap block for two shells of equal angular momentum l.
///
/// On a common center the overlap is diagonal in (l, m) and independent of m
/// (eq. 8 of concentric_spherical_overlap.pdf). In VeloxChem's unnormalized
/// solid-harmonic convention (the (l,0) normalization absorbed into the
/// contraction coefficients) the per-primitive-pair value is
/// (pi/p)^{3/2} (2l-1)!! / (2^l p^l), p = alpha + beta, so the contracted
/// block value is S_ab = sum_ij c_i c_j (pi/p)^{3/2} (2l-1)!! / (2^l p^l) on
/// every diagonal (m, m) entry and zero off the diagonal.
auto
overlap_kernel_diagonal(const CBasisFunction &bra, const CBasisFunction &ket) -> Block
{
    const auto l = bra.get_angular_momentum();  // == ket angular momentum (caller guarantees l_a == l_b)

    const auto exps_a = bra.get_exponents();

    const auto coefs_a = bra.get_normalization_factors();

    const auto exps_b = ket.get_exponents();

    const auto coefs_b = ket.get_normalization_factors();

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

            sab += coefs_a[i] * coefs_b[j] * std::pow(pi / p, 1.5) * dfact / (two_l * std::pow(static_cast<double>(p), l));
        }
    }

    // block is diagonal in m: only (m, m) entries are non-zero, all equal to S_ab
    const auto n = static_cast<std::size_t>(2 * l + 1);

    std::vector<double> data(n * n, 0.0);

    for (std::size_t m = 0; m < n; m++) data[m * n + m] = sab;

    return Block{n, n, data};
}

/// @brief STUB: two-center overlap block for two shells on different atoms.
/// Returns a zero block of the correct dimensions until the kernel is implemented.
auto
overlap_kernel(const CBasisFunction &bra,
               const CBasisFunction &ket,
               const TPoint<double> &bra_center,
               const TPoint<double> &ket_center) -> Block
{
    const std::size_t nrows = 2 * bra.get_angular_momentum() + 1;

    const std::size_t ncols = 2 * ket.get_angular_momentum() + 1;

    // TODO: evaluate two-center overlap integrals for this shell pair.
    return Block{nrows, ncols, std::vector<double>(nrows * ncols, 0.0)};
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

    // triangular loop over atom pairs with bra atom a <= ket atom b
    for (int a = 0; a < natoms; a++)
    {
        const auto bra_shells = atom_bases[basis_indices[a]].basis_functions();

        const auto bra_idx = outline.basis_function_indices(a);

        for (int b = a; b < natoms; b++)
        {
            if (a == b)
            {
                // same atom: only l == l' shell pairs contribute; triangular over shells
                for (std::size_t p = 0; p < bra_shells.size(); p++)
                {
                    for (std::size_t q = p; q < bra_shells.size(); q++)
                    {
                        if (bra_shells[p].get_angular_momentum() != bra_shells[q].get_angular_momentum()) continue;

                        matrix.add(bra_idx[p], bra_idx[q], overlap_kernel_diagonal(bra_shells[p], bra_shells[q]));
                    }
                }
            }
            else
            {
                // different atoms: all (l, l') shell pairs, screened
                const auto ket_shells = atom_bases[basis_indices[b]].basis_functions();

                const auto ket_idx = outline.basis_function_indices(b);

                for (std::size_t p = 0; p < bra_shells.size(); p++)
                {
                    for (std::size_t q = 0; q < ket_shells.size(); q++)
                    {
                        if (!overlap_screener(bra_shells[p], ket_shells[q], coords[a], coords[b], threshold)) continue;

                        matrix.add(bra_idx[p], ket_idx[q], overlap_kernel(bra_shells[p], ket_shells[q], coords[a], coords[b]));
                    }
                }
            }
        }
    }

    return matrix;
}

}  // namespace newints
