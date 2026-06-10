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

#include "ObaraSaikaFunc.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace osfunc {  // osfunc namespace

auto
compute_pa(const CScreenedBasisFunctionPair &pair) -> CDenseMatrix
{
    const auto bra_exponents = pair.bra_function().get_exponents();

    const auto ket_exponents = pair.ket_function().get_exponents();

    const auto npb = bra_exponents.size();

    const auto npk = ket_exponents.size();

    const auto npairs = static_cast<int>(pair.number_of_pairs());

    // bra and ket atom coordinates (one entry per atom pair)

    const auto *ax = pair.bra_x().data();
    const auto *ay = pair.bra_y().data();
    const auto *az = pair.bra_z().data();
    const auto *bx = pair.ket_x().data();
    const auto *by = pair.ket_y().data();
    const auto *bz = pair.ket_z().data();

    CDenseMatrix pamat(static_cast<int>(3 * npb * npk), npairs);

    auto *pvals = pamat.values();

    const auto block = npb * npk;  // number of primitive pairs per Cartesian component

    // exponents are identical for all atom pairs, so the prefactor depends only
    // on the primitive pair: loop over primitive pairs and vectorize over atoms

    for (size_t i = 0; i < npb; i++)
    {
        const auto a = bra_exponents[i];

        for (size_t j = 0; j < npk; j++)
        {
            const auto b = ket_exponents[j];

            const auto factor = -b / (a + b);

            const auto ij = i * npk + j;

            auto *pa_x = pvals + (0 * block + ij) * npairs;
            auto *pa_y = pvals + (1 * block + ij) * npairs;
            auto *pa_z = pvals + (2 * block + ij) * npairs;

#pragma omp simd
            for (int p = 0; p < npairs; p++)
            {
                pa_x[p] = factor * (ax[p] - bx[p]);
                pa_y[p] = factor * (ay[p] - by[p]);
                pa_z[p] = factor * (az[p] - bz[p]);
            }
        }
    }

    return pamat;
}

auto
compute_pb(const CScreenedBasisFunctionPair &pair) -> CDenseMatrix
{
    const auto bra_exponents = pair.bra_function().get_exponents();

    const auto ket_exponents = pair.ket_function().get_exponents();

    const auto npb = bra_exponents.size();

    const auto npk = ket_exponents.size();

    const auto npairs = static_cast<int>(pair.number_of_pairs());

    // bra and ket atom coordinates (one entry per atom pair)

    const auto *ax = pair.bra_x().data();
    const auto *ay = pair.bra_y().data();
    const auto *az = pair.bra_z().data();
    const auto *bx = pair.ket_x().data();
    const auto *by = pair.ket_y().data();
    const auto *bz = pair.ket_z().data();

    CDenseMatrix pbmat(static_cast<int>(3 * npb * npk), npairs);

    auto *pvals = pbmat.values();

    const auto block = npb * npk;  // number of primitive pairs per Cartesian component

    // exponents are identical for all atom pairs, so the prefactor depends only
    // on the primitive pair: loop over primitive pairs and vectorize over atoms

    for (size_t i = 0; i < npb; i++)
    {
        const auto a = bra_exponents[i];

        for (size_t j = 0; j < npk; j++)
        {
            const auto b = ket_exponents[j];

            const auto factor = a / (a + b);

            const auto ij = i * npk + j;

            auto *pb_x = pvals + (0 * block + ij) * npairs;
            auto *pb_y = pvals + (1 * block + ij) * npairs;
            auto *pb_z = pvals + (2 * block + ij) * npairs;

#pragma omp simd
            for (int p = 0; p < npairs; p++)
            {
                pb_x[p] = factor * (ax[p] - bx[p]);
                pb_y[p] = factor * (ay[p] - by[p]);
                pb_z[p] = factor * (az[p] - bz[p]);
            }
        }
    }

    return pbmat;
}

auto
compute_overlap(const CScreenedBasisFunctionPair &pair) -> CDenseMatrix
{
    const auto bra_exponents = pair.bra_function().get_exponents();

    const auto bra_coefficients = pair.bra_function().get_normalization_factors();

    const auto ket_exponents = pair.ket_function().get_exponents();

    const auto ket_coefficients = pair.ket_function().get_normalization_factors();

    const auto npb = bra_exponents.size();

    const auto npk = ket_exponents.size();

    const auto npairs = static_cast<int>(pair.number_of_pairs());

    // bra and ket atom coordinates (one entry per atom pair)

    const auto *ax = pair.bra_x().data();
    const auto *ay = pair.bra_y().data();
    const auto *az = pair.bra_z().data();
    const auto *bx = pair.ket_x().data();
    const auto *by = pair.ket_y().data();
    const auto *bz = pair.ket_z().data();

    CDenseMatrix smat(static_cast<int>(npb * npk), npairs);

    auto *svals = smat.values();

    const auto fpi = mathconst::pi_value();

    // exponents and coefficients are identical for all atom pairs, so the
    // prefactor and decay exponent depend only on the primitive pair: loop over
    // primitive pairs and vectorize over atoms

    for (size_t i = 0; i < npb; i++)
    {
        const auto a = bra_exponents[i];

        const auto ca = bra_coefficients[i];

        for (size_t j = 0; j < npk; j++)
        {
            const auto b = ket_exponents[j];

            const auto cb = ket_coefficients[j];

            const auto sab = a + b;

            const auto mu = a * b / sab;

            const auto prefactor = ca * cb * (fpi / sab) * std::sqrt(fpi / sab);

            auto *srow = svals + (i * npk + j) * npairs;

#pragma omp simd
            for (int p = 0; p < npairs; p++)
            {
                const auto abx = ax[p] - bx[p];
                const auto aby = ay[p] - by[p];
                const auto abz = az[p] - bz[p];

                const auto r2 = abx * abx + aby * aby + abz * abz;

                srow[p] = prefactor * std::exp(-mu * r2);
            }
        }
    }

    return smat;
}

auto
contract(CDenseMatrix &contracted, const CDenseMatrix &primitives) -> void
{
    const auto n_out = contracted.getNumberOfRows();

    const auto n_in = primitives.getNumberOfRows();

    const auto npairs = primitives.getNumberOfColumns();

    const auto factor = n_in / n_out;  // number of primitive pairs per contracted row

    auto *cvals = contracted.values();

    const auto *pvals = primitives.values();

    // accumulate the factor primitive rows of each contracted row, vectorizing
    // over the contiguous atom-pair dimension (the inner loop must be the simd
    // loop; an omp simd over an outer loop with an inner reduction does not
    // vectorize)

    for (int r = 0; r < n_out; r++)
    {
        auto *crow = cvals + static_cast<size_t>(r) * npairs;

        for (int k = 0; k < factor; k++)
        {
            const auto *prow = pvals + (static_cast<size_t>(r) * factor + k) * npairs;

#pragma omp simd
            for (int p = 0; p < npairs; p++)
            {
                crow[p] += prow[p];
            }
        }
    }
}

}  // namespace osfunc
