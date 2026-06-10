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

}  // namespace osfunc
