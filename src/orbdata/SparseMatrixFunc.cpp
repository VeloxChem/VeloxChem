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

#include "SparseMatrixFunc.hpp"

#include "Matrix.hpp"

namespace sparsematfunc {  // sparsematfunc namespace

auto
to_dense_matrix(const CSparseMatrix                           &matrix,
                const std::vector<CScreenedBasisFunctionPair> &screened_pairs,
                const CMolecularBasis                         &basis) -> CDenseMatrix
{
    const auto naos = static_cast<int>(basis.dimensions_of_basis());

    CDenseMatrix dmat(naos, naos);

    dmat.zero();

    if (naos == 0) return dmat;

    const auto max_l = basis.max_angular_momentum();

    // per angular momentum: base atomic orbital index of the l block and the
    // component stride (number of contracted functions of that momentum)

    std::vector<int> base(max_l + 1);

    std::vector<int> stride(max_l + 1);

    for (int l = 0; l <= max_l; l++)
    {
        base[l]   = static_cast<int>(basis.dimensions_of_basis(l));
        stride[l] = static_cast<int>(basis.number_of_basis_functions(l));
    }

    // first-component atomic orbital index of each (atom, local basis function),
    // assigned in molecule order with per-momentum position counters

    const auto indices = basis.basis_sets_indices();

    const auto basis_sets = basis.basis_sets();

    const auto natoms = indices.size();

    std::vector<std::vector<int>> ao_first(natoms);

    std::vector<int> counters(max_l + 1, 0);

    for (size_t i = 0; i < natoms; i++)
    {
        const auto functions = basis_sets[indices[i]].basis_functions();

        ao_first[i].resize(functions.size());

        for (size_t k = 0; k < functions.size(); k++)
        {
            const auto l = functions[k].get_angular_momentum();

            ao_first[i][k] = base[l] + counters[l];

            counters[l]++;
        }
    }

    // scatter the allocated blocks into the dense matrix

    const auto mtype = matrix.type();

    const auto sign = (mtype == mat_t::antisymmetric) ? -1.0 : 1.0;

    auto *dvals = dmat.values();

    for (const auto key : matrix.keys())
    {
        const auto &pair = screened_pairs[key];

        const auto &block = matrix.block(key);

        const auto la = pair.bra_function().get_angular_momentum();

        const auto lb = pair.ket_function().get_angular_momentum();

        const auto ka = pair.bra_bf_index();

        const auto kb = pair.ket_bf_index();

        const auto nca = 2 * la + 1;

        const auto ncb = 2 * lb + 1;

        const auto bra_stride = stride[la];

        const auto ket_stride = stride[lb];

        const auto &bra_atoms = pair.bra_atoms();

        const auto &ket_atoms = pair.ket_atoms();

        const auto npairs = block.getNumberOfColumns();

        const auto *bvals = block.values();

        for (int i = 0; i < npairs; i++)
        {
            const auto bra_first = ao_first[bra_atoms[i]][ka];

            const auto ket_first = ao_first[ket_atoms[i]][kb];

            for (int ma = 0; ma < nca; ma++)
            {
                const auto row = bra_first + ma * bra_stride;

                for (int mb = 0; mb < ncb; mb++)
                {
                    const auto col = ket_first + mb * ket_stride;

                    const auto value = bvals[(ma * ncb + mb) * npairs + i];

                    dvals[row * naos + col] = value;

                    if ((mtype != mat_t::general) && (row != col))
                    {
                        dvals[col * naos + row] = sign * value;
                    }
                }
            }
        }
    }

    return dmat;
}

}  // namespace sparsematfunc
