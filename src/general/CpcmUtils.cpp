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

#include "CpcmUtils.hpp"

#include <omp.h>

#include <cmath>
#include <cstring>

#include "MathFunc.hpp"

namespace cpcm {  // cpcm namespace

auto
form_matrix_A(const double* ptr_grid_data,
              const int     row_start,
              const int     row_end,
              const int     ncols,
              const double* ptr_sw_func) -> std::vector<double>
{
    const auto npoints = row_end - row_start;

    std::vector<double> Amat(npoints * npoints);

    const auto ptr_A = Amat.data();

    const double sqrt_2_invpi = std::sqrt(2.0 / mathconst::pi_value());

    #pragma omp parallel for schedule(static)
    for (int i = row_start; i < row_end; i++)
    {
        const double xi = ptr_grid_data[i * ncols + 0];
        const double yi = ptr_grid_data[i * ncols + 1];
        const double zi = ptr_grid_data[i * ncols + 2];
        const double wi = ptr_grid_data[i * ncols + 3];

        const double zeta_i = ptr_grid_data[i * ncols + 4];
        const double zeta_i2 = zeta_i * zeta_i;

        for (int j = 0; j < npoints; j++)
        {
            if (j == i)
            {
                ptr_A[i * npoints + i] = zeta_i * sqrt_2_invpi / ptr_sw_func[i];
            }
            else
            {
                const double xj = ptr_grid_data[j * ncols + 0];
                const double yj = ptr_grid_data[j * ncols + 1];
                const double zj = ptr_grid_data[j * ncols + 2];
                const double wj = ptr_grid_data[j * ncols + 3];

                const double zeta_j = ptr_grid_data[j * ncols + 4];
                const double zeta_j2 = zeta_j * zeta_j;

                const double zeta_ij = zeta_i * zeta_j / std::sqrt(zeta_i2 + zeta_j2);

                const double r_ij = std::sqrt((xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) + (zi - zj)*(zi - zj));

                ptr_A[i * npoints + j] = std::erf(zeta_ij * r_ij) / r_ij;
            }
        }
    }

    return Amat;
}

}  // namespace cpcm
