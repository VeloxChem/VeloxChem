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

#define PAD_SIZE 8

namespace cpcm {  // cpcm namespace

auto
local_matrix_A_diagonals(const double* ptr_grid_data,
                         const int     row_start,
                         const int     row_end,
                         const int     ncols,
                         const double* ptr_sw_func) -> std::vector<double> 
{
    const double sqrt_2_invpi = std::sqrt(2.0 / mathconst::pi_value());

    std::vector<double> Adiag(row_end - row_start, 0.0);

    auto ptr_Adiag = Adiag.data();

    #pragma omp parallel for schedule(static)
    for (int i = row_start; i < row_end; i++)
    {
        const double xi = ptr_grid_data[i * ncols + 0];
        const double yi = ptr_grid_data[i * ncols + 1];
        const double zi = ptr_grid_data[i * ncols + 2];
        const double wi = ptr_grid_data[i * ncols + 3];

        const double zeta_i = ptr_grid_data[i * ncols + 4];
        const double zeta_i2 = zeta_i * zeta_i;

        ptr_Adiag[i - row_start] = zeta_i * sqrt_2_invpi / ptr_sw_func[i];
    }

    return Adiag;
}

auto
local_matrix_A_dot_vector(const int     npoints,
                          const double* ptr_grid_data,
                          const int     row_start,
                          const int     row_end,
                          const int     ncols,
                          const double* ptr_sw_func,
                          const double* ptr_vector) -> std::vector<double>
{
    const double sqrt_2_invpi = std::sqrt(2.0 / mathconst::pi_value());

    std::vector<double> product(row_end - row_start, 0.0);

    auto ptr_product = product.data();

    #pragma omp parallel for schedule(static, PAD_SIZE)
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
                const auto Aij = zeta_i * sqrt_2_invpi / ptr_sw_func[i];

                ptr_product[i - row_start] += Aij * ptr_vector[j];
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

                const auto Aij = std::erf(zeta_ij * r_ij) / r_ij;

                ptr_product[i - row_start] += Aij * ptr_vector[j];
            }
        }
    }

    return product;
}

auto
comp_grad_Aij(const double* ptr_grid_coords,
              const double* ptr_zeta,
              const int*    ptr_atom_indices,
              const double* ptr_q,
              const int     row_start,
              const int     row_end,
              const int     npoints,
              const int     natoms) -> std::vector<double>
{
    const double two_sqrt_invpi = 2.0 / std::sqrt(mathconst::pi_value());

    auto nthreads = omp_get_max_threads();

    std::vector<double> omp_grad_Aij(nthreads * natoms * 3, 0.0);

    auto ptr_omp_grad_Aij = omp_grad_Aij.data();

    #pragma omp parallel for schedule(static)
    for (int i = row_start; i < row_end; i++)
    {
        const auto thread_id = omp_get_thread_num();

        const auto x_i = ptr_grid_coords[i * 3 + 0];
        const auto y_i = ptr_grid_coords[i * 3 + 1];
        const auto z_i = ptr_grid_coords[i * 3 + 2];

        const auto zeta_i = ptr_zeta[i];
        const auto zeta_i2 = zeta_i * zeta_i;

        const auto atomidx_i = ptr_atom_indices[i];

        for (int j = 0; j < npoints; j++)
        {
            // Note: skip i == j since delta_ij will be zero
            if (i == j) continue;

            const auto x_j = ptr_grid_coords[j * 3 + 0];
            const auto y_j = ptr_grid_coords[j * 3 + 1];
            const auto z_j = ptr_grid_coords[j * 3 + 2];

            const auto zeta_j = ptr_zeta[j];
            const auto zeta_j2 = zeta_j * zeta_j;

            const auto atomidx_j = ptr_atom_indices[j];

            const auto r_ij_2 = ((x_i - x_j) * (x_i - x_j) +
                                 (y_i - y_j) * (y_i - y_j) +
                                 (z_i - z_j) * (z_i - z_j));

            const auto r_ij = std::sqrt(r_ij_2);

            const double uvec_ij[3] = {(x_i - x_j) / r_ij,
                                       (y_i - y_j) / r_ij,
                                       (z_i - z_j) / r_ij};

            const auto zeta_ij = (zeta_i * zeta_j) / std::sqrt(zeta_i2 + zeta_j2);
            const auto erf_term = std::erf(zeta_ij * r_ij);
            const auto exp_term = std::exp(-zeta_ij * zeta_ij * r_ij_2);

            const auto dA_dr = -1.0 * (erf_term - two_sqrt_invpi * zeta_ij * r_ij * exp_term) / r_ij_2;

            for (int a = 0; a < natoms; a++)
            {
                const auto delta_ij = static_cast<double>(static_cast<int>(atomidx_i == a) -
                                                          static_cast<int>(atomidx_j == a));

                for (int c = 0; c < 3; c++)
                {
                    // np.einsum('ij,aij,ijc,i,j->ac', dA_dr, delta_ij, dr_rij, q, q)

                    ptr_omp_grad_Aij[thread_id * natoms * 3 + a * 3 + c] +=

                        ptr_q[i] * delta_ij * dA_dr * uvec_ij[c] * ptr_q[j];
                }
            }
        }
    }

    std::vector<double> grad_Aij(natoms * 3, 0.0);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int a = 0; a < natoms; a++)
        {
            for (int c = 0; c < 3; c++)
            {
                grad_Aij[a * 3 + c] += omp_grad_Aij[thread_id * natoms * 3 + a * 3 + c];
            }
        }
    }

    return grad_Aij;
}

auto
comp_grad_Aii(const double* ptr_grid_coords,
              const double* ptr_zeta,
              const double* ptr_sw_f,
              const int*    ptr_atom_indices,
              const double* ptr_q,
              const double* ptr_atom_coords,
              const double* ptr_atom_radii,
              const int     row_start,
              const int     row_end,
              const int     npoints,
              const int     natoms) -> std::vector<double>
{
    const double sqrt_2_inv_pi = std::sqrt(2.0 / mathconst::pi_value());

    const double sqrt_invpi = 1.0 / std::sqrt(mathconst::pi_value());

    auto nthreads = omp_get_max_threads();

    std::vector<double> omp_grad_Aii(nthreads * natoms * 3, 0.0);

    auto ptr_omp_grad_Aii = omp_grad_Aii.data();

    #pragma omp parallel for schedule(static)
    for (int i = row_start; i < row_end; i++)
    {
        const auto thread_id = omp_get_thread_num();

        const auto x_i = ptr_grid_coords[i * 3 + 0];
        const auto y_i = ptr_grid_coords[i * 3 + 1];
        const auto z_i = ptr_grid_coords[i * 3 + 2];

        const auto zeta_i = ptr_zeta[i];
        const auto zeta_i2 = zeta_i * zeta_i;

        const auto sw_f_i = ptr_sw_f[i];
        const auto q_i = ptr_q[i];

        const auto atomidx_i = ptr_atom_indices[i];

        const auto factor_i = sqrt_2_inv_pi * (zeta_i / sw_f_i) * (q_i * q_i);

        for (int b = 0; b < natoms; b++)
        {
            const auto x_b = ptr_atom_coords[b * 3 + 0];
            const auto y_b = ptr_atom_coords[b * 3 + 1];
            const auto z_b = ptr_atom_coords[b * 3 + 2];

            const auto R_b = ptr_atom_radii[b];

            const auto r_ib_2 = ((x_i - x_b) * (x_i - x_b) +
                                 (y_i - y_b) * (y_i - y_b) +
                                 (z_i - z_b) * (z_i - z_b));

            const auto r_ib = std::sqrt(r_ib_2);

            const double uvec_ib[3] = {(x_i - x_b) / r_ib,
                                       (y_i - y_b) / r_ib,
                                       (z_i - z_b) / r_ib};

            const auto term_m = zeta_i * (R_b - r_ib);
            const auto term_p = zeta_i * (R_b + r_ib);
            const auto fib = 1.0 - 0.5 * (std::erf(term_m) + std::erf(term_p));

            const auto ratio_fib = (sqrt_invpi * zeta_i / fib) * (
                -std::exp(-term_m * term_m) + std::exp(-term_p * term_p));

            for (int a = 0; a < natoms; a++)
            {
                const auto delta = static_cast<double>(static_cast<int>(a == atomidx_i) -
                                                       static_cast<int>(a == b));

                for (int c = 0; c < 3; c++)
                {
                    // np.einsum('aib,ib,ibc,i->ac', delta, ratio_fib, dr_iJ, factor_i)

                    ptr_omp_grad_Aii[thread_id * natoms * 3 + a * 3 + c] +=

                        delta * ratio_fib * uvec_ib[c] * factor_i;
                }
            }
        }
    }

    std::vector<double> grad_Aii(natoms * 3, 0.0);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int a = 0; a < natoms; a++)
        {
            for (int c = 0; c < 3; c++)
            {
                grad_Aii[a * 3 + c] += omp_grad_Aii[thread_id * natoms * 3 + a * 3 + c];
            }
        }
    }

    return grad_Aii;
}

}  // namespace cpcm
