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

#include "QuadrupoleIntegrals.hpp"

#include <omp.h>

#include <algorithm>
#include <unordered_map>

#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"

#define PAD_SIZE 8

#define MATH_CONST_PI 3.14159265358979323846

namespace onee {  // onee namespace

auto
computeQuadrupoleIntegrals(const CMolecule&           molecule,
                           const CMolecularBasis&     basis,
                           const std::vector<double>& origin) -> std::vector<CDenseMatrix>
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto naos = gtofunc::getNumberOfAtomicOrbitals(gto_blocks);

    std::vector<CDenseMatrix> matrices_mu(6);

    matrices_mu[0] = CDenseMatrix(naos, naos);
    matrices_mu[1] = CDenseMatrix(naos, naos);
    matrices_mu[2] = CDenseMatrix(naos, naos);
    matrices_mu[3] = CDenseMatrix(naos, naos);
    matrices_mu[4] = CDenseMatrix(naos, naos);
    matrices_mu[5] = CDenseMatrix(naos, naos);

    auto& MXX = matrices_mu[0];
    auto& MXY = matrices_mu[1];
    auto& MXZ = matrices_mu[2];
    auto& MYY = matrices_mu[3];
    auto& MYZ = matrices_mu[4];
    auto& MZZ = matrices_mu[5];

    MXX.zero();
    MXY.zero();
    MXZ.zero();
    MYY.zero();
    MYZ.zero();
    MZZ.zero();

    // gto blocks

    int s_prim_count = 0;
    int p_prim_count = 0;
    int d_prim_count = 0;
    int f_prim_count = 0;

    int ncgtos_d = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.number_of_basis_functions();
        const auto npgtos = gto_block.number_of_primitives();

        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 0)
        {
            s_prim_count += npgtos * ncgtos;
        }
        else if (gto_ang == 1)
        {
            p_prim_count += npgtos * ncgtos;
        }
        else if (gto_ang == 2)
        {
            d_prim_count += npgtos * ncgtos;

            ncgtos_d += ncgtos;
        }
        else if (gto_ang == 3)
        {
            f_prim_count += npgtos * ncgtos;
        }
        else
        {
            std::string errangmom("computeQuadrupoleIntegrals: Only implemented up to f-orbitals");

            errors::assertMsgCritical(false, errangmom);
        }
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_p;
    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_d;
    std::unordered_map<int, std::vector<std::pair<int, double>>> cart_sph_f;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 1)
        {
            auto p_map = gto_block.getCartesianToSphericalMappingForP();

            for (const auto& [cart_ind, sph_ind_coef] : p_map)
            {
                cart_sph_p[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 2)
        {
            auto d_map = gto_block.getCartesianToSphericalMappingForD();

            for (const auto& [cart_ind, sph_ind_coef] : d_map)
            {
                cart_sph_d[cart_ind] = sph_ind_coef;
            }
        }
        else if (gto_ang == 3)
        {
            auto f_map = gto_block.getCartesianToSphericalMappingForF(ncgtos_d);

            for (const auto& [cart_ind, sph_ind_coef] : f_map)
            {
                cart_sph_f[cart_ind] = sph_ind_coef;
            }
        }
    }

    // S gto block

    std::vector<double> s_prim_info(5 * s_prim_count);
    std::vector<int>    s_prim_aoinds(1 * s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);

    // P gto block

    std::vector<double> p_prim_info(5 * p_prim_count);
    std::vector<int>    p_prim_aoinds(3 * p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);

    // D gto block

    std::vector<double> d_prim_info(5 * d_prim_count);
    std::vector<int>    d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    // F gto block

    std::vector<double> f_prim_info(5 * f_prim_count);
    std::vector<int>    f_prim_aoinds(10 * f_prim_count);

    gtoinfo::updatePrimitiveInfoForF(f_prim_info.data(), f_prim_aoinds.data(), f_prim_count, gto_blocks, ncgtos_d);

    // primitive pairs

    std::vector<std::tuple<int, int>> pair_inds_ss;
    std::vector<std::tuple<int, int>> pair_inds_sp;
    std::vector<std::tuple<int, int>> pair_inds_sd;
    std::vector<std::tuple<int, int>> pair_inds_sf;
    std::vector<std::tuple<int, int>> pair_inds_pp;
    std::vector<std::tuple<int, int>> pair_inds_pd;
    std::vector<std::tuple<int, int>> pair_inds_pf;
    std::vector<std::tuple<int, int>> pair_inds_dd;
    std::vector<std::tuple<int, int>> pair_inds_df;
    std::vector<std::tuple<int, int>> pair_inds_ff;

    for (int i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int j = i; j < s_prim_count; j++)
        {
            pair_inds_ss.push_back(std::make_tuple(i, j));
        }

        // S-P gto block pair

        for (int j = 0; j < p_prim_count; j++)
        {
            for (int s = 0; s < 3; s++)
            {
                pair_inds_sp.push_back(std::make_tuple(i, j * 3 + s));
            }
        }

        // S-D gto block pair

        for (int j = 0; j < d_prim_count; j++)
        {
            for (int s = 0; s < 6; s++)
            {
                pair_inds_sd.push_back(std::make_tuple(i, j * 6 + s));
            }
        }

        // S-F gto block pair

        for (int j = 0; j < f_prim_count; j++)
        {
            for (int s = 0; s < 10; s++)
            {
                pair_inds_sf.push_back(std::make_tuple(i, j * 10 + s));
            }
        }
    }

    for (int i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int j = i; j < p_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 3; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    pair_inds_pp.push_back(std::make_tuple(i * 3 + i_cart, j * 3 + j_cart));
                }
            }
        }

        // P-D gto block pair

        for (int j = 0; j < d_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int j_cart = 0; j_cart < 6; j_cart++)
                {
                    pair_inds_pd.push_back(std::make_tuple(i * 3 + i_cart, j * 6 + j_cart));
                }
            }
        }

        // P-F gto block pair

        for (int j = 0; j < f_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int j_cart = 0; j_cart < 10; j_cart++)
                {
                    pair_inds_pf.push_back(std::make_tuple(i * 3 + i_cart, j * 10 + j_cart));
                }
            }
        }
    }

    for (int i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int j = i; j < d_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 6; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    pair_inds_dd.push_back(std::make_tuple(i * 6 + i_cart, j * 6 + j_cart));
                }
            }
        }

        // D-F gto block pair

        for (int j = 0; j < f_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 6; i_cart++)
            {
                for (int j_cart = 0; j_cart < 10; j_cart++)
                {
                    pair_inds_df.push_back(std::make_tuple(i * 6 + i_cart, j * 10 + j_cart));
                }
            }
        }
    }

    for (int i = 0; i < f_prim_count; i++)
    {
        // F-F gto block pair

        for (int j = i; j < f_prim_count; j++)
        {
            for (int i_cart = 0; i_cart < 10; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int j_cart = j_cart_start; j_cart < 10; j_cart++)
                {
                    pair_inds_ff.push_back(std::make_tuple(i * 10 + i_cart, j * 10 + j_cart));
                }
            }
        }
    }

    const auto ss_prim_pair_count = static_cast<int>(pair_inds_ss.size());
    const auto sp_prim_pair_count = static_cast<int>(pair_inds_sp.size());
    const auto sd_prim_pair_count = static_cast<int>(pair_inds_sd.size());
    const auto sf_prim_pair_count = static_cast<int>(pair_inds_sf.size());
    const auto pp_prim_pair_count = static_cast<int>(pair_inds_pp.size());
    const auto pd_prim_pair_count = static_cast<int>(pair_inds_pd.size());
    const auto pf_prim_pair_count = static_cast<int>(pair_inds_pf.size());
    const auto dd_prim_pair_count = static_cast<int>(pair_inds_dd.size());
    const auto df_prim_pair_count = static_cast<int>(pair_inds_df.size());
    const auto ff_prim_pair_count = static_cast<int>(pair_inds_ff.size());

    const auto max_prim_pair_count = std::max({
            ss_prim_pair_count,
            sp_prim_pair_count,
            sd_prim_pair_count,
            sf_prim_pair_count,
            pp_prim_pair_count,
            pd_prim_pair_count,
            pf_prim_pair_count,
            dd_prim_pair_count,
            df_prim_pair_count,
            ff_prim_pair_count
    });

    std::vector<double> mat_mu_xx(max_prim_pair_count);
    std::vector<double> mat_mu_xy(max_prim_pair_count);
    std::vector<double> mat_mu_xz(max_prim_pair_count);
    std::vector<double> mat_mu_yy(max_prim_pair_count);
    std::vector<double> mat_mu_yz(max_prim_pair_count);
    std::vector<double> mat_mu_zz(max_prim_pair_count);

    auto nthreads = omp_get_max_threads();

    const double delta[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    const int d_cart_inds[6][2] = {
        {0,0}, {0,1}, {0,2},
        {1,1}, {1,2},
        {2,2}
    };

    const int f_cart_inds[10][3] = {
        {0,0,0}, {0,0,1}, {0,0,2}, {0,1,1}, {0,1,2}, {0,2,2},
        {1,1,1}, {1,1,2}, {1,2,2},
        {2,2,2}
    };

    // auto-generated code begins here

    // S-S block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_ss[ij]);
        const auto j = std::get<1>(pair_inds_ss[ij]);

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = s_prim_info[j + s_prim_count * 0];
        const auto c_j = s_prim_info[j + s_prim_count * 1];
        const auto x_j = s_prim_info[j + s_prim_count * 2];
        const auto y_j = s_prim_info[j + s_prim_count * 3];
        const auto z_j = s_prim_info[j + s_prim_count * 4];




        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};




        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        1.0 * (
                            1.0
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_ss[ij]);
        const auto j = std::get<1>(pair_inds_ss[ij]);

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = s_prim_aoinds[j];

        // Cartesian to spherical
        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            {
                auto j_cgto_sph = j_cgto;
                double j_coef_sph = 1.0;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                if (i != j) MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                if (i != j) MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                if (i != j) MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                if (i != j) MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                if (i != j) MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                if (i != j) MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // S-P block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sp[ij]);
        const auto j = std::get<1>(pair_inds_sp[ij]);

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const auto x_j = p_prim_info[j / 3 + p_prim_count * 2];
        const auto y_j = p_prim_info[j / 3 + p_prim_count * 3];
        const auto z_j = p_prim_info[j / 3 + p_prim_count * 4];


        const auto b0 = j % 3;


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};


        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    PC[m] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][n]
                        )

                    )

                    + PC[n] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][m]
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        (
                            PB_0
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sp[ij]);
        const auto j = std::get<1>(pair_inds_sp[ij]);

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

        // Cartesian to spherical
        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // S-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < sd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sd[ij]);
        const auto j = std::get<1>(pair_inds_sd[ij]);

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];


        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};


        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m])
                        )

                    )

                    + PC[m] * (

                        0.5 / (a_i + a_j) * (
                            delta[b1][n] * (PB_0)
                            + delta[b0][n] * (PB_1)
                        )

                    )

                    + PC[n] * (

                        0.5 / (a_i + a_j) * (
                            delta[b1][m] * (PB_0)
                            + delta[b0][m] * (PB_1)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][b1]
                        )

                        + (
                            PB_0 * PB_1
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < sd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sd[ij]);
        const auto j = std::get<1>(pair_inds_sd[ij]);

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        // Cartesian to spherical
        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // S-F block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < sf_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sf[ij]);
        const auto j = std::get<1>(pair_inds_sf[ij]);

        const auto a_i = s_prim_info[i + s_prim_count * 0];
        const auto c_i = s_prim_info[i + s_prim_count * 1];
        const auto x_i = s_prim_info[i + s_prim_count * 2];
        const auto y_i = s_prim_info[i + s_prim_count * 3];
        const auto z_i = s_prim_info[i + s_prim_count * 4];

        const auto a_j = f_prim_info[j / 10 + f_prim_count * 0];
        const auto c_j = f_prim_info[j / 10 + f_prim_count * 1];
        const auto x_j = f_prim_info[j / 10 + f_prim_count * 2];
        const auto y_j = f_prim_info[j / 10 + f_prim_count * 3];
        const auto z_j = f_prim_info[j / 10 + f_prim_count * 4];


        const auto b0 = f_cart_inds[j % 10][0];
        const auto b1 = f_cart_inds[j % 10][1];
        const auto b2 = f_cart_inds[j % 10][2];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};


        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];
        const auto PB_2 = (-a_i / (a_i + a_j)) * rij[b2];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b1][m] * delta[b2][n] + delta[b1][n] * delta[b2][m]) * (PB_0)
                            + (delta[b0][m] * delta[b2][n] + delta[b0][n] * delta[b2][m]) * (PB_1)
                            + (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PB_2)
                        )

                    )

                    + PC[m] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][n] + delta[b0][b2] * delta[b1][n] + delta[b0][n] * delta[b1][b2])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][n] * (PB_0 * PB_1)
                            + delta[b1][n] * (PB_0 * PB_2)
                            + delta[b0][n] * (PB_1 * PB_2)
                        )

                    )

                    + PC[n] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][m] * (PB_0 * PB_1)
                            + delta[b1][m] * (PB_0 * PB_2)
                            + delta[b0][m] * (PB_1 * PB_2)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.5 / (a_i + a_j) * (
                            delta[b1][b2] * (PB_0)
                            + delta[b0][b2] * (PB_1)
                            + delta[b0][b1] * (PB_2)
                        )

                        + (
                            PB_0 * PB_1 * PB_2
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < sf_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_sf[ij]);
        const auto j = std::get<1>(pair_inds_sf[ij]);

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        // Cartesian to spherical
        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // P-P block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pp[ij]);
        const auto j = std::get<1>(pair_inds_pp[ij]);

        const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
        const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
        const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
        const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
        const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

        const auto a_j = p_prim_info[j / 3 + p_prim_count * 0];
        const auto c_j = p_prim_info[j / 3 + p_prim_count * 1];
        const auto x_j = p_prim_info[j / 3 + p_prim_count * 2];
        const auto y_j = p_prim_info[j / 3 + p_prim_count * 3];
        const auto z_j = p_prim_info[j / 3 + p_prim_count * 4];

        const auto a0 = i % 3;

        const auto b0 = j % 3;


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m])
                        )

                    )

                    + PC[m] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][n] * (PA_0)
                            + delta[a0][n] * (PB_0)
                        )

                    )

                    + PC[n] * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][m] * (PA_0)
                            + delta[a0][m] * (PB_0)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.5 / (a_i + a_j) * (
                            delta[a0][b0]
                        )

                        + (
                            PA_0 * PB_0
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pp[ij]);
        const auto j = std::get<1>(pair_inds_pp[ij]);

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                if (i != j) MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                if (i != j) MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                if (i != j) MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                if (i != j) MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                if (i != j) MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                if (i != j) MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // P-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < pd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pd[ij]);
        const auto j = std::get<1>(pair_inds_pd[ij]);

        const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
        const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
        const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
        const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
        const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];

        const auto a0 = i % 3;

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PA_0)
                            + (delta[a0][m] * delta[b1][n] + delta[a0][n] * delta[b1][m]) * (PB_0)
                            + (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m]) * (PB_1)
                        )

                    )

                    + PC[m] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][b0] * delta[b1][n] + delta[a0][b1] * delta[b0][n] + delta[a0][n] * delta[b0][b1])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][n] * (PA_0 * PB_0)
                            + delta[b0][n] * (PA_0 * PB_1)
                            + delta[a0][n] * (PB_0 * PB_1)
                        )

                    )

                    + PC[n] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][m] * (PA_0 * PB_0)
                            + delta[b0][m] * (PA_0 * PB_1)
                            + delta[a0][m] * (PB_0 * PB_1)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PA_0)
                            + delta[a0][b1] * (PB_0)
                            + delta[a0][b0] * (PB_1)
                        )

                        + (
                            PA_0 * PB_0 * PB_1
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < pd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pd[ij]);
        const auto j = std::get<1>(pair_inds_pd[ij]);

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // P-F block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < pf_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pf[ij]);
        const auto j = std::get<1>(pair_inds_pf[ij]);

        const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
        const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
        const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
        const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
        const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

        const auto a_j = f_prim_info[j / 10 + f_prim_count * 0];
        const auto c_j = f_prim_info[j / 10 + f_prim_count * 1];
        const auto x_j = f_prim_info[j / 10 + f_prim_count * 2];
        const auto y_j = f_prim_info[j / 10 + f_prim_count * 3];
        const auto z_j = f_prim_info[j / 10 + f_prim_count * 4];

        const auto a0 = i % 3;

        const auto b0 = f_cart_inds[j % 10][0];
        const auto b1 = f_cart_inds[j % 10][1];
        const auto b2 = f_cart_inds[j % 10][2];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];
        const auto PB_2 = (-a_i / (a_i + a_j)) * rij[b2];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][b0] * delta[b1][m] * delta[b2][n] + delta[a0][b0] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[b0][m] * delta[b2][n] + delta[a0][b1] * delta[b0][n] * delta[b2][m] + delta[a0][b2] * delta[b0][m] * delta[b1][n] + delta[a0][b2] * delta[b0][n] * delta[b1][m] + delta[a0][m] * delta[b0][b1] * delta[b2][n] + delta[a0][m] * delta[b0][b2] * delta[b1][n] + delta[a0][m] * delta[b0][n] * delta[b1][b2] + delta[a0][n] * delta[b0][b1] * delta[b2][m] + delta[a0][n] * delta[b0][b2] * delta[b1][m] + delta[a0][n] * delta[b0][m] * delta[b1][b2])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b1][m] * delta[b2][n] + delta[b1][n] * delta[b2][m]) * (PA_0 * PB_0)
                            + (delta[b0][m] * delta[b2][n] + delta[b0][n] * delta[b2][m]) * (PA_0 * PB_1)
                            + (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PA_0 * PB_2)
                            + (delta[a0][m] * delta[b2][n] + delta[a0][n] * delta[b2][m]) * (PB_0 * PB_1)
                            + (delta[a0][m] * delta[b1][n] + delta[a0][n] * delta[b1][m]) * (PB_0 * PB_2)
                            + (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m]) * (PB_1 * PB_2)
                        )

                    )

                    + PC[m] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][n] + delta[b0][b2] * delta[b1][n] + delta[b0][n] * delta[b1][b2]) * (PA_0)
                            + (delta[a0][b1] * delta[b2][n] + delta[a0][b2] * delta[b1][n] + delta[a0][n] * delta[b1][b2]) * (PB_0)
                            + (delta[a0][b0] * delta[b2][n] + delta[a0][b2] * delta[b0][n] + delta[a0][n] * delta[b0][b2]) * (PB_1)
                            + (delta[a0][b0] * delta[b1][n] + delta[a0][b1] * delta[b0][n] + delta[a0][n] * delta[b0][b1]) * (PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][n] * (PA_0 * PB_0 * PB_1)
                            + delta[b1][n] * (PA_0 * PB_0 * PB_2)
                            + delta[b0][n] * (PA_0 * PB_1 * PB_2)
                            + delta[a0][n] * (PB_0 * PB_1 * PB_2)
                        )

                    )

                    + PC[n] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][m] * (PA_0 * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PB_0 * PB_2)
                            + delta[b0][m] * (PA_0 * PB_1 * PB_2)
                            + delta[a0][m] * (PB_0 * PB_1 * PB_2)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][b2] * (PA_0 * PB_0)
                            + delta[b0][b2] * (PA_0 * PB_1)
                            + delta[b0][b1] * (PA_0 * PB_2)
                            + delta[a0][b2] * (PB_0 * PB_1)
                            + delta[a0][b1] * (PB_0 * PB_2)
                            + delta[a0][b0] * (PB_1 * PB_2)
                        )

                        + (
                            PA_0 * PB_0 * PB_1 * PB_2
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < pf_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_pf[ij]);
        const auto j = std::get<1>(pair_inds_pf[ij]);

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // D-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < dd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_dd[ij]);
        const auto j = std::get<1>(pair_inds_dd[ij]);

        const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
        const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
        const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
        const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
        const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

        const auto a_j = d_prim_info[j / 6 + d_prim_count * 0];
        const auto c_j = d_prim_info[j / 6 + d_prim_count * 1];
        const auto x_j = d_prim_info[j / 6 + d_prim_count * 2];
        const auto y_j = d_prim_info[j / 6 + d_prim_count * 3];
        const auto z_j = d_prim_info[j / 6 + d_prim_count * 4];

        const auto a0 = d_cart_inds[i % 6][0];
        const auto a1 = d_cart_inds[i % 6][1];

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[b0][m] * delta[b1][n] + delta[a0][a1] * delta[b0][n] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][n] + delta[a0][b0] * delta[a1][n] * delta[b1][m] + delta[a0][b1] * delta[a1][m] * delta[b0][n] + delta[a0][b1] * delta[a1][n] * delta[b0][m] + delta[a0][m] * delta[a1][b0] * delta[b1][n] + delta[a0][m] * delta[a1][b1] * delta[b0][n] + delta[a0][m] * delta[a1][n] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[b1][m] + delta[a0][n] * delta[a1][b1] * delta[b0][m] + delta[a0][n] * delta[a1][m] * delta[b0][b1])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PA_0 * PA_1)
                            + (delta[a1][m] * delta[b1][n] + delta[a1][n] * delta[b1][m]) * (PA_0 * PB_0)
                            + (delta[a1][m] * delta[b0][n] + delta[a1][n] * delta[b0][m]) * (PA_0 * PB_1)
                            + (delta[a0][m] * delta[b1][n] + delta[a0][n] * delta[b1][m]) * (PA_1 * PB_0)
                            + (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m]) * (PA_1 * PB_1)
                            + (delta[a0][m] * delta[a1][n] + delta[a0][n] * delta[a1][m]) * (PB_0 * PB_1)
                        )

                    )

                    + PC[m] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][b0] * delta[b1][n] + delta[a1][b1] * delta[b0][n] + delta[a1][n] * delta[b0][b1]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][n] + delta[a0][b1] * delta[b0][n] + delta[a0][n] * delta[b0][b1]) * (PA_1)
                            + (delta[a0][a1] * delta[b1][n] + delta[a0][b1] * delta[a1][n] + delta[a0][n] * delta[a1][b1]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][n] + delta[a0][b0] * delta[a1][n] + delta[a0][n] * delta[a1][b0]) * (PB_1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][n] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][n] * (PA_0 * PA_1 * PB_1)
                            + delta[a1][n] * (PA_0 * PB_0 * PB_1)
                            + delta[a0][n] * (PA_1 * PB_0 * PB_1)
                        )

                    )

                    + PC[n] * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][m] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m] * (PA_0 * PA_1 * PB_1)
                            + delta[a1][m] * (PA_0 * PB_0 * PB_1)
                            + delta[a0][m] * (PA_1 * PB_0 * PB_1)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0])
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b0][b1] * (PA_0 * PA_1)
                            + delta[a1][b1] * (PA_0 * PB_0)
                            + delta[a1][b0] * (PA_0 * PB_1)
                            + delta[a0][b1] * (PA_1 * PB_0)
                            + delta[a0][b0] * (PA_1 * PB_1)
                            + delta[a0][a1] * (PB_0 * PB_1)
                        )

                        + (
                            PA_0 * PA_1 * PB_0 * PB_1
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < dd_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_dd[ij]);
        const auto j = std::get<1>(pair_inds_dd[ij]);

        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                if (i != j) MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                if (i != j) MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                if (i != j) MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                if (i != j) MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                if (i != j) MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                if (i != j) MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // D-F block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < df_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_df[ij]);
        const auto j = std::get<1>(pair_inds_df[ij]);

        const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
        const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
        const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
        const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
        const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

        const auto a_j = f_prim_info[j / 10 + f_prim_count * 0];
        const auto c_j = f_prim_info[j / 10 + f_prim_count * 1];
        const auto x_j = f_prim_info[j / 10 + f_prim_count * 2];
        const auto y_j = f_prim_info[j / 10 + f_prim_count * 3];
        const auto z_j = f_prim_info[j / 10 + f_prim_count * 4];

        const auto a0 = d_cart_inds[i % 6][0];
        const auto a1 = d_cart_inds[i % 6][1];

        const auto b0 = f_cart_inds[j % 10][0];
        const auto b1 = f_cart_inds[j % 10][1];
        const auto b2 = f_cart_inds[j % 10][2];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];
        const auto PB_2 = (-a_i / (a_i + a_j)) * rij[b2];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][b0] * delta[b1][m] * delta[b2][n] + delta[a1][b0] * delta[b1][n] * delta[b2][m] + delta[a1][b1] * delta[b0][m] * delta[b2][n] + delta[a1][b1] * delta[b0][n] * delta[b2][m] + delta[a1][b2] * delta[b0][m] * delta[b1][n] + delta[a1][b2] * delta[b0][n] * delta[b1][m] + delta[a1][m] * delta[b0][b1] * delta[b2][n] + delta[a1][m] * delta[b0][b2] * delta[b1][n] + delta[a1][m] * delta[b0][n] * delta[b1][b2] + delta[a1][n] * delta[b0][b1] * delta[b2][m] + delta[a1][n] * delta[b0][b2] * delta[b1][m] + delta[a1][n] * delta[b0][m] * delta[b1][b2]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][m] * delta[b2][n] + delta[a0][b0] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[b0][m] * delta[b2][n] + delta[a0][b1] * delta[b0][n] * delta[b2][m] + delta[a0][b2] * delta[b0][m] * delta[b1][n] + delta[a0][b2] * delta[b0][n] * delta[b1][m] + delta[a0][m] * delta[b0][b1] * delta[b2][n] + delta[a0][m] * delta[b0][b2] * delta[b1][n] + delta[a0][m] * delta[b0][n] * delta[b1][b2] + delta[a0][n] * delta[b0][b1] * delta[b2][m] + delta[a0][n] * delta[b0][b2] * delta[b1][m] + delta[a0][n] * delta[b0][m] * delta[b1][b2]) * (PA_1)
                            + (delta[a0][a1] * delta[b1][m] * delta[b2][n] + delta[a0][a1] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[a1][m] * delta[b2][n] + delta[a0][b1] * delta[a1][n] * delta[b2][m] + delta[a0][b2] * delta[a1][m] * delta[b1][n] + delta[a0][b2] * delta[a1][n] * delta[b1][m] + delta[a0][m] * delta[a1][b1] * delta[b2][n] + delta[a0][m] * delta[a1][b2] * delta[b1][n] + delta[a0][m] * delta[a1][n] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[b2][m] + delta[a0][n] * delta[a1][b2] * delta[b1][m] + delta[a0][n] * delta[a1][m] * delta[b1][b2]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][m] * delta[b2][n] + delta[a0][a1] * delta[b0][n] * delta[b2][m] + delta[a0][b0] * delta[a1][m] * delta[b2][n] + delta[a0][b0] * delta[a1][n] * delta[b2][m] + delta[a0][b2] * delta[a1][m] * delta[b0][n] + delta[a0][b2] * delta[a1][n] * delta[b0][m] + delta[a0][m] * delta[a1][b0] * delta[b2][n] + delta[a0][m] * delta[a1][b2] * delta[b0][n] + delta[a0][m] * delta[a1][n] * delta[b0][b2] + delta[a0][n] * delta[a1][b0] * delta[b2][m] + delta[a0][n] * delta[a1][b2] * delta[b0][m] + delta[a0][n] * delta[a1][m] * delta[b0][b2]) * (PB_1)
                            + (delta[a0][a1] * delta[b0][m] * delta[b1][n] + delta[a0][a1] * delta[b0][n] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][n] + delta[a0][b0] * delta[a1][n] * delta[b1][m] + delta[a0][b1] * delta[a1][m] * delta[b0][n] + delta[a0][b1] * delta[a1][n] * delta[b0][m] + delta[a0][m] * delta[a1][b0] * delta[b1][n] + delta[a0][m] * delta[a1][b1] * delta[b0][n] + delta[a0][m] * delta[a1][n] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[b1][m] + delta[a0][n] * delta[a1][b1] * delta[b0][m] + delta[a0][n] * delta[a1][m] * delta[b0][b1]) * (PB_2)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b1][m] * delta[b2][n] + delta[b1][n] * delta[b2][m]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][m] * delta[b2][n] + delta[b0][n] * delta[b2][m]) * (PA_0 * PA_1 * PB_1)
                            + (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PA_0 * PA_1 * PB_2)
                            + (delta[a1][m] * delta[b2][n] + delta[a1][n] * delta[b2][m]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][m] * delta[b1][n] + delta[a1][n] * delta[b1][m]) * (PA_0 * PB_0 * PB_2)
                            + (delta[a1][m] * delta[b0][n] + delta[a1][n] * delta[b0][m]) * (PA_0 * PB_1 * PB_2)
                            + (delta[a0][m] * delta[b2][n] + delta[a0][n] * delta[b2][m]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][m] * delta[b1][n] + delta[a0][n] * delta[b1][m]) * (PA_1 * PB_0 * PB_2)
                            + (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m]) * (PA_1 * PB_1 * PB_2)
                            + (delta[a0][m] * delta[a1][n] + delta[a0][n] * delta[a1][m]) * (PB_0 * PB_1 * PB_2)
                        )

                    )

                    + PC[m] * (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[b0][b1] * delta[b2][n] + delta[a0][a1] * delta[b0][b2] * delta[b1][n] + delta[a0][a1] * delta[b0][n] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][n] + delta[a0][b0] * delta[a1][b2] * delta[b1][n] + delta[a0][b0] * delta[a1][n] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][n] + delta[a0][b1] * delta[a1][b2] * delta[b0][n] + delta[a0][b1] * delta[a1][n] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][n] + delta[a0][b2] * delta[a1][b1] * delta[b0][n] + delta[a0][b2] * delta[a1][n] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[b0][b2] + delta[a0][n] * delta[a1][b2] * delta[b0][b1])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][n] + delta[b0][b2] * delta[b1][n] + delta[b0][n] * delta[b1][b2]) * (PA_0 * PA_1)
                            + (delta[a1][b1] * delta[b2][n] + delta[a1][b2] * delta[b1][n] + delta[a1][n] * delta[b1][b2]) * (PA_0 * PB_0)
                            + (delta[a1][b0] * delta[b2][n] + delta[a1][b2] * delta[b0][n] + delta[a1][n] * delta[b0][b2]) * (PA_0 * PB_1)
                            + (delta[a1][b0] * delta[b1][n] + delta[a1][b1] * delta[b0][n] + delta[a1][n] * delta[b0][b1]) * (PA_0 * PB_2)
                            + (delta[a0][b1] * delta[b2][n] + delta[a0][b2] * delta[b1][n] + delta[a0][n] * delta[b1][b2]) * (PA_1 * PB_0)
                            + (delta[a0][b0] * delta[b2][n] + delta[a0][b2] * delta[b0][n] + delta[a0][n] * delta[b0][b2]) * (PA_1 * PB_1)
                            + (delta[a0][b0] * delta[b1][n] + delta[a0][b1] * delta[b0][n] + delta[a0][n] * delta[b0][b1]) * (PA_1 * PB_2)
                            + (delta[a0][a1] * delta[b2][n] + delta[a0][b2] * delta[a1][n] + delta[a0][n] * delta[a1][b2]) * (PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][n] + delta[a0][b1] * delta[a1][n] + delta[a0][n] * delta[a1][b1]) * (PB_0 * PB_2)
                            + (delta[a0][a1] * delta[b0][n] + delta[a0][b0] * delta[a1][n] + delta[a0][n] * delta[a1][b0]) * (PB_1 * PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][n] * (PA_0 * PA_1 * PB_0 * PB_1)
                            + delta[b1][n] * (PA_0 * PA_1 * PB_0 * PB_2)
                            + delta[b0][n] * (PA_0 * PA_1 * PB_1 * PB_2)
                            + delta[a1][n] * (PA_0 * PB_0 * PB_1 * PB_2)
                            + delta[a0][n] * (PA_1 * PB_0 * PB_1 * PB_2)
                        )

                    )

                    + PC[n] * (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1)
                            + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PB_0)
                            + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_2)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_2)
                            + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_2)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][m] * (PA_0 * PA_1 * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PB_2)
                            + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PB_2)
                            + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PB_2)
                            + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PB_2)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1)
                            + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PB_1)
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][b2] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][b2] * (PA_0 * PA_1 * PB_1)
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_2)
                            + delta[a1][b2] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_2)
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_2)
                            + delta[a0][b2] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_2)
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_2)
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_2)
                        )

                        + (
                            PA_0 * PA_1 * PB_0 * PB_1 * PB_2
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < df_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_df[ij]);
        const auto j = std::get<1>(pair_inds_df[ij]);

        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }


    // F-F block

    #pragma omp parallel for schedule(static, PAD_SIZE)
    for (int ij = 0; ij < ff_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_ff[ij]);
        const auto j = std::get<1>(pair_inds_ff[ij]);

        const auto a_i = f_prim_info[i / 10 + f_prim_count * 0];
        const auto c_i = f_prim_info[i / 10 + f_prim_count * 1];
        const auto x_i = f_prim_info[i / 10 + f_prim_count * 2];
        const auto y_i = f_prim_info[i / 10 + f_prim_count * 3];
        const auto z_i = f_prim_info[i / 10 + f_prim_count * 4];

        const auto a_j = f_prim_info[j / 10 + f_prim_count * 0];
        const auto c_j = f_prim_info[j / 10 + f_prim_count * 1];
        const auto x_j = f_prim_info[j / 10 + f_prim_count * 2];
        const auto y_j = f_prim_info[j / 10 + f_prim_count * 3];
        const auto z_j = f_prim_info[j / 10 + f_prim_count * 4];

        const auto a0 = f_cart_inds[i % 10][0];
        const auto a1 = f_cart_inds[i % 10][1];
        const auto a2 = f_cart_inds[i % 10][2];

        const auto b0 = f_cart_inds[j % 10][0];
        const auto b1 = f_cart_inds[j % 10][1];
        const auto b2 = f_cart_inds[j % 10][2];


        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - origin[0],
                              (a_i * y_i + a_j * y_j) / (a_i + a_j) - origin[1],
                              (a_i * z_i + a_j * z_j) / (a_i + a_j) - origin[2]};

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];
        const auto PA_2 = (a_j / (a_i + a_j)) * rij[a2];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];
        const auto PB_2 = (-a_i / (a_i + a_j)) * rij[b2];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double* mat_mu[6] = {mat_mu_xx.data(), mat_mu_xy.data(), mat_mu_xz.data(),
                             mat_mu_yy.data(), mat_mu_yz.data(), mat_mu_zz.data()};

        std::vector<std::vector<int>> m_n_inds{{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};

        for (int idx = 0; idx < 6; idx++)
        {
            auto m = m_n_inds[idx][0];
            auto n = m_n_inds[idx][1];

            // Note: minus sign from electron charge

            mat_mu[idx][ij] = (-1.0) * S_ij_00 * (

                    (

                        0.0625 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[a2][b0] * delta[b1][m] * delta[b2][n] + delta[a0][a1] * delta[a2][b0] * delta[b1][n] * delta[b2][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] * delta[b2][n] + delta[a0][a1] * delta[a2][b1] * delta[b0][n] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] * delta[b1][n] + delta[a0][a1] * delta[a2][b2] * delta[b0][n] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] * delta[b2][n] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] * delta[b1][n] + delta[a0][a1] * delta[a2][m] * delta[b0][n] * delta[b1][b2] + delta[a0][a1] * delta[a2][n] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][n] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][n] * delta[b0][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] * delta[b2][n] + delta[a0][a2] * delta[a1][b0] * delta[b1][n] * delta[b2][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] * delta[b2][n] + delta[a0][a2] * delta[a1][b1] * delta[b0][n] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] * delta[b1][n] + delta[a0][a2] * delta[a1][b2] * delta[b0][n] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] * delta[b2][n] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] * delta[b1][n] + delta[a0][a2] * delta[a1][m] * delta[b0][n] * delta[b1][b2] + delta[a0][a2] * delta[a1][n] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][n] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][n] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] * delta[b2][n] + delta[a0][b0] * delta[a1][a2] * delta[b1][n] * delta[b2][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] * delta[b2][n] + delta[a0][b0] * delta[a1][b1] * delta[a2][n] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] * delta[b1][n] + delta[a0][b0] * delta[a1][b2] * delta[a2][n] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] * delta[b2][n] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] * delta[b1][n] + delta[a0][b0] * delta[a1][m] * delta[a2][n] * delta[b1][b2] + delta[a0][b0] * delta[a1][n] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][n] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][n] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] * delta[b2][n] + delta[a0][b1] * delta[a1][a2] * delta[b0][n] * delta[b2][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] * delta[b2][n] + delta[a0][b1] * delta[a1][b0] * delta[a2][n] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] * delta[b0][n] + delta[a0][b1] * delta[a1][b2] * delta[a2][n] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] * delta[b2][n] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] * delta[b0][n] + delta[a0][b1] * delta[a1][m] * delta[a2][n] * delta[b0][b2] + delta[a0][b1] * delta[a1][n] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][n] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][n] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] * delta[b1][n] + delta[a0][b2] * delta[a1][a2] * delta[b0][n] * delta[b1][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] * delta[b1][n] + delta[a0][b2] * delta[a1][b0] * delta[a2][n] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] * delta[b0][n] + delta[a0][b2] * delta[a1][b1] * delta[a2][n] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] * delta[b1][n] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] * delta[b0][n] + delta[a0][b2] * delta[a1][m] * delta[a2][n] * delta[b0][b1] + delta[a0][b2] * delta[a1][n] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][n] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][n] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] * delta[b2][n] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] * delta[b1][n] + delta[a0][m] * delta[a1][a2] * delta[b0][n] * delta[b1][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] * delta[b2][n] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] * delta[b1][n] + delta[a0][m] * delta[a1][b0] * delta[a2][n] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b0] * delta[b2][n] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] * delta[b0][n] + delta[a0][m] * delta[a1][b1] * delta[a2][n] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0] * delta[b1][n] + delta[a0][m] * delta[a1][b2] * delta[a2][b1] * delta[b0][n] + delta[a0][m] * delta[a1][b2] * delta[a2][n] * delta[b0][b1] + delta[a0][m] * delta[a1][n] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][n] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][n] * delta[a2][b2] * delta[b0][b1] + delta[a0][n] * delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][n] * delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][n] * delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][n] * delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][n] * delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][n] * delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][n] * delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][n] * delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][n] * delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][n] * delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][n] * delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][n] * delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][n] * delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][n] * delta[a1][m] * delta[a2][b2] * delta[b0][b1])
                        )

                        + 0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a2][b0] * delta[b1][m] * delta[b2][n] + delta[a2][b0] * delta[b1][n] * delta[b2][m] + delta[a2][b1] * delta[b0][m] * delta[b2][n] + delta[a2][b1] * delta[b0][n] * delta[b2][m] + delta[a2][b2] * delta[b0][m] * delta[b1][n] + delta[a2][b2] * delta[b0][n] * delta[b1][m] + delta[a2][m] * delta[b0][b1] * delta[b2][n] + delta[a2][m] * delta[b0][b2] * delta[b1][n] + delta[a2][m] * delta[b0][n] * delta[b1][b2] + delta[a2][n] * delta[b0][b1] * delta[b2][m] + delta[a2][n] * delta[b0][b2] * delta[b1][m] + delta[a2][n] * delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1)
                            + (delta[a1][b0] * delta[b1][m] * delta[b2][n] + delta[a1][b0] * delta[b1][n] * delta[b2][m] + delta[a1][b1] * delta[b0][m] * delta[b2][n] + delta[a1][b1] * delta[b0][n] * delta[b2][m] + delta[a1][b2] * delta[b0][m] * delta[b1][n] + delta[a1][b2] * delta[b0][n] * delta[b1][m] + delta[a1][m] * delta[b0][b1] * delta[b2][n] + delta[a1][m] * delta[b0][b2] * delta[b1][n] + delta[a1][m] * delta[b0][n] * delta[b1][b2] + delta[a1][n] * delta[b0][b1] * delta[b2][m] + delta[a1][n] * delta[b0][b2] * delta[b1][m] + delta[a1][n] * delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_2)
                            + (delta[a1][a2] * delta[b1][m] * delta[b2][n] + delta[a1][a2] * delta[b1][n] * delta[b2][m] + delta[a1][b1] * delta[a2][m] * delta[b2][n] + delta[a1][b1] * delta[a2][n] * delta[b2][m] + delta[a1][b2] * delta[a2][m] * delta[b1][n] + delta[a1][b2] * delta[a2][n] * delta[b1][m] + delta[a1][m] * delta[a2][b1] * delta[b2][n] + delta[a1][m] * delta[a2][b2] * delta[b1][n] + delta[a1][m] * delta[a2][n] * delta[b1][b2] + delta[a1][n] * delta[a2][b1] * delta[b2][m] + delta[a1][n] * delta[a2][b2] * delta[b1][m] + delta[a1][n] * delta[a2][m] * delta[b1][b2]) * (PA_0 * PB_0)
                            + (delta[a1][a2] * delta[b0][m] * delta[b2][n] + delta[a1][a2] * delta[b0][n] * delta[b2][m] + delta[a1][b0] * delta[a2][m] * delta[b2][n] + delta[a1][b0] * delta[a2][n] * delta[b2][m] + delta[a1][b2] * delta[a2][m] * delta[b0][n] + delta[a1][b2] * delta[a2][n] * delta[b0][m] + delta[a1][m] * delta[a2][b0] * delta[b2][n] + delta[a1][m] * delta[a2][b2] * delta[b0][n] + delta[a1][m] * delta[a2][n] * delta[b0][b2] + delta[a1][n] * delta[a2][b0] * delta[b2][m] + delta[a1][n] * delta[a2][b2] * delta[b0][m] + delta[a1][n] * delta[a2][m] * delta[b0][b2]) * (PA_0 * PB_1)
                            + (delta[a1][a2] * delta[b0][m] * delta[b1][n] + delta[a1][a2] * delta[b0][n] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][n] + delta[a1][b0] * delta[a2][n] * delta[b1][m] + delta[a1][b1] * delta[a2][m] * delta[b0][n] + delta[a1][b1] * delta[a2][n] * delta[b0][m] + delta[a1][m] * delta[a2][b0] * delta[b1][n] + delta[a1][m] * delta[a2][b1] * delta[b0][n] + delta[a1][m] * delta[a2][n] * delta[b0][b1] + delta[a1][n] * delta[a2][b0] * delta[b1][m] + delta[a1][n] * delta[a2][b1] * delta[b0][m] + delta[a1][n] * delta[a2][m] * delta[b0][b1]) * (PA_0 * PB_2)
                            + (delta[a0][b0] * delta[b1][m] * delta[b2][n] + delta[a0][b0] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[b0][m] * delta[b2][n] + delta[a0][b1] * delta[b0][n] * delta[b2][m] + delta[a0][b2] * delta[b0][m] * delta[b1][n] + delta[a0][b2] * delta[b0][n] * delta[b1][m] + delta[a0][m] * delta[b0][b1] * delta[b2][n] + delta[a0][m] * delta[b0][b2] * delta[b1][n] + delta[a0][m] * delta[b0][n] * delta[b1][b2] + delta[a0][n] * delta[b0][b1] * delta[b2][m] + delta[a0][n] * delta[b0][b2] * delta[b1][m] + delta[a0][n] * delta[b0][m] * delta[b1][b2]) * (PA_1 * PA_2)
                            + (delta[a0][a2] * delta[b1][m] * delta[b2][n] + delta[a0][a2] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[a2][m] * delta[b2][n] + delta[a0][b1] * delta[a2][n] * delta[b2][m] + delta[a0][b2] * delta[a2][m] * delta[b1][n] + delta[a0][b2] * delta[a2][n] * delta[b1][m] + delta[a0][m] * delta[a2][b1] * delta[b2][n] + delta[a0][m] * delta[a2][b2] * delta[b1][n] + delta[a0][m] * delta[a2][n] * delta[b1][b2] + delta[a0][n] * delta[a2][b1] * delta[b2][m] + delta[a0][n] * delta[a2][b2] * delta[b1][m] + delta[a0][n] * delta[a2][m] * delta[b1][b2]) * (PA_1 * PB_0)
                            + (delta[a0][a2] * delta[b0][m] * delta[b2][n] + delta[a0][a2] * delta[b0][n] * delta[b2][m] + delta[a0][b0] * delta[a2][m] * delta[b2][n] + delta[a0][b0] * delta[a2][n] * delta[b2][m] + delta[a0][b2] * delta[a2][m] * delta[b0][n] + delta[a0][b2] * delta[a2][n] * delta[b0][m] + delta[a0][m] * delta[a2][b0] * delta[b2][n] + delta[a0][m] * delta[a2][b2] * delta[b0][n] + delta[a0][m] * delta[a2][n] * delta[b0][b2] + delta[a0][n] * delta[a2][b0] * delta[b2][m] + delta[a0][n] * delta[a2][b2] * delta[b0][m] + delta[a0][n] * delta[a2][m] * delta[b0][b2]) * (PA_1 * PB_1)
                            + (delta[a0][a2] * delta[b0][m] * delta[b1][n] + delta[a0][a2] * delta[b0][n] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][n] + delta[a0][b0] * delta[a2][n] * delta[b1][m] + delta[a0][b1] * delta[a2][m] * delta[b0][n] + delta[a0][b1] * delta[a2][n] * delta[b0][m] + delta[a0][m] * delta[a2][b0] * delta[b1][n] + delta[a0][m] * delta[a2][b1] * delta[b0][n] + delta[a0][m] * delta[a2][n] * delta[b0][b1] + delta[a0][n] * delta[a2][b0] * delta[b1][m] + delta[a0][n] * delta[a2][b1] * delta[b0][m] + delta[a0][n] * delta[a2][m] * delta[b0][b1]) * (PA_1 * PB_2)
                            + (delta[a0][a1] * delta[b1][m] * delta[b2][n] + delta[a0][a1] * delta[b1][n] * delta[b2][m] + delta[a0][b1] * delta[a1][m] * delta[b2][n] + delta[a0][b1] * delta[a1][n] * delta[b2][m] + delta[a0][b2] * delta[a1][m] * delta[b1][n] + delta[a0][b2] * delta[a1][n] * delta[b1][m] + delta[a0][m] * delta[a1][b1] * delta[b2][n] + delta[a0][m] * delta[a1][b2] * delta[b1][n] + delta[a0][m] * delta[a1][n] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[b2][m] + delta[a0][n] * delta[a1][b2] * delta[b1][m] + delta[a0][n] * delta[a1][m] * delta[b1][b2]) * (PA_2 * PB_0)
                            + (delta[a0][a1] * delta[b0][m] * delta[b2][n] + delta[a0][a1] * delta[b0][n] * delta[b2][m] + delta[a0][b0] * delta[a1][m] * delta[b2][n] + delta[a0][b0] * delta[a1][n] * delta[b2][m] + delta[a0][b2] * delta[a1][m] * delta[b0][n] + delta[a0][b2] * delta[a1][n] * delta[b0][m] + delta[a0][m] * delta[a1][b0] * delta[b2][n] + delta[a0][m] * delta[a1][b2] * delta[b0][n] + delta[a0][m] * delta[a1][n] * delta[b0][b2] + delta[a0][n] * delta[a1][b0] * delta[b2][m] + delta[a0][n] * delta[a1][b2] * delta[b0][m] + delta[a0][n] * delta[a1][m] * delta[b0][b2]) * (PA_2 * PB_1)
                            + (delta[a0][a1] * delta[b0][m] * delta[b1][n] + delta[a0][a1] * delta[b0][n] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][n] + delta[a0][b0] * delta[a1][n] * delta[b1][m] + delta[a0][b1] * delta[a1][m] * delta[b0][n] + delta[a0][b1] * delta[a1][n] * delta[b0][m] + delta[a0][m] * delta[a1][b0] * delta[b1][n] + delta[a0][m] * delta[a1][b1] * delta[b0][n] + delta[a0][m] * delta[a1][n] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[b1][m] + delta[a0][n] * delta[a1][b1] * delta[b0][m] + delta[a0][n] * delta[a1][m] * delta[b0][b1]) * (PA_2 * PB_2)
                            + (delta[a0][a1] * delta[a2][m] * delta[b2][n] + delta[a0][a1] * delta[a2][n] * delta[b2][m] + delta[a0][a2] * delta[a1][m] * delta[b2][n] + delta[a0][a2] * delta[a1][n] * delta[b2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][n] + delta[a0][b2] * delta[a1][n] * delta[a2][m] + delta[a0][m] * delta[a1][a2] * delta[b2][n] + delta[a0][m] * delta[a1][b2] * delta[a2][n] + delta[a0][m] * delta[a1][n] * delta[a2][b2] + delta[a0][n] * delta[a1][a2] * delta[b2][m] + delta[a0][n] * delta[a1][b2] * delta[a2][m] + delta[a0][n] * delta[a1][m] * delta[a2][b2]) * (PB_0 * PB_1)
                            + (delta[a0][a1] * delta[a2][m] * delta[b1][n] + delta[a0][a1] * delta[a2][n] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][n] + delta[a0][a2] * delta[a1][n] * delta[b1][m] + delta[a0][b1] * delta[a1][m] * delta[a2][n] + delta[a0][b1] * delta[a1][n] * delta[a2][m] + delta[a0][m] * delta[a1][a2] * delta[b1][n] + delta[a0][m] * delta[a1][b1] * delta[a2][n] + delta[a0][m] * delta[a1][n] * delta[a2][b1] + delta[a0][n] * delta[a1][a2] * delta[b1][m] + delta[a0][n] * delta[a1][b1] * delta[a2][m] + delta[a0][n] * delta[a1][m] * delta[a2][b1]) * (PB_0 * PB_2)
                            + (delta[a0][a1] * delta[a2][m] * delta[b0][n] + delta[a0][a1] * delta[a2][n] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][n] + delta[a0][a2] * delta[a1][n] * delta[b0][m] + delta[a0][b0] * delta[a1][m] * delta[a2][n] + delta[a0][b0] * delta[a1][n] * delta[a2][m] + delta[a0][m] * delta[a1][a2] * delta[b0][n] + delta[a0][m] * delta[a1][b0] * delta[a2][n] + delta[a0][m] * delta[a1][n] * delta[a2][b0] + delta[a0][n] * delta[a1][a2] * delta[b0][m] + delta[a0][n] * delta[a1][b0] * delta[a2][m] + delta[a0][n] * delta[a1][m] * delta[a2][b0]) * (PB_1 * PB_2)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b1][m] * delta[b2][n] + delta[b1][n] * delta[b2][m]) * (PA_0 * PA_1 * PA_2 * PB_0)
                            + (delta[b0][m] * delta[b2][n] + delta[b0][n] * delta[b2][m]) * (PA_0 * PA_1 * PA_2 * PB_1)
                            + (delta[b0][m] * delta[b1][n] + delta[b0][n] * delta[b1][m]) * (PA_0 * PA_1 * PA_2 * PB_2)
                            + (delta[a2][m] * delta[b2][n] + delta[a2][n] * delta[b2][m]) * (PA_0 * PA_1 * PB_0 * PB_1)
                            + (delta[a2][m] * delta[b1][n] + delta[a2][n] * delta[b1][m]) * (PA_0 * PA_1 * PB_0 * PB_2)
                            + (delta[a2][m] * delta[b0][n] + delta[a2][n] * delta[b0][m]) * (PA_0 * PA_1 * PB_1 * PB_2)
                            + (delta[a1][m] * delta[b2][n] + delta[a1][n] * delta[b2][m]) * (PA_0 * PA_2 * PB_0 * PB_1)
                            + (delta[a1][m] * delta[b1][n] + delta[a1][n] * delta[b1][m]) * (PA_0 * PA_2 * PB_0 * PB_2)
                            + (delta[a1][m] * delta[b0][n] + delta[a1][n] * delta[b0][m]) * (PA_0 * PA_2 * PB_1 * PB_2)
                            + (delta[a1][m] * delta[a2][n] + delta[a1][n] * delta[a2][m]) * (PA_0 * PB_0 * PB_1 * PB_2)
                            + (delta[a0][m] * delta[b2][n] + delta[a0][n] * delta[b2][m]) * (PA_1 * PA_2 * PB_0 * PB_1)
                            + (delta[a0][m] * delta[b1][n] + delta[a0][n] * delta[b1][m]) * (PA_1 * PA_2 * PB_0 * PB_2)
                            + (delta[a0][m] * delta[b0][n] + delta[a0][n] * delta[b0][m]) * (PA_1 * PA_2 * PB_1 * PB_2)
                            + (delta[a0][m] * delta[a2][n] + delta[a0][n] * delta[a2][m]) * (PA_1 * PB_0 * PB_1 * PB_2)
                            + (delta[a0][m] * delta[a1][n] + delta[a0][n] * delta[a1][m]) * (PA_2 * PB_0 * PB_1 * PB_2)
                        )

                    )

                    + PC[m] * (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][a2] * delta[b0][b1] * delta[b2][n] + delta[a1][a2] * delta[b0][b2] * delta[b1][n] + delta[a1][a2] * delta[b0][n] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][n] + delta[a1][b0] * delta[a2][b2] * delta[b1][n] + delta[a1][b0] * delta[a2][n] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][n] + delta[a1][b1] * delta[a2][b2] * delta[b0][n] + delta[a1][b1] * delta[a2][n] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][n] + delta[a1][b2] * delta[a2][b1] * delta[b0][n] + delta[a1][b2] * delta[a2][n] * delta[b0][b1] + delta[a1][n] * delta[a2][b0] * delta[b1][b2] + delta[a1][n] * delta[a2][b1] * delta[b0][b2] + delta[a1][n] * delta[a2][b2] * delta[b0][b1]) * (PA_0)
                            + (delta[a0][a2] * delta[b0][b1] * delta[b2][n] + delta[a0][a2] * delta[b0][b2] * delta[b1][n] + delta[a0][a2] * delta[b0][n] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][n] + delta[a0][b0] * delta[a2][b2] * delta[b1][n] + delta[a0][b0] * delta[a2][n] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][n] + delta[a0][b1] * delta[a2][b2] * delta[b0][n] + delta[a0][b1] * delta[a2][n] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][n] + delta[a0][b2] * delta[a2][b1] * delta[b0][n] + delta[a0][b2] * delta[a2][n] * delta[b0][b1] + delta[a0][n] * delta[a2][b0] * delta[b1][b2] + delta[a0][n] * delta[a2][b1] * delta[b0][b2] + delta[a0][n] * delta[a2][b2] * delta[b0][b1]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[b2][n] + delta[a0][a1] * delta[b0][b2] * delta[b1][n] + delta[a0][a1] * delta[b0][n] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][n] + delta[a0][b0] * delta[a1][b2] * delta[b1][n] + delta[a0][b0] * delta[a1][n] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][n] + delta[a0][b1] * delta[a1][b2] * delta[b0][n] + delta[a0][b1] * delta[a1][n] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][n] + delta[a0][b2] * delta[a1][b1] * delta[b0][n] + delta[a0][b2] * delta[a1][n] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[b0][b2] + delta[a0][n] * delta[a1][b2] * delta[b0][b1]) * (PA_2)
                            + (delta[a0][a1] * delta[a2][b1] * delta[b2][n] + delta[a0][a1] * delta[a2][b2] * delta[b1][n] + delta[a0][a1] * delta[a2][n] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][n] + delta[a0][a2] * delta[a1][b2] * delta[b1][n] + delta[a0][a2] * delta[a1][n] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][n] + delta[a0][b1] * delta[a1][b2] * delta[a2][n] + delta[a0][b1] * delta[a1][n] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][n] + delta[a0][b2] * delta[a1][b1] * delta[a2][n] + delta[a0][b2] * delta[a1][n] * delta[a2][b1] + delta[a0][n] * delta[a1][a2] * delta[b1][b2] + delta[a0][n] * delta[a1][b1] * delta[a2][b2] + delta[a0][n] * delta[a1][b2] * delta[a2][b1]) * (PB_0)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b2][n] + delta[a0][a1] * delta[a2][b2] * delta[b0][n] + delta[a0][a1] * delta[a2][n] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][n] + delta[a0][a2] * delta[a1][b2] * delta[b0][n] + delta[a0][a2] * delta[a1][n] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][n] + delta[a0][b0] * delta[a1][b2] * delta[a2][n] + delta[a0][b0] * delta[a1][n] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][n] + delta[a0][b2] * delta[a1][b0] * delta[a2][n] + delta[a0][b2] * delta[a1][n] * delta[a2][b0] + delta[a0][n] * delta[a1][a2] * delta[b0][b2] + delta[a0][n] * delta[a1][b0] * delta[a2][b2] + delta[a0][n] * delta[a1][b2] * delta[a2][b0]) * (PB_1)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b1][n] + delta[a0][a1] * delta[a2][b1] * delta[b0][n] + delta[a0][a1] * delta[a2][n] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][n] + delta[a0][a2] * delta[a1][b1] * delta[b0][n] + delta[a0][a2] * delta[a1][n] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][n] + delta[a0][b0] * delta[a1][b1] * delta[a2][n] + delta[a0][b0] * delta[a1][n] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][n] + delta[a0][b1] * delta[a1][b0] * delta[a2][n] + delta[a0][b1] * delta[a1][n] * delta[a2][b0] + delta[a0][n] * delta[a1][a2] * delta[b0][b1] + delta[a0][n] * delta[a1][b0] * delta[a2][b1] + delta[a0][n] * delta[a1][b1] * delta[a2][b0]) * (PB_2)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][n] + delta[b0][b2] * delta[b1][n] + delta[b0][n] * delta[b1][b2]) * (PA_0 * PA_1 * PA_2)
                            + (delta[a2][b1] * delta[b2][n] + delta[a2][b2] * delta[b1][n] + delta[a2][n] * delta[b1][b2]) * (PA_0 * PA_1 * PB_0)
                            + (delta[a2][b0] * delta[b2][n] + delta[a2][b2] * delta[b0][n] + delta[a2][n] * delta[b0][b2]) * (PA_0 * PA_1 * PB_1)
                            + (delta[a2][b0] * delta[b1][n] + delta[a2][b1] * delta[b0][n] + delta[a2][n] * delta[b0][b1]) * (PA_0 * PA_1 * PB_2)
                            + (delta[a1][b1] * delta[b2][n] + delta[a1][b2] * delta[b1][n] + delta[a1][n] * delta[b1][b2]) * (PA_0 * PA_2 * PB_0)
                            + (delta[a1][b0] * delta[b2][n] + delta[a1][b2] * delta[b0][n] + delta[a1][n] * delta[b0][b2]) * (PA_0 * PA_2 * PB_1)
                            + (delta[a1][b0] * delta[b1][n] + delta[a1][b1] * delta[b0][n] + delta[a1][n] * delta[b0][b1]) * (PA_0 * PA_2 * PB_2)
                            + (delta[a1][a2] * delta[b2][n] + delta[a1][b2] * delta[a2][n] + delta[a1][n] * delta[a2][b2]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][a2] * delta[b1][n] + delta[a1][b1] * delta[a2][n] + delta[a1][n] * delta[a2][b1]) * (PA_0 * PB_0 * PB_2)
                            + (delta[a1][a2] * delta[b0][n] + delta[a1][b0] * delta[a2][n] + delta[a1][n] * delta[a2][b0]) * (PA_0 * PB_1 * PB_2)
                            + (delta[a0][b1] * delta[b2][n] + delta[a0][b2] * delta[b1][n] + delta[a0][n] * delta[b1][b2]) * (PA_1 * PA_2 * PB_0)
                            + (delta[a0][b0] * delta[b2][n] + delta[a0][b2] * delta[b0][n] + delta[a0][n] * delta[b0][b2]) * (PA_1 * PA_2 * PB_1)
                            + (delta[a0][b0] * delta[b1][n] + delta[a0][b1] * delta[b0][n] + delta[a0][n] * delta[b0][b1]) * (PA_1 * PA_2 * PB_2)
                            + (delta[a0][a2] * delta[b2][n] + delta[a0][b2] * delta[a2][n] + delta[a0][n] * delta[a2][b2]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][a2] * delta[b1][n] + delta[a0][b1] * delta[a2][n] + delta[a0][n] * delta[a2][b1]) * (PA_1 * PB_0 * PB_2)
                            + (delta[a0][a2] * delta[b0][n] + delta[a0][b0] * delta[a2][n] + delta[a0][n] * delta[a2][b0]) * (PA_1 * PB_1 * PB_2)
                            + (delta[a0][a1] * delta[b2][n] + delta[a0][b2] * delta[a1][n] + delta[a0][n] * delta[a1][b2]) * (PA_2 * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][n] + delta[a0][b1] * delta[a1][n] + delta[a0][n] * delta[a1][b1]) * (PA_2 * PB_0 * PB_2)
                            + (delta[a0][a1] * delta[b0][n] + delta[a0][b0] * delta[a1][n] + delta[a0][n] * delta[a1][b0]) * (PA_2 * PB_1 * PB_2)
                            + (delta[a0][a1] * delta[a2][n] + delta[a0][a2] * delta[a1][n] + delta[a0][n] * delta[a1][a2]) * (PB_0 * PB_1 * PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][n] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_1)
                            + delta[b1][n] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_2)
                            + delta[b0][n] * (PA_0 * PA_1 * PA_2 * PB_1 * PB_2)
                            + delta[a2][n] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_2)
                            + delta[a1][n] * (PA_0 * PA_2 * PB_0 * PB_1 * PB_2)
                            + delta[a0][n] * (PA_1 * PA_2 * PB_0 * PB_1 * PB_2)
                        )

                    )

                    + PC[n] * (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (PA_0)
                            + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (PA_2)
                            + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (PB_0)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (PB_1)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (PB_2)
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (PA_0 * PA_1 * PA_2)
                            + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (PA_0 * PA_1 * PB_0)
                            + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (PA_0 * PA_1 * PB_1)
                            + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (PA_0 * PA_1 * PB_2)
                            + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (PA_0 * PA_2 * PB_0)
                            + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (PA_0 * PA_2 * PB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_2 * PB_2)
                            + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (PA_0 * PB_0 * PB_2)
                            + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (PA_0 * PB_1 * PB_2)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (PA_1 * PA_2 * PB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (PA_1 * PA_2 * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_2 * PB_2)
                            + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (PA_1 * PB_0 * PB_2)
                            + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (PA_1 * PB_1 * PB_2)
                            + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (PA_2 * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_2 * PB_0 * PB_2)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_2 * PB_1 * PB_2)
                            + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (PB_0 * PB_1 * PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b2][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_1 * PA_2 * PB_0 * PB_2)
                            + delta[b0][m] * (PA_0 * PA_1 * PA_2 * PB_1 * PB_2)
                            + delta[a2][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_2)
                            + delta[a1][m] * (PA_0 * PA_2 * PB_0 * PB_1 * PB_2)
                            + delta[a0][m] * (PA_1 * PA_2 * PB_0 * PB_1 * PB_2)
                        )

                    )

                    + (PC[n] * PC[m] + 0.5 / (a_i + a_j) * delta[n][m]) * (

                        0.125 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0])
                        )

                        + 0.25 / ( (a_i + a_j) * (a_i + a_j) ) * (
                            (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (PA_0 * PA_1)
                            + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (PA_0 * PA_2)
                            + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (PA_0 * PB_0)
                            + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (PA_0 * PB_1)
                            + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (PA_0 * PB_2)
                            + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (PA_1 * PA_2)
                            + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (PA_1 * PB_0)
                            + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (PA_1 * PB_1)
                            + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (PA_1 * PB_2)
                            + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (PA_2 * PB_0)
                            + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (PA_2 * PB_1)
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_2 * PB_2)
                            + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (PB_0 * PB_1)
                            + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (PB_0 * PB_2)
                            + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (PB_1 * PB_2)
                        )

                        + 0.5 / (a_i + a_j) * (
                            delta[b1][b2] * (PA_0 * PA_1 * PA_2 * PB_0)
                            + delta[b0][b2] * (PA_0 * PA_1 * PA_2 * PB_1)
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_2 * PB_2)
                            + delta[a2][b2] * (PA_0 * PA_1 * PB_0 * PB_1)
                            + delta[a2][b1] * (PA_0 * PA_1 * PB_0 * PB_2)
                            + delta[a2][b0] * (PA_0 * PA_1 * PB_1 * PB_2)
                            + delta[a1][b2] * (PA_0 * PA_2 * PB_0 * PB_1)
                            + delta[a1][b1] * (PA_0 * PA_2 * PB_0 * PB_2)
                            + delta[a1][b0] * (PA_0 * PA_2 * PB_1 * PB_2)
                            + delta[a1][a2] * (PA_0 * PB_0 * PB_1 * PB_2)
                            + delta[a0][b2] * (PA_1 * PA_2 * PB_0 * PB_1)
                            + delta[a0][b1] * (PA_1 * PA_2 * PB_0 * PB_2)
                            + delta[a0][b0] * (PA_1 * PA_2 * PB_1 * PB_2)
                            + delta[a0][a2] * (PA_1 * PB_0 * PB_1 * PB_2)
                            + delta[a0][a1] * (PA_2 * PB_0 * PB_1 * PB_2)
                        )

                        + (
                            PA_0 * PA_1 * PA_2 * PB_0 * PB_1 * PB_2
                        )

                    )

            );
        }
    }

    for (int ij = 0; ij < ff_prim_pair_count; ij++)
    {
        const auto i = std::get<0>(pair_inds_ff[ij]);
        const auto j = std::get<1>(pair_inds_ff[ij]);

        const auto i_cgto = f_prim_aoinds[(i / 10) + f_prim_count * (i % 10)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        // Cartesian to spherical
        for (const auto& i_cgto_sph_ind_coef : cart_sph_f[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                MXX.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                if (i != j) MXX.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xx[ij] * coef_sph;

                MXY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                if (i != j) MXY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xy[ij] * coef_sph;

                MXZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                if (i != j) MXZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_xz[ij] * coef_sph;

                MYY.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                if (i != j) MYY.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yy[ij] * coef_sph;

                MYZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                if (i != j) MYZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_yz[ij] * coef_sph;

                MZZ.row(i_cgto_sph)[j_cgto_sph] += mat_mu_zz[ij] * coef_sph;

                if (i != j) MZZ.row(j_cgto_sph)[i_cgto_sph] += mat_mu_zz[ij] * coef_sph;
            }
        }
    }

    // auto-generated code ends here

    return matrices_mu;
}

}  // namespace onee
