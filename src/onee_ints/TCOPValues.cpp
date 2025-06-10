//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "TCOPValues.hpp"

#include <omp.h>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"

#define CHUNK_SIZE 4

#define MATH_CONST_PI 3.14159265358979323846

namespace onee {  // onee namespace

auto
computeTCOPValues(const             CMolecule& molecule, 
                  const             CMolecularBasis& basis, 
                  const double*     point_coords, 
                  const int         npoints, 
                  const double*     point_exp, 
                  const double*     point_amp, 
                  const double*     point_norms,
                  const double*     D, 
                  const int         naos) -> std::vector<double>
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    errors::assertMsgCritical(
            naos == gtofunc::getNumberOfAtomicOrbitals(gto_blocks),
            std::string("computeTCOPValues: Inconsistent number of AOs"));

    auto nthreads = omp_get_max_threads();

    std::vector<std::vector<double>> f_tilde_values_omp(nthreads);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        f_tilde_values_omp[thread_id] = std::vector<double>(npoints, 0.0);
    }

    // points info

    std::vector<double> points_info(npoints * 8);

    for (int c = 0; c < npoints; c++)
    {
        points_info[c + npoints * 0] = point_coords[c * 3 + 0];
        points_info[c + npoints * 1] = point_coords[c * 3 + 1];
        points_info[c + npoints * 2] = point_coords[c * 3 + 2];
        points_info[c + npoints * 3] = point_exp[c];
        points_info[c + npoints * 4] = point_amp[c];
        points_info[c + npoints * 5] = point_norms[c * 3 + 0];
        points_info[c + npoints * 6] = point_norms[c * 3 + 1];
        points_info[c + npoints * 7] = point_norms[c * 3 + 2];
    }

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
            std::string errangmom("computeFock2GOSTcontrib: Only implemented up to f-orbitals");

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

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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



        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = s_prim_aoinds[j];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            {
                auto j_cgto_sph = j_cgto;
                double j_coef_sph = 1.0;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = ((i == j) ? Dij : (Dij + Dji));

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};



            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        1.0 * (
                            1.0
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // S-P block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        (
                            GB_0
                        )

                    )

                    + G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b0][m]
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // S-D block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < sd_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b0][b1]
                        )

                        + (
                            GB_0 * GB_1
                        )

                    )

                    + G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b1][m] * (GB_0)
                            + delta[b0][m] * (GB_1)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // S-F block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < sf_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = s_prim_aoinds[i];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        {
            auto i_cgto_sph = i_cgto;
            double i_coef_sph = 1.0;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);
            const auto GB_2 = (-a_i * rij[b2] + zeta_c * rcj[b2]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b1][b2] * (GB_0)
                            + delta[b0][b2] * (GB_1)
                            + delta[b0][b1] * (GB_2)
                        )

                        + (
                            GB_0 * GB_1 * GB_2
                        )

                    )

                    + G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2])
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b2][m] * (GB_0 * GB_1)
                            + delta[b1][m] * (GB_0 * GB_2)
                            + delta[b0][m] * (GB_1 * GB_2)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // P-P block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_p[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = ((i == j) ? Dij : (Dij + Dji));

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[a0][b0]
                        )

                        + (
                            GA_0 * GB_0
                        )

                    )

                    + G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b0][m] * (GA_0)
                            + delta[a0][m] * (GB_0)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // P-D block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < pd_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.5 / (zeta + zeta_c) * (
                            delta[b0][b1] * (GA_0)
                            + delta[a0][b1] * (GB_0)
                            + delta[a0][b0] * (GB_1)
                        )

                        + (
                            GA_0 * GB_0 * GB_1
                        )

                    )

                    + G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1])
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b1][m] * (GA_0 * GB_0)
                            + delta[b0][m] * (GA_0 * GB_1)
                            + delta[a0][m] * (GB_0 * GB_1)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // P-F block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < pf_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_p[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);
            const auto GB_2 = (-a_i * rij[b2] + zeta_c * rcj[b2]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1])
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b1][b2] * (GA_0 * GB_0)
                            + delta[b0][b2] * (GA_0 * GB_1)
                            + delta[b0][b1] * (GA_0 * GB_2)
                            + delta[a0][b2] * (GB_0 * GB_1)
                            + delta[a0][b1] * (GB_0 * GB_2)
                            + delta[a0][b0] * (GB_1 * GB_2)
                        )

                        + (
                            GA_0 * GB_0 * GB_1 * GB_2
                        )

                    )

                    + G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (GA_0)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (GB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (GB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (GB_2)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b2][m] * (GA_0 * GB_0 * GB_1)
                            + delta[b1][m] * (GA_0 * GB_0 * GB_2)
                            + delta[b0][m] * (GA_0 * GB_1 * GB_2)
                            + delta[a0][m] * (GB_0 * GB_1 * GB_2)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // D-D block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < dd_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
        const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_d[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = ((i == j) ? Dij : (Dij + Dji));

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);
            const auto GA_1 = (a_j * rij[a1] + zeta_c * rci[a1]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0])
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b0][b1] * (GA_0 * GA_1)
                            + delta[a1][b1] * (GA_0 * GB_0)
                            + delta[a1][b0] * (GA_0 * GB_1)
                            + delta[a0][b1] * (GA_1 * GB_0)
                            + delta[a0][b0] * (GA_1 * GB_1)
                            + delta[a0][a1] * (GB_0 * GB_1)
                        )

                        + (
                            GA_0 * GA_1 * GB_0 * GB_1
                        )

                    )

                    + G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (GA_0)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (GA_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (GB_0)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (GB_1)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b1][m] * (GA_0 * GA_1 * GB_0)
                            + delta[b0][m] * (GA_0 * GA_1 * GB_1)
                            + delta[a1][m] * (GA_0 * GB_0 * GB_1)
                            + delta[a0][m] * (GA_1 * GB_0 * GB_1)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // D-F block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < df_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_d[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = (Dij + Dji);

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);
            const auto GA_1 = (a_j * rij[a1] + zeta_c * rci[a1]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);
            const auto GB_2 = (-a_i * rij[b2] + zeta_c * rcj[b2]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (GA_0)
                            + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (GA_1)
                            + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (GB_0)
                            + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (GB_1)
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (GB_2)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b1][b2] * (GA_0 * GA_1 * GB_0)
                            + delta[b0][b2] * (GA_0 * GA_1 * GB_1)
                            + delta[b0][b1] * (GA_0 * GA_1 * GB_2)
                            + delta[a1][b2] * (GA_0 * GB_0 * GB_1)
                            + delta[a1][b1] * (GA_0 * GB_0 * GB_2)
                            + delta[a1][b0] * (GA_0 * GB_1 * GB_2)
                            + delta[a0][b2] * (GA_1 * GB_0 * GB_1)
                            + delta[a0][b1] * (GA_1 * GB_0 * GB_2)
                            + delta[a0][b0] * (GA_1 * GB_1 * GB_2)
                            + delta[a0][a1] * (GB_0 * GB_1 * GB_2)
                        )

                        + (
                            GA_0 * GA_1 * GB_0 * GB_1 * GB_2
                        )

                    )

                    + G_ij_00 * (

                        0.125 / ( (zeta + zeta_c) * (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1])
                        )

                        + 0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (GA_0 * GA_1)
                            + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (GA_0 * GB_0)
                            + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (GA_0 * GB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (GA_0 * GB_2)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (GA_1 * GB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (GA_1 * GB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (GA_1 * GB_2)
                            + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (GB_0 * GB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (GB_0 * GB_2)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (GB_1 * GB_2)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b2][m] * (GA_0 * GA_1 * GB_0 * GB_1)
                            + delta[b1][m] * (GA_0 * GA_1 * GB_0 * GB_2)
                            + delta[b0][m] * (GA_0 * GA_1 * GB_1 * GB_2)
                            + delta[a1][m] * (GA_0 * GB_0 * GB_1 * GB_2)
                            + delta[a0][m] * (GA_1 * GB_0 * GB_1 * GB_2)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }


    // F-F block

    #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
        
    for (int ij = 0; ij < ff_prim_pair_count; ij++)
    {
        const auto thread_id = omp_get_thread_num();

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

        const auto i_cgto = f_prim_aoinds[(i / 10) + f_prim_count * (i % 10)];
        const auto j_cgto = f_prim_aoinds[(j / 10) + f_prim_count * (j % 10)];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto zeta = a_i + a_j;

        // product of density and cart-sph transformation coefficients

        double dens_coef_prod = 0.0;

        for (const auto& i_cgto_sph_ind_coef : cart_sph_f[i_cgto])
        {
            auto i_cgto_sph = i_cgto_sph_ind_coef.first;
            auto i_coef_sph = i_cgto_sph_ind_coef.second;

            for (const auto& j_cgto_sph_ind_coef : cart_sph_f[j_cgto])
            {
                auto j_cgto_sph = j_cgto_sph_ind_coef.first;
                auto j_coef_sph = j_cgto_sph_ind_coef.second;

                auto coef_sph = i_coef_sph * j_coef_sph;

                auto Dij = D[i_cgto_sph * naos + j_cgto_sph];
                auto Dji = D[j_cgto_sph * naos + i_cgto_sph];

                double D_sym = ((i == j) ? Dij : (Dij + Dji));

                dens_coef_prod += coef_sph * D_sym;
            }
        }

        // J. Chem. Phys. 84, 3963-3974 (1986)
            
        for (int c = 0; c < npoints; c++)
        {  
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto zeta_c = points_info[c + npoints * 3];
            const auto p_c = points_info[c + npoints * 4];
            const auto nx_c = points_info[c + npoints * 5];
            const auto ny_c = points_info[c + npoints * 6];
            const auto nz_c = points_info[c + npoints * 7];

            const double n_c[3] = {nx_c, ny_c, nz_c};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            const double rci[3] = {x_c - x_i, y_c - y_i, z_c - z_i};

            const double rcj[3] = {x_c - x_j, y_c - y_j, z_c - z_j};

            const auto G_ij_00 = std::pow((a_i + a_j) / (a_i + a_j + zeta_c), 1.5) * std::exp(-(a_i + a_j) * zeta_c / (a_i + a_j + zeta_c) * r2_PC);

            const double GC[3] = {(a_i * x_i + a_j * x_j + zeta_c * x_c) / (a_i +  a_j + zeta_c) - x_c,
                            (a_i * y_i + a_j * y_j + zeta_c * y_c) / (a_i + a_j + zeta_c) - y_c,
                            (a_i * z_i + a_j * z_j + zeta_c * z_c) / (a_i + a_j + zeta_c) - z_c};
            const auto GA_0 = (a_j * rij[a0] + zeta_c * rci[a0]) / (a_i + a_j + zeta_c);
            const auto GA_1 = (a_j * rij[a1] + zeta_c * rci[a1]) / (a_i + a_j + zeta_c);
            const auto GA_2 = (a_j * rij[a2] + zeta_c * rci[a2]) / (a_i + a_j + zeta_c);

            const auto GB_0 = (-a_i * rij[b0] + zeta_c * rcj[b0]) / (a_i + a_j + zeta_c);
            const auto GB_1 = (-a_i * rij[b1] + zeta_c * rcj[b1]) / (a_i + a_j + zeta_c);
            const auto GB_2 = (-a_i * rij[b2] + zeta_c * rcj[b2]) / (a_i + a_j + zeta_c);


            double tco_p_m_val = 0.0;

            for (int m = 0; m < 3; m++)
                    {
                        tco_p_m_val -= 2 * zeta_c * n_c[m] * p_c * S_ij_00 * (


                    GC[m] * G_ij_00 * (

                        0.125 / ( (zeta + zeta_c) * (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a0][a1] * delta[a2][b0] * delta[b1][b2] + delta[a0][a1] * delta[a2][b1] * delta[b0][b2] + delta[a0][a1] * delta[a2][b2] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b0][b2] + delta[a0][a2] * delta[a1][b2] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[a2][b2] + delta[a0][b0] * delta[a1][b2] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][b2] + delta[a0][b1] * delta[a1][b0] * delta[a2][b2] + delta[a0][b1] * delta[a1][b2] * delta[a2][b0] + delta[a0][b2] * delta[a1][a2] * delta[b0][b1] + delta[a0][b2] * delta[a1][b0] * delta[a2][b1] + delta[a0][b2] * delta[a1][b1] * delta[a2][b0])
                        )

                        + 0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a2][b0] * delta[b1][b2] + delta[a2][b1] * delta[b0][b2] + delta[a2][b2] * delta[b0][b1]) * (GA_0 * GA_1)
                            + (delta[a1][b0] * delta[b1][b2] + delta[a1][b1] * delta[b0][b2] + delta[a1][b2] * delta[b0][b1]) * (GA_0 * GA_2)
                            + (delta[a1][a2] * delta[b1][b2] + delta[a1][b1] * delta[a2][b2] + delta[a1][b2] * delta[a2][b1]) * (GA_0 * GB_0)
                            + (delta[a1][a2] * delta[b0][b2] + delta[a1][b0] * delta[a2][b2] + delta[a1][b2] * delta[a2][b0]) * (GA_0 * GB_1)
                            + (delta[a1][a2] * delta[b0][b1] + delta[a1][b0] * delta[a2][b1] + delta[a1][b1] * delta[a2][b0]) * (GA_0 * GB_2)
                            + (delta[a0][b0] * delta[b1][b2] + delta[a0][b1] * delta[b0][b2] + delta[a0][b2] * delta[b0][b1]) * (GA_1 * GA_2)
                            + (delta[a0][a2] * delta[b1][b2] + delta[a0][b1] * delta[a2][b2] + delta[a0][b2] * delta[a2][b1]) * (GA_1 * GB_0)
                            + (delta[a0][a2] * delta[b0][b2] + delta[a0][b0] * delta[a2][b2] + delta[a0][b2] * delta[a2][b0]) * (GA_1 * GB_1)
                            + (delta[a0][a2] * delta[b0][b1] + delta[a0][b0] * delta[a2][b1] + delta[a0][b1] * delta[a2][b0]) * (GA_1 * GB_2)
                            + (delta[a0][a1] * delta[b1][b2] + delta[a0][b1] * delta[a1][b2] + delta[a0][b2] * delta[a1][b1]) * (GA_2 * GB_0)
                            + (delta[a0][a1] * delta[b0][b2] + delta[a0][b0] * delta[a1][b2] + delta[a0][b2] * delta[a1][b0]) * (GA_2 * GB_1)
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (GA_2 * GB_2)
                            + (delta[a0][a1] * delta[a2][b2] + delta[a0][a2] * delta[a1][b2] + delta[a0][b2] * delta[a1][a2]) * (GB_0 * GB_1)
                            + (delta[a0][a1] * delta[a2][b1] + delta[a0][a2] * delta[a1][b1] + delta[a0][b1] * delta[a1][a2]) * (GB_0 * GB_2)
                            + (delta[a0][a1] * delta[a2][b0] + delta[a0][a2] * delta[a1][b0] + delta[a0][b0] * delta[a1][a2]) * (GB_1 * GB_2)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b1][b2] * (GA_0 * GA_1 * GA_2 * GB_0)
                            + delta[b0][b2] * (GA_0 * GA_1 * GA_2 * GB_1)
                            + delta[b0][b1] * (GA_0 * GA_1 * GA_2 * GB_2)
                            + delta[a2][b2] * (GA_0 * GA_1 * GB_0 * GB_1)
                            + delta[a2][b1] * (GA_0 * GA_1 * GB_0 * GB_2)
                            + delta[a2][b0] * (GA_0 * GA_1 * GB_1 * GB_2)
                            + delta[a1][b2] * (GA_0 * GA_2 * GB_0 * GB_1)
                            + delta[a1][b1] * (GA_0 * GA_2 * GB_0 * GB_2)
                            + delta[a1][b0] * (GA_0 * GA_2 * GB_1 * GB_2)
                            + delta[a1][a2] * (GA_0 * GB_0 * GB_1 * GB_2)
                            + delta[a0][b2] * (GA_1 * GA_2 * GB_0 * GB_1)
                            + delta[a0][b1] * (GA_1 * GA_2 * GB_0 * GB_2)
                            + delta[a0][b0] * (GA_1 * GA_2 * GB_1 * GB_2)
                            + delta[a0][a2] * (GA_1 * GB_0 * GB_1 * GB_2)
                            + delta[a0][a1] * (GA_2 * GB_0 * GB_1 * GB_2)
                        )

                        + (
                            GA_0 * GA_1 * GA_2 * GB_0 * GB_1 * GB_2
                        )

                    )

                    + G_ij_00 * (

                        0.125 / ( (zeta + zeta_c) * (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[a1][a2] * delta[b0][b1] * delta[b2][m] + delta[a1][a2] * delta[b0][b2] * delta[b1][m] + delta[a1][a2] * delta[b0][m] * delta[b1][b2] + delta[a1][b0] * delta[a2][b1] * delta[b2][m] + delta[a1][b0] * delta[a2][b2] * delta[b1][m] + delta[a1][b0] * delta[a2][m] * delta[b1][b2] + delta[a1][b1] * delta[a2][b0] * delta[b2][m] + delta[a1][b1] * delta[a2][b2] * delta[b0][m] + delta[a1][b1] * delta[a2][m] * delta[b0][b2] + delta[a1][b2] * delta[a2][b0] * delta[b1][m] + delta[a1][b2] * delta[a2][b1] * delta[b0][m] + delta[a1][b2] * delta[a2][m] * delta[b0][b1] + delta[a1][m] * delta[a2][b0] * delta[b1][b2] + delta[a1][m] * delta[a2][b1] * delta[b0][b2] + delta[a1][m] * delta[a2][b2] * delta[b0][b1]) * (GA_0)
                            + (delta[a0][a2] * delta[b0][b1] * delta[b2][m] + delta[a0][a2] * delta[b0][b2] * delta[b1][m] + delta[a0][a2] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a2][b1] * delta[b2][m] + delta[a0][b0] * delta[a2][b2] * delta[b1][m] + delta[a0][b0] * delta[a2][m] * delta[b1][b2] + delta[a0][b1] * delta[a2][b0] * delta[b2][m] + delta[a0][b1] * delta[a2][b2] * delta[b0][m] + delta[a0][b1] * delta[a2][m] * delta[b0][b2] + delta[a0][b2] * delta[a2][b0] * delta[b1][m] + delta[a0][b2] * delta[a2][b1] * delta[b0][m] + delta[a0][b2] * delta[a2][m] * delta[b0][b1] + delta[a0][m] * delta[a2][b0] * delta[b1][b2] + delta[a0][m] * delta[a2][b1] * delta[b0][b2] + delta[a0][m] * delta[a2][b2] * delta[b0][b1]) * (GA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[b2][m] + delta[a0][a1] * delta[b0][b2] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][b2] + delta[a0][b0] * delta[a1][b1] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][b0] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][b2] + delta[a0][b2] * delta[a1][b0] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[b0][m] + delta[a0][b2] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[b0][b2] + delta[a0][m] * delta[a1][b2] * delta[b0][b1]) * (GA_2)
                            + (delta[a0][a1] * delta[a2][b1] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b1][m] + delta[a0][a1] * delta[a2][m] * delta[b1][b2] + delta[a0][a2] * delta[a1][b1] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b1][m] + delta[a0][a2] * delta[a1][m] * delta[b1][b2] + delta[a0][b1] * delta[a1][a2] * delta[b2][m] + delta[a0][b1] * delta[a1][b2] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b1][m] + delta[a0][b2] * delta[a1][b1] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b1] + delta[a0][m] * delta[a1][a2] * delta[b1][b2] + delta[a0][m] * delta[a1][b1] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b1]) * (GB_0)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b2][m] + delta[a0][a1] * delta[a2][b2] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b2] + delta[a0][a2] * delta[a1][b0] * delta[b2][m] + delta[a0][a2] * delta[a1][b2] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b2] + delta[a0][b0] * delta[a1][a2] * delta[b2][m] + delta[a0][b0] * delta[a1][b2] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b2] + delta[a0][b2] * delta[a1][a2] * delta[b0][m] + delta[a0][b2] * delta[a1][b0] * delta[a2][m] + delta[a0][b2] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b2] + delta[a0][m] * delta[a1][b0] * delta[a2][b2] + delta[a0][m] * delta[a1][b2] * delta[a2][b0]) * (GB_1)
                            + (delta[a0][a1] * delta[a2][b0] * delta[b1][m] + delta[a0][a1] * delta[a2][b1] * delta[b0][m] + delta[a0][a1] * delta[a2][m] * delta[b0][b1] + delta[a0][a2] * delta[a1][b0] * delta[b1][m] + delta[a0][a2] * delta[a1][b1] * delta[b0][m] + delta[a0][a2] * delta[a1][m] * delta[b0][b1] + delta[a0][b0] * delta[a1][a2] * delta[b1][m] + delta[a0][b0] * delta[a1][b1] * delta[a2][m] + delta[a0][b0] * delta[a1][m] * delta[a2][b1] + delta[a0][b1] * delta[a1][a2] * delta[b0][m] + delta[a0][b1] * delta[a1][b0] * delta[a2][m] + delta[a0][b1] * delta[a1][m] * delta[a2][b0] + delta[a0][m] * delta[a1][a2] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[a2][b1] + delta[a0][m] * delta[a1][b1] * delta[a2][b0]) * (GB_2)
                        )

                        + 0.25 / ( (zeta + zeta_c) * (zeta + zeta_c) ) * (
                            (delta[b0][b1] * delta[b2][m] + delta[b0][b2] * delta[b1][m] + delta[b0][m] * delta[b1][b2]) * (GA_0 * GA_1 * GA_2)
                            + (delta[a2][b1] * delta[b2][m] + delta[a2][b2] * delta[b1][m] + delta[a2][m] * delta[b1][b2]) * (GA_0 * GA_1 * GB_0)
                            + (delta[a2][b0] * delta[b2][m] + delta[a2][b2] * delta[b0][m] + delta[a2][m] * delta[b0][b2]) * (GA_0 * GA_1 * GB_1)
                            + (delta[a2][b0] * delta[b1][m] + delta[a2][b1] * delta[b0][m] + delta[a2][m] * delta[b0][b1]) * (GA_0 * GA_1 * GB_2)
                            + (delta[a1][b1] * delta[b2][m] + delta[a1][b2] * delta[b1][m] + delta[a1][m] * delta[b1][b2]) * (GA_0 * GA_2 * GB_0)
                            + (delta[a1][b0] * delta[b2][m] + delta[a1][b2] * delta[b0][m] + delta[a1][m] * delta[b0][b2]) * (GA_0 * GA_2 * GB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (GA_0 * GA_2 * GB_2)
                            + (delta[a1][a2] * delta[b2][m] + delta[a1][b2] * delta[a2][m] + delta[a1][m] * delta[a2][b2]) * (GA_0 * GB_0 * GB_1)
                            + (delta[a1][a2] * delta[b1][m] + delta[a1][b1] * delta[a2][m] + delta[a1][m] * delta[a2][b1]) * (GA_0 * GB_0 * GB_2)
                            + (delta[a1][a2] * delta[b0][m] + delta[a1][b0] * delta[a2][m] + delta[a1][m] * delta[a2][b0]) * (GA_0 * GB_1 * GB_2)
                            + (delta[a0][b1] * delta[b2][m] + delta[a0][b2] * delta[b1][m] + delta[a0][m] * delta[b1][b2]) * (GA_1 * GA_2 * GB_0)
                            + (delta[a0][b0] * delta[b2][m] + delta[a0][b2] * delta[b0][m] + delta[a0][m] * delta[b0][b2]) * (GA_1 * GA_2 * GB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (GA_1 * GA_2 * GB_2)
                            + (delta[a0][a2] * delta[b2][m] + delta[a0][b2] * delta[a2][m] + delta[a0][m] * delta[a2][b2]) * (GA_1 * GB_0 * GB_1)
                            + (delta[a0][a2] * delta[b1][m] + delta[a0][b1] * delta[a2][m] + delta[a0][m] * delta[a2][b1]) * (GA_1 * GB_0 * GB_2)
                            + (delta[a0][a2] * delta[b0][m] + delta[a0][b0] * delta[a2][m] + delta[a0][m] * delta[a2][b0]) * (GA_1 * GB_1 * GB_2)
                            + (delta[a0][a1] * delta[b2][m] + delta[a0][b2] * delta[a1][m] + delta[a0][m] * delta[a1][b2]) * (GA_2 * GB_0 * GB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (GA_2 * GB_0 * GB_2)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (GA_2 * GB_1 * GB_2)
                            + (delta[a0][a1] * delta[a2][m] + delta[a0][a2] * delta[a1][m] + delta[a0][m] * delta[a1][a2]) * (GB_0 * GB_1 * GB_2)
                        )

                        + 0.5 / (zeta + zeta_c) * (
                            delta[b2][m] * (GA_0 * GA_1 * GA_2 * GB_0 * GB_1)
                            + delta[b1][m] * (GA_0 * GA_1 * GA_2 * GB_0 * GB_2)
                            + delta[b0][m] * (GA_0 * GA_1 * GA_2 * GB_1 * GB_2)
                            + delta[a2][m] * (GA_0 * GA_1 * GB_0 * GB_1 * GB_2)
                            + delta[a1][m] * (GA_0 * GA_2 * GB_0 * GB_1 * GB_2)
                            + delta[a0][m] * (GA_1 * GA_2 * GB_0 * GB_1 * GB_2)
                        )

                    )

            );

            }

            f_tilde_values_omp[thread_id][c] += tco_p_m_val * dens_coef_prod;
        }
    }

    // auto-generated code ends here

    std::vector<double> f_tilde_values(npoints, 0.0);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int c = 0; c < npoints; c++)
        {
            f_tilde_values[c] += f_tilde_values_omp[thread_id][c];
        }
    }

    return f_tilde_values;
}

}  // namespace onee
