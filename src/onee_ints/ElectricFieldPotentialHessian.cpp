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

#include "ElectricFieldPotentialHessian.hpp"

#include <omp.h>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"

#define PAD_SIZE 8

#define MATH_CONST_PI 3.14159265358979323846

#define MATH_CONST_TWO_OVER_SQRT_PI 1.12837916709551255856

namespace onee {  // onee namespace

auto
computeElectricFieldPotentialHessian(const CMolecule& molecule, const CMolecularBasis& basis, const double* dipole_coords, const double* dipole_moments, const int ndipoles, const double* D, const int naos) -> CDenseMatrix
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    errors::assertMsgCritical(
            naos == gtofunc::getNumberOfAtomicOrbitals(gto_blocks),
            std::string("computeElectricFieldPotentialHessian: Inconsistent number of AOs"));

    auto natoms = molecule.number_of_atoms();

    auto nthreads = omp_get_max_threads();

    std::vector<CDenseMatrix> V_hess_omp(nthreads);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        V_hess_omp[thread_id] = CDenseMatrix(natoms * 3, natoms * 3);

        V_hess_omp[thread_id].zero();
    }

    // Boys function

    const auto boys_func_table = onee::getFullBoysFuncTable();

    const auto boys_func_ft = onee::getBoysFuncFactors();

    // dipoles info

    std::vector<double> dipoles_info(ndipoles * 6);

    for (int c = 0; c < ndipoles; c++)
    {
        dipoles_info[c + ndipoles * 0] = dipole_coords[c * 3 + 0];
        dipoles_info[c + ndipoles * 1] = dipole_coords[c * 3 + 1];
        dipoles_info[c + ndipoles * 2] = dipole_coords[c * 3 + 2];
        dipoles_info[c + ndipoles * 3] = dipole_moments[c * 3 + 0];
        dipoles_info[c + ndipoles * 4] = dipole_moments[c * 3 + 1];
        dipoles_info[c + ndipoles * 5] = dipole_moments[c * 3 + 2];
    }

    // gto blocks

    int s_prim_count = 0;
    int p_prim_count = 0;
    int d_prim_count = 0;
    int f_prim_count = 0;

    int ncgtos_d = 0;
    int naos_cart = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.number_of_basis_functions();
        const auto npgtos = gto_block.number_of_primitives();

        const auto gto_ang = gto_block.angular_momentum();

        if (gto_ang == 0)
        {
            s_prim_count += npgtos * ncgtos;

            naos_cart += ncgtos;
        }
        else if (gto_ang == 1)
        {
            p_prim_count += npgtos * ncgtos;

            naos_cart += ncgtos * 3;
        }
        else if (gto_ang == 2)
        {
            d_prim_count += npgtos * ncgtos;

            ncgtos_d += ncgtos;

            naos_cart += ncgtos * 6;
        }
        else if (gto_ang == 3)
        {
            f_prim_count += npgtos * ncgtos;

            naos_cart += ncgtos * 10;
        }
        else
        {
            std::string errangmom("computeElectricFieldPotentialHessian: Only implemented up to f-orbitals");

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

    // Cartesian AO to atom mapping

    std::vector<int> cart_ao_to_atom_ids(naos_cart);

    for (int iatom = 0; iatom < natoms; iatom++)
    {
        std::vector<int> atomidx({iatom});

        const auto atom_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, atomidx);

        for (const auto& atom_gto_block : atom_gto_blocks)
        {
            const auto gto_ao_inds = atom_gto_block.getAtomicOrbitalsIndexesForCartesian(ncgtos_d);

            for (const auto aoidx : gto_ao_inds)
            {
                cart_ao_to_atom_ids[aoidx] = iatom;
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

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);




        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F3_t[4];

            onee::computeBoysFunction(F3_t, (a_i + a_j) * r2_PC, 3, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F3_t[1] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_g)
                            + delta[g][m] * (PA_h)
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_g * PC[h] * PC[m]
                            + PA_h * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PC[g])
                            + delta[g][m] * (PC[h])
                        )

                    )

                    + F3_t[3] * (

                        8.0 * (a_i + a_j) * a_i * a_i * (
                            PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F3_t[1] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_g)
                            + delta[g][m] * (PB_h)
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_g * PC[h] * PC[m]
                            + PB_h * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[g])
                            + delta[g][m] * (PC[h])
                        )

                    )

                    + F3_t[3] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F3_t[1] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_h)
                            + delta[h][m] * (PB_g)
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_h * PC[g] * PC[m]
                            + PB_g * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[g])
                            + delta[g][m] * (PC[h])
                        )

                    )

                    + F3_t[3] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F3_t[1] * (

                        (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PB_g)
                            + delta[g][m] * (PB_h)
                        )

                    )

                    + F3_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PB_g * PC[h] * PC[m]
                            + PB_h * PC[g] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PC[g])
                            + delta[g][m] * (PC[h])
                        )

                    )

                    + F3_t[3] * (

                        8.0 * (a_i + a_j) * a_j * a_j * (
                            PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }


    // S-P block

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);


        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F4_t[5];

            onee::computeBoysFunction(F4_t, (a_i + a_j) * r2_PC, 4, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F4_t[1] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PB_0 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][m] * delta[g][h]
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][h] * (PA_g * PC[m])
                            + delta[b0][g] * (PA_h * PC[m])
                            + delta[g][h] * (PB_0 * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PB_0 * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[b0][m] * (PA_g * PA_h)
                            + delta[h][m] * (PA_g * PB_0)
                            + delta[g][m] * (PA_h * PB_0)
                        )

                    )

                    + F4_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][h] * (PA_g * PC[m] + PC[g] * PC[m])
                            + delta[b0][g] * (PA_h * PC[m] + PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PC[b0] * PC[m]
                            + PA_g * PB_0 * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PA_g * PC[b0] + PB_0 * PC[g])
                            + delta[b0][m] * (PA_g * PC[h] + PA_h * PC[g])
                            + delta[g][m] * (PA_h * PC[b0] + PB_0 * PC[h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[b0] * PC[m])
                        )

                    )

                    + F4_t[3] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[b0] * PC[m])
                            + delta[b0][h] * (PC[g] * PC[m])
                            + delta[b0][g] * (PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_g * PC[b0] * PC[h] * PC[m]
                            + PA_h * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PC[g] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PC[b0] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[h])
                            + delta[b0][m] * (PC[g] * PC[h])
                        )

                    )

                    + F4_t[4] * (

                        (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F4_t[1] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][h] * (PA_g * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][h] * delta[g][m]
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][h] * (PA_g * PC[m])
                            + delta[g][h] * (PB_0 * PC[m])
                            + delta[b0][g] * (PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_0 * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_g * PB_0)
                            + delta[b0][m] * (PA_g * PB_h)
                            + delta[g][m] * (PB_0 * PB_h)
                        )

                    )

                    + F4_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][h] * (PA_g * PC[m] + PC[g] * PC[m])
                            + delta[g][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[b0][g] * (PB_h * PC[m] + PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_0 * PC[h] * PC[m]
                            + PA_g * PB_h * PC[b0] * PC[m]
                            + PB_0 * PB_h * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_g * PC[b0] + PB_0 * PC[g])
                            + delta[b0][m] * (PA_g * PC[h] + PB_h * PC[g])
                            + delta[g][m] * (PB_0 * PC[h] + PB_h * PC[b0])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b0][h] * (PC[g] * PC[m])
                        )

                    )

                    + F4_t[3] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[b0] * PC[m])
                            + delta[b0][h] * (PC[g] * PC[m])
                            + delta[b0][g] * (PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_g * PC[b0] * PC[h] * PC[m]
                            + PB_0 * PC[g] * PC[h] * PC[m]
                            + PB_h * PC[b0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PC[b0] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[h])
                            + delta[b0][m] * (PC[g] * PC[h])
                        )

                    )

                    + F4_t[4] * (

                        (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F4_t[1] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][g] * (PA_h * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][g] * delta[h][m]
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][g] * (PA_h * PC[m])
                            + delta[g][h] * (PB_0 * PC[m])
                            + delta[b0][h] * (PB_g * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_0 * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_h * PB_0)
                            + delta[b0][m] * (PA_h * PB_g)
                            + delta[h][m] * (PB_0 * PB_g)
                        )

                    )

                    + F4_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][g] * (PA_h * PC[m] + PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[b0][h] * (PB_g * PC[m] + PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_0 * PC[g] * PC[m]
                            + PA_h * PB_g * PC[b0] * PC[m]
                            + PB_0 * PB_g * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[g][m] * (PA_h * PC[b0] + PB_0 * PC[h])
                            + delta[b0][m] * (PA_h * PC[g] + PB_g * PC[h])
                            + delta[h][m] * (PB_0 * PC[g] + PB_g * PC[b0])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b0][g] * (PC[h] * PC[m])
                        )

                    )

                    + F4_t[3] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[b0] * PC[m])
                            + delta[b0][h] * (PC[g] * PC[m])
                            + delta[b0][g] * (PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_h * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[b0] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PC[b0] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[h])
                            + delta[b0][m] * (PC[g] * PC[h])
                        )

                    )

                    + F4_t[4] * (

                        (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F4_t[1] * (

                        (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PB_0 * PC[m])
                            + delta[b0][h] * (PB_g * PC[m])
                            + delta[b0][g] * (PB_h * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PB_0 * PC[m])
                            + delta[b0][h] * (PB_g * PC[m])
                            + delta[b0][g] * (PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PB_0 * PB_g)
                            + delta[g][m] * (PB_0 * PB_h)
                            + delta[b0][m] * (PB_g * PB_h)
                        )

                    )

                    + F4_t[2] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[b0][h] * (PB_g * PC[m] + PC[g] * PC[m])
                            + delta[b0][g] * (PB_h * PC[m] + PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PB_g * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[g] * PC[m]
                            + PB_g * PB_h * PC[b0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PB_0 * PC[g] + PB_g * PC[b0])
                            + delta[g][m] * (PB_0 * PC[h] + PB_h * PC[b0])
                            + delta[b0][m] * (PB_g * PC[h] + PB_h * PC[g])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[b0] * PC[m])
                            + delta[b0][h] * (PC[g] * PC[m])
                            + delta[b0][g] * (PC[h] * PC[m])
                        )

                    )

                    + F4_t[3] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[b0] * PC[m])
                            + delta[b0][h] * (PC[g] * PC[m])
                            + delta[b0][g] * (PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[b0] * PC[h] * PC[m]
                            + PB_h * PC[b0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PC[b0] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[h])
                            + delta[b0][m] * (PC[g] * PC[h])
                        )

                    )

                    + F4_t[4] * (

                        (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }


    // S-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);


        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F5_t[6];

            onee::computeBoysFunction(F5_t, (a_i + a_j) * r2_PC, 5, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PB_0 * PB_1 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][m] * delta[g][h] * (PB_0)
                            + delta[b0][m] * delta[g][h] * (PB_1)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_g * PA_h * PC[m])
                            + delta[b1][h] * (PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_g * PB_1 * PC[m])
                            + delta[b1][g] * (PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_h * PB_1 * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PB_0 * PB_1 * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[b1][m] * (PA_g * PA_h * PB_0)
                            + delta[b0][m] * (PA_g * PA_h * PB_1)
                            + delta[h][m] * (PA_g * PB_0 * PB_1)
                            + delta[g][m] * (PA_h * PB_0 * PB_1)
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_g * PA_h * PC[m] + PA_g * PC[h] * PC[m] + PA_h * PC[g] * PC[m])
                            + delta[b1][h] * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[b0][h] * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m])
                            + delta[b1][g] * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + delta[b0][g] * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_g + PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_h + PC[h])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1 + PC[b1])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PB_0 * PC[b1] * PC[m]
                            + PA_g * PA_h * PB_1 * PC[b0] * PC[m]
                            + PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[b1][m] * (PA_g * PA_h * PC[b0] + PA_g * PB_0 * PC[h] + PA_h * PB_0 * PC[g])
                            + delta[b0][m] * (PA_g * PA_h * PC[b1] + PA_g * PB_1 * PC[h] + PA_h * PB_1 * PC[g])
                            + delta[h][m] * (PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PB_0 * PB_1 * PC[g])
                            + delta[g][m] * (PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PB_0 * PB_1 * PC[h])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][m] * delta[g][h] * (PC[b0])
                            + delta[b0][m] * delta[g][h] * (PC[b1])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b1][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m] + PC[b1] * PC[g] * PC[m])
                            + delta[b0][b1] * (PA_g * PC[h] * PC[m] + PA_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[b1][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m] + PC[b1] * PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_g * PA_h * PC[b0] * PC[b1] * PC[m]
                            + PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_g * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g])
                            + delta[b1][m] * (PA_g * PC[b0] * PC[h] + PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h])
                            + delta[b0][m] * (PA_g * PC[b1] * PC[h] + PA_h * PC[b1] * PC[g] + PB_1 * PC[g] * PC[h])
                            + delta[g][m] * (PA_h * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_i * (
                            PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_g * PB_1 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_g)
                            + delta[b1][h] * delta[g][m] * (PB_0)
                            + delta[b0][h] * delta[g][m] * (PB_1)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_g * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_g * PB_h * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m])
                            + delta[b1][g] * (PB_0 * PB_h * PC[m])
                            + delta[b0][g] * (PB_1 * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_g)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_0 * PB_1 * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_g * PB_0 * PB_1)
                            + delta[b1][m] * (PA_g * PB_0 * PB_h)
                            + delta[b0][m] * (PA_g * PB_1 * PB_h)
                            + delta[g][m] * (PB_0 * PB_1 * PB_h)
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[b0][h] * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m])
                            + delta[b0][b1] * (PA_g * PB_h * PC[m] + PA_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + delta[b1][g] * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + delta[b0][g] * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_g + PC[g])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1 + PC[b1])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PB_h + PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_g * PB_0 * PB_h * PC[b1] * PC[m]
                            + PA_g * PB_1 * PB_h * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PB_0 * PB_1 * PC[g])
                            + delta[b1][m] * (PA_g * PB_0 * PC[h] + PA_g * PB_h * PC[b0] + PB_0 * PB_h * PC[g])
                            + delta[b0][m] * (PA_g * PB_1 * PC[h] + PA_g * PB_h * PC[b1] + PB_1 * PB_h * PC[g])
                            + delta[g][m] * (PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[b0])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[b0][h] * (PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][h] * delta[g][m] * (PC[b0])
                            + delta[b0][h] * delta[g][m] * (PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m] + PC[b1] * PC[g] * PC[m])
                            + delta[b0][b1] * (PA_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                            + delta[b1][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m] + PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_g * PB_h * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_g * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g])
                            + delta[b1][m] * (PA_g * PC[b0] * PC[h] + PB_0 * PC[g] * PC[h] + PB_h * PC[b0] * PC[g])
                            + delta[b0][m] * (PA_g * PC[b1] * PC[h] + PB_1 * PC[g] * PC[h] + PB_h * PC[b1] * PC[g])
                            + delta[g][m] * (PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_h * PB_1 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_h)
                            + delta[b1][g] * delta[h][m] * (PB_0)
                            + delta[b0][g] * delta[h][m] * (PB_1)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_h * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_h * PB_g * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PB_0 * PB_g * PC[m])
                            + delta[b0][h] * (PB_1 * PB_g * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PB_g)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_0 * PB_1 * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_h * PB_0 * PB_1)
                            + delta[b1][m] * (PA_h * PB_0 * PB_g)
                            + delta[b0][m] * (PA_h * PB_1 * PB_g)
                            + delta[h][m] * (PB_0 * PB_1 * PB_g)
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + delta[b0][g] * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_h * PB_g * PC[m] + PA_h * PC[g] * PC[m] + PB_g * PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + delta[b1][h] * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + delta[b0][h] * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_h + PC[h])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1 + PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PB_g + PC[g])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                            + PA_h * PB_0 * PB_g * PC[b1] * PC[m]
                            + PA_h * PB_1 * PB_g * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[g][m] * (PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PB_0 * PB_1 * PC[h])
                            + delta[b1][m] * (PA_h * PB_0 * PC[g] + PA_h * PB_g * PC[b0] + PB_0 * PB_g * PC[h])
                            + delta[b0][m] * (PA_h * PB_1 * PC[g] + PA_h * PB_g * PC[b1] + PB_1 * PB_g * PC[h])
                            + delta[h][m] * (PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[b0])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + delta[b0][g] * (PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][g] * delta[h][m] * (PC[b0])
                            + delta[b0][g] * delta[h][m] * (PC[b1])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m] + PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_h * PC[g] * PC[m] + PB_g * PC[h] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m] + PC[b1] * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PA_h * PB_g * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_h * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h])
                            + delta[b1][m] * (PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h])
                            + delta[b0][m] * (PA_h * PC[b1] * PC[g] + PB_1 * PC[g] * PC[h] + PB_g * PC[b1] * PC[h])
                            + delta[h][m] * (PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PC[m] * (-1.0))
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m] * (-2.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PB_1 * PB_h * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1)
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PB_g)
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PB_h)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PB_1 * PB_h * PC[m])
                            + delta[b0][b1] * (PB_g * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PB_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PB_1 * PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PB_0 * PB_1 * PB_g)
                            + delta[g][m] * (PB_0 * PB_1 * PB_h)
                            + delta[b1][m] * (PB_0 * PB_g * PB_h)
                            + delta[b0][m] * (PB_1 * PB_g * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + delta[b1][h] * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + delta[b1][g] * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + delta[b0][h] * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m])
                            + delta[b0][g] * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m])
                            + delta[b0][b1] * (PB_g * PB_h * PC[m] + PB_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PB_1 + PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PB_g + PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PB_h + PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                            + PB_0 * PB_g * PB_h * PC[b1] * PC[m]
                            + PB_1 * PB_g * PB_h * PC[b0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[b0])
                            + delta[g][m] * (PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[b0])
                            + delta[b1][m] * (PB_0 * PB_g * PC[h] + PB_0 * PB_h * PC[g] + PB_g * PB_h * PC[b0])
                            + delta[b0][m] * (PB_1 * PB_g * PC[h] + PB_1 * PB_h * PC[g] + PB_g * PB_h * PC[b1])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PC[m])
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m] * 2.0)
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + delta[b1][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + delta[b1][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + delta[b0][h] * (PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m])
                            + delta[b0][g] * (PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m])
                        )

                        + 2.0 * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m] + PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m] + PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PB_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PB_g * PB_h * PC[b0] * PC[b1] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1])
                            + delta[g][m] * (PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1])
                            + delta[b1][m] * (PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h] + PB_h * PC[b0] * PC[g])
                            + delta[b0][m] * (PB_1 * PC[g] * PC[h] + PB_g * PC[b1] * PC[h] + PB_h * PC[b1] * PC[g])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_j * a_j * (
                            PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }


    // P-P block

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F5_t[6];

            onee::computeBoysFunction(F5_t, (a_i + a_j) * r2_PC, 5, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][m] * delta[g][h] * (PA_0)
                            + delta[a0][h] * delta[b0][m] * (PA_g)
                            + delta[a0][g] * delta[b0][m] * (PA_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][h] * (PA_0 * PA_g * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[a0][b0] * (PA_g * PA_h * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PA_h * PB_0 * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[b0][m] * (PA_0 * PA_g * PA_h)
                            + delta[h][m] * (PA_0 * PA_g * PB_0)
                            + delta[g][m] * (PA_0 * PA_h * PB_0)
                            + delta[a0][m] * (PA_g * PA_h * PB_0)
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][h] * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a0][b0] * (PA_g * PA_h * PC[m] + PA_g * PC[h] * PC[m] + PA_h * PC[g] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g + PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h + PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 + PC[b0])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PA_h * PC[b0] * PC[m]
                            + PA_0 * PA_g * PB_0 * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PC[g] * PC[m]
                            + PA_g * PA_h * PB_0 * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_g * PC[h] + PA_0 * PA_h * PC[g] + PA_g * PA_h * PC[a0])
                            + delta[g][m] * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + delta[a0][m] * (PA_g * PA_h * PC[b0] + PA_g * PB_0 * PC[h] + PA_h * PB_0 * PC[g])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b0][m] * delta[g][h] * (PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0])
                            + delta[a0][h] * delta[b0][m] * (PC[g])
                            + delta[a0][g] * delta[b0][m] * (PC[h])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PC[h] * PC[m] + PA_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m] + PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_h * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[g] * PC[h] * PC[m]
                            + PA_g * PA_h * PC[a0] * PC[b0] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[g] * PC[h] + PA_g * PC[a0] * PC[h] + PA_h * PC[a0] * PC[g])
                            + delta[a0][m] * (PA_g * PC[b0] * PC[h] + PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_i * (
                            PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[a0][g] * delta[b0][h] * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][g] * delta[b0][h] * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][h] * (PA_0 * PA_g * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PB_0 * PB_h * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][h] * delta[g][m] * (PA_0)
                            + delta[a0][m] * delta[b0][h] * (PA_g)
                        )

                        + (-2.0) * a_j * (
                            delta[a0][g] * delta[h][m] * (PB_0)
                            + delta[a0][g] * delta[b0][m] * (PB_h)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][h] * (PA_0 * PA_g * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PB_h * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[m])
                            + delta[a0][b0] * (PA_g * PB_h * PC[m])
                            + delta[a0][g] * (PB_0 * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PB_0 * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_g * PB_0)
                            + delta[b0][m] * (PA_0 * PA_g * PB_h)
                            + delta[g][m] * (PA_0 * PB_0 * PB_h)
                            + delta[a0][m] * (PA_g * PB_0 * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a0][g] * delta[b0][h] * (PC[m])
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][h] * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PB_h * PC[m] + PA_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g + PC[g])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_h + PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PB_0 * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_h * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[g] * PC[m]
                            + PA_g * PB_0 * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_g * PC[h] + PA_0 * PB_h * PC[g] + PA_g * PB_h * PC[a0])
                            + delta[g][m] * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + delta[a0][m] * (PA_g * PB_0 * PC[h] + PA_g * PB_h * PC[b0] + PB_0 * PB_h * PC[g])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[a0][g] * delta[b0][h] * (PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][g] * delta[b0][h] * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b0][h] * (PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b0][h] * delta[g][m] * (PC[a0])
                            + delta[a0][m] * delta[b0][h] * (PC[g])
                        )

                        + 2.0 * a_j * (
                            delta[a0][g] * delta[h][m] * (PC[b0])
                            + delta[a0][g] * delta[b0][m] * (PC[h])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m] + PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[h] * PC[m]
                            + PA_g * PB_h * PC[a0] * PC[b0] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + delta[b0][m] * (PA_0 * PC[g] * PC[h] + PA_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g])
                            + delta[a0][m] * (PA_g * PC[b0] * PC[h] + PB_0 * PC[g] * PC[h] + PB_h * PC[b0] * PC[g])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[a0][h] * delta[b0][g] * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][h] * delta[b0][g] * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][g] * (PA_0 * PA_h * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PB_0 * PB_g * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b0][g] * delta[h][m] * (PA_0)
                            + delta[a0][m] * delta[b0][g] * (PA_h)
                        )

                        + (-2.0) * a_j * (
                            delta[a0][h] * delta[g][m] * (PB_0)
                            + delta[a0][h] * delta[b0][m] * (PB_g)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][g] * (PA_0 * PA_h * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PB_g * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[m])
                            + delta[a0][b0] * (PA_h * PB_g * PC[m])
                            + delta[a0][h] * (PB_0 * PB_g * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_g)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PB_0 * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_h * PB_0)
                            + delta[b0][m] * (PA_0 * PA_h * PB_g)
                            + delta[h][m] * (PA_0 * PB_0 * PB_g)
                            + delta[a0][m] * (PA_h * PB_0 * PB_g)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a0][h] * delta[b0][g] * (PC[m])
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b0][g] * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_h * PB_g * PC[m] + PA_h * PC[g] * PC[m] + PB_g * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h + PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_g + PC[g])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PB_0 * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_g * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_g * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_g * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_h * PC[g] + PA_0 * PB_g * PC[h] + PA_h * PB_g * PC[a0])
                            + delta[h][m] * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + delta[a0][m] * (PA_h * PB_0 * PC[g] + PA_h * PB_g * PC[b0] + PB_0 * PB_g * PC[h])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[a0][h] * delta[b0][g] * (PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][h] * delta[b0][g] * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b0][g] * (PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b0][g] * delta[h][m] * (PC[a0])
                            + delta[a0][m] * delta[b0][g] * (PC[h])
                        )

                        + 2.0 * a_j * (
                            delta[a0][h] * delta[g][m] * (PC[b0])
                            + delta[a0][h] * delta[b0][m] * (PC[g])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_h * PC[g] * PC[m] + PB_g * PC[h] * PC[m] + PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m] + PC[b0] * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[g] * PC[m]
                            + PA_h * PB_g * PC[a0] * PC[b0] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[g] * PC[h] + PA_h * PC[a0] * PC[g] + PB_g * PC[a0] * PC[h])
                            + delta[a0][m] * (PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F5_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PB_h * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0)
                            + delta[a0][m] * delta[g][h] * (PB_0)
                            + delta[a0][m] * delta[b0][h] * (PB_g)
                            + delta[a0][m] * delta[b0][g] * (PB_h)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PB_h * PC[m])
                            + delta[a0][h] * (PB_0 * PB_g * PC[m])
                            + delta[a0][g] * (PB_0 * PB_h * PC[m])
                            + delta[a0][b0] * (PB_g * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_g)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PB_0 * PB_g)
                            + delta[g][m] * (PA_0 * PB_0 * PB_h)
                            + delta[b0][m] * (PA_0 * PB_g * PB_h)
                            + delta[a0][m] * (PB_0 * PB_g * PB_h)
                        )

                    )

                    + F5_t[2] * (

                        (-4.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + delta[a0][b0] * (PB_g * PB_h * PC[m] + PB_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 + PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_g + PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_h + PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PB_g * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[g] * PC[m]
                            + PA_0 * PB_g * PB_h * PC[b0] * PC[m]
                            + PB_0 * PB_g * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + delta[g][m] * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + delta[b0][m] * (PA_0 * PB_g * PC[h] + PA_0 * PB_h * PC[g] + PB_g * PB_h * PC[a0])
                            + delta[a0][m] * (PB_0 * PB_g * PC[h] + PB_0 * PB_h * PC[g] + PB_g * PB_h * PC[b0])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m])
                        )

                        + 2.0 * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0])
                            + delta[a0][m] * delta[g][h] * (PC[b0])
                            + delta[a0][m] * delta[b0][h] * (PC[g])
                            + delta[a0][m] * delta[b0][g] * (PC[h])
                        )

                    )

                    + F5_t[3] * (

                        (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PB_g * PC[h] * PC[m] + PB_h * PC[g] * PC[m] + PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[g] * PC[m]
                            + PB_g * PB_h * PC[a0] * PC[b0] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + delta[b0][m] * (PA_0 * PC[g] * PC[h] + PB_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g])
                            + delta[a0][m] * (PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h] + PB_h * PC[b0] * PC[g])
                        )

                    )

                    + F5_t[4] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[g] * PC[h])
                        )

                    )

                    + F5_t[5] * (

                        8.0 * (a_i + a_j) * a_j * a_j * (
                            PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }


    // P-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F6_t[7];

            onee::computeBoysFunction(F6_t, (a_i + a_j) * r2_PC, 6, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F6_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_g * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PB_1 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][m] * delta[g][h] * (PA_0 * PB_0)
                            + delta[b0][m] * delta[g][h] * (PA_0 * PB_1)
                            + delta[a0][h] * delta[b1][m] * (PA_g * PB_0)
                            + delta[a0][h] * delta[b0][m] * (PA_g * PB_1)
                            + delta[a0][g] * delta[b1][m] * (PA_h * PB_0)
                            + delta[a0][g] * delta[b0][m] * (PA_h * PB_1)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_0 * PA_g * PA_h * PC[m])
                            + delta[b1][h] * (PA_0 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PB_1 * PC[m])
                            + delta[b1][g] * (PA_0 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PB_1 * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_g * PA_h * PB_0 * PC[m])
                            + delta[a0][b0] * (PA_g * PA_h * PB_1 * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PB_1 * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PA_h)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PB_1)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PB_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PB_1)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PA_h * PB_0 * PB_1 * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[b1][m] * (PA_0 * PA_g * PA_h * PB_0)
                            + delta[b0][m] * (PA_0 * PA_g * PA_h * PB_1)
                            + delta[h][m] * (PA_0 * PA_g * PB_0 * PB_1)
                            + delta[g][m] * (PA_0 * PA_h * PB_0 * PB_1)
                            + delta[a0][m] * (PA_g * PA_h * PB_0 * PB_1)
                        )

                    )

                    + F6_t[2] * (

                        2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * (-2.0) + PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m] * (-2.0) + PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m] * (-2.0) + PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] * (-2.0) + PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] * (-2.0) + PC[b1] * PC[m] * (-1.0))
                        )

                        + (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_0 * PA_g * PA_h * PC[m] + PA_0 * PA_g * PC[h] * PC[m] + PA_0 * PA_h * PC[g] * PC[m] + PA_g * PA_h * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_g * PB_0 * PC[m] + PA_0 * PA_g * PC[b0] * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PB_1 * PC[m] + PA_0 * PA_g * PC[b1] * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_h * PB_0 * PC[m] + PA_0 * PA_h * PC[b0] * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PB_1 * PC[m] + PA_0 * PA_h * PC[b1] * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m] + PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a0][b1] * (PA_g * PA_h * PB_0 * PC[m] + PA_g * PA_h * PC[b0] * PC[m] + PA_g * PB_0 * PC[h] * PC[m] + PA_h * PB_0 * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PA_h * PB_1 * PC[m] + PA_g * PA_h * PC[b1] * PC[m] + PA_g * PB_1 * PC[h] * PC[m] + PA_h * PB_1 * PC[g] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PB_1 * PC[m] + PA_g * PB_0 * PC[b1] * PC[m] + PA_g * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PB_1 * PC[m] + PA_h * PB_0 * PC[b1] * PC[m] + PA_h * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_g + PA_0 * PC[g] + PA_g * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_h + PA_0 * PC[h] + PA_h * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PA_h + PA_g * PC[h] + PA_h * PC[g])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PB_0 + PA_g * PC[b0] + PB_0 * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PB_1 + PA_g * PC[b1] + PB_1 * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PB_0 + PA_h * PC[b0] + PB_0 * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PB_1 + PA_h * PC[b1] + PB_1 * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PA_h * PB_0 * PC[b1] * PC[m]
                            + PA_0 * PA_g * PA_h * PB_1 * PC[b0] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                            + PA_g * PA_h * PB_0 * PB_1 * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[b1][m] * (PA_0 * PA_g * PA_h * PC[b0] + PA_0 * PA_g * PB_0 * PC[h] + PA_0 * PA_h * PB_0 * PC[g] + PA_g * PA_h * PB_0 * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_g * PA_h * PC[b1] + PA_0 * PA_g * PB_1 * PC[h] + PA_0 * PA_h * PB_1 * PC[g] + PA_g * PA_h * PB_1 * PC[a0])
                            + delta[h][m] * (PA_0 * PA_g * PB_0 * PC[b1] + PA_0 * PA_g * PB_1 * PC[b0] + PA_0 * PB_0 * PB_1 * PC[g] + PA_g * PB_0 * PB_1 * PC[a0])
                            + delta[g][m] * (PA_0 * PA_h * PB_0 * PC[b1] + PA_0 * PA_h * PB_1 * PC[b0] + PA_0 * PB_0 * PB_1 * PC[h] + PA_h * PB_0 * PB_1 * PC[a0])
                            + delta[a0][m] * (PA_g * PA_h * PB_0 * PC[b1] + PA_g * PA_h * PB_1 * PC[b0] + PA_g * PB_0 * PB_1 * PC[h] + PA_h * PB_0 * PB_1 * PC[g])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[m] + PC[a0] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_g * PC[m] + PC[g] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_h * PC[m] + PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[b1] * PC[m] + PA_g * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[b1] * PC[m] + PA_h * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[h] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][m] * delta[g][h] * (PA_0 * PC[b0] + PB_0 * PC[a0])
                            + delta[b0][m] * delta[g][h] * (PA_0 * PC[b1] + PB_1 * PC[a0])
                            + delta[a0][h] * delta[b1][m] * (PA_g * PC[b0] + PB_0 * PC[g])
                            + delta[a0][h] * delta[b0][m] * (PA_g * PC[b1] + PB_1 * PC[g])
                            + delta[a0][g] * delta[b1][m] * (PA_h * PC[b0] + PB_0 * PC[h])
                            + delta[a0][g] * delta[b0][m] * (PA_h * PC[b1] + PB_1 * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PC[b1] + PB_1 * PC[b0])
                        )

                    )

                    + F6_t[3] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PC[g] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PC[h] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][m] * delta[g][h] * (PC[a0] * PC[b0])
                            + delta[b0][m] * delta[g][h] * (PC[a0] * PC[b1])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0] * PC[b1])
                            + delta[a0][h] * delta[b1][m] * (PC[b0] * PC[g])
                            + delta[a0][g] * delta[b1][m] * (PC[b0] * PC[h])
                            + delta[a0][h] * delta[b0][m] * (PC[b1] * PC[g])
                            + delta[a0][g] * delta[b0][m] * (PC[b1] * PC[h])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m] + PC[g] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m] + PC[h] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m] * 2.0)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b1][h] * (PA_0 * PA_g * PC[b0] * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PC[b1] * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_g * PC[h] * PC[m] + PA_0 * PA_h * PC[g] * PC[m] + PA_0 * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a0] * PC[m] + PA_g * PC[a0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_h * PC[b0] * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PC[b1] * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][b1] * (PA_g * PA_h * PC[b0] * PC[m] + PA_g * PB_0 * PC[h] * PC[m] + PA_g * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[g] * PC[m] + PA_h * PC[b0] * PC[g] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_g * PA_h * PC[b1] * PC[m] + PA_g * PB_1 * PC[h] * PC[m] + PA_g * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[g] * PC[m] + PA_h * PC[b1] * PC[g] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[b1] * PC[m] + PA_g * PB_1 * PC[b0] * PC[m] + PA_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[b1] * PC[m] + PA_h * PB_1 * PC[b0] * PC[m] + PA_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PA_g * PC[a0] + PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PA_h * PC[a0] + PC[a0] * PC[h])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PC[b0] + PB_0 * PC[g] + PC[b0] * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PC[b1] + PB_1 * PC[g] + PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PC[h] + PA_h * PC[g] + PC[g] * PC[h])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PC[b0] + PB_0 * PC[h] + PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PC[b1] + PB_1 * PC[h] + PC[b1] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PA_h * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_g * PA_h * PB_0 * PC[a0] * PC[b1] * PC[m]
                            + PA_g * PA_h * PB_1 * PC[a0] * PC[b0] * PC[m]
                            + PA_g * PB_0 * PB_1 * PC[a0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PC[a0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_0 * PA_g * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[b1] * PC[g] + PA_0 * PB_1 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[g])
                            + delta[b1][m] * (PA_0 * PA_g * PC[b0] * PC[h] + PA_0 * PA_h * PC[b0] * PC[g] + PA_0 * PB_0 * PC[g] * PC[h] + PA_g * PA_h * PC[a0] * PC[b0] + PA_g * PB_0 * PC[a0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[g])
                            + delta[b0][m] * (PA_0 * PA_g * PC[b1] * PC[h] + PA_0 * PA_h * PC[b1] * PC[g] + PA_0 * PB_1 * PC[g] * PC[h] + PA_g * PA_h * PC[a0] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[h] + PA_h * PB_1 * PC[a0] * PC[g])
                            + delta[g][m] * (PA_0 * PA_h * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[b1] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[b1] + PA_h * PB_1 * PC[a0] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[h])
                            + delta[a0][m] * (PA_g * PA_h * PC[b0] * PC[b1] + PA_g * PB_0 * PC[b1] * PC[h] + PA_g * PB_1 * PC[b0] * PC[h] + PA_h * PB_0 * PC[b1] * PC[g] + PA_h * PB_1 * PC[b0] * PC[g] + PB_0 * PB_1 * PC[g] * PC[h])
                        )

                    )

                    + F6_t[4] * (

                        (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[h] * PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_g * PC[b0] * PC[h] * PC[m] + PA_h * PC[b0] * PC[g] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_g * PC[b1] * PC[h] * PC[m] + PA_h * PC[b1] * PC[g] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[g] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_g * PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[h])
                            + delta[b1][m] * (PA_0 * PC[b0] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PA_g * PC[b0] * PC[b1] * PC[h] + PA_h * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[b1] * PC[g] * PC[h] + PB_1 * PC[b0] * PC[g] * PC[h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                    )

                    + F6_t[5] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F6_t[6] * (

                        (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F6_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_g * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PB_1 * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][g] * delta[b1][h] * (PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PB_1 * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PB_h * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PB_1 * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PB_0 * PB_1 * PB_h * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_g)
                            + delta[b1][h] * delta[g][m] * (PA_0 * PB_0)
                            + delta[b0][h] * delta[g][m] * (PA_0 * PB_1)
                            + delta[a0][m] * delta[b1][h] * (PA_g * PB_0)
                            + delta[a0][m] * delta[b0][h] * (PA_g * PB_1)
                        )

                        + (-2.0) * a_j * (
                            delta[a0][g] * delta[h][m] * (PB_0 * PB_1)
                            + delta[a0][g] * delta[b1][m] * (PB_0 * PB_h)
                            + delta[a0][g] * delta[b0][m] * (PB_1 * PB_h)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_g * PB_h * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PB_h * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_g * PB_0 * PB_h * PC[m])
                            + delta[a0][b0] * (PA_g * PB_1 * PB_h * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_g)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PB_h)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PB_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PB_h)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PB_0 * PB_1 * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_g * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_g * PB_0 * PB_h)
                            + delta[b0][m] * (PA_0 * PA_g * PB_1 * PB_h)
                            + delta[g][m] * (PA_0 * PB_0 * PB_1 * PB_h)
                            + delta[a0][m] * (PA_g * PB_0 * PB_1 * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a0][g] * delta[b1][h] * (PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PB_1 * PC[m])
                        )

                        + (
                            (delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h])
                        )

                    )

                    + F6_t[2] * (

                        2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * (-2.0) + PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m] * (-2.0) + PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] * (-2.0) + PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] * (-2.0) + PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m] * (-2.0) + PC[h] * PC[m] * (-1.0))
                        )

                        + (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_g * PB_0 * PC[m] + PA_0 * PA_g * PC[b0] * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PB_1 * PC[m] + PA_0 * PA_g * PC[b1] * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_g * PB_h * PC[m] + PA_0 * PA_g * PC[h] * PC[m] + PA_0 * PB_h * PC[g] * PC[m] + PA_g * PB_h * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m] + PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PB_h * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PB_h * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PB_1 * PC[m] + PA_g * PB_0 * PC[b1] * PC[m] + PA_g * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_g * PB_0 * PB_h * PC[m] + PA_g * PB_0 * PC[h] * PC[m] + PA_g * PB_h * PC[b0] * PC[m] + PB_0 * PB_h * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PB_1 * PB_h * PC[m] + PA_g * PB_1 * PC[h] * PC[m] + PA_g * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[g] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PB_h * PC[m] + PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_g + PA_0 * PC[g] + PA_g * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PB_h + PA_0 * PC[h] + PB_h * PC[a0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PB_0 + PA_g * PC[b0] + PB_0 * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PB_1 + PA_g * PC[b1] + PB_1 * PC[g])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PB_h + PA_g * PC[h] + PB_h * PC[g])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PB_h + PB_0 * PC[h] + PB_h * PC[b0])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PB_h + PB_1 * PC[h] + PB_h * PC[b1])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_h * PC[b1] * PC[m]
                            + PA_0 * PA_g * PB_1 * PB_h * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                            + PA_g * PB_0 * PB_1 * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_g * PB_0 * PC[b1] + PA_0 * PA_g * PB_1 * PC[b0] + PA_0 * PB_0 * PB_1 * PC[g] + PA_g * PB_0 * PB_1 * PC[a0])
                            + delta[b1][m] * (PA_0 * PA_g * PB_0 * PC[h] + PA_0 * PA_g * PB_h * PC[b0] + PA_0 * PB_0 * PB_h * PC[g] + PA_g * PB_0 * PB_h * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_g * PB_1 * PC[h] + PA_0 * PA_g * PB_h * PC[b1] + PA_0 * PB_1 * PB_h * PC[g] + PA_g * PB_1 * PB_h * PC[a0])
                            + delta[g][m] * (PA_0 * PB_0 * PB_1 * PC[h] + PA_0 * PB_0 * PB_h * PC[b1] + PA_0 * PB_1 * PB_h * PC[b0] + PB_0 * PB_1 * PB_h * PC[a0])
                            + delta[a0][m] * (PA_g * PB_0 * PB_1 * PC[h] + PA_g * PB_0 * PB_h * PC[b1] + PA_g * PB_1 * PB_h * PC[b0] + PB_0 * PB_1 * PB_h * PC[g])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            delta[a0][g] * delta[b1][h] * (PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_g * PC[m] + PC[g] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PB_1 * PC[m] + PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][g] * delta[b1][h] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PB_1 * PC[m] + PC[b1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PB_h * PC[m] + PC[h] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PA_g * PC[b0] * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PC[b1] * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[b0] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][h] * delta[g][m] * (PA_0 * PC[b0] + PB_0 * PC[a0])
                            + delta[b0][h] * delta[g][m] * (PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PA_g * PC[a0])
                            + delta[a0][m] * delta[b1][h] * (PA_g * PC[b0] + PB_0 * PC[g])
                            + delta[a0][m] * delta[b0][h] * (PA_g * PC[b1] + PB_1 * PC[g])
                        )

                        + 2.0 * a_j * (
                            delta[a0][g] * delta[h][m] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                            + delta[a0][g] * delta[b1][m] * (PB_0 * PC[h] + PB_h * PC[b0])
                            + delta[a0][g] * delta[b0][m] * (PB_1 * PC[h] + PB_h * PC[b1])
                        )

                    )

                    + F6_t[3] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][g] * delta[b1][h] * (PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[b1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PC[h] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][h] * delta[g][m] * (PC[a0] * PC[b0])
                            + delta[b0][h] * delta[g][m] * (PC[a0] * PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + delta[a0][m] * delta[b1][h] * (PC[b0] * PC[g])
                            + delta[a0][m] * delta[b0][h] * (PC[b1] * PC[g])
                        )

                        + (-2.0) * a_j * (
                            delta[a0][g] * delta[h][m] * (PC[b0] * PC[b1])
                            + delta[a0][g] * delta[b1][m] * (PC[b0] * PC[h])
                            + delta[a0][g] * delta[b0][m] * (PC[b1] * PC[h])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_g * PC[m] + PC[g] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m] + PC[h] * PC[m] * 2.0)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_g * PC[b0] * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_g * PC[b1] * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_g * PC[h] * PC[m] + PA_0 * PB_h * PC[g] * PC[m] + PA_0 * PC[g] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[m] + PA_g * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[g] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[m] + PA_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[m])
                            + delta[a0][h] * (PA_g * PB_0 * PC[b1] * PC[m] + PA_g * PB_1 * PC[b0] * PC[m] + PA_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_g * PB_0 * PC[h] * PC[m] + PA_g * PB_h * PC[b0] * PC[m] + PA_g * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[g] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PB_h * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_g * PB_1 * PC[h] * PC[m] + PA_g * PB_h * PC[b1] * PC[m] + PA_g * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[g] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PB_h * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[b0] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PA_g * PC[a0] + PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PB_h * PC[a0] + PC[a0] * PC[h])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_g * PC[b0] + PB_0 * PC[g] + PC[b0] * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_g * PC[b1] + PB_1 * PC[g] + PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_g * PC[h] + PB_h * PC[g] + PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PC[h] + PB_h * PC[b0] + PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PC[h] + PB_h * PC[b1] + PC[b1] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_h * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PA_g * PB_0 * PB_1 * PC[a0] * PC[h] * PC[m]
                            + PA_g * PB_0 * PB_h * PC[a0] * PC[b1] * PC[m]
                            + PA_g * PB_1 * PB_h * PC[a0] * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[a0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_g * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[b1] * PC[g] + PA_0 * PB_1 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[g])
                            + delta[b1][m] * (PA_0 * PA_g * PC[b0] * PC[h] + PA_0 * PB_0 * PC[g] * PC[h] + PA_0 * PB_h * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[h] + PA_g * PB_h * PC[a0] * PC[b0] + PB_0 * PB_h * PC[a0] * PC[g])
                            + delta[b0][m] * (PA_0 * PA_g * PC[b1] * PC[h] + PA_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_h * PC[b1] * PC[g] + PA_g * PB_1 * PC[a0] * PC[h] + PA_g * PB_h * PC[a0] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[g])
                            + delta[g][m] * (PA_0 * PB_0 * PC[b1] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[h] + PA_0 * PB_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[h] + PB_0 * PB_h * PC[a0] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[b0])
                            + delta[a0][m] * (PA_g * PB_0 * PC[b1] * PC[h] + PA_g * PB_1 * PC[b0] * PC[h] + PA_g * PB_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[g] * PC[h] + PB_0 * PB_h * PC[b1] * PC[g] + PB_1 * PB_h * PC[b0] * PC[g])
                        )

                    )

                    + F6_t[4] * (

                        (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[h] * PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[m] + PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_g * PC[b0] * PC[h] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PB_h * PC[b0] * PC[g] * PC[m] + PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_g * PC[b1] * PC[h] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PB_h * PC[b1] * PC[g] * PC[m] + PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[m] + PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[g] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_g * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[b0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a0] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[b1])
                            + delta[b1][m] * (PA_0 * PC[b0] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[g] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[g])
                            + delta[b0][m] * (PA_0 * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[g] * PC[h] + PB_h * PC[a0] * PC[b1] * PC[g])
                            + delta[a0][m] * (PA_g * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[b1] * PC[g] * PC[h] + PB_1 * PC[b0] * PC[g] * PC[h] + PB_h * PC[b0] * PC[b1] * PC[g])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                    )

                    + F6_t[5] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F6_t[6] * (

                        (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F6_t[1] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_h * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PB_1 * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][h] * delta[b1][g] * (PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PB_1 * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PB_g * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PB_1 * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PB_0 * PB_1 * PB_g * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_h)
                            + delta[b1][g] * delta[h][m] * (PA_0 * PB_0)
                            + delta[b0][g] * delta[h][m] * (PA_0 * PB_1)
                            + delta[a0][m] * delta[b1][g] * (PA_h * PB_0)
                            + delta[a0][m] * delta[b0][g] * (PA_h * PB_1)
                        )

                        + (-2.0) * a_j * (
                            delta[a0][h] * delta[g][m] * (PB_0 * PB_1)
                            + delta[a0][h] * delta[b1][m] * (PB_0 * PB_g)
                            + delta[a0][h] * delta[b0][m] * (PB_1 * PB_g)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_h * PB_g * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PB_g * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_h * PB_0 * PB_g * PC[m])
                            + delta[a0][b0] * (PA_h * PB_1 * PB_g * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PB_g * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PB_g)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PB_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_h * PB_g)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PB_g)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PB_g)
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PB_0 * PB_1 * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_h * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_h * PB_0 * PB_g)
                            + delta[b0][m] * (PA_0 * PA_h * PB_1 * PB_g)
                            + delta[h][m] * (PA_0 * PB_0 * PB_1 * PB_g)
                            + delta[a0][m] * (PA_h * PB_0 * PB_1 * PB_g)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a0][h] * delta[b1][g] * (PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PB_1 * PC[m])
                        )

                        + (
                            (delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g])
                        )

                    )

                    + F6_t[2] * (

                        2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * (-2.0) + PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m] * (-2.0) + PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] * (-2.0) + PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] * (-2.0) + PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m] * (-2.0) + PC[g] * PC[m] * (-1.0))
                        )

                        + (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_h * PB_0 * PC[m] + PA_0 * PA_h * PC[b0] * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PB_1 * PC[m] + PA_0 * PA_h * PC[b1] * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_h * PB_g * PC[m] + PA_0 * PA_h * PC[g] * PC[m] + PA_0 * PB_g * PC[h] * PC[m] + PA_h * PB_g * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m] + PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PB_g * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PB_g * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PB_1 * PC[m] + PA_h * PB_0 * PC[b1] * PC[m] + PA_h * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_h * PB_0 * PB_g * PC[m] + PA_h * PB_0 * PC[g] * PC[m] + PA_h * PB_g * PC[b0] * PC[m] + PB_0 * PB_g * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_h * PB_1 * PB_g * PC[m] + PA_h * PB_1 * PC[g] * PC[m] + PA_h * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PB_g * PC[m] + PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_h + PA_0 * PC[h] + PA_h * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PB_g + PA_0 * PC[g] + PB_g * PC[a0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PB_0 + PA_h * PC[b0] + PB_0 * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PB_1 + PA_h * PC[b1] + PB_1 * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_h * PB_g + PA_h * PC[g] + PB_g * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PB_g + PB_0 * PC[g] + PB_g * PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PB_g + PB_1 * PC[g] + PB_g * PC[b1])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_g * PC[b1] * PC[m]
                            + PA_0 * PA_h * PB_1 * PB_g * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PB_g * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_h * PB_0 * PC[b1] + PA_0 * PA_h * PB_1 * PC[b0] + PA_0 * PB_0 * PB_1 * PC[h] + PA_h * PB_0 * PB_1 * PC[a0])
                            + delta[b1][m] * (PA_0 * PA_h * PB_0 * PC[g] + PA_0 * PA_h * PB_g * PC[b0] + PA_0 * PB_0 * PB_g * PC[h] + PA_h * PB_0 * PB_g * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_h * PB_1 * PC[g] + PA_0 * PA_h * PB_g * PC[b1] + PA_0 * PB_1 * PB_g * PC[h] + PA_h * PB_1 * PB_g * PC[a0])
                            + delta[h][m] * (PA_0 * PB_0 * PB_1 * PC[g] + PA_0 * PB_0 * PB_g * PC[b1] + PA_0 * PB_1 * PB_g * PC[b0] + PB_0 * PB_1 * PB_g * PC[a0])
                            + delta[a0][m] * (PA_h * PB_0 * PB_1 * PC[g] + PA_h * PB_0 * PB_g * PC[b1] + PA_h * PB_1 * PB_g * PC[b0] + PB_0 * PB_1 * PB_g * PC[h])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            delta[a0][h] * delta[b1][g] * (PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_h * PC[m] + PC[h] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PB_1 * PC[m] + PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][h] * delta[b1][g] * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PB_1 * PC[m] + PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PB_g * PC[m] + PC[g] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PA_h * PC[b0] * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PC[b1] * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[b0] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][g] * delta[h][m] * (PA_0 * PC[b0] + PB_0 * PC[a0])
                            + delta[b0][g] * delta[h][m] * (PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PA_h * PC[a0])
                            + delta[a0][m] * delta[b1][g] * (PA_h * PC[b0] + PB_0 * PC[h])
                            + delta[a0][m] * delta[b0][g] * (PA_h * PC[b1] + PB_1 * PC[h])
                        )

                        + 2.0 * a_j * (
                            delta[a0][h] * delta[g][m] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                            + delta[a0][h] * delta[b1][m] * (PB_0 * PC[g] + PB_g * PC[b0])
                            + delta[a0][h] * delta[b0][m] * (PB_1 * PC[g] + PB_g * PC[b1])
                        )

                    )

                    + F6_t[3] * (

                        (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a0][h] * delta[b1][g] * (PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PC[g] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][g] * delta[h][m] * (PC[a0] * PC[b0])
                            + delta[b0][g] * delta[h][m] * (PC[a0] * PC[b1])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + delta[a0][m] * delta[b1][g] * (PC[b0] * PC[h])
                            + delta[a0][m] * delta[b0][g] * (PC[b1] * PC[h])
                        )

                        + (-2.0) * a_j * (
                            delta[a0][h] * delta[g][m] * (PC[b0] * PC[b1])
                            + delta[a0][h] * delta[b1][m] * (PC[b0] * PC[g])
                            + delta[a0][h] * delta[b0][m] * (PC[b1] * PC[g])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_h * PC[m] + PC[h] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m] + PC[g] * PC[m] * 2.0)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_h * PC[b0] * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_h * PC[b1] * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_h * PC[g] * PC[m] + PA_0 * PB_g * PC[h] * PC[m] + PA_0 * PC[g] * PC[h] * PC[m] + PA_h * PB_g * PC[a0] * PC[m] + PA_h * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[h] * PC[m])
                            + delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_g * PC[b1] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[m])
                            + delta[a0][g] * (PA_h * PB_0 * PC[b1] * PC[m] + PA_h * PB_1 * PC[b0] * PC[m] + PA_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_h * PB_0 * PC[g] * PC[m] + PA_h * PB_g * PC[b0] * PC[m] + PA_h * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[h] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_h * PB_1 * PC[g] * PC[m] + PA_h * PB_g * PC[b1] * PC[m] + PA_h * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[h] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PB_g * PC[b1] * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PB_g * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[b0] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PB_g * PC[a0] + PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PA_h * PC[a0] + PC[a0] * PC[h])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_h * PC[b0] + PB_0 * PC[h] + PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_h * PC[b1] + PB_1 * PC[h] + PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_h * PC[g] + PB_g * PC[h] + PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PC[g] + PB_g * PC[b0] + PC[b0] * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PC[g] + PB_g * PC[b1] + PC[b1] * PC[g])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_g * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PC[a0] * PC[g] * PC[m]
                            + PA_h * PB_0 * PB_g * PC[a0] * PC[b1] * PC[m]
                            + PA_h * PB_1 * PB_g * PC[a0] * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_g * PC[a0] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_h * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[b1] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[b1] + PA_h * PB_1 * PC[a0] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[h])
                            + delta[b1][m] * (PA_0 * PA_h * PC[b0] * PC[g] + PA_0 * PB_0 * PC[g] * PC[h] + PA_0 * PB_g * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[g] + PA_h * PB_g * PC[a0] * PC[b0] + PB_0 * PB_g * PC[a0] * PC[h])
                            + delta[b0][m] * (PA_0 * PA_h * PC[b1] * PC[g] + PA_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_g * PC[b1] * PC[h] + PA_h * PB_1 * PC[a0] * PC[g] + PA_h * PB_g * PC[a0] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[h])
                            + delta[h][m] * (PA_0 * PB_0 * PC[b1] * PC[g] + PA_0 * PB_1 * PC[b0] * PC[g] + PA_0 * PB_g * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[g] + PB_0 * PB_g * PC[a0] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[b0])
                            + delta[a0][m] * (PA_h * PB_0 * PC[b1] * PC[g] + PA_h * PB_1 * PC[b0] * PC[g] + PA_h * PB_g * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[g] * PC[h] + PB_0 * PB_g * PC[b1] * PC[h] + PB_1 * PB_g * PC[b0] * PC[h])
                        )

                    )

                    + F6_t[4] * (

                        (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[h] * PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[g] * PC[h] * PC[m] + PA_h * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[h] * PC[m] + PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PA_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_h * PC[b0] * PC[g] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[h] * PC[m] + PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_h * PC[b1] * PC[g] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PB_g * PC[b1] * PC[h] * PC[m] + PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[g] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PA_h * PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PB_1 * PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[b0] * PC[g] + PB_g * PC[a0] * PC[b0] * PC[b1])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[h])
                            + delta[b1][m] * (PA_0 * PC[b0] * PC[g] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b0] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[b1] * PC[g] * PC[h] + PA_h * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b1] * PC[h])
                            + delta[a0][m] * (PA_h * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[b1] * PC[g] * PC[h] + PB_1 * PC[b0] * PC[g] * PC[h] + PB_g * PC[b0] * PC[b1] * PC[h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                        )

                    )

                    + F6_t[5] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F6_t[6] * (

                        (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F6_t[1] * (

                        2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[m] * (-1.0))
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PB_g * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PB_h * PC[m] * (-1.0))
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (-1.0)
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (-2.0)
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PB_h * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1)
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PB_g)
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PB_h)
                            + delta[a0][m] * delta[g][h] * (PB_0 * PB_1)
                            + delta[a0][m] * delta[b1][h] * (PB_0 * PB_g)
                            + delta[a0][m] * delta[b1][g] * (PB_0 * PB_h)
                            + delta[a0][m] * delta[b0][h] * (PB_1 * PB_g)
                            + delta[a0][m] * delta[b0][g] * (PB_1 * PB_h)
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PB_h * PC[m])
                            + delta[b0][b1] * (PA_0 * PB_g * PB_h * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][b1] * (PB_0 * PB_g * PB_h * PC[m])
                            + delta[a0][b0] * (PB_1 * PB_g * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PB_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PB_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PB_g)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PB_h)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PB_g)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PB_h)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_g * PB_h)
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PB_1 * PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PB_0 * PB_1 * PB_g)
                            + delta[g][m] * (PA_0 * PB_0 * PB_1 * PB_h)
                            + delta[b1][m] * (PA_0 * PB_0 * PB_g * PB_h)
                            + delta[b0][m] * (PA_0 * PB_1 * PB_g * PB_h)
                            + delta[a0][m] * (PB_0 * PB_1 * PB_g * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m])
                        )

                        + (
                            (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                    )

                    + F6_t[2] * (

                        2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * (-2.0) + PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] * (-2.0) + PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] * (-2.0) + PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m] * (-2.0) + PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m] * (-2.0) + PC[h] * PC[m] * (-1.0))
                        )

                        + (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PB_1 * PC[m] + PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PB_g * PC[m] + PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PB_h * PC[m] + PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PB_g * PC[m] + PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PB_h * PC[m] + PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PB_g * PB_h * PC[m] + PA_0 * PB_g * PC[h] * PC[m] + PA_0 * PB_h * PC[g] * PC[m] + PB_g * PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PB_g * PC[m] + PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[b0] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PB_h * PC[m] + PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[b0] * PC[m])
                            + delta[a0][b1] * (PB_0 * PB_g * PB_h * PC[m] + PB_0 * PB_g * PC[h] * PC[m] + PB_0 * PB_h * PC[g] * PC[m] + PB_g * PB_h * PC[b0] * PC[m])
                            + delta[a0][b0] * (PB_1 * PB_g * PB_h * PC[m] + PB_1 * PB_g * PC[h] * PC[m] + PB_1 * PB_h * PC[g] * PC[m] + PB_g * PB_h * PC[b1] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PB_g + PA_0 * PC[g] + PB_g * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PB_h + PA_0 * PC[h] + PB_h * PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PB_g + PB_0 * PC[g] + PB_g * PC[b0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PB_h + PB_0 * PC[h] + PB_h * PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PB_g + PB_1 * PC[g] + PB_g * PC[b1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PB_h + PB_1 * PC[h] + PB_h * PC[b1])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_g * PB_h + PB_g * PC[h] + PB_h * PC[g])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                            + PA_0 * PB_0 * PB_g * PB_h * PC[b1] * PC[m]
                            + PA_0 * PB_1 * PB_g * PB_h * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_g * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PA_0 * PB_0 * PB_1 * PC[g] + PA_0 * PB_0 * PB_g * PC[b1] + PA_0 * PB_1 * PB_g * PC[b0] + PB_0 * PB_1 * PB_g * PC[a0])
                            + delta[g][m] * (PA_0 * PB_0 * PB_1 * PC[h] + PA_0 * PB_0 * PB_h * PC[b1] + PA_0 * PB_1 * PB_h * PC[b0] + PB_0 * PB_1 * PB_h * PC[a0])
                            + delta[b1][m] * (PA_0 * PB_0 * PB_g * PC[h] + PA_0 * PB_0 * PB_h * PC[g] + PA_0 * PB_g * PB_h * PC[b0] + PB_0 * PB_g * PB_h * PC[a0])
                            + delta[b0][m] * (PA_0 * PB_1 * PB_g * PC[h] + PA_0 * PB_1 * PB_h * PC[g] + PA_0 * PB_g * PB_h * PC[b1] + PB_1 * PB_g * PB_h * PC[a0])
                            + delta[a0][m] * (PB_0 * PB_1 * PB_g * PC[h] + PB_0 * PB_1 * PB_h * PC[g] + PB_0 * PB_g * PB_h * PC[b1] + PB_1 * PB_g * PB_h * PC[b0])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[m] + PC[a0] * PC[m])
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] * 2.0 + PC[a0] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PB_g * PC[m] + PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PB_h * PC[m] + PC[h] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h])
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * 2.0
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_g * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[m])
                        )

                        + 2.0 * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[b0] + PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[b1] + PB_1 * PC[a0])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PB_g * PC[a0])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PB_h * PC[a0])
                            + delta[a0][m] * delta[g][h] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                            + delta[a0][m] * delta[b1][h] * (PB_0 * PC[g] + PB_g * PC[b0])
                            + delta[a0][m] * delta[b1][g] * (PB_0 * PC[h] + PB_h * PC[b0])
                            + delta[a0][m] * delta[b0][h] * (PB_1 * PC[g] + PB_g * PC[b1])
                            + delta[a0][m] * delta[b0][g] * (PB_1 * PC[h] + PB_h * PC[b1])
                        )

                    )

                    + F6_t[3] * (

                        2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PC[a0] * PC[m] * (-1.0))
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PC[h] * PC[m] * (-1.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + delta[a0][m] * delta[g][h] * (PC[b0] * PC[b1])
                            + delta[a0][m] * delta[b1][h] * (PC[b0] * PC[g])
                            + delta[a0][m] * delta[b1][g] * (PC[b0] * PC[h])
                            + delta[a0][m] * delta[b0][h] * (PC[b1] * PC[g])
                            + delta[a0][m] * delta[b0][g] * (PC[b1] * PC[h])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[m] + PC[a0] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PB_0 * PC[m] + PC[b0] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PB_1 * PC[m] + PC[b1] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PB_g * PC[m] + PC[g] * PC[m] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PB_h * PC[m] + PC[h] * PC[m] * 2.0)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PB_0 * PC[b1] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PB_0 * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PB_0 * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[m] + PA_0 * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_g * PC[b1] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[m])
                            + delta[b0][g] * (PA_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[m] + PA_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PB_g * PC[h] * PC[m] + PA_0 * PB_h * PC[g] * PC[m] + PA_0 * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a0] * PC[m] + PB_g * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[g] * PC[m])
                            + delta[a0][h] * (PB_0 * PB_1 * PC[g] * PC[m] + PB_0 * PB_g * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[b0] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][g] * (PB_0 * PB_1 * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[m] + PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[b0] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][b1] * (PB_0 * PB_g * PC[h] * PC[m] + PB_0 * PB_h * PC[g] * PC[m] + PB_0 * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[b0] * PC[m] + PB_g * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PB_1 * PB_g * PC[h] * PC[m] + PB_1 * PB_h * PC[g] * PC[m] + PB_1 * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[b1] * PC[m] + PB_g * PC[b1] * PC[h] * PC[m] + PB_h * PC[b1] * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[g] + PB_g * PC[a0] + PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[h] + PB_h * PC[a0] + PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PB_0 * PC[g] + PB_g * PC[b0] + PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PB_0 * PC[h] + PB_h * PC[b0] + PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PB_1 * PC[g] + PB_g * PC[b1] + PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PB_1 * PC[h] + PB_h * PC[b1] + PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PB_g * PC[h] + PB_h * PC[g] + PC[g] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_g * PB_h * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PB_g * PC[a0] * PC[h] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[a0] * PC[g] * PC[m]
                            + PB_0 * PB_g * PB_h * PC[a0] * PC[b1] * PC[m]
                            + PB_1 * PB_g * PB_h * PC[a0] * PC[b0] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PB_0 * PC[b1] * PC[g] + PA_0 * PB_1 * PC[b0] * PC[g] + PA_0 * PB_g * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[g] + PB_0 * PB_g * PC[a0] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[b0])
                            + delta[g][m] * (PA_0 * PB_0 * PC[b1] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[h] + PA_0 * PB_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[h] + PB_0 * PB_h * PC[a0] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[b0])
                            + delta[b1][m] * (PA_0 * PB_0 * PC[g] * PC[h] + PA_0 * PB_g * PC[b0] * PC[h] + PA_0 * PB_h * PC[b0] * PC[g] + PB_0 * PB_g * PC[a0] * PC[h] + PB_0 * PB_h * PC[a0] * PC[g] + PB_g * PB_h * PC[a0] * PC[b0])
                            + delta[b0][m] * (PA_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_g * PC[b1] * PC[h] + PA_0 * PB_h * PC[b1] * PC[g] + PB_1 * PB_g * PC[a0] * PC[h] + PB_1 * PB_h * PC[a0] * PC[g] + PB_g * PB_h * PC[a0] * PC[b1])
                            + delta[a0][m] * (PB_0 * PB_1 * PC[g] * PC[h] + PB_0 * PB_g * PC[b1] * PC[h] + PB_0 * PB_h * PC[b1] * PC[g] + PB_1 * PB_g * PC[b0] * PC[h] + PB_1 * PB_h * PC[b0] * PC[g] + PB_g * PB_h * PC[b0] * PC[b1])
                        )

                    )

                    + F6_t[4] * (

                        (-2.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[h] * PC[m])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[m] + PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[m] + PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[g] * PC[m] + PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PB_0 * PC[b1] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PB_0 * PC[b1] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[m] + PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PB_0 * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[h] * PC[m] + PB_h * PC[b0] * PC[g] * PC[m] + PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PB_1 * PC[g] * PC[h] * PC[m] + PB_g * PC[b1] * PC[h] * PC[m] + PB_h * PC[b1] * PC[g] * PC[m] + PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[g] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PB_1 * PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PB_g * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[b0] * PC[g] + PB_g * PC[a0] * PC[b0] * PC[b1])
                            + delta[g][m] * (PA_0 * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a0] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[b1])
                            + delta[b1][m] * (PA_0 * PC[b0] * PC[g] * PC[h] + PB_0 * PC[a0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[g])
                            + delta[b0][m] * (PA_0 * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b1] * PC[h] + PB_h * PC[a0] * PC[b1] * PC[g])
                            + delta[a0][m] * (PB_0 * PC[b1] * PC[g] * PC[h] + PB_1 * PC[b0] * PC[g] * PC[h] + PB_g * PC[b0] * PC[b1] * PC[h] + PB_h * PC[b0] * PC[b1] * PC[g])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                        )

                    )

                    + F6_t[5] * (

                        4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F6_t[6] * (

                        (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }


    // D-D block

    #pragma omp parallel for schedule(static, PAD_SIZE)
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

        const auto i_atom = cart_ao_to_atom_ids[i_cgto];
        const auto j_atom = cart_ao_to_atom_ids[j_cgto];

        const double rij[3] = {x_j - x_i, y_j - y_i, z_j - z_i};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        const auto PA_0 = (a_j / (a_i + a_j)) * rij[a0];
        const auto PA_1 = (a_j / (a_i + a_j)) * rij[a1];

        const auto PB_0 = (-a_i / (a_i + a_j)) * rij[b0];
        const auto PB_1 = (-a_i / (a_i + a_j)) * rij[b1];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        double V_const = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j) * S_ij_00;

        // loop over hessian components
        for (int g = 0; g < 3; g++)
        {
            const auto PA_g = (a_j / (a_i + a_j)) * rij[g];
            const auto PB_g = (-a_i / (a_i + a_j)) * rij[g];

        for (int h = 0; h < 3; h++)
        {
            const auto PA_h = (a_j / (a_i + a_j)) * rij[h];
            const auto PB_h = (-a_i / (a_i + a_j)) * rij[h];

        for (int c = 0; c < ndipoles; c++)
        {
            const auto x_c   = dipoles_info[c + ndipoles * 0];
            const auto y_c   = dipoles_info[c + ndipoles * 1];
            const auto z_c   = dipoles_info[c + ndipoles * 2];
            const auto x_dip = dipoles_info[c + ndipoles * 3];
            const auto y_dip = dipoles_info[c + ndipoles * 4];
            const auto z_dip = dipoles_info[c + ndipoles * 5];

            const double dip[3] = {x_dip, y_dip, z_dip};

            const double PC[3] = {(a_i * x_i + a_j * x_j) / (a_i + a_j) - x_c,
                                  (a_i * y_i + a_j * y_j) / (a_i + a_j) - y_c,
                                  (a_i * z_i + a_j * z_j) / (a_i + a_j) - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F7_t[8];

            onee::computeBoysFunction(F7_t, (a_i + a_j) * r2_PC, 7, boys_func_table.data(), boys_func_ft.data());

            for (int m = 0; m < 3; m++)
            {
                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ii = (-1.0) * V_const * dip[m] * (

                    F7_t[1] * (

                        1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m] * (-2.0))
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PA_1 * PC[m] * (-1.0))
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PA_g * PC[m] * (-1.0))
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PA_h * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-1.0))
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PA_g * PC[m] * (-1.0))
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PA_h * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PB_1 * PC[m] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * (-2.0))
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PA_0 * (-1.0))
                            + (delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PA_1 * (-1.0))
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g * (-1.0))
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h * (-1.0))
                            + (delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][g] * delta[b1][m]) * (PB_0 * (-2.0))
                            + (delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][g] * delta[b0][m]) * (PB_1 * (-2.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PB_1 * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][m] * delta[g][h] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][m] * delta[g][h] * (PA_0 * PA_1 * PB_1)
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PA_g * PB_0)
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PA_g * PB_1)
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PA_h * PB_0)
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PA_h * PB_1)
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1)
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PA_g * PB_0)
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PA_g * PB_1)
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PA_h * PB_0)
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PA_h * PB_1)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PB_1)
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PB_1)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PA_h * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_0 * PA_1 * PA_g * PA_h * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PB_1 * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PB_1 * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PA_h * PB_0 * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PA_h * PB_1 * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PA_h * PB_0 * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PA_h * PB_1 * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a0][a1] * (PA_g * PA_h * PB_0 * PB_1 * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PA_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PA_h)
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PB_0)
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PB_1)
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PB_0)
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PB_1)
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PA_h)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PB_1)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PB_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PB_1)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PA_h * PB_0)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PA_h * PB_1)
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PB_1)
                        )

                        + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_1 * PA_g * PA_h * PB_0 * PB_1 * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[b1][m] * (PA_0 * PA_1 * PA_g * PA_h * PB_0)
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PA_h * PB_1)
                            + delta[h][m] * (PA_0 * PA_1 * PA_g * PB_0 * PB_1)
                            + delta[g][m] * (PA_0 * PA_1 * PA_h * PB_0 * PB_1)
                            + delta[a1][m] * (PA_0 * PA_g * PA_h * PB_0 * PB_1)
                            + delta[a0][m] * (PA_1 * PA_g * PA_h * PB_0 * PB_1)
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                        )

                        + (
                            (delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][g] * delta[b1][m]) * (PB_0)
                            + (delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][g] * delta[b0][m]) * (PB_1)
                        )

                    )

                    + F7_t[2] * (

                        (-3.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * (-2.0) + PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m] * (-2.0) + PA_0 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m] * (-2.0) + PA_0 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-2.0) + PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-2.0) + PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m] * (-2.0) + PA_1 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m] * (-2.0) + PA_1 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-2.0) + PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-2.0) + PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PA_h * PC[m] * (-2.0) + PA_g * PC[h] * PC[m] * (-1.0) + PA_h * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] * (-2.0) + PA_g * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] * (-2.0) + PA_g * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] * (-2.0) + PA_h * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] * (-2.0) + PA_h * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * (-2.0) + PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0))
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g * (-2.0) + PC[g] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h * (-2.0) + PC[h] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b0][b1] * (PA_0 * PA_1 * PA_g * PA_h * PC[m] + PA_0 * PA_1 * PA_g * PC[h] * PC[m] + PA_0 * PA_1 * PA_h * PC[g] * PC[m] + PA_0 * PA_g * PA_h * PC[a1] * PC[m] + PA_1 * PA_g * PA_h * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PA_g * PB_0 * PC[m] + PA_0 * PA_1 * PA_g * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_g * PB_0 * PC[a1] * PC[m] + PA_1 * PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PB_1 * PC[m] + PA_0 * PA_1 * PA_g * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_g * PB_1 * PC[a1] * PC[m] + PA_1 * PA_g * PB_1 * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PA_h * PB_0 * PC[m] + PA_0 * PA_1 * PA_h * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[a1] * PC[m] + PA_1 * PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PB_1 * PC[m] + PA_0 * PA_1 * PA_h * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[a1] * PC[m] + PA_1 * PA_h * PB_1 * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PA_h * PB_0 * PC[m] + PA_0 * PA_g * PA_h * PC[b0] * PC[m] + PA_0 * PA_g * PB_0 * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[g] * PC[m] + PA_g * PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PA_h * PB_1 * PC[m] + PA_0 * PA_g * PA_h * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[g] * PC[m] + PA_g * PA_h * PB_1 * PC[a0] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PB_1 * PC[m] + PA_0 * PA_g * PB_0 * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PB_1 * PC[m] + PA_0 * PA_h * PB_0 * PC[b1] * PC[m] + PA_0 * PA_h * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PA_h * PB_0 * PC[m] + PA_1 * PA_g * PA_h * PC[b0] * PC[m] + PA_1 * PA_g * PB_0 * PC[h] * PC[m] + PA_1 * PA_h * PB_0 * PC[g] * PC[m] + PA_g * PA_h * PB_0 * PC[a1] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PA_h * PB_1 * PC[m] + PA_1 * PA_g * PA_h * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[h] * PC[m] + PA_1 * PA_h * PB_1 * PC[g] * PC[m] + PA_g * PA_h * PB_1 * PC[a1] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PB_1 * PC[m] + PA_1 * PA_g * PB_0 * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a1] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PB_1 * PC[m] + PA_1 * PA_h * PB_0 * PC[b1] * PC[m] + PA_1 * PA_h * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a1] * PC[m])
                            + delta[a0][a1] * (PA_g * PA_h * PB_0 * PB_1 * PC[m] + PA_g * PA_h * PB_0 * PC[b1] * PC[m] + PA_g * PA_h * PB_1 * PC[b0] * PC[m] + PA_g * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PA_g + PA_0 * PA_1 * PC[g] + PA_0 * PA_g * PC[a1] + PA_1 * PA_g * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PA_h + PA_0 * PA_1 * PC[h] + PA_0 * PA_h * PC[a1] + PA_1 * PA_h * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PA_h + PA_0 * PA_g * PC[h] + PA_0 * PA_h * PC[g] + PA_g * PA_h * PC[a0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PB_0 + PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PB_1 + PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_g * PB_1 * PC[a0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PB_0 + PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PB_1 + PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_h * PB_1 * PC[a0])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PA_h + PA_1 * PA_g * PC[h] + PA_1 * PA_h * PC[g] + PA_g * PA_h * PC[a1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PB_0 + PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_g * PB_0 * PC[a1])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PB_1 + PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_g * PB_1 * PC[a1])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PB_0 + PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_h * PB_0 * PC[a1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PB_1 + PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_h * PB_1 * PC[a1])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PA_h * PB_0 + PA_g * PA_h * PC[b0] + PA_g * PB_0 * PC[h] + PA_h * PB_0 * PC[g])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PA_h * PB_1 + PA_g * PA_h * PC[b1] + PA_g * PB_1 * PC[h] + PA_h * PB_1 * PC[g])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PB_1 + PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PB_0 * PB_1 * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PB_1 + PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PB_0 * PB_1 * PC[h])
                        )

                        + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_1 * PA_g * PA_h * PB_0 * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PA_g * PA_h * PB_1 * PC[b0] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                            + PA_0 * PA_g * PA_h * PB_0 * PB_1 * PC[a1] * PC[m]
                            + PA_1 * PA_g * PA_h * PB_0 * PB_1 * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[b1][m] * (PA_0 * PA_1 * PA_g * PA_h * PC[b0] + PA_0 * PA_1 * PA_g * PB_0 * PC[h] + PA_0 * PA_1 * PA_h * PB_0 * PC[g] + PA_0 * PA_g * PA_h * PB_0 * PC[a1] + PA_1 * PA_g * PA_h * PB_0 * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PA_h * PC[b1] + PA_0 * PA_1 * PA_g * PB_1 * PC[h] + PA_0 * PA_1 * PA_h * PB_1 * PC[g] + PA_0 * PA_g * PA_h * PB_1 * PC[a1] + PA_1 * PA_g * PA_h * PB_1 * PC[a0])
                            + delta[h][m] * (PA_0 * PA_1 * PA_g * PB_0 * PC[b1] + PA_0 * PA_1 * PA_g * PB_1 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_1 * PC[g] + PA_0 * PA_g * PB_0 * PB_1 * PC[a1] + PA_1 * PA_g * PB_0 * PB_1 * PC[a0])
                            + delta[g][m] * (PA_0 * PA_1 * PA_h * PB_0 * PC[b1] + PA_0 * PA_1 * PA_h * PB_1 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_1 * PC[h] + PA_0 * PA_h * PB_0 * PB_1 * PC[a1] + PA_1 * PA_h * PB_0 * PB_1 * PC[a0])
                            + delta[a1][m] * (PA_0 * PA_g * PA_h * PB_0 * PC[b1] + PA_0 * PA_g * PA_h * PB_1 * PC[b0] + PA_0 * PA_g * PB_0 * PB_1 * PC[h] + PA_0 * PA_h * PB_0 * PB_1 * PC[g] + PA_g * PA_h * PB_0 * PB_1 * PC[a0])
                            + delta[a0][m] * (PA_1 * PA_g * PA_h * PB_0 * PC[b1] + PA_1 * PA_g * PA_h * PB_1 * PC[b0] + PA_1 * PA_g * PB_0 * PB_1 * PC[h] + PA_1 * PA_h * PB_0 * PB_1 * PC[g] + PA_g * PA_h * PB_0 * PB_1 * PC[a1])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                        )

                        + (-1.0) * (
                            (delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][g] * delta[b1][m]) * (PC[b0])
                            + (delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][g] * delta[b0][m]) * (PC[b1])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * 2.0)
                            + (delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m] * 4.0)
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PA_g * PC[m] + PA_1 * PC[g] * PC[m] + PA_g * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PA_h * PC[m] + PA_1 * PC[h] * PC[m] + PA_h * PC[a1] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m])
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * 2.0 + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0)
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PA_1 + PC[a1])
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g + PC[g])
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h + PC[h])
                            + (delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0])
                            + (delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][g] * delta[b1][m]) * (PB_0 * 2.0 + PC[b0] * 2.0)
                            + (delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1])
                            + (delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][g] * delta[b0][m]) * (PB_1 * 2.0 + PC[b1] * 2.0)
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PC[b1] * PC[m] + PA_0 * PA_h * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a1] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PC[b1] * PC[m] + PA_1 * PA_h * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a1] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][m] * delta[g][h] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + delta[b0][m] * delta[g][h] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_g * PB_1 * PC[a0])
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_h * PB_1 * PC[a0])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_g * PB_0 * PC[a1])
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_g * PB_1 * PC[a1])
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_h * PB_0 * PC[a1])
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_h * PB_1 * PC[a1])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PB_0 * PB_1 * PC[g])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PB_0 * PB_1 * PC[h])
                        )

                    )

                    + F7_t[3] * (

                        1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m] * (-2.0))
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-1.0))
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-1.0))
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-1.0))
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-1.0))
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[g] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[g] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[h] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[h] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] * (-2.0) + PB_1 * PC[b0] * PC[m] * (-2.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PC[a0] * (-1.0))
                            + (delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PC[a1] * (-1.0))
                            + (delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][g] * delta[b1][m]) * (PC[b0] * (-2.0))
                            + (delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1] * (-1.0))
                            + (delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][g] * delta[b0][m]) * (PC[b1] * (-2.0))
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g] * (-1.0))
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h] * (-1.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][m] * delta[g][h] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + delta[b0][m] * delta[g][h] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g])
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PC[b0] * PC[g] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g])
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PC[b0] * PC[h] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h])
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PC[b1] * PC[g] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g])
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PC[b1] * PC[h] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h])
                        )

                        + 3.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] * 2.0 + PA_g * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] * 2.0 + PA_h * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[h] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m] + PA_1 * PC[g] * PC[m] * 2.0 + PA_g * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m] + PA_1 * PC[h] * PC[m] * 2.0 + PA_h * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PA_h * PC[m] + PA_g * PC[h] * PC[m] * 2.0 + PA_h * PC[g] * PC[m] * 2.0 + PC[g] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] * 2.0 + PB_0 * PC[g] * PC[m] * 2.0 + PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] * 2.0 + PB_1 * PC[g] * PC[m] * 2.0 + PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] * 2.0 + PB_0 * PC[h] * PC[m] * 2.0 + PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] * 2.0 + PB_1 * PC[h] * PC[m] * 2.0 + PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b1] * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g + PC[g] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h + PC[h] * 2.0)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0] * 2.0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1] * 2.0)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PA_g * PB_0 * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PB_0 * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PA_g * PB_1 * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PB_1 * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_g * PC[h] * PC[m] + PA_0 * PA_1 * PA_h * PC[g] * PC[m] + PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_g * PA_h * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PA_h * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[g] * PC[m] + PA_g * PA_h * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[a1] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PB_0 * PC[a0] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[a1] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PB_1 * PC[a0] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PA_h * PC[b0] * PC[m] + PA_0 * PA_g * PB_0 * PC[h] * PC[m] + PA_0 * PA_g * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[g] * PC[m] + PA_0 * PA_h * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PB_0 * PC[a0] * PC[m] + PA_g * PA_h * PC[a0] * PC[b0] * PC[m] + PA_g * PB_0 * PC[a0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[g] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PA_h * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[h] * PC[m] + PA_0 * PA_g * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[g] * PC[m] + PA_0 * PA_h * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PB_1 * PC[a0] * PC[m] + PA_g * PA_h * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[b0] * PC[m] + PA_0 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a0] * PC[m] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PC[b1] * PC[m] + PA_0 * PA_h * PB_1 * PC[b0] * PC[m] + PA_0 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a0] * PC[m] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PA_h * PC[b0] * PC[m] + PA_1 * PA_g * PB_0 * PC[h] * PC[m] + PA_1 * PA_g * PC[b0] * PC[h] * PC[m] + PA_1 * PA_h * PB_0 * PC[g] * PC[m] + PA_1 * PA_h * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PB_0 * PC[a1] * PC[m] + PA_g * PA_h * PC[a1] * PC[b0] * PC[m] + PA_g * PB_0 * PC[a1] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PA_h * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[h] * PC[m] + PA_1 * PA_g * PC[b1] * PC[h] * PC[m] + PA_1 * PA_h * PB_1 * PC[g] * PC[m] + PA_1 * PA_h * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PB_1 * PC[a1] * PC[m] + PA_g * PA_h * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a1] * PC[g] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[b0] * PC[m] + PA_1 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a1] * PC[m] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PC[b1] * PC[m] + PA_1 * PA_h * PB_1 * PC[b0] * PC[m] + PA_1 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a1] * PC[m] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m])
                            + delta[a0][a1] * (PA_g * PA_h * PB_0 * PC[b1] * PC[m] + PA_g * PA_h * PB_1 * PC[b0] * PC[m] + PA_g * PA_h * PC[b0] * PC[b1] * PC[m] + PA_g * PB_0 * PB_1 * PC[h] * PC[m] + PA_g * PB_0 * PC[b1] * PC[h] * PC[m] + PA_g * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[g] * PC[m] + PA_h * PB_0 * PC[b1] * PC[g] * PC[m] + PA_h * PB_1 * PC[b0] * PC[g] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PA_g * PC[a1] + PA_0 * PC[a1] * PC[g] + PA_1 * PA_g * PC[a0] + PA_1 * PC[a0] * PC[g] + PA_g * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PA_h * PC[a1] + PA_0 * PC[a1] * PC[h] + PA_1 * PA_h * PC[a0] + PA_1 * PC[a0] * PC[h] + PA_h * PC[a0] * PC[a1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_0 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_0 * PC[b1] * PC[g] + PA_g * PB_1 * PC[a0] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PC[h] + PA_0 * PA_h * PC[g] + PA_0 * PC[g] * PC[h] + PA_g * PA_h * PC[a0] + PA_g * PC[a0] * PC[h] + PA_h * PC[a0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_0 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_0 * PC[b1] * PC[h] + PA_h * PB_1 * PC[a0] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_1 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a1] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_1 * PC[b1] * PC[g] + PA_g * PB_1 * PC[a1] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PC[h] + PA_1 * PA_h * PC[g] + PA_1 * PC[g] * PC[h] + PA_g * PA_h * PC[a1] + PA_g * PC[a1] * PC[h] + PA_h * PC[a1] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_1 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a1] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_1 * PC[b1] * PC[h] + PA_h * PB_1 * PC[a1] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PA_h * PC[b0] + PA_g * PB_0 * PC[h] + PA_g * PC[b0] * PC[h] + PA_h * PB_0 * PC[g] + PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PA_h * PC[b1] + PA_g * PB_1 * PC[h] + PA_g * PC[b1] * PC[h] + PA_h * PB_1 * PC[g] + PA_h * PC[b1] * PC[g] + PB_1 * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PA_g * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[g] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PA_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[h] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_1 * PA_g * PA_h * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_g * PA_h * PB_0 * PC[a1] * PC[b1] * PC[m]
                            + PA_0 * PA_g * PA_h * PB_1 * PC[a1] * PC[b0] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_1 * PC[a1] * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_1 * PC[a1] * PC[g] * PC[m]
                            + PA_1 * PA_g * PA_h * PB_0 * PC[a0] * PC[b1] * PC[m]
                            + PA_1 * PA_g * PA_h * PB_1 * PC[a0] * PC[b0] * PC[m]
                            + PA_1 * PA_g * PB_0 * PB_1 * PC[a0] * PC[h] * PC[m]
                            + PA_1 * PA_h * PB_0 * PB_1 * PC[a0] * PC[g] * PC[m]
                            + PA_g * PA_h * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[b1] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] + PA_0 * PA_g * PB_0 * PC[a1] * PC[b1] + PA_0 * PA_g * PB_1 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] + PA_1 * PA_g * PB_0 * PC[a0] * PC[b1] + PA_1 * PA_g * PB_1 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] + PA_g * PB_0 * PB_1 * PC[a0] * PC[a1])
                            + delta[b1][m] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[h] + PA_0 * PA_1 * PA_h * PC[b0] * PC[g] + PA_0 * PA_1 * PB_0 * PC[g] * PC[h] + PA_0 * PA_g * PA_h * PC[a1] * PC[b0] + PA_0 * PA_g * PB_0 * PC[a1] * PC[h] + PA_0 * PA_h * PB_0 * PC[a1] * PC[g] + PA_1 * PA_g * PA_h * PC[a0] * PC[b0] + PA_1 * PA_g * PB_0 * PC[a0] * PC[h] + PA_1 * PA_h * PB_0 * PC[a0] * PC[g] + PA_g * PA_h * PB_0 * PC[a0] * PC[a1])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PC[b1] * PC[h] + PA_0 * PA_1 * PA_h * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[g] * PC[h] + PA_0 * PA_g * PA_h * PC[a1] * PC[b1] + PA_0 * PA_g * PB_1 * PC[a1] * PC[h] + PA_0 * PA_h * PB_1 * PC[a1] * PC[g] + PA_1 * PA_g * PA_h * PC[a0] * PC[b1] + PA_1 * PA_g * PB_1 * PC[a0] * PC[h] + PA_1 * PA_h * PB_1 * PC[a0] * PC[g] + PA_g * PA_h * PB_1 * PC[a0] * PC[a1])
                            + delta[g][m] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[b1] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[h] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[h] + PA_0 * PA_h * PB_0 * PC[a1] * PC[b1] + PA_0 * PA_h * PB_1 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[h] + PA_1 * PA_h * PB_0 * PC[a0] * PC[b1] + PA_1 * PA_h * PB_1 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[h] + PA_h * PB_0 * PB_1 * PC[a0] * PC[a1])
                            + delta[a1][m] * (PA_0 * PA_g * PA_h * PC[b0] * PC[b1] + PA_0 * PA_g * PB_0 * PC[b1] * PC[h] + PA_0 * PA_g * PB_1 * PC[b0] * PC[h] + PA_0 * PA_h * PB_0 * PC[b1] * PC[g] + PA_0 * PA_h * PB_1 * PC[b0] * PC[g] + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] + PA_g * PA_h * PB_0 * PC[a0] * PC[b1] + PA_g * PA_h * PB_1 * PC[a0] * PC[b0] + PA_g * PB_0 * PB_1 * PC[a0] * PC[h] + PA_h * PB_0 * PB_1 * PC[a0] * PC[g])
                            + delta[a0][m] * (PA_1 * PA_g * PA_h * PC[b0] * PC[b1] + PA_1 * PA_g * PB_0 * PC[b1] * PC[h] + PA_1 * PA_g * PB_1 * PC[b0] * PC[h] + PA_1 * PA_h * PB_0 * PC[b1] * PC[g] + PA_1 * PA_h * PB_1 * PC[b0] * PC[g] + PA_1 * PB_0 * PB_1 * PC[g] * PC[h] + PA_g * PA_h * PB_0 * PC[a1] * PC[b1] + PA_g * PA_h * PB_1 * PC[a1] * PC[b0] + PA_g * PB_0 * PB_1 * PC[a1] * PC[h] + PA_h * PB_0 * PB_1 * PC[a1] * PC[g])
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                        )

                    )

                    + F7_t[4] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[g] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[g] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PC[h] * PC[m] * (-1.0) + PA_h * PC[g] * PC[m] * (-1.0) + PC[g] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[h] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[h] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g])
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_g * PC[a1] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PA_g * PC[a0] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a0] * PC[b0] * PC[m] + PA_g * PB_0 * PC[a0] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[g] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[g] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PC[b0] * PC[h] * PC[m] + PA_1 * PA_h * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a1] * PC[b0] * PC[m] + PA_g * PB_0 * PC[a1] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[g] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PC[b1] * PC[h] * PC[m] + PA_1 * PA_h * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PA_h * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a1] * PC[g] * PC[m] + PA_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[a0][a1] * (PA_g * PA_h * PC[b0] * PC[b1] * PC[m] + PA_g * PB_0 * PC[b1] * PC[h] * PC[m] + PA_g * PB_1 * PC[b0] * PC[h] * PC[m] + PA_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_0 * PC[b1] * PC[g] * PC[m] + PA_h * PB_1 * PC[b0] * PC[g] * PC[m] + PA_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PA_g * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PA_h * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g] + PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h] + PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g] + PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h] + PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[g] * PC[h] + PA_g * PC[a0] * PC[h] + PA_h * PC[a0] * PC[g] + PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PC[b0] * PC[g] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g] + PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PC[b0] * PC[h] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h] + PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PC[b1] * PC[g] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g] + PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PC[b1] * PC[h] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h] + PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[g] * PC[h] + PA_g * PC[a1] * PC[h] + PA_h * PC[a1] * PC[g] + PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PC[b0] * PC[h] + PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h] + PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PC[b1] * PC[h] + PA_h * PC[b1] * PC[g] + PB_1 * PC[g] * PC[h] + PC[b1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PC[b0] * PC[b1] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_1 * PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_g * PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PA_g * PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PA_g * PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PA_1 * PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_1 * PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PA_g * PA_h * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m]
                            + PA_g * PA_h * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m]
                            + PA_g * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g])
                            + delta[g][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[b1] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h])
                            + delta[b1][m] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[g] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[h] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[g] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[h] + PA_g * PA_h * PC[a0] * PC[a1] * PC[b0] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[h] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[g])
                            + delta[b0][m] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[h] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[h] + PA_g * PA_h * PC[a0] * PC[a1] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[h] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[g])
                            + delta[a1][m] * (PA_0 * PA_g * PC[b0] * PC[b1] * PC[h] + PA_0 * PA_h * PC[b0] * PC[b1] * PC[g] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_g * PA_h * PC[a0] * PC[b0] * PC[b1] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h])
                            + delta[a0][m] * (PA_1 * PA_g * PC[b0] * PC[b1] * PC[h] + PA_1 * PA_h * PC[b0] * PC[b1] * PC[g] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_g * PA_h * PC[a1] * PC[b0] * PC[b1] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[g] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[h])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            delta[b0][b1] * delta[g][h] * (PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PC[a0] * PC[g] * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PC[a1] * PC[g] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PC[a1] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m] * 2.0)
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][m] * delta[g][h] * (PC[a0] * PC[a1] * PC[b0])
                            + delta[b0][m] * delta[g][h] * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PC[a0] * PC[b0] * PC[b1])
                            + delta[a1][h] * delta[b1][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[a1][g] * delta[b1][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[a1][h] * delta[b0][m] * (PC[a0] * PC[b1] * PC[g])
                            + delta[a1][g] * delta[b0][m] * (PC[a0] * PC[b1] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[a1] * PC[b0] * PC[b1])
                            + delta[a0][h] * delta[b1][m] * (PC[a1] * PC[b0] * PC[g])
                            + delta[a0][g] * delta[b1][m] * (PC[a1] * PC[b0] * PC[h])
                            + delta[a0][h] * delta[b0][m] * (PC[a1] * PC[b1] * PC[g])
                            + delta[a0][g] * delta[b0][m] * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                        )

                    )

                    + F7_t[5] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_i * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[g] * PC[h] * PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PA_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_i * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_g * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_h * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_g * PA_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_i * (
                            delta[h][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h])
                            + delta[b1][m] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h])
                            + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h])
                            + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h])
                        )

                    )

                    + F7_t[6] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_i * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_i * (
                            PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_i * (
                            delta[h][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h])
                            + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F7_t[7] * (

                        8.0 * (a_i + a_j) * a_i * a_i * (
                            PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ij = (-1.0) * V_const * dip[m] * (

                    F7_t[1] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h]) * (PC[m])
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PA_g * PC[m])
                            + delta[a1][g] * delta[b1][h] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PA_g * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PB_1 * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PA_g * PB_0 * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PA_g * PB_1 * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][g] * delta[b1][h] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PB_1 * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PB_h * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PB_1 * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PB_h * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PA_g)
                            + (delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PB_1)
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h]) * (PA_0)
                            + (delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h]) * (PA_1)
                            + (delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PB_0)
                            + (delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PB_1)
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h)
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PB_1 * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][g] * (PA_0 * PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PB_h * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PA_g)
                            + delta[b1][h] * delta[g][m] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][h] * delta[g][m] * (PA_0 * PA_1 * PB_1)
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PA_g * PB_0)
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PA_g * PB_1)
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PA_g * PB_0)
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PA_g * PB_1)
                        )

                        + (-2.0) * a_j * (
                            delta[a1][g] * delta[h][m] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PB_0 * PB_h)
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PB_1 * PB_h)
                            + delta[a0][g] * delta[h][m] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PB_0 * PB_h)
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PB_1 * PB_h)
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PB_h)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PB_h * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PB_0 * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_g * PB_h * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PB_h * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PB_h * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PB_0 * PB_h * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PB_1 * PB_h * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PB_0 * PB_h * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PB_1 * PB_h * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][a1] * (PA_g * PB_0 * PB_1 * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PA_g)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PB_h)
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PB_0)
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PB_h)
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PB_h)
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PB_h)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PB_0)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PB_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PB_h)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PB_h)
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PB_0 * PB_h)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PB_1 * PB_h)
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PB_h)
                        )

                        + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h]) * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_g * PB_0 * PB_1 * PB_h * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PA_g * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_1 * PA_g * PB_0 * PB_h)
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PB_1 * PB_h)
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_h)
                            + delta[a1][m] * (PA_0 * PA_g * PB_0 * PB_1 * PB_h)
                            + delta[a0][m] * (PA_1 * PA_g * PB_0 * PB_1 * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a1][g] * delta[b1][h] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PB_1 * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PB_1 * PC[m])
                        )

                        + (
                            (delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h]) * (PA_0)
                            + (delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h]) * (PA_1)
                            + (delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PB_0)
                            + (delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PB_1)
                        )

                    )

                    + F7_t[2] * (

                        (-3.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * (-2.0) + PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m] * (-2.0) + PA_0 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-2.0) + PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-2.0) + PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m] * (-2.0) + PA_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m] * (-2.0) + PA_1 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-2.0) + PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-2.0) + PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m] * (-2.0) + PA_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] * (-2.0) + PA_g * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] * (-2.0) + PA_g * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PB_h * PC[m] * (-2.0) + PA_g * PC[h] * PC[m] * (-1.0) + PB_h * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * (-2.0) + PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m] * (-2.0) + PB_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m] * (-2.0) + PB_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b1] * PC[m] * (-1.0))
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g * (-2.0) + PC[g] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h * (-2.0) + PC[h] * (-1.0))
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PB_0 * PC[m] + PA_0 * PA_1 * PA_g * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_g * PB_0 * PC[a1] * PC[m] + PA_1 * PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PB_1 * PC[m] + PA_0 * PA_1 * PA_g * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_g * PB_1 * PC[a1] * PC[m] + PA_1 * PA_g * PB_1 * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_g * PB_h * PC[m] + PA_0 * PA_1 * PA_g * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[g] * PC[m] + PA_0 * PA_g * PB_h * PC[a1] * PC[m] + PA_1 * PA_g * PB_h * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PB_h * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b0] * PC[m] + PA_0 * PB_0 * PB_h * PC[a1] * PC[m] + PA_1 * PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PB_h * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[a1] * PC[m] + PA_1 * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PB_1 * PC[m] + PA_0 * PA_g * PB_0 * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PB_0 * PB_h * PC[m] + PA_0 * PA_g * PB_0 * PC[h] * PC[m] + PA_0 * PA_g * PB_h * PC[b0] * PC[m] + PA_0 * PB_0 * PB_h * PC[g] * PC[m] + PA_g * PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PB_1 * PB_h * PC[m] + PA_0 * PA_g * PB_1 * PC[h] * PC[m] + PA_0 * PA_g * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[g] * PC[m] + PA_g * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PB_h * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PB_1 * PC[m] + PA_1 * PA_g * PB_0 * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a1] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PB_0 * PB_h * PC[m] + PA_1 * PA_g * PB_0 * PC[h] * PC[m] + PA_1 * PA_g * PB_h * PC[b0] * PC[m] + PA_1 * PB_0 * PB_h * PC[g] * PC[m] + PA_g * PB_0 * PB_h * PC[a1] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PB_1 * PB_h * PC[m] + PA_1 * PA_g * PB_1 * PC[h] * PC[m] + PA_1 * PA_g * PB_h * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[g] * PC[m] + PA_g * PB_1 * PB_h * PC[a1] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PB_h * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a1] * PC[m])
                            + delta[a0][a1] * (PA_g * PB_0 * PB_1 * PB_h * PC[m] + PA_g * PB_0 * PB_1 * PC[h] * PC[m] + PA_g * PB_0 * PB_h * PC[b1] * PC[m] + PA_g * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PA_g + PA_0 * PA_1 * PC[g] + PA_0 * PA_g * PC[a1] + PA_1 * PA_g * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PB_h + PA_0 * PA_1 * PC[h] + PA_0 * PB_h * PC[a1] + PA_1 * PB_h * PC[a0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PB_0 + PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PB_1 + PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_g * PB_1 * PC[a0])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PB_h + PA_0 * PA_g * PC[h] + PA_0 * PB_h * PC[g] + PA_g * PB_h * PC[a0])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PB_h + PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PB_h + PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[a0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PB_0 + PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_g * PB_0 * PC[a1])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PB_1 + PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_g * PB_1 * PC[a1])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PB_h + PA_1 * PA_g * PC[h] + PA_1 * PB_h * PC[g] + PA_g * PB_h * PC[a1])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PB_h + PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PB_0 * PB_h * PC[a1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PB_h + PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PB_1 * PB_h * PC[a1])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PB_1 + PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PB_0 * PB_1 * PC[g])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PB_0 * PB_h + PA_g * PB_0 * PC[h] + PA_g * PB_h * PC[b0] + PB_0 * PB_h * PC[g])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PB_1 * PB_h + PA_g * PB_1 * PC[h] + PA_g * PB_h * PC[b1] + PB_1 * PB_h * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PB_h + PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[b0])
                        )

                        + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h]) * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_g * PB_0 * PB_1 * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_0 * PB_h * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_1 * PB_h * PC[b0] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_1 * PB_h * PC[a1] * PC[m]
                            + PA_1 * PA_g * PB_0 * PB_1 * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PA_g * PB_0 * PC[b1] + PA_0 * PA_1 * PA_g * PB_1 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_1 * PC[g] + PA_0 * PA_g * PB_0 * PB_1 * PC[a1] + PA_1 * PA_g * PB_0 * PB_1 * PC[a0])
                            + delta[b1][m] * (PA_0 * PA_1 * PA_g * PB_0 * PC[h] + PA_0 * PA_1 * PA_g * PB_h * PC[b0] + PA_0 * PA_1 * PB_0 * PB_h * PC[g] + PA_0 * PA_g * PB_0 * PB_h * PC[a1] + PA_1 * PA_g * PB_0 * PB_h * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PB_1 * PC[h] + PA_0 * PA_1 * PA_g * PB_h * PC[b1] + PA_0 * PA_1 * PB_1 * PB_h * PC[g] + PA_0 * PA_g * PB_1 * PB_h * PC[a1] + PA_1 * PA_g * PB_1 * PB_h * PC[a0])
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[h] + PA_0 * PA_1 * PB_0 * PB_h * PC[b1] + PA_0 * PA_1 * PB_1 * PB_h * PC[b0] + PA_0 * PB_0 * PB_1 * PB_h * PC[a1] + PA_1 * PB_0 * PB_1 * PB_h * PC[a0])
                            + delta[a1][m] * (PA_0 * PA_g * PB_0 * PB_1 * PC[h] + PA_0 * PA_g * PB_0 * PB_h * PC[b1] + PA_0 * PA_g * PB_1 * PB_h * PC[b0] + PA_0 * PB_0 * PB_1 * PB_h * PC[g] + PA_g * PB_0 * PB_1 * PB_h * PC[a0])
                            + delta[a0][m] * (PA_1 * PA_g * PB_0 * PB_1 * PC[h] + PA_1 * PA_g * PB_0 * PB_h * PC[b1] + PA_1 * PA_g * PB_1 * PB_h * PC[b0] + PA_1 * PB_0 * PB_1 * PB_h * PC[g] + PA_g * PB_0 * PB_1 * PB_h * PC[a1])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            delta[a1][g] * delta[b1][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                        )

                        + (-1.0) * (
                            (delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h]) * (PC[a0])
                            + (delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h]) * (PC[a1])
                            + (delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PC[b0])
                            + (delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PC[b1])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b1][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PA_g * PC[m] + PA_1 * PC[g] * PC[m] + PA_g * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][g] * delta[b1][h] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PB_h * PC[m] + PA_1 * PC[h] * PC[m] + PB_h * PC[a1] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1])
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PA_g + PC[g])
                            + (delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PB_0 + PC[b0])
                            + (delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PB_1 + PC[b1])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h]) * (PA_0 + PC[a0])
                            + (delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h]) * (PA_1 + PC[a1])
                            + (delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PB_0 + PC[b0])
                            + (delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PB_1 + PC[b1])
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h + PC[h])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_g * PB_0 * PC[a1] * PC[m] + PA_1 * PA_g * PB_0 * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_g * PB_1 * PC[a1] * PC[m] + PA_1 * PA_g * PB_1 * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a1][g] * (PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a1] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][h] * delta[g][m] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + delta[b0][h] * delta[g][m] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PA_g * PC[a1] + PA_1 * PA_g * PC[a0])
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_g * PB_0 * PC[a0])
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_g * PB_1 * PC[a0])
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_g * PB_0 * PC[a1])
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_g * PB_1 * PC[a1])
                        )

                        + 2.0 * a_j * (
                            delta[a1][g] * delta[h][m] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[a0])
                            + delta[a0][g] * delta[h][m] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PB_0 * PB_h * PC[a1])
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PB_1 * PB_h * PC[a1])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[b0])
                        )

                    )

                    + F7_t[3] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h]) * (PC[m])
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m] + PC[a0] * PC[a1] * PC[m])
                            + delta[a1][g] * delta[b1][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PC[g] * PC[m] + PA_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m] + PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m] + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PC[g] * PC[m] + PA_g * PC[a1] * PC[m] + PC[a1] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PA_g * PC[b0] * PC[m] + PB_0 * PC[g] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PA_g * PC[b1] * PC[m] + PB_1 * PC[g] * PC[m] + PC[b1] * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][g] * delta[b1][h] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m] + PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m] + PC[a1] * PC[b1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PA_1 * PC[h] * PC[m] + PB_h * PC[a1] * PC[m] + PC[a1] * PC[h] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m] + PC[b1] * PC[h] * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PC[g])
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h]) * (PC[a0])
                            + (delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h]) * (PC[a1])
                            + (delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][m] * delta[a1][g] * delta[b1][h]) * (PC[b0])
                            + (delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][m] * delta[a1][g] * delta[b0][h]) * (PC[b1])
                            + (delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][g] * (PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][h] * delta[g][m] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + delta[b0][h] * delta[g][m] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PA_g * PC[a0] * PC[a1])
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g])
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PC[b0] * PC[g] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g])
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PC[b1] * PC[g] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g])
                        )

                        + (-2.0) * a_j * (
                            delta[a1][g] * delta[h][m] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + delta[a1][g] * delta[b1][m] * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + delta[a1][g] * delta[b0][m] * (PA_0 * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1])
                            + delta[a0][g] * delta[h][m] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + delta[a0][g] * delta[b1][m] * (PA_1 * PC[b0] * PC[h] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0])
                            + delta[a0][g] * delta[b0][m] * (PA_1 * PC[b1] * PC[h] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1])
                        )

                        + 3.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PA_g * PC[m] + PA_0 * PC[g] * PC[m] * 2.0 + PA_g * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] * 2.0 + PB_h * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PA_g * PC[m] + PA_1 * PC[g] * PC[m] * 2.0 + PA_g * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m] + PA_1 * PC[h] * PC[m] * 2.0 + PB_h * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PB_0 * PC[m] + PA_g * PC[b0] * PC[m] * 2.0 + PB_0 * PC[g] * PC[m] * 2.0 + PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PB_1 * PC[m] + PA_g * PC[b1] * PC[m] * 2.0 + PB_1 * PC[g] * PC[m] * 2.0 + PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PB_h * PC[m] + PA_g * PC[h] * PC[m] * 2.0 + PB_h * PC[g] * PC[m] * 2.0 + PC[g] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] * 2.0 + PB_h * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] * 2.0 + PB_h * PC[b1] * PC[m] * 2.0 + PC[b1] * PC[h] * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PA_g + PC[g] * 2.0)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0] * 2.0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h + PC[h] * 2.0)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][h] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PA_g * PB_0 * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PB_0 * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PA_g * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PA_g * PB_1 * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PB_1 * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_g * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[g] * PC[m] + PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_g * PB_h * PC[a1] * PC[m] + PA_0 * PA_g * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[g] * PC[m] + PA_1 * PA_g * PB_h * PC[a0] * PC[m] + PA_1 * PA_g * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[g] * PC[m] + PA_g * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_h * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b1] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[a1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[a0] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PB_0 * PC[b1] * PC[m] + PA_0 * PA_g * PB_1 * PC[b0] * PC[m] + PA_0 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a0] * PC[m] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PB_0 * PC[h] * PC[m] + PA_0 * PA_g * PB_h * PC[b0] * PC[m] + PA_0 * PA_g * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[g] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_h * PC[a0] * PC[m] + PA_g * PB_0 * PC[a0] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[g] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PB_1 * PC[h] * PC[m] + PA_0 * PA_g * PB_h * PC[b1] * PC[m] + PA_0 * PA_g * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[g] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PB_h * PC[a0] * PC[m] + PA_g * PB_1 * PC[a0] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[b0] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_h * PC[a0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PB_0 * PC[b1] * PC[m] + PA_1 * PA_g * PB_1 * PC[b0] * PC[m] + PA_1 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_1 * PC[a1] * PC[m] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PB_0 * PC[h] * PC[m] + PA_1 * PA_g * PB_h * PC[b0] * PC[m] + PA_1 * PA_g * PC[b0] * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[g] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PB_h * PC[a1] * PC[m] + PA_g * PB_0 * PC[a1] * PC[h] * PC[m] + PA_g * PB_h * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_h * PC[a1] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PB_1 * PC[h] * PC[m] + PA_1 * PA_g * PB_h * PC[b1] * PC[m] + PA_1 * PA_g * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PB_h * PC[g] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_h * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PB_h * PC[a1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[h] * PC[m] + PA_g * PB_h * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PB_h * PC[b0] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_h * PC[a1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[m])
                            + delta[a0][a1] * (PA_g * PB_0 * PB_1 * PC[h] * PC[m] + PA_g * PB_0 * PB_h * PC[b1] * PC[m] + PA_g * PB_0 * PC[b1] * PC[h] * PC[m] + PA_g * PB_1 * PB_h * PC[b0] * PC[m] + PA_g * PB_1 * PC[b0] * PC[h] * PC[m] + PA_g * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_h * PC[g] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[g] * PC[m] + PB_1 * PB_h * PC[b0] * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PA_g * PC[a1] + PA_0 * PC[a1] * PC[g] + PA_1 * PA_g * PC[a0] + PA_1 * PC[a0] * PC[g] + PA_g * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PB_h * PC[a1] + PA_0 * PC[a1] * PC[h] + PA_1 * PB_h * PC[a0] + PA_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[a1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PA_g * PC[b0] + PA_0 * PB_0 * PC[g] + PA_0 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PA_g * PC[b1] + PA_0 * PB_1 * PC[g] + PA_0 * PC[b1] * PC[g] + PA_g * PB_1 * PC[a0] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_g * PC[h] + PA_0 * PB_h * PC[g] + PA_0 * PC[g] * PC[h] + PA_g * PB_h * PC[a0] + PA_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PA_0 * PC[b0] * PC[h] + PB_0 * PB_h * PC[a0] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PA_0 * PC[b1] * PC[h] + PB_1 * PB_h * PC[a0] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PA_g * PC[b0] + PA_1 * PB_0 * PC[g] + PA_1 * PC[b0] * PC[g] + PA_g * PB_0 * PC[a1] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PA_g * PC[b1] + PA_1 * PB_1 * PC[g] + PA_1 * PC[b1] * PC[g] + PA_g * PB_1 * PC[a1] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_g * PC[h] + PA_1 * PB_h * PC[g] + PA_1 * PC[g] * PC[h] + PA_g * PB_h * PC[a1] + PA_g * PC[a1] * PC[h] + PB_h * PC[a1] * PC[g])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PA_1 * PC[b0] * PC[h] + PB_0 * PB_h * PC[a1] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PA_1 * PC[b1] * PC[h] + PB_1 * PB_h * PC[a1] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PB_0 * PC[b1] + PA_g * PB_1 * PC[b0] + PA_g * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[g] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PB_0 * PC[h] + PA_g * PB_h * PC[b0] + PA_g * PC[b0] * PC[h] + PB_0 * PB_h * PC[g] + PB_0 * PC[g] * PC[h] + PB_h * PC[b0] * PC[g])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PB_1 * PC[h] + PA_g * PB_h * PC[b1] + PA_g * PC[b1] * PC[h] + PB_1 * PB_h * PC[g] + PB_1 * PC[g] * PC[h] + PB_h * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PB_h * PC[b0] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_g * PB_0 * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_1 * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PA_g * PB_h * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_1 * PC[a1] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_0 * PB_h * PC[a1] * PC[b1] * PC[m]
                            + PA_0 * PA_g * PB_1 * PB_h * PC[a1] * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_h * PC[a1] * PC[g] * PC[m]
                            + PA_1 * PA_g * PB_0 * PB_1 * PC[a0] * PC[h] * PC[m]
                            + PA_1 * PA_g * PB_0 * PB_h * PC[a0] * PC[b1] * PC[m]
                            + PA_1 * PA_g * PB_1 * PB_h * PC[a0] * PC[b0] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PB_h * PC[a0] * PC[g] * PC[m]
                            + PA_g * PB_0 * PB_1 * PB_h * PC[a0] * PC[a1] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[b1] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] + PA_0 * PA_g * PB_0 * PC[a1] * PC[b1] + PA_0 * PA_g * PB_1 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] + PA_1 * PA_g * PB_0 * PC[a0] * PC[b1] + PA_1 * PA_g * PB_1 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] + PA_g * PB_0 * PB_1 * PC[a0] * PC[a1])
                            + delta[b1][m] * (PA_0 * PA_1 * PA_g * PC[b0] * PC[h] + PA_0 * PA_1 * PB_0 * PC[g] * PC[h] + PA_0 * PA_1 * PB_h * PC[b0] * PC[g] + PA_0 * PA_g * PB_0 * PC[a1] * PC[h] + PA_0 * PA_g * PB_h * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_h * PC[a1] * PC[g] + PA_1 * PA_g * PB_0 * PC[a0] * PC[h] + PA_1 * PA_g * PB_h * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_h * PC[a0] * PC[g] + PA_g * PB_0 * PB_h * PC[a0] * PC[a1])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_g * PC[b1] * PC[h] + PA_0 * PA_1 * PB_1 * PC[g] * PC[h] + PA_0 * PA_1 * PB_h * PC[b1] * PC[g] + PA_0 * PA_g * PB_1 * PC[a1] * PC[h] + PA_0 * PA_g * PB_h * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_h * PC[a1] * PC[g] + PA_1 * PA_g * PB_1 * PC[a0] * PC[h] + PA_1 * PA_g * PB_h * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_h * PC[a0] * PC[g] + PA_g * PB_1 * PB_h * PC[a0] * PC[a1])
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[h] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[h] + PA_0 * PA_1 * PB_h * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[h] + PA_0 * PB_0 * PB_h * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_h * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[h] + PA_1 * PB_0 * PB_h * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_h * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_h * PC[a0] * PC[a1])
                            + delta[a1][m] * (PA_0 * PA_g * PB_0 * PC[b1] * PC[h] + PA_0 * PA_g * PB_1 * PC[b0] * PC[h] + PA_0 * PA_g * PB_h * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_0 * PB_h * PC[b1] * PC[g] + PA_0 * PB_1 * PB_h * PC[b0] * PC[g] + PA_g * PB_0 * PB_1 * PC[a0] * PC[h] + PA_g * PB_0 * PB_h * PC[a0] * PC[b1] + PA_g * PB_1 * PB_h * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_h * PC[a0] * PC[g])
                            + delta[a0][m] * (PA_1 * PA_g * PB_0 * PC[b1] * PC[h] + PA_1 * PA_g * PB_1 * PC[b0] * PC[h] + PA_1 * PA_g * PB_h * PC[b0] * PC[b1] + PA_1 * PB_0 * PB_1 * PC[g] * PC[h] + PA_1 * PB_0 * PB_h * PC[b1] * PC[g] + PA_1 * PB_1 * PB_h * PC[b0] * PC[g] + PA_g * PB_0 * PB_1 * PC[a1] * PC[h] + PA_g * PB_0 * PB_h * PC[a1] * PC[b1] + PA_g * PB_1 * PB_h * PC[a1] * PC[b0] + PB_0 * PB_1 * PB_h * PC[a1] * PC[g])
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a1][g] * delta[b1][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PC[a0] * PC[b1] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[a1] * PC[b1] * PC[m])
                        )

                    )

                    + F7_t[4] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PC[g] * PC[m] * (-1.0) + PA_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PA_g * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[g] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PA_g * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[g] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_g * PC[h] * PC[m] * (-1.0) + PB_h * PC[g] * PC[m] * (-1.0) + PC[g] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-2.0))
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g])
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_g * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PA_g * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[a1] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[a1][h] * (PA_0 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_g * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[g] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PB_0 * PC[a0] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[b0] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_g * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[g] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PB_1 * PC[a0] * PC[h] * PC[m] + PA_g * PB_h * PC[a0] * PC[b1] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][h] * (PA_1 * PA_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_g * PC[b0] * PC[h] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[g] * PC[m] + PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PB_0 * PC[a1] * PC[h] * PC[m] + PA_g * PB_h * PC[a1] * PC[b0] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_g * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_h * PC[b1] * PC[g] * PC[m] + PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PB_1 * PC[a1] * PC[h] * PC[m] + PA_g * PB_h * PC[a1] * PC[b1] * PC[m] + PA_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][a1] * (PA_g * PB_0 * PC[b1] * PC[h] * PC[m] + PA_g * PB_1 * PC[b0] * PC[h] * PC[m] + PA_g * PB_h * PC[b0] * PC[b1] * PC[m] + PA_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PB_h * PC[b0] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PA_g * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PC[b0] * PC[g] + PA_g * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[g] + PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PC[b1] * PC[g] + PA_g * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[g] + PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[g] * PC[h] + PA_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g] + PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PC[b0] * PC[g] + PA_g * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[g] + PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PC[b0] * PC[h] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PC[b1] * PC[g] + PA_g * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[g] + PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PC[b1] * PC[h] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[g] * PC[h] + PA_g * PC[a1] * PC[h] + PB_h * PC[a1] * PC[g] + PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PA_g * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_g * PC[b0] * PC[h] + PB_0 * PC[g] * PC[h] + PB_h * PC[b0] * PC[g] + PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_g * PC[b1] * PC[h] + PB_1 * PC[g] * PC[h] + PB_h * PC[b1] * PC[g] + PC[b1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_g * PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_g * PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_1 * PB_h * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PA_1 * PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_1 * PA_g * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PB_1 * PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PA_g * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m]
                            + PA_g * PB_0 * PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m]
                            + PA_g * PB_1 * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[a0] * PC[a1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[b1] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g])
                            + delta[g][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[h] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[b0])
                            + delta[b1][m] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PA_g * PC[a1] * PC[b0] * PC[h] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[g] + PA_1 * PA_g * PC[a0] * PC[b0] * PC[h] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[g] + PA_g * PB_0 * PC[a0] * PC[a1] * PC[h] + PA_g * PB_h * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[g])
                            + delta[b0][m] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[h] + PA_0 * PA_g * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[g] + PA_1 * PA_g * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[g] + PA_g * PB_1 * PC[a0] * PC[a1] * PC[h] + PA_g * PB_h * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[g])
                            + delta[a1][m] * (PA_0 * PA_g * PC[b0] * PC[b1] * PC[h] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[g] + PA_g * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_g * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_g * PB_h * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[g])
                            + delta[a0][m] * (PA_1 * PA_g * PC[b0] * PC[b1] * PC[h] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[g] + PA_g * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_g * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_g * PB_h * PC[a1] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[g] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[g])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + delta[a1][g] * delta[b1][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PC[a0] * PC[g] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PC[a1] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][g] * delta[b1][h] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * delta[b0][h] * (PC[a0] * PC[b1] * PC[m])
                            + delta[a1][g] * delta[b0][b1] * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][g] * delta[b1][h] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * delta[b0][h] * (PC[a1] * PC[b1] * PC[m])
                            + delta[a0][g] * delta[b0][b1] * (PC[a1] * PC[h] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][h] * delta[g][m] * (PC[a0] * PC[a1] * PC[b0])
                            + delta[b0][h] * delta[g][m] * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + delta[a1][m] * delta[b1][h] * (PC[a0] * PC[b0] * PC[g])
                            + delta[a1][m] * delta[b0][h] * (PC[a0] * PC[b1] * PC[g])
                            + delta[a0][m] * delta[b1][h] * (PC[a1] * PC[b0] * PC[g])
                            + delta[a0][m] * delta[b0][h] * (PC[a1] * PC[b1] * PC[g])
                        )

                        + 2.0 * a_j * (
                            delta[a1][g] * delta[h][m] * (PC[a0] * PC[b0] * PC[b1])
                            + delta[a1][g] * delta[b1][m] * (PC[a0] * PC[b0] * PC[h])
                            + delta[a1][g] * delta[b0][m] * (PC[a0] * PC[b1] * PC[h])
                            + delta[a0][g] * delta[h][m] * (PC[a1] * PC[b0] * PC[b1])
                            + delta[a0][g] * delta[b1][m] * (PC[a1] * PC[b0] * PC[h])
                            + delta[a0][g] * delta[b0][m] * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                        )

                    )

                    + F7_t[5] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[g] * PC[h] * PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PA_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_g * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_g * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_g * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_g * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g])
                            + delta[g][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1])
                            + delta[b1][m] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[g])
                            + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[g])
                            + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a0] * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_g * PC[a1] * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[g])
                        )

                    )

                    + F7_t[6] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h])
                            + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F7_t[7] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_ji = (-1.0) * V_const * dip[m] * (

                    F7_t[1] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m])
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PA_h * PC[m])
                            + delta[a1][h] * delta[b1][g] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PA_h * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PB_1 * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PA_h * PB_0 * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PA_h * PB_1 * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][h] * delta[b1][g] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PB_1 * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PB_g * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PB_1 * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PB_g * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PA_h)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g]) * (PA_0)
                            + (delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g]) * (PA_1)
                            + (delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g)
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PA_1 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PB_1 * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][h] * (PA_0 * PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PB_g * PC[m])
                        )

                        + (-2.0) * a_i * (
                            (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PA_h)
                            + delta[b1][g] * delta[h][m] * (PA_0 * PA_1 * PB_0)
                            + delta[b0][g] * delta[h][m] * (PA_0 * PA_1 * PB_1)
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PA_h * PB_0)
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PA_h * PB_1)
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PA_h * PB_0)
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PA_h * PB_1)
                        )

                        + (-2.0) * a_j * (
                            delta[a1][h] * delta[g][m] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PB_0 * PB_g)
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PB_1 * PB_g)
                            + delta[a0][h] * delta[g][m] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PB_0 * PB_g)
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PB_1 * PB_g)
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PB_g)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_h * PB_g * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_1 * PA_h * PB_0 * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PB_1 * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_h * PB_g * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PB_g * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PB_g * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_h * PB_0 * PB_g * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_h * PB_1 * PB_g * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PB_1 * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_h * PB_0 * PB_g * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_h * PB_1 * PB_g * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][a1] * (PA_h * PB_0 * PB_1 * PB_g * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PA_h)
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PB_g)
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PB_0)
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PB_1)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_h * PB_g)
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PB_g)
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PB_g)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PB_0)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PB_1)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_h * PB_g)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PB_g)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PB_g)
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PB_1)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_h * PB_0 * PB_g)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_h * PB_1 * PB_g)
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PB_g)
                        )

                        + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_h * PB_0 * PB_1 * PB_g * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_1 * PA_h * PB_0 * PB_1)
                            + delta[b1][m] * (PA_0 * PA_1 * PA_h * PB_0 * PB_g)
                            + delta[b0][m] * (PA_0 * PA_1 * PA_h * PB_1 * PB_g)
                            + delta[h][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_g)
                            + delta[a1][m] * (PA_0 * PA_h * PB_0 * PB_1 * PB_g)
                            + delta[a0][m] * (PA_1 * PA_h * PB_0 * PB_1 * PB_g)
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a1][h] * delta[b1][g] * (PA_0 * PB_0 * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PB_1 * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PB_0 * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PB_1 * PC[m])
                        )

                        + (
                            (delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g]) * (PA_0)
                            + (delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g]) * (PA_1)
                            + (delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                        )

                    )

                    + F7_t[2] * (

                        (-3.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * (-2.0) + PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m] * (-2.0) + PA_0 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-2.0) + PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-2.0) + PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m] * (-2.0) + PA_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m] * (-2.0) + PA_1 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-2.0) + PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-2.0) + PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m] * (-2.0) + PA_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] * (-2.0) + PA_h * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] * (-2.0) + PA_h * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_h * PB_g * PC[m] * (-2.0) + PA_h * PC[g] * PC[m] * (-1.0) + PB_g * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * (-2.0) + PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m] * (-2.0) + PB_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m] * (-2.0) + PB_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b1] * PC[m] * (-1.0))
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h * (-2.0) + PC[h] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g * (-2.0) + PC[g] * (-1.0))
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_1 * PA_h * PB_0 * PC[m] + PA_0 * PA_1 * PA_h * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[a1] * PC[m] + PA_1 * PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PB_1 * PC[m] + PA_0 * PA_1 * PA_h * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[a1] * PC[m] + PA_1 * PA_h * PB_1 * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_h * PB_g * PC[m] + PA_0 * PA_1 * PA_h * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[h] * PC[m] + PA_0 * PA_h * PB_g * PC[a1] * PC[m] + PA_1 * PA_h * PB_g * PC[a0] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PB_g * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b0] * PC[m] + PA_0 * PB_0 * PB_g * PC[a1] * PC[m] + PA_1 * PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PB_g * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[a1] * PC[m] + PA_1 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PB_1 * PC[m] + PA_0 * PA_h * PB_0 * PC[b1] * PC[m] + PA_0 * PA_h * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_h * PB_0 * PB_g * PC[m] + PA_0 * PA_h * PB_0 * PC[g] * PC[m] + PA_0 * PA_h * PB_g * PC[b0] * PC[m] + PA_0 * PB_0 * PB_g * PC[h] * PC[m] + PA_h * PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_h * PB_1 * PB_g * PC[m] + PA_0 * PA_h * PB_1 * PC[g] * PC[m] + PA_0 * PA_h * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[h] * PC[m] + PA_h * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PB_g * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PB_1 * PC[m] + PA_1 * PA_h * PB_0 * PC[b1] * PC[m] + PA_1 * PA_h * PB_1 * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a1] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_h * PB_0 * PB_g * PC[m] + PA_1 * PA_h * PB_0 * PC[g] * PC[m] + PA_1 * PA_h * PB_g * PC[b0] * PC[m] + PA_1 * PB_0 * PB_g * PC[h] * PC[m] + PA_h * PB_0 * PB_g * PC[a1] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_h * PB_1 * PB_g * PC[m] + PA_1 * PA_h * PB_1 * PC[g] * PC[m] + PA_1 * PA_h * PB_g * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[h] * PC[m] + PA_h * PB_1 * PB_g * PC[a1] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PB_g * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a1] * PC[m])
                            + delta[a0][a1] * (PA_h * PB_0 * PB_1 * PB_g * PC[m] + PA_h * PB_0 * PB_1 * PC[g] * PC[m] + PA_h * PB_0 * PB_g * PC[b1] * PC[m] + PA_h * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PA_h + PA_0 * PA_1 * PC[h] + PA_0 * PA_h * PC[a1] + PA_1 * PA_h * PC[a0])
                            + (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PB_g + PA_0 * PA_1 * PC[g] + PA_0 * PB_g * PC[a1] + PA_1 * PB_g * PC[a0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PB_0 + PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PB_1 + PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_h * PB_1 * PC[a0])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_h * PB_g + PA_0 * PA_h * PC[g] + PA_0 * PB_g * PC[h] + PA_h * PB_g * PC[a0])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PB_g + PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PB_g + PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[a0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PB_0 + PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_h * PB_0 * PC[a1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PB_1 + PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_h * PB_1 * PC[a1])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_h * PB_g + PA_1 * PA_h * PC[g] + PA_1 * PB_g * PC[h] + PA_h * PB_g * PC[a1])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PB_g + PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PB_0 * PB_g * PC[a1])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PB_g + PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PB_1 * PB_g * PC[a1])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PB_1 + PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PB_0 * PB_1 * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_h * PB_0 * PB_g + PA_h * PB_0 * PC[g] + PA_h * PB_g * PC[b0] + PB_0 * PB_g * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_h * PB_1 * PB_g + PA_h * PB_1 * PC[g] + PA_h * PB_g * PC[b1] + PB_1 * PB_g * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PB_g + PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[b0])
                        )

                        + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_h * PB_0 * PB_1 * PC[g] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_0 * PB_g * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_1 * PB_g * PC[b0] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_1 * PB_g * PC[a1] * PC[m]
                            + PA_1 * PA_h * PB_0 * PB_1 * PB_g * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_1 * PA_h * PB_0 * PC[b1] + PA_0 * PA_1 * PA_h * PB_1 * PC[b0] + PA_0 * PA_1 * PB_0 * PB_1 * PC[h] + PA_0 * PA_h * PB_0 * PB_1 * PC[a1] + PA_1 * PA_h * PB_0 * PB_1 * PC[a0])
                            + delta[b1][m] * (PA_0 * PA_1 * PA_h * PB_0 * PC[g] + PA_0 * PA_1 * PA_h * PB_g * PC[b0] + PA_0 * PA_1 * PB_0 * PB_g * PC[h] + PA_0 * PA_h * PB_0 * PB_g * PC[a1] + PA_1 * PA_h * PB_0 * PB_g * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_h * PB_1 * PC[g] + PA_0 * PA_1 * PA_h * PB_g * PC[b1] + PA_0 * PA_1 * PB_1 * PB_g * PC[h] + PA_0 * PA_h * PB_1 * PB_g * PC[a1] + PA_1 * PA_h * PB_1 * PB_g * PC[a0])
                            + delta[h][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[g] + PA_0 * PA_1 * PB_0 * PB_g * PC[b1] + PA_0 * PA_1 * PB_1 * PB_g * PC[b0] + PA_0 * PB_0 * PB_1 * PB_g * PC[a1] + PA_1 * PB_0 * PB_1 * PB_g * PC[a0])
                            + delta[a1][m] * (PA_0 * PA_h * PB_0 * PB_1 * PC[g] + PA_0 * PA_h * PB_0 * PB_g * PC[b1] + PA_0 * PA_h * PB_1 * PB_g * PC[b0] + PA_0 * PB_0 * PB_1 * PB_g * PC[h] + PA_h * PB_0 * PB_1 * PB_g * PC[a0])
                            + delta[a0][m] * (PA_1 * PA_h * PB_0 * PB_1 * PC[g] + PA_1 * PA_h * PB_0 * PB_g * PC[b1] + PA_1 * PA_h * PB_1 * PB_g * PC[b0] + PA_1 * PB_0 * PB_1 * PB_g * PC[h] + PA_h * PB_0 * PB_1 * PB_g * PC[a1])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            delta[a1][h] * delta[b1][g] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                        )

                        + (-1.0) * (
                            (delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b1][g] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PA_h * PC[m] + PA_1 * PC[h] * PC[m] + PA_h * PC[a1] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][h] * delta[b1][g] * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PB_g * PC[m] + PA_1 * PC[g] * PC[m] + PB_g * PC[a1] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1])
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PA_h + PC[h])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g]) * (PA_0 + PC[a0])
                            + (delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g]) * (PA_1 + PC[a1])
                            + (delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0])
                            + (delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1])
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g + PC[g])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[a1] * PC[m] + PA_1 * PA_h * PB_0 * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[a1] * PC[m] + PA_1 * PA_h * PB_1 * PC[a0] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a1][h] * (PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a1] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][g] * delta[h][m] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + delta[b0][g] * delta[h][m] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PA_h * PC[a1] + PA_1 * PA_h * PC[a0])
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_h * PB_0 * PC[a0])
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_h * PB_1 * PC[a0])
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_h * PB_0 * PC[a1])
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_h * PB_1 * PC[a1])
                        )

                        + 2.0 * a_j * (
                            delta[a1][h] * delta[g][m] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[a0])
                            + delta[a0][h] * delta[g][m] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PB_0 * PB_g * PC[a1])
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PB_1 * PB_g * PC[a1])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[b0])
                        )

                    )

                    + F7_t[3] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m])
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m] + PC[a0] * PC[a1] * PC[m])
                            + delta[a1][h] * delta[b1][g] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PC[h] * PC[m] + PA_h * PC[a0] * PC[m] + PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m] + PC[a1] * PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m] + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PC[h] * PC[m] + PA_h * PC[a1] * PC[m] + PC[a1] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PA_h * PC[b0] * PC[m] + PB_0 * PC[h] * PC[m] + PC[b0] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PA_h * PC[b1] * PC[m] + PB_1 * PC[h] * PC[m] + PC[b1] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][h] * delta[b1][g] * (PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m] + PC[a0] * PC[b0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m] + PC[a0] * PC[b1] * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m] + PC[a0] * PC[g] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m] + PC[a1] * PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m] + PC[a1] * PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PA_1 * PC[g] * PC[m] + PB_g * PC[a1] * PC[m] + PC[a1] * PC[g] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m] + PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m] + PC[b0] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m] + PC[b1] * PC[g] * PC[m])
                        )

                        + (-1.0) / (a_i + a_j) * a_i * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PC[h])
                        )

                        + (-1.0) / (a_i + a_j) * a_j * (
                            (delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g])
                        )

                        + (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][h] * (PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[m])
                        )

                        + (-2.0) * a_i * (
                            delta[b1][g] * delta[h][m] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + delta[b0][g] * delta[h][m] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PA_h * PC[a0] * PC[a1])
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h])
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PC[b0] * PC[h] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h])
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PC[b1] * PC[h] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h])
                        )

                        + (-2.0) * a_j * (
                            delta[a1][h] * delta[g][m] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + delta[a1][h] * delta[b1][m] * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + delta[a1][h] * delta[b0][m] * (PA_0 * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1])
                            + delta[a0][h] * delta[g][m] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + delta[a0][h] * delta[b1][m] * (PA_1 * PC[b0] * PC[g] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0])
                            + delta[a0][h] * delta[b0][m] * (PA_1 * PC[b1] * PC[g] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1])
                        )

                        + 3.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PA_h * PC[m] + PA_0 * PC[h] * PC[m] * 2.0 + PA_h * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[h] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] * 2.0 + PB_g * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PA_h * PC[m] + PA_1 * PC[h] * PC[m] * 2.0 + PA_h * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m] + PA_1 * PC[g] * PC[m] * 2.0 + PB_g * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PB_0 * PC[m] + PA_h * PC[b0] * PC[m] * 2.0 + PB_0 * PC[h] * PC[m] * 2.0 + PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PB_1 * PC[m] + PA_h * PC[b1] * PC[m] * 2.0 + PB_1 * PC[h] * PC[m] * 2.0 + PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_h * PB_g * PC[m] + PA_h * PC[g] * PC[m] * 2.0 + PB_g * PC[h] * PC[m] * 2.0 + PC[g] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] * 2.0 + PB_g * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] * 2.0 + PB_g * PC[b1] * PC[m] * 2.0 + PC[b1] * PC[g] * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PA_h + PC[h] * 2.0)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0] * 2.0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g + PC[g] * 2.0)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[b1][g] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PB_0 * PC[a1] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PB_0 * PC[a0] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PA_h * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PB_1 * PC[a1] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PB_1 * PC[a0] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PA_h * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[h] * PC[m] + PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_h * PB_g * PC[a1] * PC[m] + PA_0 * PA_h * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[h] * PC[m] + PA_1 * PA_h * PB_g * PC[a0] * PC[m] + PA_1 * PA_h * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[h] * PC[m] + PA_h * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_g * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b1] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PB_g * PC[a1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[a0] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PB_0 * PC[b1] * PC[m] + PA_0 * PA_h * PB_1 * PC[b0] * PC[m] + PA_0 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a0] * PC[m] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_h * PB_0 * PC[g] * PC[m] + PA_0 * PA_h * PB_g * PC[b0] * PC[m] + PA_0 * PA_h * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[h] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_g * PC[a0] * PC[m] + PA_h * PB_0 * PC[a0] * PC[g] * PC[m] + PA_h * PB_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_h * PB_1 * PC[g] * PC[m] + PA_0 * PA_h * PB_g * PC[b1] * PC[m] + PA_0 * PA_h * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PB_g * PC[h] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PB_g * PC[a0] * PC[m] + PA_h * PB_1 * PC[a0] * PC[g] * PC[m] + PA_h * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PB_g * PC[b0] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_g * PC[a0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PB_0 * PC[b1] * PC[m] + PA_1 * PA_h * PB_1 * PC[b0] * PC[m] + PA_1 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_1 * PC[a1] * PC[m] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_h * PB_0 * PC[g] * PC[m] + PA_1 * PA_h * PB_g * PC[b0] * PC[m] + PA_1 * PA_h * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[h] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PB_g * PC[a1] * PC[m] + PA_h * PB_0 * PC[a1] * PC[g] * PC[m] + PA_h * PB_g * PC[a1] * PC[b0] * PC[m] + PB_0 * PB_g * PC[a1] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_h * PB_1 * PC[g] * PC[m] + PA_1 * PA_h * PB_g * PC[b1] * PC[m] + PA_1 * PA_h * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PB_g * PC[h] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PB_g * PC[a1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[g] * PC[m] + PA_h * PB_g * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a1] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PB_g * PC[b0] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_g * PC[a1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[m])
                            + delta[a0][a1] * (PA_h * PB_0 * PB_1 * PC[g] * PC[m] + PA_h * PB_0 * PB_g * PC[b1] * PC[m] + PA_h * PB_0 * PC[b1] * PC[g] * PC[m] + PA_h * PB_1 * PB_g * PC[b0] * PC[m] + PA_h * PB_1 * PC[b0] * PC[g] * PC[m] + PA_h * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_g * PC[h] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PC[b1] * PC[h] * PC[m] + PB_1 * PB_g * PC[b0] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PB_g * PC[a1] + PA_0 * PC[a1] * PC[g] + PA_1 * PB_g * PC[a0] + PA_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PA_h * PC[a1] + PA_0 * PC[a1] * PC[h] + PA_1 * PA_h * PC[a0] + PA_1 * PC[a0] * PC[h] + PA_h * PC[a0] * PC[a1])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PA_h * PC[b0] + PA_0 * PB_0 * PC[h] + PA_0 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PA_h * PC[b1] + PA_0 * PB_1 * PC[h] + PA_0 * PC[b1] * PC[h] + PA_h * PB_1 * PC[a0] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PA_h * PC[g] + PA_0 * PB_g * PC[h] + PA_0 * PC[g] * PC[h] + PA_h * PB_g * PC[a0] + PA_h * PC[a0] * PC[g] + PB_g * PC[a0] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PA_0 * PC[b0] * PC[g] + PB_0 * PB_g * PC[a0] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PA_0 * PC[b1] * PC[g] + PB_1 * PB_g * PC[a0] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PA_h * PC[b0] + PA_1 * PB_0 * PC[h] + PA_1 * PC[b0] * PC[h] + PA_h * PB_0 * PC[a1] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PA_h * PC[b1] + PA_1 * PB_1 * PC[h] + PA_1 * PC[b1] * PC[h] + PA_h * PB_1 * PC[a1] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PA_h * PC[g] + PA_1 * PB_g * PC[h] + PA_1 * PC[g] * PC[h] + PA_h * PB_g * PC[a1] + PA_h * PC[a1] * PC[g] + PB_g * PC[a1] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PA_1 * PC[b0] * PC[g] + PB_0 * PB_g * PC[a1] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PA_1 * PC[b1] * PC[g] + PB_1 * PB_g * PC[a1] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PB_0 * PC[b1] + PA_h * PB_1 * PC[b0] + PA_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[h] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_h * PB_0 * PC[g] + PA_h * PB_g * PC[b0] + PA_h * PC[b0] * PC[g] + PB_0 * PB_g * PC[h] + PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_h * PB_1 * PC[g] + PA_h * PB_g * PC[b1] + PA_h * PC[b1] * PC[g] + PB_1 * PB_g * PC[h] + PB_1 * PC[g] * PC[h] + PB_g * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PB_g * PC[b0] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_h * PB_0 * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_1 * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PA_h * PB_g * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_1 * PC[a1] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_0 * PB_g * PC[a1] * PC[b1] * PC[m]
                            + PA_0 * PA_h * PB_1 * PB_g * PC[a1] * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_g * PC[a1] * PC[h] * PC[m]
                            + PA_1 * PA_h * PB_0 * PB_1 * PC[a0] * PC[g] * PC[m]
                            + PA_1 * PA_h * PB_0 * PB_g * PC[a0] * PC[b1] * PC[m]
                            + PA_1 * PA_h * PB_1 * PB_g * PC[a0] * PC[b0] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PB_g * PC[a0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PB_g * PC[a0] * PC[a1] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[g][m] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[b1] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[h] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[h] + PA_0 * PA_h * PB_0 * PC[a1] * PC[b1] + PA_0 * PA_h * PB_1 * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[h] + PA_1 * PA_h * PB_0 * PC[a0] * PC[b1] + PA_1 * PA_h * PB_1 * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[h] + PA_h * PB_0 * PB_1 * PC[a0] * PC[a1])
                            + delta[b1][m] * (PA_0 * PA_1 * PA_h * PC[b0] * PC[g] + PA_0 * PA_1 * PB_0 * PC[g] * PC[h] + PA_0 * PA_1 * PB_g * PC[b0] * PC[h] + PA_0 * PA_h * PB_0 * PC[a1] * PC[g] + PA_0 * PA_h * PB_g * PC[a1] * PC[b0] + PA_0 * PB_0 * PB_g * PC[a1] * PC[h] + PA_1 * PA_h * PB_0 * PC[a0] * PC[g] + PA_1 * PA_h * PB_g * PC[a0] * PC[b0] + PA_1 * PB_0 * PB_g * PC[a0] * PC[h] + PA_h * PB_0 * PB_g * PC[a0] * PC[a1])
                            + delta[b0][m] * (PA_0 * PA_1 * PA_h * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[g] * PC[h] + PA_0 * PA_1 * PB_g * PC[b1] * PC[h] + PA_0 * PA_h * PB_1 * PC[a1] * PC[g] + PA_0 * PA_h * PB_g * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_g * PC[a1] * PC[h] + PA_1 * PA_h * PB_1 * PC[a0] * PC[g] + PA_1 * PA_h * PB_g * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_g * PC[a0] * PC[h] + PA_h * PB_1 * PB_g * PC[a0] * PC[a1])
                            + delta[h][m] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] + PA_0 * PA_1 * PB_g * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] + PA_0 * PB_0 * PB_g * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_g * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] + PA_1 * PB_0 * PB_g * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_g * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_g * PC[a0] * PC[a1])
                            + delta[a1][m] * (PA_0 * PA_h * PB_0 * PC[b1] * PC[g] + PA_0 * PA_h * PB_1 * PC[b0] * PC[g] + PA_0 * PA_h * PB_g * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_0 * PB_g * PC[b1] * PC[h] + PA_0 * PB_1 * PB_g * PC[b0] * PC[h] + PA_h * PB_0 * PB_1 * PC[a0] * PC[g] + PA_h * PB_0 * PB_g * PC[a0] * PC[b1] + PA_h * PB_1 * PB_g * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_g * PC[a0] * PC[h])
                            + delta[a0][m] * (PA_1 * PA_h * PB_0 * PC[b1] * PC[g] + PA_1 * PA_h * PB_1 * PC[b0] * PC[g] + PA_1 * PA_h * PB_g * PC[b0] * PC[b1] + PA_1 * PB_0 * PB_1 * PC[g] * PC[h] + PA_1 * PB_0 * PB_g * PC[b1] * PC[h] + PA_1 * PB_1 * PB_g * PC[b0] * PC[h] + PA_h * PB_0 * PB_1 * PC[a1] * PC[g] + PA_h * PB_0 * PB_g * PC[a1] * PC[b1] + PA_h * PB_1 * PB_g * PC[a1] * PC[b0] + PB_0 * PB_1 * PB_g * PC[a1] * PC[h])
                        )

                        + 2.0 * (a_i + a_j) * (
                            delta[a1][h] * delta[b1][g] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PC[a0] * PC[b1] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[a1] * PC[b1] * PC[m])
                        )

                    )

                    + F7_t[4] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PC[h] * PC[m] * (-1.0) + PA_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PA_h * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[h] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PA_h * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[h] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_h * PC[g] * PC[m] * (-1.0) + PB_g * PC[h] * PC[m] * (-1.0) + PC[g] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-2.0))
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g])
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PA_h * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[h] * PC[m] + PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PA_h * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_h * PB_g * PC[a0] * PC[a1] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[a1][g] * (PA_0 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PA_h * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b0] * PC[h] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_h * PB_0 * PC[a0] * PC[g] * PC[m] + PA_h * PB_g * PC[a0] * PC[b0] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PA_h * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b1] * PC[h] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_h * PB_1 * PC[a0] * PC[g] * PC[m] + PA_h * PB_g * PC[a0] * PC[b1] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][g] * (PA_1 * PA_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PA_h * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b0] * PC[h] * PC[m] + PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_h * PB_0 * PC[a1] * PC[g] * PC[m] + PA_h * PB_g * PC[a1] * PC[b0] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PA_h * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b1] * PC[h] * PC[m] + PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_h * PB_1 * PC[a1] * PC[g] * PC[m] + PA_h * PB_g * PC[a1] * PC[b1] * PC[m] + PA_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][a1] * (PA_h * PB_0 * PC[b1] * PC[g] * PC[m] + PA_h * PB_1 * PC[b0] * PC[g] * PC[m] + PA_h * PB_g * PC[b0] * PC[b1] * PC[m] + PA_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PC[b1] * PC[h] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PB_g * PC[b0] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PA_h * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PC[b0] * PC[h] + PA_h * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[h] + PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PC[b1] * PC[h] + PA_h * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[h] + PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[g] * PC[h] + PA_h * PC[a0] * PC[g] + PB_g * PC[a0] * PC[h] + PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PC[b0] * PC[g] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PC[b0] * PC[h] + PA_h * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[h] + PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PC[b1] * PC[g] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PC[b1] * PC[h] + PA_h * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[h] + PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[g] * PC[h] + PA_h * PC[a1] * PC[g] + PB_g * PC[a1] * PC[h] + PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PA_h * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PA_h * PC[b0] * PC[g] + PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h] + PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PA_h * PC[b1] * PC[g] + PB_1 * PC[g] * PC[h] + PB_g * PC[b1] * PC[h] + PC[b1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[g])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PA_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_h * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_h * PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_g * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PB_g * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_1 * PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PA_1 * PA_h * PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_h * PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m]
                            + PA_h * PB_0 * PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m]
                            + PA_h * PB_1 * PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_g * PC[a0] * PC[a1] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[b0])
                            + delta[g][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[b1] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[b1] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h])
                            + delta[b1][m] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b0] * PC[g] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b0] * PC[g] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[h] + PA_h * PB_0 * PC[a0] * PC[a1] * PC[g] + PA_h * PB_g * PC[a0] * PC[a1] * PC[b0] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[h])
                            + delta[b0][m] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[h] + PA_0 * PA_h * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[h] + PA_1 * PA_h * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[h] + PA_h * PB_1 * PC[a0] * PC[a1] * PC[g] + PA_h * PB_g * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[h])
                            + delta[a1][m] * (PA_0 * PA_h * PC[b0] * PC[b1] * PC[g] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[h] + PA_h * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_h * PB_1 * PC[a0] * PC[b0] * PC[g] + PA_h * PB_g * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[h])
                            + delta[a0][m] * (PA_1 * PA_h * PC[b0] * PC[b1] * PC[g] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[h] + PA_h * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_h * PB_1 * PC[a1] * PC[b0] * PC[g] + PA_h * PB_g * PC[a1] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[h] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[h])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_i * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + delta[a1][h] * delta[b1][g] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PC[a0] * PC[h] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PC[a1] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[a1][h] * delta[b1][g] * (PC[a0] * PC[b0] * PC[m])
                            + delta[a1][h] * delta[b0][g] * (PC[a0] * PC[b1] * PC[m])
                            + delta[a1][h] * delta[b0][b1] * (PC[a0] * PC[g] * PC[m])
                            + delta[a0][h] * delta[b1][g] * (PC[a1] * PC[b0] * PC[m])
                            + delta[a0][h] * delta[b0][g] * (PC[a1] * PC[b1] * PC[m])
                            + delta[a0][h] * delta[b0][b1] * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 * a_i * (
                            delta[b1][g] * delta[h][m] * (PC[a0] * PC[a1] * PC[b0])
                            + delta[b0][g] * delta[h][m] * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + delta[a1][m] * delta[b1][g] * (PC[a0] * PC[b0] * PC[h])
                            + delta[a1][m] * delta[b0][g] * (PC[a0] * PC[b1] * PC[h])
                            + delta[a0][m] * delta[b1][g] * (PC[a1] * PC[b0] * PC[h])
                            + delta[a0][m] * delta[b0][g] * (PC[a1] * PC[b1] * PC[h])
                        )

                        + 2.0 * a_j * (
                            delta[a1][h] * delta[g][m] * (PC[a0] * PC[b0] * PC[b1])
                            + delta[a1][h] * delta[b1][m] * (PC[a0] * PC[b0] * PC[g])
                            + delta[a1][h] * delta[b0][m] * (PC[a0] * PC[b1] * PC[g])
                            + delta[a0][h] * delta[g][m] * (PC[a1] * PC[b0] * PC[b1])
                            + delta[a0][h] * delta[b1][m] * (PC[a1] * PC[b0] * PC[g])
                            + delta[a0][h] * delta[b0][m] * (PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                        )

                    )

                    + F7_t[5] * (

                        (-4.0) * (a_i + a_j) * a_i * (
                            delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_i * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[g] * PC[h] * PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PA_h * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PA_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PA_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_i * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_h * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_h * PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_h * PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PA_h * PB_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PB_1 * PB_g * PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m]
                        )

                        + 4.0 * a_i * a_j * (
                            delta[h][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[b1])
                            + delta[g][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h])
                            + delta[b1][m] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[h])
                            + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[h] + PA_h * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[h])
                            + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_h * PC[a0] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h])
                            + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] + PA_h * PC[a1] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[h])
                        )

                    )

                    + F7_t[6] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_i * a_j * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_i * a_j * (
                            PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                        )

                        + (-4.0) * a_i * a_j * (
                            delta[h][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h])
                            + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F7_t[7] * (

                        8.0 * (a_i + a_j) * a_i * a_j * (
                            PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

                // Note: minus sign from electric field - electric dipole interaction

                double V_hess_jj = (-1.0) * V_const * dip[m] * (

                    F7_t[1] * (

                        1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g]) * (PC[m] * (-2.0))
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PA_1 * PC[m] * (-1.0))
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PB_g * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PB_h * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PB_g * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PB_h * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PB_1 * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b1][h] * (PB_0 * PB_g * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b1][g] * (PB_0 * PB_h * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b0][h] * (PB_1 * PB_g * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b0][g] * (PB_1 * PB_h * PC[m] * (-1.0))
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PA_0 * (-1.0))
                            + (delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PA_1 * (-1.0))
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PB_g * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PB_h * (-1.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PB_h * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1)
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PB_g)
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PB_h)
                            + delta[a1][m] * delta[g][h] * (PA_0 * PB_0 * PB_1)
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PB_0 * PB_g)
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PB_0 * PB_h)
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PB_1 * PB_g)
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PB_1 * PB_h)
                            + delta[a0][m] * delta[g][h] * (PA_1 * PB_0 * PB_1)
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PB_0 * PB_g)
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PB_0 * PB_h)
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PB_1 * PB_g)
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PB_1 * PB_h)
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_g * PB_h * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PB_g * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PB_h * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PB_g * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PB_h * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_g * PB_h * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_g * PB_h * PC[m])
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_g * PB_h * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PB_g * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PB_h * PC[m])
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_g * PB_h * PC[m])
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_g * PB_h * PC[m])
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_g * PB_h * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0)
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1)
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PB_g)
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PB_h)
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1)
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PB_g)
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PB_h)
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PB_g)
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PB_h)
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_g * PB_h)
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1)
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PB_g)
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PB_h)
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PB_g)
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PB_h)
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_g * PB_h)
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PB_g)
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PB_h)
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_g * PB_h)
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_g * PB_h)
                        )

                        + 1.0 / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PA_1 * PB_0 * PB_1 * PB_g * PB_h * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_g)
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PB_h)
                            + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PB_g * PB_h)
                            + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PB_g * PB_h)
                            + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PB_g * PB_h)
                            + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PB_g * PB_h)
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m])
                        )

                        + (
                            (delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0)
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1)
                        )

                    )

                    + F7_t[2] * (

                        (-3.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * (-2.0) + PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] * (-2.0) + PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] * (-2.0) + PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m] * (-2.0) + PA_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m] * (-2.0) + PA_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a0] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] * (-2.0) + PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] * (-2.0) + PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m] * (-2.0) + PA_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m] * (-2.0) + PA_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] * (-2.0) + PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m] * (-2.0) + PB_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m] * (-2.0) + PB_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m] * (-2.0) + PB_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m] * (-2.0) + PB_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_g * PB_h * PC[m] * (-2.0) + PB_g * PC[h] * PC[m] * (-1.0) + PB_h * PC[g] * PC[m] * (-1.0))
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g * (-2.0) + PC[g] * (-1.0))
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h * (-2.0) + PC[h] * (-1.0))
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[m] + PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PB_g * PC[m] + PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b0] * PC[m] + PA_0 * PB_0 * PB_g * PC[a1] * PC[m] + PA_1 * PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PB_h * PC[m] + PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b0] * PC[m] + PA_0 * PB_0 * PB_h * PC[a1] * PC[m] + PA_1 * PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PB_g * PC[m] + PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[a1] * PC[m] + PA_1 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PB_h * PC[m] + PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[a1] * PC[m] + PA_1 * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_g * PB_h * PC[m] + PA_0 * PA_1 * PB_g * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[g] * PC[m] + PA_0 * PB_g * PB_h * PC[a1] * PC[m] + PA_1 * PB_g * PB_h * PC[a0] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PB_g * PC[m] + PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PB_h * PC[m] + PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a0] * PC[m])
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_g * PB_h * PC[m] + PA_0 * PB_0 * PB_g * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[g] * PC[m] + PA_0 * PB_g * PB_h * PC[b0] * PC[m] + PB_0 * PB_g * PB_h * PC[a0] * PC[m])
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_g * PB_h * PC[m] + PA_0 * PB_1 * PB_g * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[g] * PC[m] + PA_0 * PB_g * PB_h * PC[b1] * PC[m] + PB_1 * PB_g * PB_h * PC[a0] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PB_g * PC[m] + PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[b0] * PC[m] + PB_0 * PB_1 * PB_g * PC[a1] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PB_h * PC[m] + PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[b0] * PC[m] + PB_0 * PB_1 * PB_h * PC[a1] * PC[m])
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_g * PB_h * PC[m] + PA_1 * PB_0 * PB_g * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[g] * PC[m] + PA_1 * PB_g * PB_h * PC[b0] * PC[m] + PB_0 * PB_g * PB_h * PC[a1] * PC[m])
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_g * PB_h * PC[m] + PA_1 * PB_1 * PB_g * PC[h] * PC[m] + PA_1 * PB_1 * PB_h * PC[g] * PC[m] + PA_1 * PB_g * PB_h * PC[b1] * PC[m] + PB_1 * PB_g * PB_h * PC[a1] * PC[m])
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_g * PB_h * PC[m] + PB_0 * PB_1 * PB_g * PC[h] * PC[m] + PB_0 * PB_1 * PB_h * PC[g] * PC[m] + PB_0 * PB_g * PB_h * PC[b1] * PC[m] + PB_1 * PB_g * PB_h * PC[b0] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PB_g + PA_0 * PA_1 * PC[g] + PA_0 * PB_g * PC[a1] + PA_1 * PB_g * PC[a0])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PB_h + PA_0 * PA_1 * PC[h] + PA_0 * PB_h * PC[a1] + PA_1 * PB_h * PC[a0])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PB_g + PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PB_h + PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PB_g + PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[a0])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PB_h + PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[a0])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_g * PB_h + PA_0 * PB_g * PC[h] + PA_0 * PB_h * PC[g] + PB_g * PB_h * PC[a0])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PB_g + PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PB_0 * PB_g * PC[a1])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PB_h + PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PB_0 * PB_h * PC[a1])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PB_g + PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PB_1 * PB_g * PC[a1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PB_h + PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PB_1 * PB_h * PC[a1])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_g * PB_h + PA_1 * PB_g * PC[h] + PA_1 * PB_h * PC[g] + PB_g * PB_h * PC[a1])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PB_g + PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[b0])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PB_h + PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[b0])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_g * PB_h + PB_0 * PB_g * PC[h] + PB_0 * PB_h * PC[g] + PB_g * PB_h * PC[b0])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_g * PB_h + PB_1 * PB_g * PC[h] + PB_1 * PB_h * PC[g] + PB_g * PB_h * PC[b1])
                        )

                        + (-1.0) / (a_i + a_j) * (a_i + a_j) * (
                            (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g]) * (PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PA_1 * PB_0 * PB_1 * PB_g * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_1 * PB_h * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_g * PB_h * PC[b1] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PB_g * PB_h * PC[b0] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_g * PB_h * PC[a1] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PB_g * PB_h * PC[a0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[g] + PA_0 * PA_1 * PB_0 * PB_g * PC[b1] + PA_0 * PA_1 * PB_1 * PB_g * PC[b0] + PA_0 * PB_0 * PB_1 * PB_g * PC[a1] + PA_1 * PB_0 * PB_1 * PB_g * PC[a0])
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PB_1 * PC[h] + PA_0 * PA_1 * PB_0 * PB_h * PC[b1] + PA_0 * PA_1 * PB_1 * PB_h * PC[b0] + PA_0 * PB_0 * PB_1 * PB_h * PC[a1] + PA_1 * PB_0 * PB_1 * PB_h * PC[a0])
                            + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PB_g * PC[h] + PA_0 * PA_1 * PB_0 * PB_h * PC[g] + PA_0 * PA_1 * PB_g * PB_h * PC[b0] + PA_0 * PB_0 * PB_g * PB_h * PC[a1] + PA_1 * PB_0 * PB_g * PB_h * PC[a0])
                            + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PB_g * PC[h] + PA_0 * PA_1 * PB_1 * PB_h * PC[g] + PA_0 * PA_1 * PB_g * PB_h * PC[b1] + PA_0 * PB_1 * PB_g * PB_h * PC[a1] + PA_1 * PB_1 * PB_g * PB_h * PC[a0])
                            + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PB_g * PC[h] + PA_0 * PB_0 * PB_1 * PB_h * PC[g] + PA_0 * PB_0 * PB_g * PB_h * PC[b1] + PA_0 * PB_1 * PB_g * PB_h * PC[b0] + PB_0 * PB_1 * PB_g * PB_h * PC[a0])
                            + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PB_g * PC[h] + PA_1 * PB_0 * PB_1 * PB_h * PC[g] + PA_1 * PB_0 * PB_g * PB_h * PC[b1] + PA_1 * PB_1 * PB_g * PB_h * PC[b0] + PB_0 * PB_1 * PB_g * PB_h * PC[a1])
                        )

                        + (-2.0) * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                        )

                        + (-1.0) * (
                            (delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * 2.0)
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g]) * (PC[m] * 4.0)
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] + PA_1 * PC[a0] * PC[m])
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] * 2.0 + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0)
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] + PB_0 * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] + PB_g * PC[a0] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] + PB_h * PC[a0] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] + PB_0 * PC[a1] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] + PB_1 * PC[a1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PB_g * PC[m] + PA_1 * PC[g] * PC[m] + PB_g * PC[a1] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PB_h * PC[m] + PA_1 * PC[h] * PC[m] + PB_h * PC[a1] * PC[m])
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] + PB_1 * PC[b0] * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] + PB_g * PC[b0] * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] + PB_h * PC[b0] * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] + PB_g * PC[b1] * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] + PB_h * PC[b1] * PC[m])
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PA_0 + PC[a0])
                            + (delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 * 2.0 + PC[a0] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PA_1 + PC[a1])
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 * 2.0 + PC[a1] * 2.0)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1])
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PB_g + PC[g])
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PB_h + PC[h])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b0] * PC[m] + PA_0 * PB_0 * PB_g * PC[a1] * PC[m] + PA_1 * PB_0 * PB_g * PC[a0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b0] * PC[m] + PA_0 * PB_0 * PB_h * PC[a1] * PC[m] + PA_1 * PB_0 * PB_h * PC[a0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b1] * PC[m] + PA_0 * PB_1 * PB_g * PC[a1] * PC[m] + PA_1 * PB_1 * PB_g * PC[a0] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b1] * PC[m] + PA_0 * PB_1 * PB_h * PC[a1] * PC[m] + PA_1 * PB_1 * PB_h * PC[a0] * PC[m])
                        )

                        + 2.0 * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PB_g * PC[a1] + PA_1 * PB_g * PC[a0])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PB_h * PC[a1] + PA_1 * PB_h * PC[a0])
                            + delta[a1][m] * delta[g][h] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PB_0 * PB_g * PC[a0])
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PB_0 * PB_h * PC[a0])
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PB_1 * PB_g * PC[a0])
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PB_1 * PB_h * PC[a0])
                            + delta[a0][m] * delta[g][h] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PB_0 * PB_g * PC[a1])
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PB_0 * PB_h * PC[a1])
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PB_1 * PB_g * PC[a1])
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PB_1 * PB_h * PC[a1])
                        )

                    )

                    + F7_t[3] * (

                        1.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g]) * (PC[m] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g]) * (PC[m] * (-2.0))
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-1.0))
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] * (-2.0) + PA_1 * PC[a0] * PC[m] * (-2.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PA_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-1.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PA_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-1.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PA_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-1.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PA_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[g][h] * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b1][h] * (PB_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b1][g] * (PB_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b0][h] * (PB_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-1.0))
                            + delta[a0][a1] * delta[b0][g] * (PB_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-1.0))
                        )

                        + 1.0 / (a_i + a_j) * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h]) * (PC[a0] * (-1.0))
                            + (delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0] * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h]) * (PC[a1] * (-1.0))
                            + (delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1] * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1] * (-1.0))
                            + (delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h]) * (PC[g] * (-1.0))
                            + (delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g]) * (PC[h] * (-1.0))
                        )

                        + (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[m])
                        )

                        + (-2.0) * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[a1])
                            + delta[a1][m] * delta[g][h] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + delta[a1][m] * delta[b1][h] * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + delta[a1][m] * delta[b1][g] * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + delta[a1][m] * delta[b0][h] * (PA_0 * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1])
                            + delta[a1][m] * delta[b0][g] * (PA_0 * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1])
                            + delta[a0][m] * delta[g][h] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + delta[a0][m] * delta[b1][h] * (PA_1 * PC[b0] * PC[g] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0])
                            + delta[a0][m] * delta[b1][g] * (PA_1 * PC[b0] * PC[h] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0])
                            + delta[a0][m] * delta[b0][h] * (PA_1 * PC[b1] * PC[g] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1])
                            + delta[a0][m] * delta[b0][g] * (PA_1 * PC[b1] * PC[h] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1])
                        )

                        + 3.0 / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PA_1 * PC[m] + PA_0 * PC[a1] * PC[m] * 2.0 + PA_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PB_0 * PC[m] + PA_0 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PB_1 * PC[m] + PA_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PB_g * PC[m] + PA_0 * PC[g] * PC[m] * 2.0 + PB_g * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PB_h * PC[m] + PA_0 * PC[h] * PC[m] * 2.0 + PB_h * PC[a0] * PC[m] * 2.0 + PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PB_0 * PC[m] + PA_1 * PC[b0] * PC[m] * 2.0 + PB_0 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PB_1 * PC[m] + PA_1 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PB_g * PC[m] + PA_1 * PC[g] * PC[m] * 2.0 + PB_g * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PB_h * PC[m] + PA_1 * PC[h] * PC[m] * 2.0 + PB_h * PC[a1] * PC[m] * 2.0 + PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PB_1 * PC[m] + PB_0 * PC[b1] * PC[m] * 2.0 + PB_1 * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PB_g * PC[m] + PB_0 * PC[g] * PC[m] * 2.0 + PB_g * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PB_h * PC[m] + PB_0 * PC[h] * PC[m] * 2.0 + PB_h * PC[b0] * PC[m] * 2.0 + PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PB_g * PC[m] + PB_1 * PC[g] * PC[m] * 2.0 + PB_g * PC[b1] * PC[m] * 2.0 + PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PB_h * PC[m] + PB_1 * PC[h] * PC[m] * 2.0 + PB_h * PC[b1] * PC[m] * 2.0 + PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_g * PB_h * PC[m] + PB_g * PC[h] * PC[m] * 2.0 + PB_h * PC[g] * PC[m] * 2.0 + PC[g] * PC[h] * PC[m])
                        )

                        + 1.0 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PA_0 + PC[a0] * 2.0)
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PA_1 + PC[a1] * 2.0)
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PB_0 + PC[b0] * 2.0)
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PB_1 + PC[b1] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PB_g + PC[g] * 2.0)
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PB_h + PC[h] * 2.0)
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[m] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PB_0 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_g * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PB_0 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b0] * PC[m] + PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[a1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[m] + PA_1 * PB_0 * PB_h * PC[a0] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PB_1 * PC[g] * PC[m] + PA_0 * PA_1 * PB_g * PC[b1] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PB_g * PC[a1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PB_g * PC[a0] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PB_1 * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[b1] * PC[m] + PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[a1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PA_1 * PB_1 * PB_h * PC[a0] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PB_g * PC[h] * PC[m] + PA_0 * PA_1 * PB_h * PC[g] * PC[m] + PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PB_h * PC[a1] * PC[m] + PA_0 * PB_g * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[g] * PC[m] + PA_1 * PB_g * PB_h * PC[a0] * PC[m] + PA_1 * PB_g * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[g] * PC[m] + PB_g * PB_h * PC[a0] * PC[a1] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PB_1 * PC[g] * PC[m] + PA_0 * PB_0 * PB_g * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PB_g * PC[b0] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_g * PC[a0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PB_1 * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[b1] * PC[m] + PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[b0] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_h * PC[a0] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[a1][b1] * (PA_0 * PB_0 * PB_g * PC[h] * PC[m] + PA_0 * PB_0 * PB_h * PC[g] * PC[m] + PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PB_h * PC[b0] * PC[m] + PA_0 * PB_g * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PB_h * PC[a0] * PC[m] + PB_0 * PB_g * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[g] * PC[m] + PB_g * PB_h * PC[a0] * PC[b0] * PC[m])
                            + delta[a1][b0] * (PA_0 * PB_1 * PB_g * PC[h] * PC[m] + PA_0 * PB_1 * PB_h * PC[g] * PC[m] + PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PB_h * PC[b1] * PC[m] + PA_0 * PB_g * PC[b1] * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PB_h * PC[a0] * PC[m] + PB_1 * PB_g * PC[a0] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[g] * PC[m] + PB_g * PB_h * PC[a0] * PC[b1] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PB_1 * PC[g] * PC[m] + PA_1 * PB_0 * PB_g * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PB_g * PC[b0] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_g * PC[a1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PB_1 * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[b1] * PC[m] + PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PB_h * PC[b0] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PB_h * PC[a1] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[m])
                            + delta[a0][b1] * (PA_1 * PB_0 * PB_g * PC[h] * PC[m] + PA_1 * PB_0 * PB_h * PC[g] * PC[m] + PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PB_h * PC[b0] * PC[m] + PA_1 * PB_g * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PB_h * PC[a1] * PC[m] + PB_0 * PB_g * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[g] * PC[m] + PB_g * PB_h * PC[a1] * PC[b0] * PC[m])
                            + delta[a0][b0] * (PA_1 * PB_1 * PB_g * PC[h] * PC[m] + PA_1 * PB_1 * PB_h * PC[g] * PC[m] + PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PB_h * PC[b1] * PC[m] + PA_1 * PB_g * PC[b1] * PC[h] * PC[m] + PA_1 * PB_h * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PB_h * PC[a1] * PC[m] + PB_1 * PB_g * PC[a1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a1] * PC[g] * PC[m] + PB_g * PB_h * PC[a1] * PC[b1] * PC[m])
                            + delta[a0][a1] * (PB_0 * PB_1 * PB_g * PC[h] * PC[m] + PB_0 * PB_1 * PB_h * PC[g] * PC[m] + PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PB_h * PC[b1] * PC[m] + PB_0 * PB_g * PC[b1] * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PB_h * PC[b0] * PC[m] + PB_1 * PB_g * PC[b0] * PC[h] * PC[m] + PB_1 * PB_h * PC[b0] * PC[g] * PC[m] + PB_g * PB_h * PC[b0] * PC[b1] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PA_1 * PC[g] + PA_0 * PB_g * PC[a1] + PA_0 * PC[a1] * PC[g] + PA_1 * PB_g * PC[a0] + PA_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[a1])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PA_1 * PC[h] + PA_0 * PB_h * PC[a1] + PA_0 * PC[a1] * PC[h] + PA_1 * PB_h * PC[a0] + PA_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[a1])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PB_0 * PC[g] + PA_0 * PB_g * PC[b0] + PA_0 * PC[b0] * PC[g] + PB_0 * PB_g * PC[a0] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PB_0 * PC[h] + PA_0 * PB_h * PC[b0] + PA_0 * PC[b0] * PC[h] + PB_0 * PB_h * PC[a0] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PB_1 * PC[g] + PA_0 * PB_g * PC[b1] + PA_0 * PC[b1] * PC[g] + PB_1 * PB_g * PC[a0] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PB_1 * PC[h] + PA_0 * PB_h * PC[b1] + PA_0 * PC[b1] * PC[h] + PB_1 * PB_h * PC[a0] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PB_g * PC[h] + PA_0 * PB_h * PC[g] + PA_0 * PC[g] * PC[h] + PB_g * PB_h * PC[a0] + PB_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PB_0 * PC[g] + PA_1 * PB_g * PC[b0] + PA_1 * PC[b0] * PC[g] + PB_0 * PB_g * PC[a1] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PB_0 * PC[h] + PA_1 * PB_h * PC[b0] + PA_1 * PC[b0] * PC[h] + PB_0 * PB_h * PC[a1] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PB_1 * PC[g] + PA_1 * PB_g * PC[b1] + PA_1 * PC[b1] * PC[g] + PB_1 * PB_g * PC[a1] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PB_1 * PC[h] + PA_1 * PB_h * PC[b1] + PA_1 * PC[b1] * PC[h] + PB_1 * PB_h * PC[a1] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PB_g * PC[h] + PA_1 * PB_h * PC[g] + PA_1 * PC[g] * PC[h] + PB_g * PB_h * PC[a1] + PB_g * PC[a1] * PC[h] + PB_h * PC[a1] * PC[g])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PB_1 * PC[g] + PB_0 * PB_g * PC[b1] + PB_0 * PC[b1] * PC[g] + PB_1 * PB_g * PC[b0] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PB_1 * PC[h] + PB_0 * PB_h * PC[b1] + PB_0 * PC[b1] * PC[h] + PB_1 * PB_h * PC[b0] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PB_g * PC[h] + PB_0 * PB_h * PC[g] + PB_0 * PC[g] * PC[h] + PB_g * PB_h * PC[b0] + PB_g * PC[b0] * PC[h] + PB_h * PC[b0] * PC[g])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PB_g * PC[h] + PB_1 * PB_h * PC[g] + PB_1 * PC[g] * PC[h] + PB_g * PB_h * PC[b1] + PB_g * PC[b1] * PC[h] + PB_h * PC[b1] * PC[g])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PA_1 * PB_0 * PB_1 * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_g * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_0 * PB_h * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PB_g * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PB_h * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PA_1 * PB_g * PB_h * PC[b0] * PC[b1] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_g * PC[a1] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PB_h * PC[a1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PB_g * PB_h * PC[a1] * PC[b1] * PC[m]
                            + PA_0 * PB_1 * PB_g * PB_h * PC[a1] * PC[b0] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PB_g * PC[a0] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PB_h * PC[a0] * PC[g] * PC[m]
                            + PA_1 * PB_0 * PB_g * PB_h * PC[a0] * PC[b1] * PC[m]
                            + PA_1 * PB_1 * PB_g * PB_h * PC[a0] * PC[b0] * PC[m]
                            + PB_0 * PB_1 * PB_g * PB_h * PC[a0] * PC[a1] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] + PA_0 * PA_1 * PB_g * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] + PA_0 * PB_0 * PB_g * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_g * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] + PA_1 * PB_0 * PB_g * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_g * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_g * PC[a0] * PC[a1])
                            + delta[g][m] * (PA_0 * PA_1 * PB_0 * PC[b1] * PC[h] + PA_0 * PA_1 * PB_1 * PC[b0] * PC[h] + PA_0 * PA_1 * PB_h * PC[b0] * PC[b1] + PA_0 * PB_0 * PB_1 * PC[a1] * PC[h] + PA_0 * PB_0 * PB_h * PC[a1] * PC[b1] + PA_0 * PB_1 * PB_h * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_1 * PC[a0] * PC[h] + PA_1 * PB_0 * PB_h * PC[a0] * PC[b1] + PA_1 * PB_1 * PB_h * PC[a0] * PC[b0] + PB_0 * PB_1 * PB_h * PC[a0] * PC[a1])
                            + delta[b1][m] * (PA_0 * PA_1 * PB_0 * PC[g] * PC[h] + PA_0 * PA_1 * PB_g * PC[b0] * PC[h] + PA_0 * PA_1 * PB_h * PC[b0] * PC[g] + PA_0 * PB_0 * PB_g * PC[a1] * PC[h] + PA_0 * PB_0 * PB_h * PC[a1] * PC[g] + PA_0 * PB_g * PB_h * PC[a1] * PC[b0] + PA_1 * PB_0 * PB_g * PC[a0] * PC[h] + PA_1 * PB_0 * PB_h * PC[a0] * PC[g] + PA_1 * PB_g * PB_h * PC[a0] * PC[b0] + PB_0 * PB_g * PB_h * PC[a0] * PC[a1])
                            + delta[b0][m] * (PA_0 * PA_1 * PB_1 * PC[g] * PC[h] + PA_0 * PA_1 * PB_g * PC[b1] * PC[h] + PA_0 * PA_1 * PB_h * PC[b1] * PC[g] + PA_0 * PB_1 * PB_g * PC[a1] * PC[h] + PA_0 * PB_1 * PB_h * PC[a1] * PC[g] + PA_0 * PB_g * PB_h * PC[a1] * PC[b1] + PA_1 * PB_1 * PB_g * PC[a0] * PC[h] + PA_1 * PB_1 * PB_h * PC[a0] * PC[g] + PA_1 * PB_g * PB_h * PC[a0] * PC[b1] + PB_1 * PB_g * PB_h * PC[a0] * PC[a1])
                            + delta[a1][m] * (PA_0 * PB_0 * PB_1 * PC[g] * PC[h] + PA_0 * PB_0 * PB_g * PC[b1] * PC[h] + PA_0 * PB_0 * PB_h * PC[b1] * PC[g] + PA_0 * PB_1 * PB_g * PC[b0] * PC[h] + PA_0 * PB_1 * PB_h * PC[b0] * PC[g] + PA_0 * PB_g * PB_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PB_g * PC[a0] * PC[h] + PB_0 * PB_1 * PB_h * PC[a0] * PC[g] + PB_0 * PB_g * PB_h * PC[a0] * PC[b1] + PB_1 * PB_g * PB_h * PC[a0] * PC[b0])
                            + delta[a0][m] * (PA_1 * PB_0 * PB_1 * PC[g] * PC[h] + PA_1 * PB_0 * PB_g * PC[b1] * PC[h] + PA_1 * PB_0 * PB_h * PC[b1] * PC[g] + PA_1 * PB_1 * PB_g * PC[b0] * PC[h] + PA_1 * PB_1 * PB_h * PC[b0] * PC[g] + PA_1 * PB_g * PB_h * PC[b0] * PC[b1] + PB_0 * PB_1 * PB_g * PC[a1] * PC[h] + PB_0 * PB_1 * PB_h * PC[a1] * PC[g] + PB_0 * PB_g * PB_h * PC[a1] * PC[b1] + PB_1 * PB_g * PB_h * PC[a1] * PC[b0])
                        )

                        + 2.0 * (a_i + a_j) * (
                            (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                        )

                    )

                    + F7_t[4] * (

                        (-1.0) / ( (a_i + a_j) * (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[a0][a1] * delta[b0][b1] * delta[g][h] + delta[a0][a1] * delta[b0][g] * delta[b1][h] + delta[a0][a1] * delta[b0][h] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[b1][h] + delta[a0][b0] * delta[a1][h] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[b0][h] + delta[a0][b1] * delta[a1][h] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][h] + delta[a0][g] * delta[a1][b1] * delta[b0][h] + delta[a0][g] * delta[a1][h] * delta[b0][b1] + delta[a0][h] * delta[a1][b0] * delta[b1][g] + delta[a0][h] * delta[a1][b1] * delta[b0][g] + delta[a0][h] * delta[a1][g] * delta[b0][b1]) * (PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PA_0 * PC[a1] * PC[m] * (-1.0) + PA_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[a1] * PC[m] * (-2.0))
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PA_0 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PA_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PA_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PA_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a0] * PC[m] * (-1.0) + PC[a0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PA_1 * PC[b0] * PC[m] * (-1.0) + PB_0 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b0] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PA_1 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PA_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PA_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[a1] * PC[m] * (-1.0) + PC[a1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PB_0 * PC[b1] * PC[m] * (-1.0) + PB_1 * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[b1] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PB_0 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PB_0 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b0] * PC[m] * (-1.0) + PC[b0] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PB_1 * PC[g] * PC[m] * (-1.0) + PB_g * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[g] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PB_1 * PC[h] * PC[m] * (-1.0) + PB_h * PC[b1] * PC[m] * (-1.0) + PC[b1] * PC[h] * PC[m] * (-2.0))
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_g * PC[h] * PC[m] * (-1.0) + PB_h * PC[g] * PC[m] * (-1.0) + PC[g] * PC[h] * PC[m] * (-2.0))
                        )

                        + (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_j * a_j * (
                            (delta[a1][b0] * delta[b1][g] * delta[h][m] + delta[a1][b0] * delta[b1][h] * delta[g][m] + delta[a1][b0] * delta[b1][m] * delta[g][h] + delta[a1][b1] * delta[b0][g] * delta[h][m] + delta[a1][b1] * delta[b0][h] * delta[g][m] + delta[a1][b1] * delta[b0][m] * delta[g][h] + delta[a1][g] * delta[b0][b1] * delta[h][m] + delta[a1][g] * delta[b0][h] * delta[b1][m] + delta[a1][g] * delta[b0][m] * delta[b1][h] + delta[a1][h] * delta[b0][b1] * delta[g][m] + delta[a1][h] * delta[b0][g] * delta[b1][m] + delta[a1][h] * delta[b0][m] * delta[b1][g] + delta[a1][m] * delta[b0][b1] * delta[g][h] + delta[a1][m] * delta[b0][g] * delta[b1][h] + delta[a1][m] * delta[b0][h] * delta[b1][g]) * (PC[a0])
                            + (delta[a0][b0] * delta[b1][g] * delta[h][m] + delta[a0][b0] * delta[b1][h] * delta[g][m] + delta[a0][b0] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[b0][g] * delta[h][m] + delta[a0][b1] * delta[b0][h] * delta[g][m] + delta[a0][b1] * delta[b0][m] * delta[g][h] + delta[a0][g] * delta[b0][b1] * delta[h][m] + delta[a0][g] * delta[b0][h] * delta[b1][m] + delta[a0][g] * delta[b0][m] * delta[b1][h] + delta[a0][h] * delta[b0][b1] * delta[g][m] + delta[a0][h] * delta[b0][g] * delta[b1][m] + delta[a0][h] * delta[b0][m] * delta[b1][g] + delta[a0][m] * delta[b0][b1] * delta[g][h] + delta[a0][m] * delta[b0][g] * delta[b1][h] + delta[a0][m] * delta[b0][h] * delta[b1][g]) * (PC[a1])
                            + (delta[a0][a1] * delta[b1][g] * delta[h][m] + delta[a0][a1] * delta[b1][h] * delta[g][m] + delta[a0][a1] * delta[b1][m] * delta[g][h] + delta[a0][b1] * delta[a1][g] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[g][m] + delta[a0][b1] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b1] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b1][m] + delta[a0][g] * delta[a1][m] * delta[b1][h] + delta[a0][h] * delta[a1][b1] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b1][m] + delta[a0][h] * delta[a1][m] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b1][h] + delta[a0][m] * delta[a1][h] * delta[b1][g]) * (PC[b0])
                            + (delta[a0][a1] * delta[b0][g] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[g][m] + delta[a0][a1] * delta[b0][m] * delta[g][h] + delta[a0][b0] * delta[a1][g] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[g][m] + delta[a0][b0] * delta[a1][m] * delta[g][h] + delta[a0][g] * delta[a1][b0] * delta[h][m] + delta[a0][g] * delta[a1][h] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[g][m] + delta[a0][h] * delta[a1][g] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][g] + delta[a0][m] * delta[a1][b0] * delta[g][h] + delta[a0][m] * delta[a1][g] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][g]) * (PC[b1])
                            + (delta[a0][a1] * delta[b0][b1] * delta[h][m] + delta[a0][a1] * delta[b0][h] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][h] + delta[a0][b0] * delta[a1][b1] * delta[h][m] + delta[a0][b0] * delta[a1][h] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][h] + delta[a0][b1] * delta[a1][b0] * delta[h][m] + delta[a0][b1] * delta[a1][h] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][h] + delta[a0][h] * delta[a1][b0] * delta[b1][m] + delta[a0][h] * delta[a1][b1] * delta[b0][m] + delta[a0][h] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][h] + delta[a0][m] * delta[a1][b1] * delta[b0][h] + delta[a0][m] * delta[a1][h] * delta[b0][b1]) * (PC[g])
                            + (delta[a0][a1] * delta[b0][b1] * delta[g][m] + delta[a0][a1] * delta[b0][g] * delta[b1][m] + delta[a0][a1] * delta[b0][m] * delta[b1][g] + delta[a0][b0] * delta[a1][b1] * delta[g][m] + delta[a0][b0] * delta[a1][g] * delta[b1][m] + delta[a0][b0] * delta[a1][m] * delta[b1][g] + delta[a0][b1] * delta[a1][b0] * delta[g][m] + delta[a0][b1] * delta[a1][g] * delta[b0][m] + delta[a0][b1] * delta[a1][m] * delta[b0][g] + delta[a0][g] * delta[a1][b0] * delta[b1][m] + delta[a0][g] * delta[a1][b1] * delta[b0][m] + delta[a0][g] * delta[a1][m] * delta[b0][b1] + delta[a0][m] * delta[a1][b0] * delta[b1][g] + delta[a0][m] * delta[a1][b1] * delta[b0][g] + delta[a0][m] * delta[a1][g] * delta[b0][b1]) * (PC[h])
                        )

                        + (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PA_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_0 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[m] + PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PB_0 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[m] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m])
                            + delta[b0][g] * (PA_0 * PA_1 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m])
                            + delta[b0][b1] * (PA_0 * PA_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[a1] * PC[h] * PC[m] + PA_0 * PB_h * PC[a1] * PC[g] * PC[m] + PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[a0] * PC[h] * PC[m] + PA_1 * PB_h * PC[a0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a0] * PC[a1] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[g] * PC[m])
                            + delta[a1][h] * (PA_0 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[m] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a1][g] * (PA_0 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_0 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[m] + PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PB_1 * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m])
                            + delta[a1][b1] * (PA_0 * PB_0 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b0] * PC[h] * PC[m] + PA_0 * PB_h * PC[b0] * PC[g] * PC[m] + PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PC[a0] * PC[h] * PC[m] + PB_0 * PB_h * PC[a0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a0] * PC[b0] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m])
                            + delta[a1][b0] * (PA_0 * PB_1 * PC[g] * PC[h] * PC[m] + PA_0 * PB_g * PC[b1] * PC[h] * PC[m] + PA_0 * PB_h * PC[b1] * PC[g] * PC[m] + PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PB_g * PC[a0] * PC[h] * PC[m] + PB_1 * PB_h * PC[a0] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a0] * PC[b1] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][h] * (PA_1 * PB_0 * PC[b1] * PC[g] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[m] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[m] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][g] * (PA_1 * PB_0 * PC[b1] * PC[h] * PC[m] + PA_1 * PB_1 * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PB_1 * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[a0][b1] * (PA_1 * PB_0 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b0] * PC[h] * PC[m] + PA_1 * PB_h * PC[b0] * PC[g] * PC[m] + PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PC[a1] * PC[h] * PC[m] + PB_0 * PB_h * PC[a1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a1] * PC[b0] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[a0][b0] * (PA_1 * PB_1 * PC[g] * PC[h] * PC[m] + PA_1 * PB_g * PC[b1] * PC[h] * PC[m] + PA_1 * PB_h * PC[b1] * PC[g] * PC[m] + PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PB_g * PC[a1] * PC[h] * PC[m] + PB_1 * PB_h * PC[a1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[a1] * PC[b1] * PC[m] + PB_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][a1] * (PB_0 * PB_1 * PC[g] * PC[h] * PC[m] + PB_0 * PB_g * PC[b1] * PC[h] * PC[m] + PB_0 * PB_h * PC[b1] * PC[g] * PC[m] + PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PB_g * PC[b0] * PC[h] * PC[m] + PB_1 * PB_h * PC[b0] * PC[g] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_g * PB_h * PC[b0] * PC[b1] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m])
                        )

                        + (-2.0) / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PA_0 * PC[a1] * PC[g] + PA_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PA_0 * PC[a1] * PC[h] + PA_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PA_0 * PC[b0] * PC[g] + PB_0 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PA_0 * PC[b0] * PC[h] + PB_0 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PA_0 * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[g] + PB_g * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PA_0 * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[h] + PB_h * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PA_0 * PC[g] * PC[h] + PB_g * PC[a0] * PC[h] + PB_h * PC[a0] * PC[g] + PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PA_1 * PC[b0] * PC[g] + PB_0 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PA_1 * PC[b0] * PC[h] + PB_0 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PA_1 * PC[b1] * PC[g] + PB_1 * PC[a1] * PC[g] + PB_g * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PA_1 * PC[b1] * PC[h] + PB_1 * PC[a1] * PC[h] + PB_h * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PA_1 * PC[g] * PC[h] + PB_g * PC[a1] * PC[h] + PB_h * PC[a1] * PC[g] + PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PB_0 * PC[b1] * PC[g] + PB_1 * PC[b0] * PC[g] + PB_g * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PB_0 * PC[b1] * PC[h] + PB_1 * PC[b0] * PC[h] + PB_h * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PB_0 * PC[g] * PC[h] + PB_g * PC[b0] * PC[h] + PB_h * PC[b0] * PC[g] + PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PB_1 * PC[g] * PC[h] + PB_g * PC[b1] * PC[h] + PB_h * PC[b1] * PC[g] + PC[b1] * PC[g] * PC[h])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_g * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PA_1 * PB_h * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_g * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PB_h * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PA_0 * PB_1 * PB_g * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PB_h * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PA_0 * PB_g * PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m]
                            + PA_1 * PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PB_1 * PB_g * PC[a0] * PC[b0] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PB_h * PC[a0] * PC[b0] * PC[g] * PC[m]
                            + PA_1 * PB_g * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m]
                            + PB_0 * PB_1 * PB_g * PC[a0] * PC[a1] * PC[h] * PC[m]
                            + PB_0 * PB_1 * PB_h * PC[a0] * PC[a1] * PC[g] * PC[m]
                            + PB_0 * PB_g * PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m]
                            + PB_1 * PB_g * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[b0])
                            + delta[g][m] * (PA_0 * PA_1 * PC[b0] * PC[b1] * PC[h] + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[b1] + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[h] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[b1] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[b0])
                            + delta[b1][m] * (PA_0 * PA_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PB_0 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_g * PC[a1] * PC[b0] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b0] * PC[g] + PA_1 * PB_0 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_g * PC[a0] * PC[b0] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b0] * PC[g] + PB_0 * PB_g * PC[a0] * PC[a1] * PC[h] + PB_0 * PB_h * PC[a0] * PC[a1] * PC[g] + PB_g * PB_h * PC[a0] * PC[a1] * PC[b0])
                            + delta[b0][m] * (PA_0 * PA_1 * PC[b1] * PC[g] * PC[h] + PA_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PA_0 * PB_g * PC[a1] * PC[b1] * PC[h] + PA_0 * PB_h * PC[a1] * PC[b1] * PC[g] + PA_1 * PB_1 * PC[a0] * PC[g] * PC[h] + PA_1 * PB_g * PC[a0] * PC[b1] * PC[h] + PA_1 * PB_h * PC[a0] * PC[b1] * PC[g] + PB_1 * PB_g * PC[a0] * PC[a1] * PC[h] + PB_1 * PB_h * PC[a0] * PC[a1] * PC[g] + PB_g * PB_h * PC[a0] * PC[a1] * PC[b1])
                            + delta[a1][m] * (PA_0 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_0 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_0 * PB_g * PC[b0] * PC[b1] * PC[h] + PA_0 * PB_h * PC[b0] * PC[b1] * PC[g] + PB_0 * PB_1 * PC[a0] * PC[g] * PC[h] + PB_0 * PB_g * PC[a0] * PC[b1] * PC[h] + PB_0 * PB_h * PC[a0] * PC[b1] * PC[g] + PB_1 * PB_g * PC[a0] * PC[b0] * PC[h] + PB_1 * PB_h * PC[a0] * PC[b0] * PC[g] + PB_g * PB_h * PC[a0] * PC[b0] * PC[b1])
                            + delta[a0][m] * (PA_1 * PB_0 * PC[b1] * PC[g] * PC[h] + PA_1 * PB_1 * PC[b0] * PC[g] * PC[h] + PA_1 * PB_g * PC[b0] * PC[b1] * PC[h] + PA_1 * PB_h * PC[b0] * PC[b1] * PC[g] + PB_0 * PB_1 * PC[a1] * PC[g] * PC[h] + PB_0 * PB_g * PC[a1] * PC[b1] * PC[h] + PB_0 * PB_h * PC[a1] * PC[b1] * PC[g] + PB_1 * PB_g * PC[a1] * PC[b0] * PC[h] + PB_1 * PB_h * PC[a1] * PC[b0] * PC[g] + PB_g * PB_h * PC[a1] * PC[b0] * PC[b1])
                        )

                        + 2.0 / (a_i + a_j) * (a_i + a_j) * a_j * (
                            delta[b0][b1] * delta[g][h] * (PC[a0] * PC[a1] * PC[m])
                            + (delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m] * 2.0)
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h]) * (PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g]) * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h]) * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g]) * (PC[a1] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[g][h] * (PC[b0] * PC[b1] * PC[m])
                            + delta[a0][a1] * delta[b1][h] * (PC[b0] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b1][g] * (PC[b0] * PC[h] * PC[m])
                            + delta[a0][a1] * delta[b0][h] * (PC[b1] * PC[g] * PC[m])
                            + delta[a0][a1] * delta[b0][g] * (PC[b1] * PC[h] * PC[m])
                        )

                        + 4.0 * (a_i + a_j) * a_j * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m])
                        )

                        + 2.0 * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + delta[a1][m] * delta[g][h] * (PC[a0] * PC[b0] * PC[b1])
                            + delta[a1][m] * delta[b1][h] * (PC[a0] * PC[b0] * PC[g])
                            + delta[a1][m] * delta[b1][g] * (PC[a0] * PC[b0] * PC[h])
                            + delta[a1][m] * delta[b0][h] * (PC[a0] * PC[b1] * PC[g])
                            + delta[a1][m] * delta[b0][g] * (PC[a0] * PC[b1] * PC[h])
                            + delta[a0][m] * delta[g][h] * (PC[a1] * PC[b0] * PC[b1])
                            + delta[a0][m] * delta[b1][h] * (PC[a1] * PC[b0] * PC[g])
                            + delta[a0][m] * delta[b1][g] * (PC[a1] * PC[b0] * PC[h])
                            + delta[a0][m] * delta[b0][h] * (PC[a1] * PC[b1] * PC[g])
                            + delta[a0][m] * delta[b0][g] * (PC[a1] * PC[b1] * PC[h])
                        )

                    )

                    + F7_t[5] * (

                        (-4.0) * (a_i + a_j) * a_j * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                        )

                        + 2.0 / ( (a_i + a_j) * (a_i + a_j) ) * (a_i + a_j) * a_j * a_j * (
                            (delta[b0][b1] * delta[g][h] + delta[b0][g] * delta[b1][h] + delta[b0][h] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[m])
                            + (delta[a1][b1] * delta[g][h] + delta[a1][g] * delta[b1][h] + delta[a1][h] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[m])
                            + (delta[a1][b0] * delta[g][h] + delta[a1][g] * delta[b0][h] + delta[a1][h] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[m])
                            + (delta[a1][b0] * delta[b1][h] + delta[a1][b1] * delta[b0][h] + delta[a1][h] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[m])
                            + (delta[a1][b0] * delta[b1][g] + delta[a1][b1] * delta[b0][g] + delta[a1][g] * delta[b0][b1]) * (PC[a0] * PC[h] * PC[m])
                            + (delta[a0][b1] * delta[g][h] + delta[a0][g] * delta[b1][h] + delta[a0][h] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[m])
                            + (delta[a0][b0] * delta[g][h] + delta[a0][g] * delta[b0][h] + delta[a0][h] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[m])
                            + (delta[a0][b0] * delta[b1][h] + delta[a0][b1] * delta[b0][h] + delta[a0][h] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[m])
                            + (delta[a0][b0] * delta[b1][g] + delta[a0][b1] * delta[b0][g] + delta[a0][g] * delta[b0][b1]) * (PC[a1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[g][h] + delta[a0][g] * delta[a1][h] + delta[a0][h] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[m])
                            + (delta[a0][a1] * delta[b1][h] + delta[a0][b1] * delta[a1][h] + delta[a0][h] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b1][g] + delta[a0][b1] * delta[a1][g] + delta[a0][g] * delta[a1][b1]) * (PC[b0] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][h] + delta[a0][b0] * delta[a1][h] + delta[a0][h] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[m])
                            + (delta[a0][a1] * delta[b0][g] + delta[a0][b0] * delta[a1][g] + delta[a0][g] * delta[a1][b0]) * (PC[b1] * PC[h] * PC[m])
                            + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[g] * PC[h] * PC[m])
                        )

                        + 4.0 / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PA_0 * PC[a1] * PC[b0] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[m] + PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PA_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[m] + PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PA_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PA_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[a1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[a1] * PC[g] * PC[m] + PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PA_0 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[m] + PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PA_0 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_0 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b0] * PC[g] * PC[m] + PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PA_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[a0] * PC[g] * PC[h] * PC[m] + PB_g * PC[a0] * PC[b1] * PC[h] * PC[m] + PB_h * PC[a0] * PC[b1] * PC[g] * PC[m] + PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PA_1 * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[m] + PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PA_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_0 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b0] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b0] * PC[g] * PC[m] + PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PA_1 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[a1] * PC[g] * PC[h] * PC[m] + PB_g * PC[a1] * PC[b1] * PC[h] * PC[m] + PB_h * PC[a1] * PC[b1] * PC[g] * PC[m] + PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PB_0 * PC[b1] * PC[g] * PC[h] * PC[m] + PB_1 * PC[b0] * PC[g] * PC[h] * PC[m] + PB_g * PC[b0] * PC[b1] * PC[h] * PC[m] + PB_h * PC[b0] * PC[b1] * PC[g] * PC[m] + PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + 2.0 / (a_i + a_j) * a_j * a_j * (
                            (delta[b1][g] * delta[h][m] + delta[b1][h] * delta[g][m] + delta[b1][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b0])
                            + (delta[b0][g] * delta[h][m] + delta[b0][h] * delta[g][m] + delta[b0][m] * delta[g][h]) * (PC[a0] * PC[a1] * PC[b1])
                            + (delta[b0][b1] * delta[h][m] + delta[b0][h] * delta[b1][m] + delta[b0][m] * delta[b1][h]) * (PC[a0] * PC[a1] * PC[g])
                            + (delta[b0][b1] * delta[g][m] + delta[b0][g] * delta[b1][m] + delta[b0][m] * delta[b1][g]) * (PC[a0] * PC[a1] * PC[h])
                            + (delta[a1][g] * delta[h][m] + delta[a1][h] * delta[g][m] + delta[a1][m] * delta[g][h]) * (PC[a0] * PC[b0] * PC[b1])
                            + (delta[a1][b1] * delta[h][m] + delta[a1][h] * delta[b1][m] + delta[a1][m] * delta[b1][h]) * (PC[a0] * PC[b0] * PC[g])
                            + (delta[a1][b1] * delta[g][m] + delta[a1][g] * delta[b1][m] + delta[a1][m] * delta[b1][g]) * (PC[a0] * PC[b0] * PC[h])
                            + (delta[a1][b0] * delta[h][m] + delta[a1][h] * delta[b0][m] + delta[a1][m] * delta[b0][h]) * (PC[a0] * PC[b1] * PC[g])
                            + (delta[a1][b0] * delta[g][m] + delta[a1][g] * delta[b0][m] + delta[a1][m] * delta[b0][g]) * (PC[a0] * PC[b1] * PC[h])
                            + (delta[a1][b0] * delta[b1][m] + delta[a1][b1] * delta[b0][m] + delta[a1][m] * delta[b0][b1]) * (PC[a0] * PC[g] * PC[h])
                            + (delta[a0][g] * delta[h][m] + delta[a0][h] * delta[g][m] + delta[a0][m] * delta[g][h]) * (PC[a1] * PC[b0] * PC[b1])
                            + (delta[a0][b1] * delta[h][m] + delta[a0][h] * delta[b1][m] + delta[a0][m] * delta[b1][h]) * (PC[a1] * PC[b0] * PC[g])
                            + (delta[a0][b1] * delta[g][m] + delta[a0][g] * delta[b1][m] + delta[a0][m] * delta[b1][g]) * (PC[a1] * PC[b0] * PC[h])
                            + (delta[a0][b0] * delta[h][m] + delta[a0][h] * delta[b0][m] + delta[a0][m] * delta[b0][h]) * (PC[a1] * PC[b1] * PC[g])
                            + (delta[a0][b0] * delta[g][m] + delta[a0][g] * delta[b0][m] + delta[a0][m] * delta[b0][g]) * (PC[a1] * PC[b1] * PC[h])
                            + (delta[a0][b0] * delta[b1][m] + delta[a0][b1] * delta[b0][m] + delta[a0][m] * delta[b0][b1]) * (PC[a1] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[h][m] + delta[a0][h] * delta[a1][m] + delta[a0][m] * delta[a1][h]) * (PC[b0] * PC[b1] * PC[g])
                            + (delta[a0][a1] * delta[g][m] + delta[a0][g] * delta[a1][m] + delta[a0][m] * delta[a1][g]) * (PC[b0] * PC[b1] * PC[h])
                            + (delta[a0][a1] * delta[b1][m] + delta[a0][b1] * delta[a1][m] + delta[a0][m] * delta[a1][b1]) * (PC[b0] * PC[g] * PC[h])
                            + (delta[a0][a1] * delta[b0][m] + delta[a0][b0] * delta[a1][m] + delta[a0][m] * delta[a1][b0]) * (PC[b1] * PC[g] * PC[h])
                        )

                        + 8.0 * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_0 * PB_g * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_0 * PB_h * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PA_1 * PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m]
                            + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PB_g * PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m]
                            + PB_0 * PB_h * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m]
                            + PB_1 * PB_g * PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m]
                            + PB_1 * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m]
                            + PB_g * PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m]
                        )

                        + 4.0 * a_j * a_j * (
                            delta[h][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[b1])
                            + delta[g][m] * (PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[h] + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1])
                            + delta[b1][m] * (PA_0 * PC[a1] * PC[b0] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PB_0 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[g])
                            + delta[b0][m] * (PA_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PA_1 * PC[a0] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[a1] * PC[g] * PC[h] + PB_g * PC[a0] * PC[a1] * PC[b1] * PC[h] + PB_h * PC[a0] * PC[a1] * PC[b1] * PC[g])
                            + delta[a1][m] * (PA_0 * PC[b0] * PC[b1] * PC[g] * PC[h] + PB_0 * PC[a0] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a0] * PC[b0] * PC[g] * PC[h] + PB_g * PC[a0] * PC[b0] * PC[b1] * PC[h] + PB_h * PC[a0] * PC[b0] * PC[b1] * PC[g])
                            + delta[a0][m] * (PA_1 * PC[b0] * PC[b1] * PC[g] * PC[h] + PB_0 * PC[a1] * PC[b1] * PC[g] * PC[h] + PB_1 * PC[a1] * PC[b0] * PC[g] * PC[h] + PB_g * PC[a1] * PC[b0] * PC[b1] * PC[h] + PB_h * PC[a1] * PC[b0] * PC[b1] * PC[g])
                        )

                    )

                    + F7_t[6] * (

                        (-4.0) / (a_i + a_j) * (a_i + a_j) * a_j * a_j * (
                            delta[g][h] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[m])
                            + delta[b1][h] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[m])
                            + delta[b1][g] * (PC[a0] * PC[a1] * PC[b0] * PC[h] * PC[m])
                            + delta[b0][h] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[m])
                            + delta[b0][g] * (PC[a0] * PC[a1] * PC[b1] * PC[h] * PC[m])
                            + delta[b0][b1] * (PC[a0] * PC[a1] * PC[g] * PC[h] * PC[m])
                            + delta[a1][h] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a1][g] * (PC[a0] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a1][b1] * (PC[a0] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a1][b0] * (PC[a0] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][h] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m])
                            + delta[a0][g] * (PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m])
                            + delta[a0][b1] * (PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m])
                            + delta[a0][b0] * (PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m])
                            + delta[a0][a1] * (PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m])
                        )

                        + (-8.0) * (a_i + a_j) * a_j * a_j * (
                            PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h] * PC[m]
                            + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h] * PC[m]
                            + PB_g * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h] * PC[m]
                            + PB_h * PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[m]
                        )

                        + (-4.0) * a_j * a_j * (
                            delta[h][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g])
                            + delta[g][m] * (PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[h])
                            + delta[b1][m] * (PC[a0] * PC[a1] * PC[b0] * PC[g] * PC[h])
                            + delta[b0][m] * (PC[a0] * PC[a1] * PC[b1] * PC[g] * PC[h])
                            + delta[a1][m] * (PC[a0] * PC[b0] * PC[b1] * PC[g] * PC[h])
                            + delta[a0][m] * (PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h])
                        )

                    )

                    + F7_t[7] * (

                        8.0 * (a_i + a_j) * a_j * a_j * (
                            PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[g] * PC[h] * PC[m]
                        )

                    )

                );

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

                        double hess_ii = V_hess_ii * coef_sph * D_sym;
                        double hess_ij = V_hess_ij * coef_sph * D_sym;
                        double hess_ji = V_hess_ji * coef_sph * D_sym;
                        double hess_jj = V_hess_jj * coef_sph * D_sym;

                        V_hess_omp[thread_id].row(i_atom * 3 + g)[i_atom * 3 + h] += hess_ii;
                        V_hess_omp[thread_id].row(i_atom * 3 + g)[j_atom * 3 + h] += hess_ij;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[i_atom * 3 + h] += hess_ji;
                        V_hess_omp[thread_id].row(j_atom * 3 + g)[j_atom * 3 + h] += hess_jj;
                    }
                }
            }
        }
        }
        }
    }

    // auto-generated code ends here

    CDenseMatrix V_hess(natoms * 3, natoms * 3);

    V_hess.zero();

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int a = 0; a < natoms * 3; a++)
        {
            for (int b = 0; b < natoms * 3; b++)
            {
                V_hess.row(a)[b] += V_hess_omp[thread_id].row(a)[b];
            }
        }
    }

    return V_hess;
}

}  // namespace onee
