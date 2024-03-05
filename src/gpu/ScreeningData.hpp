//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#ifndef ScreeningData_hpp
#define ScreeningData_hpp

#include <cstdint>
#include <vector>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class CScreeningData stores screening data for ERI evaluation on GPUs.

 @author X. Li
 */
class CScreeningData
{
    int64_t _num_gpus_per_node{8};

    CDenseMatrix _Q_matrix_ss;
    CDenseMatrix _Q_matrix_sp;
    CDenseMatrix _Q_matrix_sd;
    CDenseMatrix _Q_matrix_pp;
    CDenseMatrix _Q_matrix_pd;

    std::vector<std::vector<uint32_t>> _ss_first_inds_local;
    std::vector<std::vector<uint32_t>> _sp_first_inds_local;
    std::vector<std::vector<uint32_t>> _sd_first_inds_local;
    std::vector<std::vector<uint32_t>> _pp_first_inds_local;
    std::vector<std::vector<uint32_t>> _pd_first_inds_local;

    std::vector<std::vector<uint32_t>> _ss_second_inds_local;
    std::vector<std::vector<uint32_t>> _sp_second_inds_local;
    std::vector<std::vector<uint32_t>> _sd_second_inds_local;
    std::vector<std::vector<uint32_t>> _pp_second_inds_local;
    std::vector<std::vector<uint32_t>> _pd_second_inds_local;

    std::vector<std::vector<double>> _ss_mat_Q_local;
    std::vector<std::vector<double>> _sp_mat_Q_local;
    std::vector<std::vector<double>> _sd_mat_Q_local;
    std::vector<std::vector<double>> _pp_mat_Q_local;
    std::vector<std::vector<double>> _pd_mat_Q_local;

    std::vector<uint32_t> _ss_first_inds;
    std::vector<uint32_t> _sp_first_inds;
    std::vector<uint32_t> _pp_first_inds;

    std::vector<uint32_t> _ss_second_inds;
    std::vector<uint32_t> _sp_second_inds;
    std::vector<uint32_t> _pp_second_inds;

    std::vector<double> _ss_mat_Q;
    std::vector<double> _sp_mat_Q;
    std::vector<double> _pp_mat_Q;

    std::vector<double> _ss_mat_D;
    std::vector<double> _sp_mat_D;
    std::vector<double> _pp_mat_D;

    double _ss_max_D;
    double _sp_max_D;
    double _pp_max_D;

    std::vector<double>   _mat_Q_for_K_ss;
    std::vector<double>   _mat_Q_for_K_sp;
    std::vector<double>   _mat_Q_for_K_ps;
    std::vector<double>   _mat_Q_for_K_pp;

    std::vector<uint32_t> _density_inds_for_K_ss;
    std::vector<uint32_t> _density_inds_for_K_sp;
    std::vector<uint32_t> _density_inds_for_K_ps;
    std::vector<uint32_t> _density_inds_for_K_pp;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_ss;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_ss;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_sp;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_sp;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_pp;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_pp;

    auto _computeQMatrices(const CMolecule& molecule, const CMolecularBasis& basis) -> void;

    auto _sortQ(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) -> void;

   public:
    CScreeningData(const CMolecule& molecule, const CMolecularBasis& basis, const int64_t num_gpus_per_node);

    auto getNumGpusPerNode() const -> const int64_t;

    auto getQMatrixSS() const -> const CDenseMatrix&;
    auto getQMatrixSP() const -> const CDenseMatrix&;
    auto getQMatrixSD() const -> const CDenseMatrix&;
    auto getQMatrixPP() const -> const CDenseMatrix&;

    auto sortQD(const int64_t s_prim_count, const int64_t p_prim_count, const std::vector<uint32_t>& s_prim_aoinds, const std::vector<uint32_t>& p_prim_aoinds, const int64_t naos, const double* dens_ptr) -> void;

    auto get_ss_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_ss_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_ss_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_sp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_pp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;

    auto get_ss_first_inds() const -> const std::vector<uint32_t>&;
    auto get_sp_first_inds() const -> const std::vector<uint32_t>&;
    auto get_pp_first_inds() const -> const std::vector<uint32_t>&;

    auto get_ss_second_inds() const -> const std::vector<uint32_t>&;
    auto get_sp_second_inds() const -> const std::vector<uint32_t>&;
    auto get_pp_second_inds() const -> const std::vector<uint32_t>&;

    auto get_ss_mat_Q() const -> const std::vector<double>&;
    auto get_sp_mat_Q() const -> const std::vector<double>&;
    auto get_pp_mat_Q() const -> const std::vector<double>&;

    auto get_ss_mat_D() const -> const std::vector<double>&;
    auto get_sp_mat_D() const -> const std::vector<double>&;
    auto get_pp_mat_D() const -> const std::vector<double>&;

    auto get_ss_max_D() const -> double;
    auto get_sp_max_D() const -> double;
    auto get_pp_max_D() const -> double;

    auto get_mat_Q_for_K_ss() const -> const std::vector<double>&;
    auto get_mat_Q_for_K_sp() const -> const std::vector<double>&;
    auto get_mat_Q_for_K_ps() const -> const std::vector<double>&;
    auto get_mat_Q_for_K_pp() const -> const std::vector<double>&;

    auto get_density_inds_for_K_ss() const -> const std::vector<uint32_t>&;
    auto get_density_inds_for_K_sp() const -> const std::vector<uint32_t>&;
    auto get_density_inds_for_K_ps() const -> const std::vector<uint32_t>&;
    auto get_density_inds_for_K_pp() const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_mat_Q_full(const int64_t s_prim_count, const int64_t p_prim_count) const -> CDenseMatrix;
    auto get_mat_D_abs_full(const int64_t s_prim_count, const int64_t p_prim_count, const std::vector<uint32_t>& s_prim_aoinds, const std::vector<uint32_t>& p_prim_aoinds, const int64_t naos, const double* dens_ptr) const -> CDenseMatrix;

    auto form_mat_Q_and_density_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count) -> void;
    auto form_pair_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count, const CDenseMatrix& Q_prime, const double Q_prime_thresh) -> void;
};

#endif /* ScreeningData_hpp */
