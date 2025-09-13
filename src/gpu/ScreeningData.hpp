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

#ifndef ScreeningData_hpp
#define ScreeningData_hpp

#include <cstdint>
#include <vector>
#include <string>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class CScreeningData stores screening data for ERI evaluation on GPUs.
 */
class CScreeningData
{
    int _rank{0};

    int _nnodes{1};

    int64_t _num_gpus_per_node{8};

    double _pair_threshold{0.0};

    double _density_threshold{0.0};

    std::string _timer_summary;
    std::vector<std::string> _gpu_timer_summary;

    CDenseMatrix _Q_matrix_ss;
    CDenseMatrix _Q_matrix_sp;
    CDenseMatrix _Q_matrix_sd;
    CDenseMatrix _Q_matrix_pp;
    CDenseMatrix _Q_matrix_pd;
    CDenseMatrix _Q_matrix_dd;

    std::vector<std::vector<uint32_t>> _ss_first_inds_local;
    std::vector<std::vector<uint32_t>> _sp_first_inds_local;
    std::vector<std::vector<uint32_t>> _sd_first_inds_local;
    std::vector<std::vector<uint32_t>> _pp_first_inds_local;
    std::vector<std::vector<uint32_t>> _pd_first_inds_local;
    std::vector<std::vector<uint32_t>> _dd_first_inds_local;

    std::vector<std::vector<uint32_t>> _ss_second_inds_local;
    std::vector<std::vector<uint32_t>> _sp_second_inds_local;
    std::vector<std::vector<uint32_t>> _sd_second_inds_local;
    std::vector<std::vector<uint32_t>> _pp_second_inds_local;
    std::vector<std::vector<uint32_t>> _pd_second_inds_local;
    std::vector<std::vector<uint32_t>> _dd_second_inds_local;

    std::vector<std::vector<double>> _ss_mat_Q_local;
    std::vector<std::vector<double>> _sp_mat_Q_local;
    std::vector<std::vector<double>> _sd_mat_Q_local;
    std::vector<std::vector<double>> _pp_mat_Q_local;
    std::vector<std::vector<double>> _pd_mat_Q_local;
    std::vector<std::vector<double>> _dd_mat_Q_local;

    std::vector<std::vector<double>> _ss_pair_data_local;
    std::vector<std::vector<double>> _sp_pair_data_local;
    std::vector<std::vector<double>> _sd_pair_data_local;
    std::vector<std::vector<double>> _pp_pair_data_local;
    std::vector<std::vector<double>> _pd_pair_data_local;
    std::vector<std::vector<double>> _dd_pair_data_local;

    std::vector<uint32_t> _ss_first_inds;
    std::vector<uint32_t> _sp_first_inds;
    std::vector<uint32_t> _sd_first_inds;
    std::vector<uint32_t> _pp_first_inds;
    std::vector<uint32_t> _pd_first_inds;
    std::vector<uint32_t> _dd_first_inds;

    std::vector<uint32_t> _ss_second_inds;
    std::vector<uint32_t> _sp_second_inds;
    std::vector<uint32_t> _sd_second_inds;
    std::vector<uint32_t> _pp_second_inds;
    std::vector<uint32_t> _pd_second_inds;
    std::vector<uint32_t> _dd_second_inds;

    std::vector<double> _ss_mat_Q;
    std::vector<double> _sp_mat_Q;
    std::vector<double> _sd_mat_Q;
    std::vector<double> _pp_mat_Q;
    std::vector<double> _pd_mat_Q;
    std::vector<double> _dd_mat_Q;

    std::vector<double> _ss_mat_D;
    std::vector<double> _sp_mat_D;
    std::vector<double> _sd_mat_D;
    std::vector<double> _pp_mat_D;
    std::vector<double> _pd_mat_D;
    std::vector<double> _dd_mat_D;

    std::vector<double> _ss_pair_data;
    std::vector<double> _sp_pair_data;
    std::vector<double> _sd_pair_data;
    std::vector<double> _pp_pair_data;
    std::vector<double> _pd_pair_data;
    std::vector<double> _dd_pair_data;

    double _ss_max_D;
    double _sp_max_D;
    double _sd_max_D;
    double _pp_max_D;
    double _pd_max_D;
    double _dd_max_D;

    std::vector<double>   _Q_K_ss;
    std::vector<double>   _Q_K_sp;
    std::vector<double>   _Q_K_ps;
    std::vector<double>   _Q_K_sd;
    std::vector<double>   _Q_K_ds;
    std::vector<double>   _Q_K_pp;
    std::vector<double>   _Q_K_pd;
    std::vector<double>   _Q_K_dp;
    std::vector<double>   _Q_K_dd;

    std::vector<uint32_t> _D_inds_K_ss;
    std::vector<uint32_t> _D_inds_K_sp;
    std::vector<uint32_t> _D_inds_K_ps;
    std::vector<uint32_t> _D_inds_K_sd;
    std::vector<uint32_t> _D_inds_K_ds;
    std::vector<uint32_t> _D_inds_K_pp;
    std::vector<uint32_t> _D_inds_K_pd;
    std::vector<uint32_t> _D_inds_K_dp;
    std::vector<uint32_t> _D_inds_K_dd;

    std::vector<uint32_t> _pair_displs_K_ss;
    std::vector<uint32_t> _pair_displs_K_sp;
    std::vector<uint32_t> _pair_displs_K_ps;
    std::vector<uint32_t> _pair_displs_K_sd;
    std::vector<uint32_t> _pair_displs_K_ds;
    std::vector<uint32_t> _pair_displs_K_pp;
    std::vector<uint32_t> _pair_displs_K_pd;
    std::vector<uint32_t> _pair_displs_K_dp;
    std::vector<uint32_t> _pair_displs_K_dd;

    std::vector<uint32_t> _pair_counts_K_ss;
    std::vector<uint32_t> _pair_counts_K_sp;
    std::vector<uint32_t> _pair_counts_K_ps;
    std::vector<uint32_t> _pair_counts_K_sd;
    std::vector<uint32_t> _pair_counts_K_ds;
    std::vector<uint32_t> _pair_counts_K_pp;
    std::vector<uint32_t> _pair_counts_K_pd;
    std::vector<uint32_t> _pair_counts_K_dp;
    std::vector<uint32_t> _pair_counts_K_dd;

    std::vector<double>   _pair_data_K_ss;
    std::vector<double>   _pair_data_K_sp;
    std::vector<double>   _pair_data_K_ps;
    std::vector<double>   _pair_data_K_sd;
    std::vector<double>   _pair_data_K_ds;
    std::vector<double>   _pair_data_K_pp;
    std::vector<double>   _pair_data_K_pd;
    std::vector<double>   _pair_data_K_dp;
    std::vector<double>   _pair_data_K_dd;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_ss;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_ss;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_sp;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_sp;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_sd;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_sd;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_pp;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_pp;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_pd;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_pd;

    std::vector<std::vector<uint32_t>> _local_pair_inds_i_for_K_dd;
    std::vector<std::vector<uint32_t>> _local_pair_inds_k_for_K_dd;

    auto _computeQMatrices(const CMolecule& molecule, const CMolecularBasis& basis) -> void;

    auto _sortQ(const int64_t s_prim_count,
                const int64_t p_prim_count,
                const int64_t d_prim_count,
                const std::vector<double>& s_prim_info,
                const std::vector<double>& p_prim_info,
                const std::vector<double>& d_prim_info) -> void;

   public:
    CScreeningData(const CMolecule& molecule, const CMolecularBasis& basis, const int64_t num_gpus_per_node, const double pair_threshold, const double density_threshold, const int rank, const int nnodes);

    auto getNumGpusPerNode() const -> const int64_t;

    auto setTimerSummary(const std::string& timer_summary) -> void;
    auto getTimerSummary() const -> const std::string;

    auto initGpuTimers(const int64_t num_gpus_per_node) -> void;

    auto setGpuTimerSummary(const int64_t gpu_id, const std::string& timer_summary) -> void;
    auto getGpuTimerSummary() const -> const std::vector<std::string>;

    auto getQMatrixSS() const -> const CDenseMatrix&;
    auto getQMatrixSP() const -> const CDenseMatrix&;
    auto getQMatrixSD() const -> const CDenseMatrix&;
    auto getQMatrixPP() const -> const CDenseMatrix&;
    auto getQMatrixPD() const -> const CDenseMatrix&;
    auto getQMatrixDD() const -> const CDenseMatrix&;

    auto sortQD(const int64_t s_prim_count,
                const int64_t p_prim_count,
                const int64_t d_prim_count,
                const std::vector<uint32_t>& s_prim_aoinds,
                const std::vector<uint32_t>& p_prim_aoinds,
                const std::vector<uint32_t>& d_prim_aoinds,
                const std::vector<double>& s_prim_info,
                const std::vector<double>& p_prim_info,
                const std::vector<double>& d_prim_info,
                const int64_t naos,
                const double* dens_ptr,
                const double eri_threshold) -> void;

    auto findMaxDensities(const int64_t s_prim_count,
                          const int64_t p_prim_count,
                          const int64_t d_prim_count,
                          const std::vector<uint32_t>& s_prim_aoinds,
                          const std::vector<uint32_t>& p_prim_aoinds,
                          const std::vector<uint32_t>& d_prim_aoinds,
                          const int64_t naos,
                          const double* dens_ptr) -> void;

    auto get_ss_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_dd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_ss_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_sd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_pd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_dd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_ss_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_sp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_sd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_pp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_pd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_dd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>&;

    auto get_ss_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_sp_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_sd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_pp_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_pd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;
    auto get_dd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>&;

    auto get_ss_first_inds() const -> const std::vector<uint32_t>&;
    auto get_sp_first_inds() const -> const std::vector<uint32_t>&;
    auto get_sd_first_inds() const -> const std::vector<uint32_t>&;
    auto get_pp_first_inds() const -> const std::vector<uint32_t>&;
    auto get_pd_first_inds() const -> const std::vector<uint32_t>&;
    auto get_dd_first_inds() const -> const std::vector<uint32_t>&;

    auto get_ss_second_inds() const -> const std::vector<uint32_t>&;
    auto get_sp_second_inds() const -> const std::vector<uint32_t>&;
    auto get_sd_second_inds() const -> const std::vector<uint32_t>&;
    auto get_pp_second_inds() const -> const std::vector<uint32_t>&;
    auto get_pd_second_inds() const -> const std::vector<uint32_t>&;
    auto get_dd_second_inds() const -> const std::vector<uint32_t>&;

    auto get_ss_mat_Q() const -> const std::vector<double>&;
    auto get_sp_mat_Q() const -> const std::vector<double>&;
    auto get_sd_mat_Q() const -> const std::vector<double>&;
    auto get_pp_mat_Q() const -> const std::vector<double>&;
    auto get_pd_mat_Q() const -> const std::vector<double>&;
    auto get_dd_mat_Q() const -> const std::vector<double>&;

    auto get_ss_mat_D() const -> const std::vector<double>&;
    auto get_sp_mat_D() const -> const std::vector<double>&;
    auto get_sd_mat_D() const -> const std::vector<double>&;
    auto get_pp_mat_D() const -> const std::vector<double>&;
    auto get_pd_mat_D() const -> const std::vector<double>&;
    auto get_dd_mat_D() const -> const std::vector<double>&;

    auto get_ss_max_D() const -> double;
    auto get_sp_max_D() const -> double;
    auto get_sd_max_D() const -> double;
    auto get_pp_max_D() const -> double;
    auto get_pd_max_D() const -> double;
    auto get_dd_max_D() const -> double;

    auto get_ss_pair_data() const -> const std::vector<double>&;
    auto get_sp_pair_data() const -> const std::vector<double>&;
    auto get_sd_pair_data() const -> const std::vector<double>&;
    auto get_pp_pair_data() const -> const std::vector<double>&;
    auto get_pd_pair_data() const -> const std::vector<double>&;
    auto get_dd_pair_data() const -> const std::vector<double>&;

    auto get_Q_K_ss() const -> const std::vector<double>&;
    auto get_Q_K_sp() const -> const std::vector<double>&;
    auto get_Q_K_ps() const -> const std::vector<double>&;
    auto get_Q_K_sd() const -> const std::vector<double>&;
    auto get_Q_K_ds() const -> const std::vector<double>&;
    auto get_Q_K_pp() const -> const std::vector<double>&;
    auto get_Q_K_pd() const -> const std::vector<double>&;
    auto get_Q_K_dp() const -> const std::vector<double>&;
    auto get_Q_K_dd() const -> const std::vector<double>&;

    auto get_D_inds_K_ss() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_sp() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_ps() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_sd() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_ds() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_pp() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_pd() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_dp() const -> const std::vector<uint32_t>&;
    auto get_D_inds_K_dd() const -> const std::vector<uint32_t>&;

    auto get_pair_displs_K_ss() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_sp() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_ps() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_sd() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_ds() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_pp() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_pd() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_dp() const -> const std::vector<uint32_t>&;
    auto get_pair_displs_K_dd() const -> const std::vector<uint32_t>&;

    auto get_pair_counts_K_ss() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_sp() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_ps() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_sd() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_ds() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_pp() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_pd() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_dp() const -> const std::vector<uint32_t>&;
    auto get_pair_counts_K_dd() const -> const std::vector<uint32_t>&;

    auto get_pair_data_K_ss() const -> const std::vector<double>&;
    auto get_pair_data_K_sp() const -> const std::vector<double>&;
    auto get_pair_data_K_ps() const -> const std::vector<double>&;
    auto get_pair_data_K_sd() const -> const std::vector<double>&;
    auto get_pair_data_K_ds() const -> const std::vector<double>&;
    auto get_pair_data_K_pp() const -> const std::vector<double>&;
    auto get_pair_data_K_pd() const -> const std::vector<double>&;
    auto get_pair_data_K_dp() const -> const std::vector<double>&;
    auto get_pair_data_K_dd() const -> const std::vector<double>&;

    auto get_local_pair_inds_i_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_local_pair_inds_i_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;
    auto get_local_pair_inds_k_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>&;

    auto get_mat_Q_full(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) const -> CDenseMatrix;
    auto get_mat_D_abs_full(const int64_t s_prim_count,
                            const int64_t p_prim_count,
                            const int64_t d_prim_count,
                            const std::vector<uint32_t>& s_prim_aoinds,
                            const std::vector<uint32_t>& p_prim_aoinds,
                            const std::vector<uint32_t>& d_prim_aoinds,
                            const int64_t naos,
                            const double* dens_ptr) const -> CDenseMatrix;

    auto form_Q_and_D_inds_for_K(const int64_t                s_prim_count,
                                 const int64_t                p_prim_count,
                                 const int64_t                d_prim_count,
                                 const std::vector<uint32_t>& s_prim_aoinds,
                                 const std::vector<uint32_t>& p_prim_aoinds,
                                 const std::vector<uint32_t>& d_prim_aoinds,
                                 const std::vector<double>&   s_prim_info,
                                 const std::vector<double>&   p_prim_info,
                                 const std::vector<double>&   d_prim_info) -> void;

    auto form_pair_inds_for_K(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count, const CDenseMatrix& Q_prime, const double Q_prime_thresh) -> void;
};

#endif /* ScreeningData_hpp */
