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

#include "GpuRuntime.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include "BoysFuncTable.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuDevices.hpp"
#include "GpuSafeChecks.hpp"
#include "GpuWrapper.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathFunc.hpp"
#include "OneElectronIntegrals.hpp"
#include "ScreeningData.hpp"
#include "StringFormat.hpp"

#define MATH_CONST_PI 3.14159265358979323846

CScreeningData::CScreeningData(const CMolecule& molecule,
                               const CMolecularBasis& basis,
                               const int64_t num_gpus_per_node,
                               const double pair_threshold,
                               const double density_threshold,
                               const int rank,
                               const int nnodes)
{
    _rank = rank;

    _nnodes = nnodes;

    _num_gpus_per_node = num_gpus_per_node;

    _pair_threshold = pair_threshold;

    _density_threshold = density_threshold;

    _computeQMatricesOnGPU(molecule, basis);

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    _s_prim_count = 0;
    _p_prim_count = 0;
    _d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0)
        {
            _s_prim_count += npgtos * ncgtos;
        }
        else if (gto_ang == 1)
        {
            _p_prim_count += npgtos * ncgtos;
        }
        else if (gto_ang == 2)
        {
            _d_prim_count += npgtos * ncgtos;
        }
    }

    // S gto block

    _s_prim_info.resize(5 * _s_prim_count);
    _s_prim_aoinds.resize(1 * _s_prim_count);

    gtoinfo::updatePrimitiveInfoForS(_s_prim_info.data(), _s_prim_aoinds.data(), _s_prim_count, gto_blocks);

    // P gto block

    _p_prim_info.resize(5 * _p_prim_count);
    _p_prim_aoinds.resize(3 * _p_prim_count);

    gtoinfo::updatePrimitiveInfoForP(_p_prim_info.data(), _p_prim_aoinds.data(), _p_prim_count, gto_blocks);

    // D gto block

    _d_prim_info.resize(5 * _d_prim_count);
    _d_prim_aoinds.resize(6 * _d_prim_count);

    gtoinfo::updatePrimitiveInfoForD(_d_prim_info.data(), _d_prim_aoinds.data(), _d_prim_count, gto_blocks);

    _sortQ();

    form_Q_and_D_inds_for_K();

    const auto cart_naos = _s_prim_count + _p_prim_count * 3 + _d_prim_count * 6;

    // allocate device pointers

    const auto boys_func_table_size = boysfunc::getFullBoysFuncTableSize();
    const auto boys_func_ft_size    = boysfunc::getBoysFuncFactorsSize();

    const auto s_prim_info_size = 5 * _s_prim_count;
    const auto p_prim_info_size = 5 * _p_prim_count;
    const auto d_prim_info_size = 5 * _d_prim_count;

    const auto s_prim_aoinds_size = 1 * _s_prim_count;
    const auto p_prim_aoinds_size = 3 * _p_prim_count;
    const auto d_prim_aoinds_size = 6 * _d_prim_count;

    // Note: these counts are larger than those used in J computation
    const auto ss_prim_pair_count = s_prim_aoinds_size * (s_prim_aoinds_size + 1) / 2;
    const auto sp_prim_pair_count = s_prim_aoinds_size * p_prim_aoinds_size;
    const auto sd_prim_pair_count = s_prim_aoinds_size * d_prim_aoinds_size;
    const auto pp_prim_pair_count = p_prim_aoinds_size * (p_prim_aoinds_size + 1) / 2;
    const auto pd_prim_pair_count = p_prim_aoinds_size * d_prim_aoinds_size;
    const auto dd_prim_pair_count = d_prim_aoinds_size * (d_prim_aoinds_size + 1) / 2;

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, sd_prim_pair_count,
                                               pp_prim_pair_count, pd_prim_pair_count, dd_prim_pair_count});

    const auto all_prim_pair_count = (ss_prim_pair_count + sp_prim_pair_count + sd_prim_pair_count +
                                      pp_prim_pair_count + pd_prim_pair_count + dd_prim_pair_count);

    // Note: these counts are larger than those used in K computation
    const auto pair_inds_count_for_K_ss = ss_prim_pair_count;
    const auto pair_inds_count_for_K_sp = sp_prim_pair_count;
    const auto pair_inds_count_for_K_sd = sd_prim_pair_count;
    const auto pair_inds_count_for_K_pp = pp_prim_pair_count;
    const auto pair_inds_count_for_K_pd = pd_prim_pair_count;
    const auto pair_inds_count_for_K_dd = dd_prim_pair_count;

    const auto max_pair_inds_count = std::max({pair_inds_count_for_K_ss, pair_inds_count_for_K_sp, pair_inds_count_for_K_pp,
                                               pair_inds_count_for_K_sd, pair_inds_count_for_K_pd, pair_inds_count_for_K_dd});

    const auto all_pair_inds_count = (pair_inds_count_for_K_ss + pair_inds_count_for_K_sp + pair_inds_count_for_K_pp +
                                      pair_inds_count_for_K_sd + pair_inds_count_for_K_pd + pair_inds_count_for_K_dd);

    const auto all_prim_aoinds_size = (s_prim_aoinds_size + p_prim_aoinds_size + d_prim_aoinds_size);

    _d_data_boys_func       = std::vector<double*>(num_gpus_per_node);
    _d_data_spd_prim_info   = std::vector<double*>(num_gpus_per_node);
    _d_data_spd_prim_aoinds = std::vector<uint32_t*>(num_gpus_per_node);

    _d_data_matrices_ABC = std::vector<double*>(num_gpus_per_node);

    _d_data_mat_D_J                 = std::vector<double*>(num_gpus_per_node);
    _d_data_mat_Q                   = std::vector<double*>(num_gpus_per_node);
    _d_data_first_second_inds       = std::vector<uint32_t*>(num_gpus_per_node);
    _d_data_pair_data               = std::vector<double*>(num_gpus_per_node);
    _d_data_mat_Q_local             = std::vector<double*>(num_gpus_per_node);
    _d_data_first_second_inds_local = std::vector<uint32_t*>(num_gpus_per_node);
    _d_data_pair_data_local         = std::vector<double*>(num_gpus_per_node);

    _d_data_mat_K                = std::vector<double*>(num_gpus_per_node);
    _d_data_pair_inds_for_K      = std::vector<uint32_t*>(num_gpus_per_node);
    _d_data_mat_D_full_AO        = std::vector<double*>(num_gpus_per_node);
    _d_data_Q_K                  = std::vector<double*>(num_gpus_per_node);
    _d_data_D_inds_K             = std::vector<uint32_t*>(num_gpus_per_node);
    _d_data_pair_counts_displs_K = std::vector<uint32_t*>(num_gpus_per_node);
    _d_data_pair_data_K          = std::vector<double*>(num_gpus_per_node);

    _devptr_double = std::vector<double*>(num_gpus_per_node);
    _devptr_uint32 = std::vector<uint32_t*>(num_gpus_per_node);

    for (int64_t gpu_id = 0; gpu_id < num_gpus_per_node; gpu_id++)
    {
        gpuSafe(gpuSetDevice(gpu_id));

        const auto ss_prim_pair_count_local = ss_prim_pair_count;
        const auto sp_prim_pair_count_local = sp_prim_pair_count;
        const auto sd_prim_pair_count_local = sd_prim_pair_count;
        const auto pp_prim_pair_count_local = pp_prim_pair_count;
        const auto pd_prim_pair_count_local = pd_prim_pair_count;
        const auto dd_prim_pair_count_local = dd_prim_pair_count;

        const auto max_prim_pair_count_local = std::max({ss_prim_pair_count_local, sp_prim_pair_count_local,
                                                         sd_prim_pair_count_local, pp_prim_pair_count_local,
                                                         pd_prim_pair_count_local, dd_prim_pair_count_local});

        const auto all_prim_pair_count_local = (ss_prim_pair_count_local + sp_prim_pair_count_local +
                                                sd_prim_pair_count_local + pp_prim_pair_count_local +
                                                pd_prim_pair_count_local + dd_prim_pair_count_local);

        // calculate general size

        const auto size_general_double = ((boys_func_table_size + boys_func_ft_size) +
                                          (s_prim_info_size + p_prim_info_size + d_prim_info_size));
        const auto size_general_uint32 = (s_prim_aoinds_size + p_prim_aoinds_size + d_prim_aoinds_size);

        // calculate preLinK size

        const auto size_prelink_double = 3 * cart_naos * cart_naos;

        // calculate J size

        // d_data_mat_D_J, (max_prim_pair_count + max_prim_pair_count_local) * sizeof(double)
        //
        // d_data_mat_Q, (ss_prim_pair_count + sp_prim_pair_count + sd_prim_pair_count +
        //                pp_prim_pair_count + pd_prim_pair_count + dd_prim_pair_count) * sizeof(double)
        //
        // d_data_first_second_inds, (ss_prim_pair_count + ss_prim_pair_count +
        //                            sp_prim_pair_count + sp_prim_pair_count +
        //                            sd_prim_pair_count + sd_prim_pair_count +
        //                            pp_prim_pair_count + pp_prim_pair_count +
        //                            pd_prim_pair_count + pd_prim_pair_count +
        //                            dd_prim_pair_count + dd_prim_pair_count) * sizeof(uint32_t)
        //
        // d_data_pair_data, (ss_pair_data.size() + sp_pair_data.size() + sd_pair_data.size() +
        //                    pp_pair_data.size() + pd_pair_data.size() + dd_pair_data.size()) * sizeof(double)
        //
        // d_data_mat_Q_local, (ss_prim_pair_count_local + sp_prim_pair_count_local + sd_prim_pair_count_local +
        //                      pp_prim_pair_count_local + pd_prim_pair_count_local + dd_prim_pair_count_local) * sizeof(double)
        //
        // d_data_first_second_inds_local, (ss_prim_pair_count_local + ss_prim_pair_count_local +
        //                                  sp_prim_pair_count_local + sp_prim_pair_count_local +
        //                                  sd_prim_pair_count_local + sd_prim_pair_count_local +
        //                                  pp_prim_pair_count_local + pp_prim_pair_count_local +
        //                                  pd_prim_pair_count_local + pd_prim_pair_count_local +
        //                                  dd_prim_pair_count_local + dd_prim_pair_count_local) * sizeof(uint32_t)
        //
        // d_data_pair_data_local, (ss_pair_data_local.size() + sp_pair_data_local.size() + sd_pair_data_local.size() +
        //                          pp_pair_data_local.size() + pd_pair_data_local.size() + dd_pair_data_local.size()) * sizeof(double)

        // size for double* pointers
        // _d_data_mat_D_J
        // _d_data_mat_Q
        // _d_data_pair_data
        // _d_data_mat_Q_local
        // _d_data_pair_data_local
        const auto size_J_double = ((max_prim_pair_count + max_prim_pair_count_local) +
                                    all_prim_pair_count +
                                    all_prim_pair_count +
                                    all_prim_pair_count_local +
                                    all_prim_pair_count_local);

        // size for uint32_t* pointers
        // _d_data_first_second_inds
        // _d_data_first_second_inds_local
        const auto size_J_uint32 = (all_prim_pair_count * 2 +
                                    all_prim_pair_count_local * 2);

        // calculate K size

        // d_mat_K, max_pair_inds_count * sizeof(double)
        //
        // d_data_pair_inds_for_K, (pair_inds_count_for_K_ss + pair_inds_count_for_K_ss +
        //                          pair_inds_count_for_K_sp + pair_inds_count_for_K_sp +
        //                          pair_inds_count_for_K_sd + pair_inds_count_for_K_sd +
        //                          pair_inds_count_for_K_pp + pair_inds_count_for_K_pp +
        //                          pair_inds_count_for_K_pd + pair_inds_count_for_K_pd +
        //                          pair_inds_count_for_K_dd + pair_inds_count_for_K_dd) * sizeof(uint32_t)
        //
        // d_mat_D_full_AO, cart_naos * cart_naos * sizeof(double)
        //
        // d_data_Q_K, (Q_K_ss.size() + Q_K_sp.size() + Q_K_sd.size() +
        //              Q_K_ps.size() + Q_K_pp.size() + Q_K_pd.size() +
        //              Q_K_ds.size() + Q_K_dp.size() + Q_K_dd.size()) * sizeof(double)
        //
        // d_data_D_inds_K, (D_inds_K_ss.size() + D_inds_K_sp.size() + D_inds_K_sd.size() +
        //                   D_inds_K_ps.size() + D_inds_K_pp.size() + D_inds_K_pd.size() +
        //                   D_inds_K_ds.size() + D_inds_K_dp.size() + D_inds_K_dd.size()) * sizeof(uint32_t)
        //
        // d_data_pair_counts_displs_K, (pair_counts_K_ss.size() + pair_counts_K_sp.size() + pair_counts_K_sd.size() +
        //                               pair_counts_K_ps.size() + pair_counts_K_pp.size() + pair_counts_K_pd.size() +
        //                               pair_counts_K_ds.size() + pair_counts_K_dp.size() + pair_counts_K_dd.size() +
        //                               pair_displs_K_ss.size() + pair_displs_K_sp.size() + pair_displs_K_sd.size() +
        //                               pair_displs_K_ps.size() + pair_displs_K_pp.size() + pair_displs_K_pd.size() +
        //                               pair_displs_K_ds.size() + pair_displs_K_dp.size() + pair_displs_K_dd.size()) * sizeof(uint32_t)
        //
        // d_data_pair_data_K, (pair_data_K_ss.size() + pair_data_K_sp.size() + pair_data_K_sd.size() +
        //                      pair_data_K_ps.size() + pair_data_K_pp.size() + pair_data_K_pd.size() +
        //                      pair_data_K_ds.size() + pair_data_K_dp.size() + pair_data_K_dd.size()) * sizeof(double)

        // size for double* pointers
        // _d_data_mat_K
        // _d_data_mat_D_full_AO
        // _d_data_Q_K
        // _d_data_pair_data_K
        const auto size_K_double = (max_pair_inds_count +
                                    cart_naos * cart_naos +
                                    cart_naos * cart_naos +
                                    cart_naos * cart_naos);

        // size for uint32_t* pointers
        // _d_data_pair_inds_for_K
        // _d_data_D_inds_K
        // _d_data_pair_counts_displs_K
        const auto size_K_uint32 = (all_pair_inds_count * 2 +
                                    cart_naos * cart_naos +
                                    all_prim_aoinds_size * 3 * 2);

        // allocate data

        const auto size_JK_double = std::max({size_prelink_double, size_J_double, size_K_double});
        const auto size_JK_uint32 = std::max({size_J_uint32, size_K_uint32});

        gpuSafe(gpuMalloc(&_devptr_double[gpu_id], (size_general_double + size_JK_double) * sizeof(double)));
        gpuSafe(gpuMalloc(&_devptr_uint32[gpu_id], (size_general_uint32 + size_JK_uint32) * sizeof(uint32_t)));

        // prepare genral data (Boys function and primitive AOs)

        // _d_data_boys_func, (boys_func_table_size + boys_func_ft_size) * sizeof(double)
        // _d_data_spd_prim_info, (s_prim_info_size + p_prim_info_size + d_prim_info_size) * sizeof(double)
        _d_data_boys_func[gpu_id] = _devptr_double[gpu_id];
        _d_data_spd_prim_info[gpu_id] = _d_data_boys_func[gpu_id] + (boys_func_table_size + boys_func_ft_size);

        // _d_data_spd_prim_aoinds, (s_prim_aoinds_size + p_prim_aoinds_size + d_prim_aoinds_size) * sizeof(uint32_t)
        _d_data_spd_prim_aoinds[gpu_id] = _devptr_uint32[gpu_id];

        // prepare preLinK pointers (the same memory is also used by J and K)

        _d_data_matrices_ABC[gpu_id] = _d_data_spd_prim_info[gpu_id] + (s_prim_info_size + p_prim_info_size + d_prim_info_size);

        // prepare J pointers (the same memory is also used by preLinK and K)

        // double* pointers for J
        _d_data_mat_D_J[gpu_id]         = _d_data_spd_prim_info[gpu_id] + (s_prim_info_size + p_prim_info_size + d_prim_info_size);
        _d_data_mat_Q[gpu_id]           = _d_data_mat_D_J[gpu_id] + (max_prim_pair_count + max_prim_pair_count_local);
        _d_data_pair_data[gpu_id]       = _d_data_mat_Q[gpu_id] + all_prim_pair_count;
        _d_data_mat_Q_local[gpu_id]     = _d_data_pair_data[gpu_id] + all_prim_pair_count;
        _d_data_pair_data_local[gpu_id] = _d_data_mat_Q_local[gpu_id] + all_prim_pair_count_local;

        // uint32_t* pointers for J
        _d_data_first_second_inds[gpu_id]       = _d_data_spd_prim_aoinds[gpu_id] + (s_prim_aoinds_size + p_prim_aoinds_size + d_prim_aoinds_size);
        _d_data_first_second_inds_local[gpu_id] = _d_data_first_second_inds[gpu_id] + all_prim_pair_count * 2;

        // prepare K pointers (the same memory is also used by preLinK and J)

        // double* pointers for K
        _d_data_mat_K[gpu_id]         = _d_data_spd_prim_info[gpu_id] + (s_prim_info_size + p_prim_info_size + d_prim_info_size);
        _d_data_mat_D_full_AO[gpu_id] = _d_data_mat_K[gpu_id] + max_pair_inds_count;
        _d_data_Q_K[gpu_id]           = _d_data_mat_D_full_AO[gpu_id] + cart_naos * cart_naos;
        _d_data_pair_data_K[gpu_id]   = _d_data_Q_K[gpu_id] + cart_naos * cart_naos;

        // uint32_t* pointers for K
        _d_data_pair_inds_for_K[gpu_id]      = _d_data_spd_prim_aoinds[gpu_id] + (s_prim_aoinds_size + p_prim_aoinds_size + d_prim_aoinds_size);
        _d_data_D_inds_K[gpu_id]             = _d_data_pair_inds_for_K[gpu_id] + all_pair_inds_count * 2;
        _d_data_pair_counts_displs_K[gpu_id] = _d_data_D_inds_K[gpu_id] + cart_naos * cart_naos;
    }

    gpuSafe(gpuSetDevice(0));
}

CScreeningData::~CScreeningData()
{
    errors::assertMsgCritical(
        !omp_in_parallel(),
        std::string("CScreeningData destructor: should never be called in omp parallel reigion"));

    // free device pointers
    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        gpuSafe(gpuSetDevice(gpu_id));

        gpuSafe(gpuFree(_devptr_double[gpu_id]));
        gpuSafe(gpuFree(_devptr_uint32[gpu_id]));
    }

    gpuSafe(gpuSetDevice(0));
}

auto
CScreeningData::get_devptr_data_boys_func(const int64_t gpu_id) const -> double*
{
    return _d_data_boys_func[gpu_id];
}

auto
CScreeningData::get_devptr_data_spd_prim_info(const int64_t gpu_id) const -> double*
{
    return _d_data_spd_prim_info[gpu_id];
}

auto
CScreeningData::get_devptr_data_spd_prim_aoinds(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_spd_prim_aoinds[gpu_id];
}

auto
CScreeningData::get_devptr_data_matrices_ABC(const int64_t gpu_id) const -> double*
{
    return _d_data_matrices_ABC[gpu_id];
}

auto
CScreeningData::get_devptr_data_mat_D_J(const int64_t gpu_id) const -> double*
{
    return _d_data_mat_D_J[gpu_id];
}

auto
CScreeningData::get_devptr_data_mat_Q(const int64_t gpu_id) const -> double*
{
    return _d_data_mat_Q[gpu_id];
}

auto
CScreeningData::get_devptr_data_first_second_inds(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_first_second_inds[gpu_id];
}

auto
CScreeningData::get_devptr_data_pair_data(const int64_t gpu_id) const -> double*
{
    return _d_data_pair_data[gpu_id];
}

auto
CScreeningData::get_devptr_data_mat_Q_local(const int64_t gpu_id) const -> double*
{
    return _d_data_mat_Q_local[gpu_id];
}

auto
CScreeningData::get_devptr_data_first_second_inds_local(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_first_second_inds_local[gpu_id];
}

auto
CScreeningData::get_devptr_data_pair_data_local(const int64_t gpu_id) const -> double*
{
    return _d_data_pair_data_local[gpu_id];
}

auto
CScreeningData::get_devptr_data_mat_K(const int64_t gpu_id) const -> double*
{
    return _d_data_mat_K[gpu_id];
}

auto
CScreeningData::get_devptr_data_pair_inds_for_K(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_pair_inds_for_K[gpu_id];
}

auto
CScreeningData::get_devptr_data_mat_D_full_AO(const int64_t gpu_id) const -> double*
{
    return _d_data_mat_D_full_AO[gpu_id];
}

auto
CScreeningData::get_devptr_data_Q_K(const int64_t gpu_id) const -> double*
{
    return _d_data_Q_K[gpu_id];
}

auto
CScreeningData::get_devptr_data_D_inds_K(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_D_inds_K[gpu_id];
}

auto
CScreeningData::get_devptr_data_pair_counts_displs_K(const int64_t gpu_id) const -> uint32_t*
{
    return _d_data_pair_counts_displs_K[gpu_id];
}

auto
CScreeningData::get_devptr_data_pair_data_K(const int64_t gpu_id) const -> double*
{
    return _d_data_pair_data_K[gpu_id];
}

auto
CScreeningData::_computeQMatricesOnGPU(const CMolecule& molecule, const CMolecularBasis& basis) -> void
{
    gpuSafe(gpuSetDevice(0));

    // GTOs blocks and number of AOs

    const auto gto_blocks = gtofunc::makeGtoBlocks(basis, molecule);

    // gto blocks

    int64_t s_prim_count = 0;
    int64_t p_prim_count = 0;
    int64_t d_prim_count = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();
        const auto npgtos = gto_block.getNumberOfPrimitives();

        const auto gto_ang = gto_block.getAngularMomentum();

        if (gto_ang == 0) s_prim_count += npgtos * ncgtos;
        if (gto_ang == 1) p_prim_count += npgtos * ncgtos;
        if (gto_ang == 2) d_prim_count += npgtos * ncgtos;
    }

    // Cartesian to spherical index mapping for P and D

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;
    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    for (const auto& gto_block : gto_blocks)
    {
        const auto gto_ang = gto_block.getAngularMomentum();

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
    }

    // S, P, D gto block

    std::vector<double>   s_prim_info(5 * s_prim_count);
    std::vector<double>   p_prim_info(5 * p_prim_count);
    std::vector<double>   d_prim_info(5 * d_prim_count);

    std::vector<uint32_t> s_prim_aoinds(1 * s_prim_count);
    std::vector<uint32_t> p_prim_aoinds(3 * p_prim_count);
    std::vector<uint32_t> d_prim_aoinds(6 * d_prim_count);

    gtoinfo::updatePrimitiveInfoForS(s_prim_info.data(), s_prim_aoinds.data(), s_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForP(p_prim_info.data(), p_prim_aoinds.data(), p_prim_count, gto_blocks);
    gtoinfo::updatePrimitiveInfoForD(d_prim_info.data(), d_prim_aoinds.data(), d_prim_count, gto_blocks);

    // GTO block pairs

    std::vector<uint32_t> ss_first_inds;
    std::vector<uint32_t> sp_first_inds;
    std::vector<uint32_t> sd_first_inds;
    std::vector<uint32_t> pp_first_inds;
    std::vector<uint32_t> pd_first_inds;
    std::vector<uint32_t> dd_first_inds;

    std::vector<uint32_t> ss_second_inds;
    std::vector<uint32_t> sp_second_inds;
    std::vector<uint32_t> sd_second_inds;
    std::vector<uint32_t> pp_second_inds;
    std::vector<uint32_t> pd_second_inds;
    std::vector<uint32_t> dd_second_inds;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++)
        {
            ss_first_inds.push_back(i);
            ss_second_inds.push_back(j);
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                sp_first_inds.push_back(i);
                sp_second_inds.push_back(j * 3 + s);
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t s = 0; s < 6; s++)
            {
                sd_first_inds.push_back(i);
                sd_second_inds.push_back(j * 6 + s);
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    pp_first_inds.push_back(i * 3 + i_cart);
                    pp_second_inds.push_back(j * 3 + j_cart);
                }
            }
        }

        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    pd_first_inds.push_back(i * 3 + i_cart);
                    pd_second_inds.push_back(j * 6 + j_cart);
                }
            }
        }
    }

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    dd_first_inds.push_back(i * 6 + i_cart);
                    dd_second_inds.push_back(j * 6 + j_cart);
                }
            }
        }
    }

    const auto ss_prim_pair_count = static_cast<int64_t>(ss_first_inds.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sp_first_inds.size());
    const auto sd_prim_pair_count = static_cast<int64_t>(sd_first_inds.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(pp_first_inds.size());
    const auto pd_prim_pair_count = static_cast<int64_t>(pd_first_inds.size());
    const auto dd_prim_pair_count = static_cast<int64_t>(dd_first_inds.size());

    const auto max_prim_pair_count = std::max({ss_prim_pair_count, sp_prim_pair_count, sd_prim_pair_count,
                                               pp_prim_pair_count, pd_prim_pair_count, dd_prim_pair_count});

    // Q on host

    std::vector<double> h_mat_Q(max_prim_pair_count);

    _Q_matrix_ss = CDenseMatrix(s_prim_count, s_prim_count);
    _Q_matrix_sp = CDenseMatrix(s_prim_count, p_prim_count * 3);
    _Q_matrix_sd = CDenseMatrix(s_prim_count, d_prim_count * 6);
    _Q_matrix_pp = CDenseMatrix(p_prim_count * 3, p_prim_count * 3);
    _Q_matrix_pd = CDenseMatrix(p_prim_count * 3, d_prim_count * 6);
    _Q_matrix_dd = CDenseMatrix(d_prim_count * 6, d_prim_count * 6);

    _Q_matrix_ss.zero();
    _Q_matrix_sp.zero();
    _Q_matrix_sd.zero();
    _Q_matrix_pp.zero();
    _Q_matrix_pd.zero();
    _Q_matrix_dd.zero();

    // Boys function (tabulated for order 0-28)

    const auto boys_func_table = boysfunc::getFullBoysFuncTable();
    const auto boys_func_ft = boysfunc::getBoysFuncFactors();

    // size of memory allocation

    const auto count_double = (max_prim_pair_count +
                               static_cast<int64_t>(boys_func_table.size()) +
                               static_cast<int64_t>(boys_func_ft.size()) +
                               static_cast<int64_t>(s_prim_info.size()) +
                               static_cast<int64_t>(p_prim_info.size()) +
                               static_cast<int64_t>(d_prim_info.size()));

    const auto count_uint32 = (ss_prim_pair_count + sp_prim_pair_count + sd_prim_pair_count +
                               pp_prim_pair_count + pd_prim_pair_count + dd_prim_pair_count) * 2;

    // memory allocation

    double* d_data_double;
    uint32_t* d_data_uint32;

    gpuSafe(gpuMalloc(&d_data_double, count_double * sizeof(double)));
    gpuSafe(gpuMalloc(&d_data_uint32, count_uint32 * sizeof(uint32_t)));

    double* d_mat_Q           = d_data_double;
    double* d_boys_func_table = d_mat_Q + max_prim_pair_count;
    double* d_boys_func_ft    = d_boys_func_table + static_cast<int64_t>(boys_func_table.size());
    double* d_s_prim_info     = d_boys_func_ft + static_cast<int64_t>(boys_func_ft.size());
    double* d_p_prim_info     = d_s_prim_info + static_cast<int64_t>(s_prim_info.size());
    double* d_d_prim_info     = d_p_prim_info + static_cast<int64_t>(p_prim_info.size());

    uint32_t* d_ss_first_inds  = d_data_uint32;
    uint32_t* d_sp_first_inds  = d_ss_first_inds + ss_prim_pair_count;
    uint32_t* d_sd_first_inds  = d_sp_first_inds + sp_prim_pair_count;
    uint32_t* d_pp_first_inds  = d_sd_first_inds + sd_prim_pair_count;
    uint32_t* d_pd_first_inds  = d_pp_first_inds + pp_prim_pair_count;
    uint32_t* d_dd_first_inds  = d_pd_first_inds + pd_prim_pair_count;
    uint32_t* d_ss_second_inds = d_dd_first_inds + dd_prim_pair_count; 
    uint32_t* d_sp_second_inds = d_ss_second_inds + ss_prim_pair_count;
    uint32_t* d_sd_second_inds = d_sp_second_inds + sp_prim_pair_count;
    uint32_t* d_pp_second_inds = d_sd_second_inds + sd_prim_pair_count;
    uint32_t* d_pd_second_inds = d_pp_second_inds + pp_prim_pair_count;
    uint32_t* d_dd_second_inds = d_pd_second_inds + pd_prim_pair_count;

    gpuSafe(gpuMemcpy(d_boys_func_table, boys_func_table.data(), boys_func_table.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_boys_func_ft, boys_func_ft.data(), boys_func_ft.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_s_prim_info, s_prim_info.data(), s_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_p_prim_info, p_prim_info.data(), p_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_d_prim_info, d_prim_info.data(), d_prim_info.size() * sizeof(double), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_ss_first_inds, ss_first_inds.data(), ss_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_ss_second_inds, ss_second_inds.data(), ss_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sp_first_inds, sp_first_inds.data(), sp_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sp_second_inds, sp_second_inds.data(), sp_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_sd_first_inds, sd_first_inds.data(), sd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_sd_second_inds, sd_second_inds.data(), sd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pp_first_inds, pp_first_inds.data(), pp_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pp_second_inds, pp_second_inds.data(), pp_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_pd_first_inds, pd_first_inds.data(), pd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_pd_second_inds, pd_second_inds.data(), pd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuMemcpy(d_dd_first_inds, dd_first_inds.data(), dd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));
    gpuSafe(gpuMemcpy(d_dd_second_inds, dd_second_inds.data(), dd_prim_pair_count * sizeof(uint32_t), gpuMemcpyHostToDevice));

    gpuSafe(gpuDeviceSynchronize());

    // compute Q

    // Q: SS

    if (ss_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((ss_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSS<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_ss_first_inds,
                           d_ss_second_inds,
                           static_cast<uint32_t>(ss_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, ss_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
        {
            const auto i = ss_first_inds[ij];
            const auto j = ss_second_inds[ij];

            _Q_matrix_ss.row(i)[j] = h_mat_Q[ij];

            if (i != j) _Q_matrix_ss.row(j)[i] = h_mat_Q[ij];
        }
    }

    // Q: SP

    if (sp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSP<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_sp_first_inds,
                           d_sp_second_inds,
                           static_cast<uint32_t>(sp_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, sp_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
        {
            const auto i = sp_first_inds[ij];
            const auto j = sp_second_inds[ij];

            _Q_matrix_sp.row(i)[j] = h_mat_Q[ij];
        }
    }

    // Q: SD

    if (sd_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((sd_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixSD<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_s_prim_info,
                           static_cast<uint32_t>(s_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_sd_first_inds,
                           d_sd_second_inds,
                           static_cast<uint32_t>(sd_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, sd_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < sd_prim_pair_count; ij++)
        {
            const auto i = sd_first_inds[ij];
            const auto j = sd_second_inds[ij];

            _Q_matrix_sd.row(i)[j] = h_mat_Q[ij];
        }
    }

    // Q: PP

    if (pp_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pp_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixPP<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_pp_first_inds,
                           d_pp_second_inds,
                           static_cast<uint32_t>(pp_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, pp_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
        {
            const auto i = pp_first_inds[ij];
            const auto j = pp_second_inds[ij];

            _Q_matrix_pp.row(i)[j] = h_mat_Q[ij];

            if (i != j) _Q_matrix_pp.row(j)[i] = h_mat_Q[ij];
        }
    }

    // Q: PD

    if (pd_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((pd_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixPD<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_p_prim_info,
                           static_cast<uint32_t>(p_prim_count),
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_pd_first_inds,
                           d_pd_second_inds,
                           static_cast<uint32_t>(pd_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, pd_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < pd_prim_pair_count; ij++)
        {
            const auto i = pd_first_inds[ij];
            const auto j = pd_second_inds[ij];

            _Q_matrix_pd.row(i)[j] = h_mat_Q[ij];
        }
    }

    // Q: DD

    if (dd_prim_pair_count > 0)
    {
        dim3 threads_per_block(TILE_DIM);

        dim3 num_blocks((dd_prim_pair_count + threads_per_block.x - 1) / threads_per_block.x);

        gpu::computeQMatrixDD<<<num_blocks, threads_per_block>>>(
                           d_mat_Q,
                           d_d_prim_info,
                           static_cast<uint32_t>(d_prim_count),
                           d_dd_first_inds,
                           d_dd_second_inds,
                           static_cast<uint32_t>(dd_prim_pair_count),
                           d_boys_func_table,
                           d_boys_func_ft);

        gpuSafe(gpuMemcpy(h_mat_Q.data(), d_mat_Q, dd_prim_pair_count * sizeof(double), gpuMemcpyDeviceToHost));

        for (int64_t ij = 0; ij < dd_prim_pair_count; ij++)
        {
            const auto i = dd_first_inds[ij];
            const auto j = dd_second_inds[ij];

            _Q_matrix_dd.row(i)[j] = h_mat_Q[ij];

            if (i != j) _Q_matrix_dd.row(j)[i] = h_mat_Q[ij];
        }
    }

    gpuSafe(gpuDeviceSynchronize());

    gpuSafe(gpuFree(d_data_double));
    gpuSafe(gpuFree(d_data_uint32));
}

auto
CScreeningData::getNumGpusPerNode() const -> const int64_t
{
    return _num_gpus_per_node;
}

auto
CScreeningData::setTimerSummary(const std::string& timer_summary) -> void
{
    _timer_summary = timer_summary;
}

auto
CScreeningData::getTimerSummary() const -> const std::string
{
    return _timer_summary;
}

auto
CScreeningData::initGpuTimers(const int64_t num_gpus_per_node) -> void
{
    _gpu_timer_summary = std::vector<std::string>(num_gpus_per_node);
}

auto
CScreeningData::setGpuTimerSummary(const int64_t gpu_id, const std::string& timer_summary) -> void
{
    _gpu_timer_summary[gpu_id] = timer_summary;
}

auto
CScreeningData::getGpuTimerSummary() const -> const std::vector<std::string>
{
    return _gpu_timer_summary;
}

auto
CScreeningData::_sortQ() -> void
{
    const auto s_prim_count = _s_prim_count;
    const auto p_prim_count = _p_prim_count;
    const auto d_prim_count = _d_prim_count;

    const auto& s_prim_info = _s_prim_info;
    const auto& p_prim_info = _p_prim_info;
    const auto& d_prim_info = _d_prim_info;

    std::vector<std::tuple<double, int64_t, int64_t>> sorted_ss_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_sp_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_sd_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_pp_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_pd_mat_Q;
    std::vector<std::tuple<double, int64_t, int64_t>> sorted_dd_mat_Q;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            if (Q_ij > _pair_threshold) sorted_ss_mat_Q.push_back(std::make_tuple(Q_ij, i, j));
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t s = 0; s < 3; s++)
            {
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + s];
                if (Q_ij > _pair_threshold) sorted_sp_mat_Q.push_back(std::make_tuple(Q_ij, i, j * 3 + s));
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t s = 0; s < 6; s++)
            {
                const auto Q_ij = _Q_matrix_sd.row(i)[j * 6 + s];
                if (Q_ij > _pair_threshold) sorted_sd_mat_Q.push_back(std::make_tuple(Q_ij, i, j * 6 + s));
            }
        }
    }

    // P-P gto block pair and P-D gto block pair

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    if (Q_ij > _pair_threshold) sorted_pp_mat_Q.push_back(std::make_tuple(Q_ij, i * 3 + i_cart, j * 3 + j_cart));
                }
            }
        }

        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_pd.row(i * 3 + i_cart)[j * 6 + j_cart];
                    if (Q_ij > _pair_threshold) sorted_pd_mat_Q.push_back(std::make_tuple(Q_ij, i * 3 + i_cart, j * 6 + j_cart));
                }
            }
        }
    }

    // D-D gto block pair

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto Q_ij = _Q_matrix_dd.row(i * 6 + i_cart)[j * 6 + j_cart];
                    if (Q_ij > _pair_threshold) sorted_dd_mat_Q.push_back(std::make_tuple(Q_ij, i * 6 + i_cart, j * 6 + j_cart));
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q.begin(), sorted_ss_mat_Q.end());
    std::sort(sorted_sp_mat_Q.begin(), sorted_sp_mat_Q.end());
    std::sort(sorted_sd_mat_Q.begin(), sorted_sd_mat_Q.end());
    std::sort(sorted_pp_mat_Q.begin(), sorted_pp_mat_Q.end());
    std::sort(sorted_pd_mat_Q.begin(), sorted_pd_mat_Q.end());
    std::sort(sorted_dd_mat_Q.begin(), sorted_dd_mat_Q.end());

    const auto ss_prim_pair_count = static_cast<int64_t>(sorted_ss_mat_Q.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sorted_sp_mat_Q.size());
    const auto sd_prim_pair_count = static_cast<int64_t>(sorted_sd_mat_Q.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(sorted_pp_mat_Q.size());
    const auto pd_prim_pair_count = static_cast<int64_t>(sorted_pd_mat_Q.size());
    const auto dd_prim_pair_count = static_cast<int64_t>(sorted_dd_mat_Q.size());

    // std::stringstream ss;
    // ss << "Pair screening\n";
    // ss << "  SS pair: " << static_cast<double>(ss_prim_pair_count) / (s_prim_count * (s_prim_count + 1) / 2)<< "\n";
    // ss << "  SP pair: " << static_cast<double>(sp_prim_pair_count) / (s_prim_count * p_prim_count * 3)<< "\n";
    // ss << "  SD pair: " << static_cast<double>(sd_prim_pair_count) / (s_prim_count * d_prim_count * 6)<< "\n";
    // ss << "  PP pair: " << static_cast<double>(pp_prim_pair_count) / (p_prim_count * 3 * (p_prim_count * 3 + 1) / 2)<< "\n";
    // ss << "  PD pair: " << static_cast<double>(pd_prim_pair_count) / (p_prim_count * 3 * d_prim_count * 6)<< "\n";
    // ss << "  DD pair: " << static_cast<double>(dd_prim_pair_count) / (d_prim_count * 6 * (d_prim_count * 6 + 1) / 2)<< "\n";
    // std::cout << ss.str() << "\n";

    // form local vectors

    _ss_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _ss_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _ss_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _sp_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sp_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sp_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _sd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _sd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _pp_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pp_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pp_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _pd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _pd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _dd_first_inds_local  = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _dd_second_inds_local = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _dd_mat_Q_local       = std::vector<std::vector<double>>(_num_gpus_per_node);

    _ss_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);
    _sp_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);
    _sd_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);
    _pp_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);
    _pd_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);
    _dd_pair_data_local = std::vector<std::vector<double>>(_num_gpus_per_node);

    auto rank = _rank;
    auto nnodes = _nnodes;

    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
        auto gpu_count = nnodes * _num_gpus_per_node;

        auto ss_batch_size = mathfunc::batch_size(ss_prim_pair_count, gpu_rank, gpu_count);

        _ss_first_inds_local[gpu_id]  = std::vector<uint32_t>(ss_batch_size);
        _ss_second_inds_local[gpu_id] = std::vector<uint32_t>(ss_batch_size);
        _ss_mat_Q_local[gpu_id]       = std::vector<double>(ss_batch_size);

        _ss_pair_data_local[gpu_id]   = std::vector<double>(ss_batch_size);

        auto sp_batch_size = mathfunc::batch_size(sp_prim_pair_count, gpu_rank, gpu_count);

        _sp_first_inds_local[gpu_id]  = std::vector<uint32_t>(sp_batch_size);
        _sp_second_inds_local[gpu_id] = std::vector<uint32_t>(sp_batch_size);
        _sp_mat_Q_local[gpu_id]       = std::vector<double>(sp_batch_size);

        _sp_pair_data_local[gpu_id]   = std::vector<double>(sp_batch_size);

        auto sd_batch_size = mathfunc::batch_size(sd_prim_pair_count, gpu_rank, gpu_count);

        _sd_first_inds_local[gpu_id]  = std::vector<uint32_t>(sd_batch_size);
        _sd_second_inds_local[gpu_id] = std::vector<uint32_t>(sd_batch_size);
        _sd_mat_Q_local[gpu_id]       = std::vector<double>(sd_batch_size);

        _sd_pair_data_local[gpu_id]   = std::vector<double>(sd_batch_size);

        auto pp_batch_size = mathfunc::batch_size(pp_prim_pair_count, gpu_rank, gpu_count);

        _pp_first_inds_local[gpu_id]  = std::vector<uint32_t>(pp_batch_size);
        _pp_second_inds_local[gpu_id] = std::vector<uint32_t>(pp_batch_size);
        _pp_mat_Q_local[gpu_id]       = std::vector<double>(pp_batch_size);

        _pp_pair_data_local[gpu_id]   = std::vector<double>(pp_batch_size);

        auto pd_batch_size = mathfunc::batch_size(pd_prim_pair_count, gpu_rank, gpu_count);

        _pd_first_inds_local[gpu_id]  = std::vector<uint32_t>(pd_batch_size);
        _pd_second_inds_local[gpu_id] = std::vector<uint32_t>(pd_batch_size);
        _pd_mat_Q_local[gpu_id]       = std::vector<double>(pd_batch_size);

        _pd_pair_data_local[gpu_id]   = std::vector<double>(pd_batch_size);

        auto dd_batch_size = mathfunc::batch_size(dd_prim_pair_count, gpu_rank, gpu_count);

        _dd_first_inds_local[gpu_id]  = std::vector<uint32_t>(dd_batch_size);
        _dd_second_inds_local[gpu_id] = std::vector<uint32_t>(dd_batch_size);
        _dd_mat_Q_local[gpu_id]       = std::vector<double>(dd_batch_size);

        _dd_pair_data_local[gpu_id]   = std::vector<double>(dd_batch_size);

        {
            auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
            auto gpu_count = nnodes * _num_gpus_per_node;

            for (int64_t ij = gpu_rank, idx = 0; ij < ss_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_ss_mat_Q[ss_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _ss_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _ss_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _ss_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _ss_pair_data_local[gpu_id][idx] = S_ij_00;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < sp_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_sp_mat_Q[sp_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _sp_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _sp_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _sp_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _sp_pair_data_local[gpu_id][idx] = S_ij_00;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < sd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_sd_mat_Q[sd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _sd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _sd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _sd_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _sd_pair_data_local[gpu_id][idx] = S_ij_00;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < pp_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_pp_mat_Q[pp_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _pp_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _pp_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _pp_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _pp_pair_data_local[gpu_id][idx] = S_ij_00;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < pd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_pd_mat_Q[pd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _pd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _pd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _pd_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _pd_pair_data_local[gpu_id][idx] = S_ij_00;
            }

            for (int64_t ij = gpu_rank, idx = 0; ij < dd_prim_pair_count; ij+=gpu_count, idx++)
            {
                const auto& vals = sorted_dd_mat_Q[dd_prim_pair_count - 1 - ij];

                auto Q_ij = std::get<0>(vals);
                auto i    = std::get<1>(vals);
                auto j    = std::get<2>(vals);

                _dd_first_inds_local[gpu_id][idx]  = static_cast<uint32_t>(i);
                _dd_second_inds_local[gpu_id][idx] = static_cast<uint32_t>(j);
                _dd_mat_Q_local[gpu_id][idx]       = Q_ij;

                // ij pair data:

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

                const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

                const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

                _dd_pair_data_local[gpu_id][idx] = S_ij_00;
            }
        }
    }
}

auto
CScreeningData::sortQD(const int64_t s_prim_count,
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
                       const double eri_threshold) -> void
{
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_ss_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_sp_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_sd_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_pp_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_pd_mat_Q_D;
    std::vector<std::tuple<double, int64_t, int64_t, double, double>> sorted_dd_mat_Q_D;

    auto max_Q_ss_ptr = std::max_element(_Q_matrix_ss.values(), _Q_matrix_ss.values() + _Q_matrix_ss.getNumberOfElements());
    double max_Q = *max_Q_ss_ptr;

    if (p_prim_count > 0)
    {
        auto max_Q_sp_ptr = std::max_element(_Q_matrix_sp.values(), _Q_matrix_sp.values() + _Q_matrix_sp.getNumberOfElements());
        auto max_Q_pp_ptr = std::max_element(_Q_matrix_pp.values(), _Q_matrix_pp.values() + _Q_matrix_pp.getNumberOfElements());
        max_Q = std::max({max_Q, *max_Q_sp_ptr, *max_Q_pp_ptr});
    }

    if (d_prim_count > 0)
    {
        auto max_Q_sd_ptr = std::max_element(_Q_matrix_sd.values(), _Q_matrix_sd.values() + _Q_matrix_sd.getNumberOfElements());
        auto max_Q_dd_ptr = std::max_element(_Q_matrix_dd.values(), _Q_matrix_dd.values() + _Q_matrix_dd.getNumberOfElements());
        max_Q = std::max({max_Q, *max_Q_sd_ptr, *max_Q_dd_ptr});
    }

    if ((p_prim_count > 0) && (d_prim_count > 0))
    {
        auto max_Q_pd_ptr = std::max_element(_Q_matrix_pd.values(), _Q_matrix_pd.values() + _Q_matrix_pd.getNumberOfElements());
        max_Q = std::max({max_Q, *max_Q_pd_ptr});
    }

    // Coulomb: Q_bra*Q_ket*D_ket > ERI_threshold
    const double QD_threshold = eri_threshold / max_Q;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
            const auto Q_ij = _Q_matrix_ss.row(i)[j];
            const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

            if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
            {
                sorted_ss_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j, Q_ij, D_ij));
            }
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 3; j_cart++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                const auto Q_ij = _Q_matrix_sp.row(i)[j * 3 + j_cart];
                const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
                {
                    sorted_sp_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j * 3 + j_cart, Q_ij, D_ij));
                }
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 6; j_cart++)
            {
                const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                const auto Q_ij = _Q_matrix_sd.row(i)[j * 6 + j_cart];
                const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
                {
                    sorted_sd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i, j * 6 + j_cart, Q_ij, D_ij));
                }
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_pp.row(i * 3 + i_cart)[j * 3 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
                    {
                        sorted_pp_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 3 + i_cart, j * 3 + j_cart, Q_ij, D_ij));
                    }
                }
            }
        }
    
        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_pd.row(i * 3 + i_cart)[j * 6 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
                    {
                        sorted_pd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 3 + i_cart, j * 6 + j_cart, Q_ij, D_ij));
                    }
                }
            }
        }    
    }

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                const auto i_cgto = d_prim_aoinds[i + d_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];
                    const auto Q_ij = _Q_matrix_dd.row(i * 6 + i_cart)[j * 6 + j_cart];
                    const auto QD_abs_ij = Q_ij * std::fabs(D_ij);

                    if ((QD_abs_ij > QD_threshold) && (Q_ij > _pair_threshold) && (std::fabs(D_ij) > _density_threshold))
                    {
                        sorted_dd_mat_Q_D.push_back(std::make_tuple(QD_abs_ij, i * 6 + i_cart, j * 6 + j_cart, Q_ij, D_ij));
                    }
                }
            }
        }
    }

    std::sort(sorted_ss_mat_Q_D.begin(), sorted_ss_mat_Q_D.end());
    std::sort(sorted_sp_mat_Q_D.begin(), sorted_sp_mat_Q_D.end());
    std::sort(sorted_sd_mat_Q_D.begin(), sorted_sd_mat_Q_D.end());
    std::sort(sorted_pp_mat_Q_D.begin(), sorted_pp_mat_Q_D.end());
    std::sort(sorted_pd_mat_Q_D.begin(), sorted_pd_mat_Q_D.end());
    std::sort(sorted_dd_mat_Q_D.begin(), sorted_dd_mat_Q_D.end());

    const auto ss_prim_pair_count = static_cast<int64_t>(sorted_ss_mat_Q_D.size());
    const auto sp_prim_pair_count = static_cast<int64_t>(sorted_sp_mat_Q_D.size());
    const auto sd_prim_pair_count = static_cast<int64_t>(sorted_sd_mat_Q_D.size());
    const auto pp_prim_pair_count = static_cast<int64_t>(sorted_pp_mat_Q_D.size());
    const auto pd_prim_pair_count = static_cast<int64_t>(sorted_pd_mat_Q_D.size());
    const auto dd_prim_pair_count = static_cast<int64_t>(sorted_dd_mat_Q_D.size());

    _ss_mat_Q       = std::vector<double>  (ss_prim_pair_count);
    _ss_mat_D       = std::vector<double>  (ss_prim_pair_count);
    _ss_first_inds  = std::vector<uint32_t>(ss_prim_pair_count);
    _ss_second_inds = std::vector<uint32_t>(ss_prim_pair_count);

    _ss_pair_data = std::vector<double>(ss_prim_pair_count);

    for (int64_t ij = 0; ij < ss_prim_pair_count; ij++)
    {
        const auto& vals = sorted_ss_mat_Q_D[ss_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _ss_mat_Q[ij] = Q_ij;
        _ss_mat_D[ij] = D_ij;

        _ss_first_inds[ij]  = static_cast<uint32_t>(i);
        _ss_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _ss_pair_data[ij] = S_ij_00;
    }

    _sp_mat_Q       = std::vector<double>  (sp_prim_pair_count);
    _sp_mat_D       = std::vector<double>  (sp_prim_pair_count);
    _sp_first_inds  = std::vector<uint32_t>(sp_prim_pair_count);
    _sp_second_inds = std::vector<uint32_t>(sp_prim_pair_count);

    _sp_pair_data = std::vector<double>(sp_prim_pair_count);

    for (int64_t ij = 0; ij < sp_prim_pair_count; ij++)
    {
        const auto& vals = sorted_sp_mat_Q_D[sp_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _sp_mat_Q[ij] = Q_ij;
        _sp_mat_D[ij] = D_ij;

        _sp_first_inds[ij]  = static_cast<uint32_t>(i);
        _sp_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _sp_pair_data[ij] = S_ij_00;
    }

    _sd_mat_Q       = std::vector<double>  (sd_prim_pair_count);
    _sd_mat_D       = std::vector<double>  (sd_prim_pair_count);
    _sd_first_inds  = std::vector<uint32_t>(sd_prim_pair_count);
    _sd_second_inds = std::vector<uint32_t>(sd_prim_pair_count);

    _sd_pair_data = std::vector<double>(sd_prim_pair_count);

    for (int64_t ij = 0; ij < sd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_sd_mat_Q_D[sd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _sd_mat_Q[ij] = Q_ij;
        _sd_mat_D[ij] = D_ij;

        _sd_first_inds[ij]  = static_cast<uint32_t>(i);
        _sd_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _sd_pair_data[ij] = S_ij_00;
    }

    _pp_mat_Q       = std::vector<double>  (pp_prim_pair_count);
    _pp_mat_D       = std::vector<double>  (pp_prim_pair_count);
    _pp_first_inds  = std::vector<uint32_t>(pp_prim_pair_count);
    _pp_second_inds = std::vector<uint32_t>(pp_prim_pair_count);

    _pp_pair_data = std::vector<double>(pp_prim_pair_count);

    for (int64_t ij = 0; ij < pp_prim_pair_count; ij++)
    {
        const auto& vals = sorted_pp_mat_Q_D[pp_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _pp_mat_Q[ij] = Q_ij;
        _pp_mat_D[ij] = D_ij;

        _pp_first_inds[ij]  = static_cast<uint32_t>(i);
        _pp_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _pp_pair_data[ij] = S_ij_00;
    }

    _pd_mat_Q       = std::vector<double>  (pd_prim_pair_count);
    _pd_mat_D       = std::vector<double>  (pd_prim_pair_count);
    _pd_first_inds  = std::vector<uint32_t>(pd_prim_pair_count);
    _pd_second_inds = std::vector<uint32_t>(pd_prim_pair_count);

    _pd_pair_data = std::vector<double>(pd_prim_pair_count);

    for (int64_t ij = 0; ij < pd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_pd_mat_Q_D[pd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _pd_mat_Q[ij] = Q_ij;
        _pd_mat_D[ij] = D_ij;

        _pd_first_inds[ij]  = static_cast<uint32_t>(i);
        _pd_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _pd_pair_data[ij] = S_ij_00;
    }

    _dd_mat_Q       = std::vector<double>  (dd_prim_pair_count);
    _dd_mat_D       = std::vector<double>  (dd_prim_pair_count);
    _dd_first_inds  = std::vector<uint32_t>(dd_prim_pair_count);
    _dd_second_inds = std::vector<uint32_t>(dd_prim_pair_count);

    _dd_pair_data = std::vector<double>(dd_prim_pair_count);

    for (int64_t ij = 0; ij < dd_prim_pair_count; ij++)
    {
        const auto& vals = sorted_dd_mat_Q_D[dd_prim_pair_count - 1 - ij];

        auto i    = std::get<1>(vals);
        auto j    = std::get<2>(vals);
        auto Q_ij = std::get<3>(vals);
        auto D_ij = std::get<4>(vals);

        _dd_mat_Q[ij] = Q_ij;
        _dd_mat_D[ij] = D_ij;

        _dd_first_inds[ij]  = static_cast<uint32_t>(i);
        _dd_second_inds[ij] = static_cast<uint32_t>(j);

        // ij pair data:

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

        const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

        const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

        _dd_pair_data[ij] = S_ij_00;
    }
}

auto
CScreeningData::findMaxDensities(const int64_t s_prim_count,
                                 const int64_t p_prim_count,
                                 const int64_t d_prim_count,
                                 const std::vector<uint32_t>& s_prim_aoinds,
                                 const std::vector<uint32_t>& p_prim_aoinds,
                                 const std::vector<uint32_t>& d_prim_aoinds,
                                 const int64_t naos,
                                 const double* dens_ptr) -> void
{
    _ss_max_D = 0.0;
    _sp_max_D = 0.0;
    _sd_max_D = 0.0;
    _pp_max_D = 0.0;
    _pd_max_D = 0.0;
    _dd_max_D = 0.0;

    // S-S gto block pair and S-P gto block pair

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        // S-S gto block pair

        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = i; j < s_prim_count; j++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

            if (std::fabs(D_ij) > _ss_max_D) _ss_max_D = std::fabs(D_ij);
        }

        // S-P gto block pair

        for (int64_t j = 0; j < p_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 3; j_cart++)
            {
                const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                if (std::fabs(D_ij) > _sp_max_D) _sp_max_D = std::fabs(D_ij);
            }
        }

        // S-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t j_cart = 0; j_cart < 6; j_cart++)
            {
                const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                if (std::fabs(D_ij) > _sd_max_D) _sd_max_D = std::fabs(D_ij);
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count; i++)
    {
        // P-P gto block pair

        for (int64_t j = i; j < p_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 3; j_cart++)
                {
                    const auto j_cgto = p_prim_aoinds[j + p_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                    if (std::fabs(D_ij) > _pp_max_D) _pp_max_D = std::fabs(D_ij);
                }
            }
        }
    
        // P-D gto block pair

        for (int64_t j = 0; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 3; i_cart++)
            {
                const auto i_cgto = p_prim_aoinds[i + p_prim_count * i_cart];

                for (int64_t j_cart = 0; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                    if (std::fabs(D_ij) > _pd_max_D) _pd_max_D = std::fabs(D_ij);
                }
            }
        }    
    }

    for (int64_t i = 0; i < d_prim_count; i++)
    {
        // D-D gto block pair

        for (int64_t j = i; j < d_prim_count; j++)
        {
            for (int64_t i_cart = 0; i_cart < 6; i_cart++)
            {
                const auto i_cgto = d_prim_aoinds[i + d_prim_count * i_cart];

                const auto j_cart_start = (j == i ? i_cart : 0);

                for (int64_t j_cart = j_cart_start; j_cart < 6; j_cart++)
                {
                    const auto j_cgto = d_prim_aoinds[j + d_prim_count * j_cart];

                    const auto D_ij = dens_ptr[i_cgto * naos + j_cgto];

                    if (std::fabs(D_ij) > _dd_max_D) _dd_max_D = std::fabs(D_ij);
                }
            }
        }
    }
}

auto CScreeningData::get_ss_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _ss_first_inds_local[gpu_id]; }
auto CScreeningData::get_sp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sp_first_inds_local[gpu_id]; }
auto CScreeningData::get_sd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sd_first_inds_local[gpu_id]; }
auto CScreeningData::get_pp_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pp_first_inds_local[gpu_id]; }
auto CScreeningData::get_pd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pd_first_inds_local[gpu_id]; }
auto CScreeningData::get_dd_first_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _dd_first_inds_local[gpu_id]; }

auto CScreeningData::get_ss_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _ss_second_inds_local[gpu_id]; }
auto CScreeningData::get_sp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sp_second_inds_local[gpu_id]; }
auto CScreeningData::get_sd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _sd_second_inds_local[gpu_id]; }
auto CScreeningData::get_pp_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pp_second_inds_local[gpu_id]; }
auto CScreeningData::get_pd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _pd_second_inds_local[gpu_id]; }
auto CScreeningData::get_dd_second_inds_local(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _dd_second_inds_local[gpu_id]; }

auto CScreeningData::get_ss_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _ss_mat_Q_local[gpu_id]; }
auto CScreeningData::get_sp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sp_mat_Q_local[gpu_id]; }
auto CScreeningData::get_sd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sd_mat_Q_local[gpu_id]; }
auto CScreeningData::get_pp_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pp_mat_Q_local[gpu_id]; }
auto CScreeningData::get_pd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pd_mat_Q_local[gpu_id]; }
auto CScreeningData::get_dd_mat_Q_local(const int64_t gpu_id) const -> const std::vector<double>& { return _dd_mat_Q_local[gpu_id]; }

auto CScreeningData::get_ss_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _ss_pair_data_local[gpu_id]; }
auto CScreeningData::get_sp_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sp_pair_data_local[gpu_id]; }
auto CScreeningData::get_sd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _sd_pair_data_local[gpu_id]; }
auto CScreeningData::get_pp_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pp_pair_data_local[gpu_id]; }
auto CScreeningData::get_pd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _pd_pair_data_local[gpu_id]; }
auto CScreeningData::get_dd_pair_data_local(const int64_t gpu_id) const -> const std::vector<double>& { return _dd_pair_data_local[gpu_id]; }

auto CScreeningData::get_ss_first_inds() const -> const std::vector<uint32_t>& { return _ss_first_inds; }
auto CScreeningData::get_sp_first_inds() const -> const std::vector<uint32_t>& { return _sp_first_inds; }
auto CScreeningData::get_sd_first_inds() const -> const std::vector<uint32_t>& { return _sd_first_inds; }
auto CScreeningData::get_pp_first_inds() const -> const std::vector<uint32_t>& { return _pp_first_inds; }
auto CScreeningData::get_pd_first_inds() const -> const std::vector<uint32_t>& { return _pd_first_inds; }
auto CScreeningData::get_dd_first_inds() const -> const std::vector<uint32_t>& { return _dd_first_inds; }

auto CScreeningData::get_ss_second_inds() const -> const std::vector<uint32_t>& { return _ss_second_inds; }
auto CScreeningData::get_sp_second_inds() const -> const std::vector<uint32_t>& { return _sp_second_inds; }
auto CScreeningData::get_sd_second_inds() const -> const std::vector<uint32_t>& { return _sd_second_inds; }
auto CScreeningData::get_pp_second_inds() const -> const std::vector<uint32_t>& { return _pp_second_inds; }
auto CScreeningData::get_pd_second_inds() const -> const std::vector<uint32_t>& { return _pd_second_inds; }
auto CScreeningData::get_dd_second_inds() const -> const std::vector<uint32_t>& { return _dd_second_inds; }

auto CScreeningData::get_ss_mat_Q() const -> const std::vector<double>& { return _ss_mat_Q; }
auto CScreeningData::get_sp_mat_Q() const -> const std::vector<double>& { return _sp_mat_Q; }
auto CScreeningData::get_sd_mat_Q() const -> const std::vector<double>& { return _sd_mat_Q; }
auto CScreeningData::get_pp_mat_Q() const -> const std::vector<double>& { return _pp_mat_Q; }
auto CScreeningData::get_pd_mat_Q() const -> const std::vector<double>& { return _pd_mat_Q; }
auto CScreeningData::get_dd_mat_Q() const -> const std::vector<double>& { return _dd_mat_Q; }

auto CScreeningData::get_ss_mat_D() const -> const std::vector<double>& { return _ss_mat_D; }
auto CScreeningData::get_sp_mat_D() const -> const std::vector<double>& { return _sp_mat_D; }
auto CScreeningData::get_sd_mat_D() const -> const std::vector<double>& { return _sd_mat_D; }
auto CScreeningData::get_pp_mat_D() const -> const std::vector<double>& { return _pp_mat_D; }
auto CScreeningData::get_pd_mat_D() const -> const std::vector<double>& { return _pd_mat_D; }
auto CScreeningData::get_dd_mat_D() const -> const std::vector<double>& { return _dd_mat_D; }

auto CScreeningData::get_ss_max_D() const -> double { return _ss_max_D; }
auto CScreeningData::get_sp_max_D() const -> double { return _sp_max_D; }
auto CScreeningData::get_sd_max_D() const -> double { return _sd_max_D; }
auto CScreeningData::get_pp_max_D() const -> double { return _pp_max_D; }
auto CScreeningData::get_pd_max_D() const -> double { return _pd_max_D; }
auto CScreeningData::get_dd_max_D() const -> double { return _dd_max_D; }

auto CScreeningData::get_ss_pair_data() const -> const std::vector<double>& { return _ss_pair_data; }
auto CScreeningData::get_sp_pair_data() const -> const std::vector<double>& { return _sp_pair_data; }
auto CScreeningData::get_sd_pair_data() const -> const std::vector<double>& { return _sd_pair_data; }
auto CScreeningData::get_pp_pair_data() const -> const std::vector<double>& { return _pp_pair_data; }
auto CScreeningData::get_pd_pair_data() const -> const std::vector<double>& { return _pd_pair_data; }
auto CScreeningData::get_dd_pair_data() const -> const std::vector<double>& { return _dd_pair_data; }

auto CScreeningData::get_Q_K_ss() const -> const std::vector<double>& { return _Q_K_ss; }
auto CScreeningData::get_Q_K_sp() const -> const std::vector<double>& { return _Q_K_sp; }
auto CScreeningData::get_Q_K_ps() const -> const std::vector<double>& { return _Q_K_ps; }
auto CScreeningData::get_Q_K_sd() const -> const std::vector<double>& { return _Q_K_sd; }
auto CScreeningData::get_Q_K_ds() const -> const std::vector<double>& { return _Q_K_ds; }
auto CScreeningData::get_Q_K_pp() const -> const std::vector<double>& { return _Q_K_pp; }
auto CScreeningData::get_Q_K_pd() const -> const std::vector<double>& { return _Q_K_pd; }
auto CScreeningData::get_Q_K_dp() const -> const std::vector<double>& { return _Q_K_dp; }
auto CScreeningData::get_Q_K_dd() const -> const std::vector<double>& { return _Q_K_dd; }

auto CScreeningData::get_D_inds_K_ss() const -> const std::vector<uint32_t>& { return _D_inds_K_ss; }
auto CScreeningData::get_D_inds_K_sp() const -> const std::vector<uint32_t>& { return _D_inds_K_sp; }
auto CScreeningData::get_D_inds_K_ps() const -> const std::vector<uint32_t>& { return _D_inds_K_ps; }
auto CScreeningData::get_D_inds_K_sd() const -> const std::vector<uint32_t>& { return _D_inds_K_sd; }
auto CScreeningData::get_D_inds_K_ds() const -> const std::vector<uint32_t>& { return _D_inds_K_ds; }
auto CScreeningData::get_D_inds_K_pp() const -> const std::vector<uint32_t>& { return _D_inds_K_pp; }
auto CScreeningData::get_D_inds_K_pd() const -> const std::vector<uint32_t>& { return _D_inds_K_pd; }
auto CScreeningData::get_D_inds_K_dp() const -> const std::vector<uint32_t>& { return _D_inds_K_dp; }
auto CScreeningData::get_D_inds_K_dd() const -> const std::vector<uint32_t>& { return _D_inds_K_dd; }

auto CScreeningData::get_pair_displs_K_ss() const -> const std::vector<uint32_t>& { return _pair_displs_K_ss; }
auto CScreeningData::get_pair_displs_K_sp() const -> const std::vector<uint32_t>& { return _pair_displs_K_sp; }
auto CScreeningData::get_pair_displs_K_ps() const -> const std::vector<uint32_t>& { return _pair_displs_K_ps; }
auto CScreeningData::get_pair_displs_K_sd() const -> const std::vector<uint32_t>& { return _pair_displs_K_sd; }
auto CScreeningData::get_pair_displs_K_ds() const -> const std::vector<uint32_t>& { return _pair_displs_K_ds; }
auto CScreeningData::get_pair_displs_K_pp() const -> const std::vector<uint32_t>& { return _pair_displs_K_pp; }
auto CScreeningData::get_pair_displs_K_pd() const -> const std::vector<uint32_t>& { return _pair_displs_K_pd; }
auto CScreeningData::get_pair_displs_K_dp() const -> const std::vector<uint32_t>& { return _pair_displs_K_dp; }
auto CScreeningData::get_pair_displs_K_dd() const -> const std::vector<uint32_t>& { return _pair_displs_K_dd; }

auto CScreeningData::get_pair_counts_K_ss() const -> const std::vector<uint32_t>& { return _pair_counts_K_ss; }
auto CScreeningData::get_pair_counts_K_sp() const -> const std::vector<uint32_t>& { return _pair_counts_K_sp; }
auto CScreeningData::get_pair_counts_K_ps() const -> const std::vector<uint32_t>& { return _pair_counts_K_ps; }
auto CScreeningData::get_pair_counts_K_sd() const -> const std::vector<uint32_t>& { return _pair_counts_K_sd; }
auto CScreeningData::get_pair_counts_K_ds() const -> const std::vector<uint32_t>& { return _pair_counts_K_ds; }
auto CScreeningData::get_pair_counts_K_pp() const -> const std::vector<uint32_t>& { return _pair_counts_K_pp; }
auto CScreeningData::get_pair_counts_K_pd() const -> const std::vector<uint32_t>& { return _pair_counts_K_pd; }
auto CScreeningData::get_pair_counts_K_dp() const -> const std::vector<uint32_t>& { return _pair_counts_K_dp; }
auto CScreeningData::get_pair_counts_K_dd() const -> const std::vector<uint32_t>& { return _pair_counts_K_dd; }

auto CScreeningData::get_pair_data_K_ss() const -> const std::vector<double>& { return _pair_data_K_ss; }
auto CScreeningData::get_pair_data_K_sp() const -> const std::vector<double>& { return _pair_data_K_sp; }
auto CScreeningData::get_pair_data_K_ps() const -> const std::vector<double>& { return _pair_data_K_ps; }
auto CScreeningData::get_pair_data_K_sd() const -> const std::vector<double>& { return _pair_data_K_sd; }
auto CScreeningData::get_pair_data_K_ds() const -> const std::vector<double>& { return _pair_data_K_ds; }
auto CScreeningData::get_pair_data_K_pp() const -> const std::vector<double>& { return _pair_data_K_pp; }
auto CScreeningData::get_pair_data_K_pd() const -> const std::vector<double>& { return _pair_data_K_pd; }
auto CScreeningData::get_pair_data_K_dp() const -> const std::vector<double>& { return _pair_data_K_dp; }
auto CScreeningData::get_pair_data_K_dd() const -> const std::vector<double>& { return _pair_data_K_dd; }

auto CScreeningData::get_local_pair_inds_i_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_ss[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_ss(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_ss[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_sp[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_sp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_sp[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_sd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_sd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_sd[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_pp[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_pp(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_pp[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_pd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_pd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_pd[gpu_id]; }

auto CScreeningData::get_local_pair_inds_i_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_i_for_K_dd[gpu_id]; }
auto CScreeningData::get_local_pair_inds_k_for_K_dd(const int64_t gpu_id) const -> const std::vector<uint32_t>& { return _local_pair_inds_k_for_K_dd[gpu_id]; }

auto CScreeningData::get_mat_Q_full(const int64_t s_prim_count, const int64_t p_prim_count, const int64_t d_prim_count) const -> CDenseMatrix
{
    const auto cart_naos = s_prim_count + p_prim_count * 3 + d_prim_count * 6;

    CDenseMatrix mat_Q_full(cart_naos, cart_naos);

    mat_Q_full.zero();

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (_Q_matrix_ss.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i)[j] = _Q_matrix_ss.row(i)[j];
            }
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            auto j_full = s_prim_count + j;

            if (_Q_matrix_sp.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i)[j_full] = _Q_matrix_sp.row(i)[j];
                mat_Q_full.row(j_full)[i] = _Q_matrix_sp.row(i)[j];
            }
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (_Q_matrix_sd.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i)[j_full] = _Q_matrix_sd.row(i)[j];
                mat_Q_full.row(j_full)[i] = _Q_matrix_sd.row(i)[j];
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        auto i_full = s_prim_count + i;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            auto j_full = s_prim_count + j;

            if (_Q_matrix_pp.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i_full)[j_full] = _Q_matrix_pp.row(i)[j];
            }
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (_Q_matrix_pd.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i_full)[j_full] = _Q_matrix_pd.row(i)[j];
                mat_Q_full.row(j_full)[i_full] = _Q_matrix_pd.row(i)[j];
            }
        }
    }

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        auto i_full = s_prim_count + p_prim_count * 3 + i;

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (_Q_matrix_dd.row(i)[j] > _pair_threshold)
            {
                mat_Q_full.row(i_full)[j_full] = _Q_matrix_dd.row(i)[j];
            }
        }
    }

    return mat_Q_full;
}

auto CScreeningData::get_mat_D_abs_full(const int64_t s_prim_count,
                                        const int64_t p_prim_count,
                                        const int64_t d_prim_count,
                                        const std::vector<uint32_t>& s_prim_aoinds,
                                        const std::vector<uint32_t>& p_prim_aoinds,
                                        const std::vector<uint32_t>& d_prim_aoinds,
                                        const int64_t naos,
                                        const double* dens_ptr) const -> CDenseMatrix
{
    const auto cart_naos = s_prim_count + p_prim_count * 3 + d_prim_count * 6;

    CDenseMatrix mat_D_abs_full(cart_naos, cart_naos);

    mat_D_abs_full.zero();

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto i_cgto = s_prim_aoinds[i];

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            const auto j_cgto = s_prim_aoinds[j];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i)[j] = D_ij;
            }
        }

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + j;

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i)[j_full] = D_ij;
                mat_D_abs_full.row(j_full)[i] = D_ij;
            }
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i)[j_full] = D_ij;
                mat_D_abs_full.row(j_full)[i] = D_ij;
            }
        }
    }

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        const auto i_cgto = p_prim_aoinds[(i / 3) + p_prim_count * (i % 3)];

        auto i_full = s_prim_count + i;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            const auto j_cgto = p_prim_aoinds[(j / 3) + p_prim_count * (j % 3)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + j;

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i_full)[j_full] = D_ij;
            }
        }

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i_full)[j_full] = D_ij;
                mat_D_abs_full.row(j_full)[i_full] = D_ij;
            }
        }
    }

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        const auto i_cgto = d_prim_aoinds[(i / 6) + d_prim_count * (i % 6)];

        auto i_full = s_prim_count + p_prim_count * 3 + i;

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            const auto j_cgto = d_prim_aoinds[(j / 6) + d_prim_count * (j % 6)];

            const auto D_ij = std::fabs(dens_ptr[i_cgto * naos + j_cgto]);

            auto j_full = s_prim_count + p_prim_count * 3 + j;

            if (D_ij > _density_threshold)
            {
                mat_D_abs_full.row(i_full)[j_full] = D_ij;
            }
        }
    }

    return mat_D_abs_full;
}

auto CScreeningData::form_Q_and_D_inds_for_K() -> void
{
    const auto s_prim_count = _s_prim_count;
    const auto p_prim_count = _p_prim_count;
    const auto d_prim_count = _d_prim_count;

    const auto& s_prim_aoinds = _s_prim_aoinds;
    const auto& p_prim_aoinds = _p_prim_aoinds;
    const auto& d_prim_aoinds = _d_prim_aoinds;

    const auto& s_prim_info = _s_prim_info;
    const auto& p_prim_info = _p_prim_info;
    const auto& d_prim_info = _d_prim_info;

    // Q_ss and D_ss for K

    int64_t ss_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto Q_i = _Q_matrix_ss.row(i);

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (Q_i[j] > _pair_threshold) ss_prim_pair_count_for_K++;
        }
    }

    _Q_K_ss      = std::vector<double>(ss_prim_pair_count_for_K);
    _D_inds_K_ss = std::vector<uint32_t>(ss_prim_pair_count_for_K);

    _pair_data_K_ss = std::vector<double>(ss_prim_pair_count_for_K);

    _pair_displs_K_ss = std::vector<uint32_t>(s_prim_count);
    _pair_counts_K_ss = std::vector<uint32_t>(s_prim_count);

    for (int64_t i = 0, displ = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_ss.row(i);

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_ss[displ + j]      = q_val;
            _D_inds_K_ss[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = s_prim_info[i + s_prim_count * 0];
            const auto c_i = s_prim_info[i + s_prim_count * 1];
            const auto x_i = s_prim_info[i + s_prim_count * 2];
            const auto y_i = s_prim_info[i + s_prim_count * 3];
            const auto z_i = s_prim_info[i + s_prim_count * 4];

            const auto a_j = s_prim_info[j_idx + s_prim_count * 0];
            const auto c_j = s_prim_info[j_idx + s_prim_count * 1];
            const auto x_j = s_prim_info[j_idx + s_prim_count * 2];
            const auto y_j = s_prim_info[j_idx + s_prim_count * 3];
            const auto z_j = s_prim_info[j_idx + s_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_ss[displ + j] = S_ij_00;
        }

        _pair_displs_K_ss[i] = displ;
        _pair_counts_K_ss[i] = count;

        displ += count;
    }

    // Q_sp and D_sp for K

    int64_t sp_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto Q_i = _Q_matrix_sp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (Q_i[j] > _pair_threshold) sp_prim_pair_count_for_K++;
        }
    }

    _Q_K_sp      = std::vector<double>(sp_prim_pair_count_for_K);
    _D_inds_K_sp = std::vector<uint32_t>(sp_prim_pair_count_for_K);

    _pair_data_K_sp = std::vector<double>(sp_prim_pair_count_for_K);

    _pair_displs_K_sp = std::vector<uint32_t>(s_prim_count);
    _pair_counts_K_sp = std::vector<uint32_t>(s_prim_count);

    for (int64_t i = 0, displ = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_sp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_sp[displ + j]      = q_val;
            _D_inds_K_sp[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = s_prim_info[i + s_prim_count * 0];
            const auto c_i = s_prim_info[i + s_prim_count * 1];
            const auto x_i = s_prim_info[i + s_prim_count * 2];
            const auto y_i = s_prim_info[i + s_prim_count * 3];
            const auto z_i = s_prim_info[i + s_prim_count * 4];

            const auto a_j = p_prim_info[j_idx / 3 + p_prim_count * 0];
            const auto c_j = p_prim_info[j_idx / 3 + p_prim_count * 1];
            const auto x_j = p_prim_info[j_idx / 3 + p_prim_count * 2];
            const auto y_j = p_prim_info[j_idx / 3 + p_prim_count * 3];
            const auto z_j = p_prim_info[j_idx / 3 + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_sp[displ + j] = S_ij_00;
        }

        _pair_displs_K_sp[i] = displ;
        _pair_counts_K_sp[i] = count;

        displ += count;
    }

    // Q_ps and D_ps for K

    int64_t ps_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (_Q_matrix_sp.row(j)[i] > _pair_threshold) ps_prim_pair_count_for_K++;
        }
    }

    _Q_K_ps      = std::vector<double>(ps_prim_pair_count_for_K);
    _D_inds_K_ps = std::vector<uint32_t>(ps_prim_pair_count_for_K);

    _pair_data_K_ps = std::vector<double>(ps_prim_pair_count_for_K);

    _pair_displs_K_ps = std::vector<uint32_t>(p_prim_count * 3);
    _pair_counts_K_ps = std::vector<uint32_t>(p_prim_count * 3);

    for (int64_t i = 0, displ = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (_Q_matrix_sp.row(j)[i] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_sp.row(j)[i], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_ps[displ + j]      = q_val;
            _D_inds_K_ps[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
            const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
            const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
            const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
            const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

            const auto a_j = s_prim_info[j_idx + s_prim_count * 0];
            const auto c_j = s_prim_info[j_idx + s_prim_count * 1];
            const auto x_j = s_prim_info[j_idx + s_prim_count * 2];
            const auto y_j = s_prim_info[j_idx + s_prim_count * 3];
            const auto z_j = s_prim_info[j_idx + s_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_ps[displ + j] = S_ij_00;
        }

        _pair_displs_K_ps[i] = displ;
        _pair_counts_K_ps[i] = count;

        displ += count;
    }

    // Q_sd and D_sd for K

    int64_t sd_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        const auto Q_i = _Q_matrix_sd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold) sd_prim_pair_count_for_K++;
        }
    }

    _Q_K_sd      = std::vector<double>(sd_prim_pair_count_for_K);
    _D_inds_K_sd = std::vector<uint32_t>(sd_prim_pair_count_for_K);

    _pair_data_K_sd = std::vector<double>(sd_prim_pair_count_for_K);

    _pair_displs_K_sd = std::vector<uint32_t>(s_prim_count);
    _pair_counts_K_sd = std::vector<uint32_t>(s_prim_count);

    for (int64_t i = 0, displ = 0; i < s_prim_count; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_sd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_sd[displ + j]      = q_val;
            _D_inds_K_sd[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = s_prim_info[i + s_prim_count * 0];
            const auto c_i = s_prim_info[i + s_prim_count * 1];
            const auto x_i = s_prim_info[i + s_prim_count * 2];
            const auto y_i = s_prim_info[i + s_prim_count * 3];
            const auto z_i = s_prim_info[i + s_prim_count * 4];

            const auto a_j = d_prim_info[j_idx / 6 + d_prim_count * 0];
            const auto c_j = d_prim_info[j_idx / 6 + d_prim_count * 1];
            const auto x_j = d_prim_info[j_idx / 6 + d_prim_count * 2];
            const auto y_j = d_prim_info[j_idx / 6 + d_prim_count * 3];
            const auto z_j = d_prim_info[j_idx / 6 + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_sd[displ + j] = S_ij_00;
        }

        _pair_displs_K_sd[i] = displ;
        _pair_counts_K_sd[i] = count;

        displ += count;
    }

    // Q_ds and D_ds for K

    int64_t ds_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (_Q_matrix_sd.row(j)[i] > _pair_threshold) ds_prim_pair_count_for_K++;
        }
    }

    _Q_K_ds      = std::vector<double>(ds_prim_pair_count_for_K);
    _D_inds_K_ds = std::vector<uint32_t>(ds_prim_pair_count_for_K);

    _pair_data_K_ds = std::vector<double>(ds_prim_pair_count_for_K);

    _pair_displs_K_ds = std::vector<uint32_t>(d_prim_count * 6);
    _pair_counts_K_ds = std::vector<uint32_t>(d_prim_count * 6);

    for (int64_t i = 0, displ = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < s_prim_count; j++)
        {
            if (_Q_matrix_sd.row(j)[i] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_sd.row(j)[i], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_ds[displ + j]      = q_val;
            _D_inds_K_ds[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
            const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
            const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
            const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
            const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

            const auto a_j = s_prim_info[j_idx + s_prim_count * 0];
            const auto c_j = s_prim_info[j_idx + s_prim_count * 1];
            const auto x_j = s_prim_info[j_idx + s_prim_count * 2];
            const auto y_j = s_prim_info[j_idx + s_prim_count * 3];
            const auto z_j = s_prim_info[j_idx + s_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_ds[displ + j] = S_ij_00;
        }

        _pair_displs_K_ds[i] = displ;
        _pair_counts_K_ds[i] = count;

        displ += count;
    }

    // Q_pp and D_pp for K

    int64_t pp_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        const auto Q_i = _Q_matrix_pp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (Q_i[j] > _pair_threshold) pp_prim_pair_count_for_K++;
        }
    }

    _Q_K_pp      = std::vector<double>(pp_prim_pair_count_for_K);
    _D_inds_K_pp = std::vector<uint32_t>(pp_prim_pair_count_for_K);

    _pair_data_K_pp = std::vector<double>(pp_prim_pair_count_for_K);

    _pair_displs_K_pp = std::vector<uint32_t>(p_prim_count * 3);
    _pair_counts_K_pp = std::vector<uint32_t>(p_prim_count * 3);

    for (int64_t i = 0, displ = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_pp.row(i);

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_pp[displ + j]      = q_val;
            _D_inds_K_pp[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
            const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
            const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
            const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
            const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

            const auto a_j = p_prim_info[j_idx / 3 + p_prim_count * 0];
            const auto c_j = p_prim_info[j_idx / 3 + p_prim_count * 1];
            const auto x_j = p_prim_info[j_idx / 3 + p_prim_count * 2];
            const auto y_j = p_prim_info[j_idx / 3 + p_prim_count * 3];
            const auto z_j = p_prim_info[j_idx / 3 + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_pp[displ + j] = S_ij_00;
        }

        _pair_displs_K_pp[i] = displ;
        _pair_counts_K_pp[i] = count;

        displ += count;
    }

    // Q_pd and D_pd for K

    int64_t pd_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        const auto Q_i = _Q_matrix_pd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold) pd_prim_pair_count_for_K++;
        }
    }

    _Q_K_pd      = std::vector<double>(pd_prim_pair_count_for_K);
    _D_inds_K_pd = std::vector<uint32_t>(pd_prim_pair_count_for_K);

    _pair_data_K_pd = std::vector<double>(pd_prim_pair_count_for_K);

    _pair_displs_K_pd = std::vector<uint32_t>(p_prim_count * 3);
    _pair_counts_K_pd = std::vector<uint32_t>(p_prim_count * 3);

    for (int64_t i = 0, displ = 0; i < p_prim_count * 3; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_pd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_pd[displ + j]      = q_val;
            _D_inds_K_pd[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = p_prim_info[i / 3 + p_prim_count * 0];
            const auto c_i = p_prim_info[i / 3 + p_prim_count * 1];
            const auto x_i = p_prim_info[i / 3 + p_prim_count * 2];
            const auto y_i = p_prim_info[i / 3 + p_prim_count * 3];
            const auto z_i = p_prim_info[i / 3 + p_prim_count * 4];

            const auto a_j = d_prim_info[j_idx / 6 + d_prim_count * 0];
            const auto c_j = d_prim_info[j_idx / 6 + d_prim_count * 1];
            const auto x_j = d_prim_info[j_idx / 6 + d_prim_count * 2];
            const auto y_j = d_prim_info[j_idx / 6 + d_prim_count * 3];
            const auto z_j = d_prim_info[j_idx / 6 + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_pd[displ + j] = S_ij_00;
        }

        _pair_displs_K_pd[i] = displ;
        _pair_counts_K_pd[i] = count;

        displ += count;
    }

    // Q_dp and D_dp for K

    int64_t dp_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (_Q_matrix_pd.row(j)[i] > _pair_threshold) dp_prim_pair_count_for_K++;
        }
    }

    _Q_K_dp      = std::vector<double>(dp_prim_pair_count_for_K);
    _D_inds_K_dp = std::vector<uint32_t>(dp_prim_pair_count_for_K);

    _pair_data_K_dp = std::vector<double>(dp_prim_pair_count_for_K);

    _pair_displs_K_dp = std::vector<uint32_t>(d_prim_count * 6);
    _pair_counts_K_dp = std::vector<uint32_t>(d_prim_count * 6);

    for (int64_t i = 0, displ = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        for (int64_t j = 0; j < p_prim_count * 3; j++)
        {
            if (_Q_matrix_pd.row(j)[i] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(_Q_matrix_pd.row(j)[i], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_dp[displ + j]      = q_val;
            _D_inds_K_dp[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
            const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
            const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
            const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
            const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

            const auto a_j = p_prim_info[j_idx / 3 + p_prim_count * 0];
            const auto c_j = p_prim_info[j_idx / 3 + p_prim_count * 1];
            const auto x_j = p_prim_info[j_idx / 3 + p_prim_count * 2];
            const auto y_j = p_prim_info[j_idx / 3 + p_prim_count * 3];
            const auto z_j = p_prim_info[j_idx / 3 + p_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_dp[displ + j] = S_ij_00;
        }

        _pair_displs_K_dp[i] = displ;
        _pair_counts_K_dp[i] = count;

        displ += count;
    }

    // Q_dd and D_dd for K

    int64_t dd_prim_pair_count_for_K = 0;

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        const auto Q_i = _Q_matrix_dd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold) dd_prim_pair_count_for_K++;
        }
    }

    _Q_K_dd      = std::vector<double>(dd_prim_pair_count_for_K);
    _D_inds_K_dd = std::vector<uint32_t>(dd_prim_pair_count_for_K);

    _pair_data_K_dd = std::vector<double>(dd_prim_pair_count_for_K);

    _pair_displs_K_dd = std::vector<uint32_t>(d_prim_count * 6);
    _pair_counts_K_dd = std::vector<uint32_t>(d_prim_count * 6);

    for (int64_t i = 0, displ = 0; i < d_prim_count * 6; i++)
    {
        std::vector<std::tuple<double, int64_t>> Q_vec_sorted;

        const auto Q_i = _Q_matrix_dd.row(i);

        for (int64_t j = 0; j < d_prim_count * 6; j++)
        {
            if (Q_i[j] > _pair_threshold)
            {
                Q_vec_sorted.push_back(std::make_tuple(Q_i[j], j));
            }
        }

        std::sort(Q_vec_sorted.begin(), Q_vec_sorted.end());
        std::reverse(Q_vec_sorted.begin(), Q_vec_sorted.end());

        const auto count = static_cast<int64_t>(Q_vec_sorted.size());

        for (int64_t j = 0; j < count; j++)
        {
            const auto& q_ij = Q_vec_sorted[j];

            auto q_val = std::get<0>(q_ij);
            auto j_idx = std::get<1>(q_ij);

            _Q_K_dd[displ + j]      = q_val;
            _D_inds_K_dd[displ + j] = j_idx;

            // ij pair data:

            const auto a_i = d_prim_info[i / 6 + d_prim_count * 0];
            const auto c_i = d_prim_info[i / 6 + d_prim_count * 1];
            const auto x_i = d_prim_info[i / 6 + d_prim_count * 2];
            const auto y_i = d_prim_info[i / 6 + d_prim_count * 3];
            const auto z_i = d_prim_info[i / 6 + d_prim_count * 4];

            const auto a_j = d_prim_info[j_idx / 6 + d_prim_count * 0];
            const auto c_j = d_prim_info[j_idx / 6 + d_prim_count * 1];
            const auto x_j = d_prim_info[j_idx / 6 + d_prim_count * 2];
            const auto y_j = d_prim_info[j_idx / 6 + d_prim_count * 3];
            const auto z_j = d_prim_info[j_idx / 6 + d_prim_count * 4];

            const auto r2_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) + (z_j - z_i) * (z_j - z_i);

            const auto S_ij_00 = c_i * c_j * std::pow(MATH_CONST_PI / (a_i + a_j), 1.5) * std::exp(-a_i * a_j / (a_i + a_j) * r2_ij);

            _pair_data_K_dd[displ + j] = S_ij_00;
        }

        _pair_displs_K_dd[i] = displ;
        _pair_counts_K_dd[i] = count;

        displ += count;
    }
}

auto CScreeningData::form_pair_inds_for_K(const int64_t s_prim_count,
                                          const int64_t p_prim_count,
                                          const int64_t d_prim_count,
                                          const CDenseMatrix& Q_prime,
                                          const double Q_prime_thresh) -> void
{
    // TODO consider determining the maximum density associated
    // with the ik pair (i.e. max_D_jl for ik)

    std::vector<std::tuple<double, int64_t, int64_t>> ss_Qp_ik;
    std::vector<std::tuple<double, int64_t, int64_t>> sp_Qp_ik;
    std::vector<std::tuple<double, int64_t, int64_t>> sd_Qp_ik;
    std::vector<std::tuple<double, int64_t, int64_t>> pp_Qp_ik;
    std::vector<std::tuple<double, int64_t, int64_t>> pd_Qp_ik;
    std::vector<std::tuple<double, int64_t, int64_t>> dd_Qp_ik;

    // ss, sp, sd blocks

    for (int64_t i = 0; i < s_prim_count; i++)
    {
        for (int64_t k = i; k < s_prim_count; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(i)[k]);

            if (Qp_ik > Q_prime_thresh) ss_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }

        for (int64_t k = 0; k < p_prim_count * 3; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(i)[s_prim_count + k]);

            if (Qp_ik > Q_prime_thresh) sp_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }

        for (int64_t k = 0; k < d_prim_count * 6; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(i)[s_prim_count + p_prim_count * 3 + k]);

            if (Qp_ik > Q_prime_thresh) sd_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }
    }

    // pp, pd blocks

    for (int64_t i = 0; i < p_prim_count * 3; i++)
    {
        for (int64_t k = i; k < p_prim_count * 3; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(s_prim_count + i)[s_prim_count + k]);

            if (Qp_ik > Q_prime_thresh) pp_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }

        for (int64_t k = 0; k < d_prim_count * 6; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(s_prim_count + i)[s_prim_count + p_prim_count * 3 + k]);

            if (Qp_ik > Q_prime_thresh) pd_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }
    }

    // dd block

    for (int64_t i = 0; i < d_prim_count * 6; i++)
    {
        for (int64_t k = i; k < d_prim_count * 6; k++)
        {
            const auto Qp_ik = std::fabs(Q_prime.row(s_prim_count + p_prim_count * 3 + i)[s_prim_count + p_prim_count * 3 + k]);

            if (Qp_ik > Q_prime_thresh) dd_Qp_ik.push_back(std::make_tuple(Qp_ik, i, k));
        }
    }

    const auto ss_pair_count_for_K = static_cast<int64_t>(ss_Qp_ik.size());
    const auto sp_pair_count_for_K = static_cast<int64_t>(sp_Qp_ik.size());
    const auto sd_pair_count_for_K = static_cast<int64_t>(sd_Qp_ik.size());
    const auto pp_pair_count_for_K = static_cast<int64_t>(pp_Qp_ik.size());
    const auto pd_pair_count_for_K = static_cast<int64_t>(pd_Qp_ik.size());
    const auto dd_pair_count_for_K = static_cast<int64_t>(dd_Qp_ik.size());

    // std::stringstream ss;
    // ss << "preLinK screening\n";
    // ss << "  SS pair: " << static_cast<double>(ss_pair_count_for_K) / (s_prim_count * (s_prim_count + 1) / 2)<< "\n";
    // ss << "  SP pair: " << static_cast<double>(sp_pair_count_for_K) / (s_prim_count * p_prim_count * 3)<< "\n";
    // ss << "  SD pair: " << static_cast<double>(sd_pair_count_for_K) / (s_prim_count * d_prim_count * 6)<< "\n";
    // ss << "  PP pair: " << static_cast<double>(pp_pair_count_for_K) / (p_prim_count * 3 * (p_prim_count * 3 + 1) / 2)<< "\n";
    // ss << "  PD pair: " << static_cast<double>(pd_pair_count_for_K) / (p_prim_count * 3 * d_prim_count * 6)<< "\n";
    // ss << "  DD pair: " << static_cast<double>(dd_pair_count_for_K) / (d_prim_count * 6 * (d_prim_count * 6 + 1) / 2)<< "\n";
    // std::cout << ss.str() << "\n";

    // form local vectors

    _local_pair_inds_i_for_K_ss = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_ss = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_sp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_sp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_sd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_sd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_pp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_pp = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_pd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_pd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    _local_pair_inds_i_for_K_dd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);
    _local_pair_inds_k_for_K_dd = std::vector<std::vector<uint32_t>>(_num_gpus_per_node);

    auto rank = _rank;
    auto nnodes = _nnodes;

    for (int64_t gpu_id = 0; gpu_id < _num_gpus_per_node; gpu_id++)
    {
        auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
        auto gpu_count = nnodes * _num_gpus_per_node;

        auto ss_batch_size = mathfunc::batch_size(ss_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_ss[gpu_id] = std::vector<uint32_t>(ss_batch_size);
        _local_pair_inds_k_for_K_ss[gpu_id] = std::vector<uint32_t>(ss_batch_size);

        auto sp_batch_size = mathfunc::batch_size(sp_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_sp[gpu_id] = std::vector<uint32_t>(sp_batch_size);
        _local_pair_inds_k_for_K_sp[gpu_id] = std::vector<uint32_t>(sp_batch_size);

        auto sd_batch_size = mathfunc::batch_size(sd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_sd[gpu_id] = std::vector<uint32_t>(sd_batch_size);
        _local_pair_inds_k_for_K_sd[gpu_id] = std::vector<uint32_t>(sd_batch_size);

        auto pp_batch_size = mathfunc::batch_size(pp_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_pp[gpu_id] = std::vector<uint32_t>(pp_batch_size);
        _local_pair_inds_k_for_K_pp[gpu_id] = std::vector<uint32_t>(pp_batch_size);

        auto pd_batch_size = mathfunc::batch_size(pd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_pd[gpu_id] = std::vector<uint32_t>(pd_batch_size);
        _local_pair_inds_k_for_K_pd[gpu_id] = std::vector<uint32_t>(pd_batch_size);

        auto dd_batch_size = mathfunc::batch_size(dd_pair_count_for_K, gpu_rank, gpu_count);

        _local_pair_inds_i_for_K_dd[gpu_id] = std::vector<uint32_t>(dd_batch_size);
        _local_pair_inds_k_for_K_dd[gpu_id] = std::vector<uint32_t>(dd_batch_size);

        {
            auto gpu_rank = gpu_id + rank * _num_gpus_per_node;
            auto gpu_count = nnodes * _num_gpus_per_node;

            for (int64_t ik = gpu_rank, idx = 0; ik < ss_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = ss_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_ss[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_ss[gpu_id][idx] = k;
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < sp_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = sp_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_sp[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_sp[gpu_id][idx] = k;
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < sd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = sd_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_sd[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_sd[gpu_id][idx] = k;
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < pp_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = pp_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_pp[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_pp[gpu_id][idx] = k;
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < pd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = pd_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_pd[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_pd[gpu_id][idx] = k;
            }

            for (int64_t ik = gpu_rank, idx = 0; ik < dd_pair_count_for_K; ik+=gpu_count, idx++)
            {
                const auto& vals = dd_Qp_ik[ik];

                // auto Qp_ik = std::get<0>(vals);
                auto i = std::get<1>(vals);
                auto k = std::get<2>(vals);

                _local_pair_inds_i_for_K_dd[gpu_id][idx] = i;
                _local_pair_inds_k_for_K_dd[gpu_id][idx] = k;
            }
        }
    }
}
