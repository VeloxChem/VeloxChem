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

#include "RIFockDriver.hpp"

#include "ThreeCenterElectronRepulsionDriver.hpp"

CRIFockDriver::CRIFockDriver()

    : _j_metric(nullptr)

    , _k_metric(nullptr)

    , _eri_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::CRIFockDriver(const CSubMatrix& j_metric)

    : _j_metric(new CSubMatrix(j_metric))

    , _k_metric(nullptr)

    , _eri_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::CRIFockDriver(const CSubMatrix& j_metric,
                             const CSubMatrix& k_metric)

    : _j_metric(new CSubMatrix(j_metric))

    , _k_metric(new CSubMatrix(k_metric))

    , _eri_buffer(CT3FlatBuffer<double>())
{
    
}

CRIFockDriver::~CRIFockDriver()
{
    if (_j_metric != nullptr)
    {
        delete _j_metric;
    }
    
    if (_k_metric != nullptr)
    {
        delete _k_metric;
    }
}

auto
CRIFockDriver::prepare_buffers(const CMolecule       &molecule,
                               const CMolecularBasis &basis,
                               const CMolecularBasis &aux_basis) -> void
{
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    _eri_buffer = eri_drv.compute(basis, aux_basis, molecule);
}

auto
CRIFockDriver::prepare_buffers(const CMolecule&        molecule,
                               const CMolecularBasis&  basis,
                               const CMolecularBasis&  aux_basis,
                               const std::vector<int>& atoms) -> void
{
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    _eri_buffer = eri_drv.compute(basis, aux_basis, molecule, atoms);
}

auto
CRIFockDriver::compute(const CMatrix     &density,
                       const std::string &label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        auto t_vec = _trafo_gamma_vector(_comp_gamma_vector(density));
        
        fmat.assign_flat_values(_comp_j_vector(t_vec));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0); 
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::compute(const CMatrix     &density,
                       const std::vector<double>& gvector,
                       const std::string &label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        fmat.assign_flat_values(_comp_j_vector(gvector));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::local_compute(const CMatrix&             density,
                             const std::vector<double>& gvector,
                             const std::string&         label) const -> CMatrix
{
    CMatrix fmat(density);
    
    // compute Coulomb contribution to Fock matrix
    
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        fmat.assign_flat_values(_comp_local_j_vector(gvector));
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
    }
    
    return fmat;
}

auto
CRIFockDriver::compute_bq_vector(const CMatrix &density) const -> std::vector<double>
{
    return _trafo_gamma_vector(_comp_gamma_vector(density));
}

auto
CRIFockDriver::compute_local_bq_vector(const CMatrix &density) const -> std::vector<double>
{
    return _trafo_local_gamma_vector(_comp_gamma_vector(density));
}

auto
CRIFockDriver::compute_bq_vector(const CSubMatrix& lambda_p, const CSubMatrix& lambda_h) const -> std::vector<double>
{
    const auto mask = _eri_buffer.mask_indices();
    
    if (!mask.empty())
    {
        //std::cout << "mask size:" << std::endl;
    }
    
    return std::vector<double>();
}

auto
CRIFockDriver::_comp_gamma_vector(const CMatrix &density) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> gvec(ndim, 0.0);
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    auto buff_ptr = &_eri_buffer;
    
    auto gvec_ptr = gvec.data();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
         
        double fsum = 0.0;
        
#pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < nelems; j++)
        {
            fsum += dvec_ptr[j] * tint_ptr[j];
        }
        
        gvec_ptr[i] = fsum;
    }
    
    return gvec;
}

auto
CRIFockDriver::_trafo_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    std::vector<double> tvec(ndim, 0.0);
    
    auto gvec_ptr = gvector.data();
    
    auto tvec_ptr = tvec.data();
    
    auto tmat_ptr = _j_metric->data();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto row_ptr = &tmat_ptr[i * ndim];
        
        double fsum = 0.0;
        
#pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < ndim; j++)
        {
            fsum += row_ptr[j] * gvec_ptr[j];
        }
        
        tvec_ptr[i] = fsum;
    }
    
    return tvec;
}

auto
CRIFockDriver::_trafo_local_gamma_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _j_metric->number_of_rows();
    
    std::vector<double> tvec(ndim, 0.0);
    
    auto gvec_ptr = gvector.data();
    
    auto tvec_ptr = tvec.data();
    
    auto tmat_ptr = _j_metric->data();
    
    const auto mask_indices = _eri_buffer.mask_indices();
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto row_ptr = &tmat_ptr[i * ndim];
        
        double fsum = 0.0;
        
        for (const auto [gidx, lidx] : mask_indices)
        {
            fsum += row_ptr[gidx] * gvec_ptr[lidx];
        }
        
        tvec_ptr[i] = fsum;
    }
    
    return tvec;
}

auto
CRIFockDriver::_comp_j_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto ndim = _eri_buffer.aux_width();
    
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> jvec(nelems, 0.0);
    
    auto jvec_ptr = jvec.data();
    
    auto gvec_ptr = gvector.data();
    
    auto buff_ptr = &_eri_buffer;
    
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = buff_ptr->data(i);
        
        const auto fact = gvec_ptr[i];
         
#pragma omp simd
        for (size_t j = 0; j < nelems; j++)
        {
            jvec_ptr[j] += tint_ptr[j] * fact;
        }
    }
        
    return jvec;
}

auto
CRIFockDriver::_comp_local_j_vector(const std::vector<double>& gvector) const -> std::vector<double>
{
    const auto nrows = _eri_buffer.width();
    
    const auto nelems = nrows * (nrows + 1) / 2;
    
    std::vector<double> jvec(nelems, 0.0);
    
    auto jvec_ptr = jvec.data();
    
    auto gvec_ptr = gvector.data();
    
    auto buff_ptr = &_eri_buffer;
    
    const auto mask_indices = _eri_buffer.mask_indices();
    
    for (const auto [gidx, lidx] : mask_indices)
    {
        auto tint_ptr = buff_ptr->data(lidx);
        
        const auto fact = gvec_ptr[gidx];
         
#pragma omp simd
        for (size_t j = 0; j < nelems; j++)
        {
            jvec_ptr[j] += tint_ptr[j] * fact;
        }
    }
        
    return jvec;
}
