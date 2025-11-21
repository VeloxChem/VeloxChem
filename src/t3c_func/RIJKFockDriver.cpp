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

#include "RIJKFockDriver.hpp"

#include <numeric>
#include <cmath>

#include "OpenMPFunc.hpp"
#include "BatchFunc.hpp"
#include "ThreeCenterElectronRepulsionDriver.hpp"
#include "SerialDenseLinearAlgebra.hpp"

#include <iostream>

auto
CRIJKFockDriver::compute_bq_vectors(const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const CMolecularBasis&  aux_basis,
                                    const CSubMatrix&       metric,
                                    const int               rank,
                                    const int               nodes) -> void
{
    // set up basis dimensions
    
    const auto naos = basis.dimensions_of_basis();
    
    const auto naux = aux_basis.dimensions_of_basis();
    
    const auto nelems = naos * (naos + 1) / 2; 
    
    // set up active auxilary indices
    
    auto gindices = std::vector<size_t>(naux);
    
    std::iota(gindices.begin(), gindices.end(), 0);
    
    auto lindices = omp::partition_tasks(gindices, rank, nodes);
    
    // allocate B^Q vectors
    
    _bq_vectors = CT3FlatBuffer<double>(lindices, naos);
    
    // set up atomic batching
    
    const auto natoms = molecule.number_of_atoms();
    
    const auto nbatches = ((natoms % 10) == 0) ? natoms / 10 : natoms / 10 + 1;
  
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto ptr_metric = &metric;
    
    // compute atomic batches contributions to B^Q vectors
    
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    for (int i = 0; i < nbatches; i++)
    {
        const auto atoms = omp::partition_atoms(natoms, i, nbatches);
        
        const auto buffer = eri_drv.compute(basis, aux_basis, molecule, atoms);
        
        auto ptr_buffer = &buffer;
        
        const auto nlocaux = lindices.size();
     
#pragma omp parallel for shared(ptr_bq_vectors, ptr_metric, ptr_buffer) schedule(dynamic)
        for (size_t j = 0; j < nlocaux; j++)
        {
            auto ptr_bq_vec = ptr_bq_vectors->data(j);
            
            for (const auto [gidx, lidx] : ptr_buffer->mask_indices())
            {
                if (const auto fact = ptr_metric->at({gidx, lindices[j]}); std::fabs(fact) > 1.0e-15)
                {
                    auto ptr_tints = ptr_buffer->data(lidx);
                    
                    #pragma omp simd
                    for (size_t k = 0; k < nelems; k++)
                    {
                        ptr_bq_vec[k] += ptr_tints[k] * fact;
                    }
                }
            }
        }
    }
}

auto
CRIJKFockDriver::compute_screened_bq_vectors(const CT4CScreener&     screener,
                                             const CMolecule&        molecule,
                                             const CMolecularBasis&  aux_basis,
                                             const CSubMatrix&       metric,
                                             const int               ithreshold,
                                             const int               rank,
                                             const int               nodes) -> void
{
    // set up auxilary basis dimensions
    
    const auto naux = aux_basis.dimensions_of_basis();
    
    // set up active auxilary indices
    
    auto gindices = std::vector<size_t>(naux);
    
    std::iota(gindices.begin(), gindices.end(), 0);
    
    auto lindices = omp::partition_tasks(gindices, rank, nodes);
    
    // allocate B^Q vectors
    
    const auto bra_indices = omp::partition_flat_buffer(screener.gto_pair_blocks(), ithreshold);
    
    _bq_mask = omp::generate_flat_buffer_mask(screener.gto_pair_blocks(), bra_indices);
    
    _bq_vectors = CT3FlatBuffer<double>(lindices, bra_indices.back());
    
    const auto nelems = _bq_vectors.elements(); 
    
    // set up atomic batching
    
    const auto natoms = molecule.number_of_atoms();
    
    const auto nbatches = ((natoms % 10) == 0) ? natoms / 10 : natoms / 10 + 1;
    
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto ptr_metric = &metric;
    
    // compute atomic batches contributions to B^Q vectors
    
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    for (int i = 0; i < nbatches; i++)
    {
        const auto atoms = omp::partition_atoms(natoms, i, nbatches);
        
        const auto buffer = eri_drv.compute(screener, aux_basis, molecule, atoms, ithreshold);
        
        auto ptr_buffer = &buffer;
        
        const auto nlocaux = lindices.size();
     
#pragma omp parallel for shared(ptr_bq_vectors, ptr_metric, ptr_buffer) schedule(dynamic)
        for (size_t j = 0; j < nlocaux; j++)
        {
            auto ptr_bq_vec = ptr_bq_vectors->data(j);
            
            for (const auto [gidx, lidx] : ptr_buffer->mask_indices())
            {
                if (const auto fact = ptr_metric->at({gidx, lindices[j]}); std::fabs(fact) > 1.0e-15)
                {
                    auto ptr_tints = ptr_buffer->data(lidx);
                    
                    #pragma omp simd
                    for (size_t k = 0; k < nelems; k++)
                    {
                        ptr_bq_vec[k] += ptr_tints[k] * fact;
                    }
                }
            }
        }
    }
}

auto
CRIJKFockDriver::local_compute_screened_bq_vectors(const CT4CScreener&        screener,
                                                   const CMolecule&           molecule,
                                                   const CMolecularBasis&     aux_basis,
                                                   const CSubMatrix&          metric,
                                                   const int                  ithreshold,
                                                   const std::vector<size_t>& indices) -> void
{
    // allocate B^Q vectors
    
    const auto bra_indices = omp::partition_flat_buffer(screener.gto_pair_blocks(), ithreshold);
    
    _bq_mask = omp::generate_flat_buffer_mask(screener.gto_pair_blocks(), bra_indices);
    
    _bq_vectors = CT3FlatBuffer<double>(indices, bra_indices.back());
    
    const auto nelems = _bq_vectors.elements();
    
    // set up atomic batching
    
    const auto natoms = molecule.number_of_atoms();
    
    const auto nbatches = ((natoms % 10) == 0) ? natoms / 10 : natoms / 10 + 1;
    
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto ptr_metric = &metric;
    
    // compute atomic batches contributions to B^Q vectors
    
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    for (int i = 0; i < nbatches; i++)
    {
        const auto atoms = omp::partition_atoms(natoms, i, nbatches);
        
        const auto buffer = eri_drv.compute(screener, aux_basis, molecule, atoms, ithreshold);
        
        auto ptr_buffer = &buffer;
        
        const auto nlocaux = indices.size();
     
#pragma omp parallel for shared(ptr_bq_vectors, ptr_metric, ptr_buffer) schedule(dynamic)
        for (size_t j = 0; j < nlocaux; j++)
        {
            auto ptr_bq_vec = ptr_bq_vectors->data(j);
            
            for (const auto [gidx, lidx] : ptr_buffer->mask_indices())
            {
                if (const auto fact = ptr_metric->at({gidx, indices[j]}); std::fabs(fact) > 1.0e-15)
                {
                    auto ptr_tints = ptr_buffer->data(lidx);
                    
                    #pragma omp simd
                    for (size_t k = 0; k < nelems; k++)
                    {
                        ptr_bq_vec[k] += ptr_tints[k] * fact;
                    }
                }
            }
        }
    }
}

auto
CRIJKFockDriver::compute_j_fock(const CMatrix     &density,
                                const std::string &label) const -> CMatrix
{
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        CMatrix fmat(density);
    
        fmat.assign_flat_values(_comp_j_vector(_comp_m_vector(density)));
        
        // rescale Fock matrix for closed shell case
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
        
        return fmat;
    }
    else
    {
        return CMatrix();
    }
}

auto
CRIJKFockDriver::compute_screened_j_fock(const CMatrix     &density,
                                         const std::string &label) const -> CMatrix
{
    if ((label == "2jk") || (label == "2jkx") || (label == "j") || (label == "j_rs"))
    {
        CMatrix fmat(density);
        
        fmat.zero(); 
        
        const auto naos = fmat.number_of_rows();
    
        fmat.assign_values(CSubMatrix(_comp_j_vector(_comp_m_vector(density)), _bq_mask, {0, 0, naos, naos}));
        
        // rescale Fock matrix for closed shell case
        
        if ((label == "2jk") || (label == "2jkx"))
        {
            fmat.scale(2.0);
        }
        
        return fmat;
    }
    else
    {
        return CMatrix();
    }
}

auto
CRIJKFockDriver::compute_k_fock(const CMatrix &density, const CSubMatrix &molorbs) const -> CMatrix
{
    CMatrix fmat(density);
    
    // set up dimensions
    
    const auto naos = molorbs.number_of_rows();
    
    const auto nmos = molorbs.number_of_columns();
    
    const auto naux = _bq_vectors.aux_width();
    
    auto kmat = CSubMatrix({0, 0, naos, naos}, 0.0);
    
    // set up full B^Q matrices
    
    auto bqao = CSubMatrix({0, 0, naos, naos});
    
    auto bqmo = CSubMatrix({0, 0, nmos, naos});
    
    for (size_t i = 0; i < naux; i++)
    {
        _bq_vectors.unpack_data(bqao, i);
        
        bqmo.zero();
        
        sdenblas::serialMultAtB(bqmo, molorbs, bqao);
        
        sdenblas::serialMultAtB(kmat, bqmo, bqmo);
    }
    
    fmat.assign_values(kmat);
    
    return fmat;
}

auto
CRIJKFockDriver::compute_screened_k_fock(const CMatrix &density, const CSubMatrix &molorbs) const -> CMatrix
{
    CMatrix fmat(density);
    
    // set up dimensions
    
    const auto naos = molorbs.number_of_rows();
    
    const auto nmos = molorbs.number_of_columns();
    
    const auto naux = _bq_vectors.aux_width();
    
    auto kmat = CSubMatrix({0, 0, naos, naos}, 0.0);
    
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto ptr_molorbs = &molorbs;
    
    auto ptr_kmat = kmat.data();
    
    auto ptr_bq_mask = &_bq_mask;
    
    // set up batch size for OMP region
    
    const size_t bsize = 100;
    
    #pragma omp parallel shared(ptr_bq_vectors, ptr_bq_mask, ptr_molorbs, ptr_kmat, naux, naos, nmos, bsize)
    {
        #pragma omp single nowait
        {
            const auto nbatches = batch::number_of_batches(naux, bsize);
        
            for (size_t i = 0; i < nbatches; i++)
            {
                const auto bpair = batch::batch_range(i, naux, bsize);
            
                #pragma omp task firstprivate(bpair)
                {
                    // set up full B^Q matrices
                    
                    auto bqao = CSubMatrix({0, 0, naos, naos});
                    
                    auto bqmo = CSubMatrix({0, 0, nmos, naos});
                    
                    auto lmat = CSubMatrix({0, 0, naos, naos}, 0.0);
                    
                    auto lmask = *ptr_bq_mask;
                    
                    auto lmos = *ptr_molorbs;
                    
                    for (auto j = bpair.first; j < bpair.second; j++)
                    {
                        ptr_bq_vectors->reduced_unpack_data(bqao, lmask, j);
                        
                        bqmo.zero();
                        
                        sdenblas::serialMultAtB(bqmo, lmos, bqao);
                        
                        sdenblas::serialMultAtB(lmat, bqmo, bqmo);
                    }
                    
                    #pragma omp critical
                    {
                        const auto nelems = naos * naos;
                        
                        auto ptr_lmat = lmat.data();
                        
                        #pragma omp simd
                        for (size_t j = 0; j < nelems; j++)
                        {
                            ptr_kmat[j] += ptr_lmat[j];
                        }
                    }
                }
            }
        }
    }
   
    fmat.assign_values(kmat);
    
    return fmat;
}

auto
CRIJKFockDriver::compute_mo_bq_vectors(const CSubMatrix& lambda_p, const CSubMatrix& lambda_h, const size_t bstart, const size_t bend) const -> std::vector<CSubMatrix>
{
    if (const auto bdim = bend - bstart; bdim > 0)
    {
        std::vector<CSubMatrix> movecs(bend - bstart, CSubMatrix());
        
        // set up AOs, MOs indices
           
        const auto nmos = lambda_p.number_of_columns();
               
        const auto naos = lambda_p.number_of_rows();
        
        // set up pointers to OMP data
        
        auto ptr_matp = &lambda_p;
        
        auto ptr_math = &lambda_h;
        
        auto ptr_movecs = movecs.data();
        
        auto ptr_bq_vectors = &_bq_vectors;
        
        auto ptr_bq_mask = &_bq_mask;
        
        #pragma omp parallel for shared(ptr_matp, ptr_math, ptr_movecs, ptr_bq_vectors, ptr_bq_mask)
        for (size_t i = bstart; i < bend; i++)
        {
            auto bqao = CSubMatrix({0, 0, naos, naos});
            
            ptr_bq_vectors->reduced_unpack_data(bqao, *ptr_bq_mask, i);
            
            auto tmat = CSubMatrix({0, 0, naos, nmos}, 0.0);
            
            sdenblas::serialMultAB(tmat, bqao, *ptr_math);
            
            auto lmat = CSubMatrix({0, 0, nmos, nmos}, 0.0);
            
            sdenblas::serialMultAtB(lmat, *ptr_matp, tmat);
            
            ptr_movecs[i - bstart] = lmat;
        }
        
        return movecs;
    }
    else
    {
        return std::vector<CSubMatrix>();
    }
}

auto
CRIJKFockDriver::_comp_m_vector(const CMatrix &density) const -> std::vector<double>
{
    const auto ndim = _bq_vectors.aux_width();
    
    const auto nelems = _bq_vectors.elements();
    
    std::vector<double> mvec(ndim, 0.0);
    
    const auto dvec = _bq_vectors.is_reduced() ? density.reduced_flat_values(_bq_mask) : density.flat_values();
    
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto dvec_ptr = dvec.data();
    
    auto mvec_ptr = mvec.data();
    
    #pragma omp parallel for shared(ptr_bq_vectors, mvec_ptr, dvec_ptr) schedule(dynamic)
    for (size_t i = 0; i < ndim; i++)
    {
        auto tint_ptr = ptr_bq_vectors->data(i);
         
        double fsum = 0.0;
        
        #pragma omp simd reduction(+:fsum)
        for (size_t j = 0; j < nelems; j++)
        {
            fsum += dvec_ptr[j] * tint_ptr[j];
        }
        
        mvec_ptr[i] = fsum;
    }
        
    return mvec;
}

auto
CRIJKFockDriver::_comp_j_vector(const std::vector<double>& mvector) const -> std::vector<double>
{
    const auto nelems = _bq_vectors.elements();
    
    std::vector<double> jvec(nelems, 0.0);
    
    // set up pointers for OMP region
    
    auto ptr_bq_vectors = &_bq_vectors;
    
    auto mvec_ptr = mvector.data();
    
    auto jvec_ptr = jvec.data();
    
    // set up batch size for OMP region
    
    const size_t bsize = 500;
    
    #pragma omp parallel shared(ptr_bq_vectors, mvec_ptr, jvec_ptr, nelems, bsize)
    {
        #pragma omp single nowait
        {
            const auto nbatches = batch::number_of_batches(nelems, bsize);
            
            for (size_t i = 0; i < nbatches; i++)
            {
                const auto bpair = batch::batch_range(i, nelems, bsize);
                
                #pragma omp task firstprivate(bpair)
                {
                    const size_t bstart = bpair.first;
                    
                    const size_t bend = bpair.second;
                    
                    for (size_t j = 0; j < ptr_bq_vectors->aux_width(); j++)
                    {
                        auto tint_ptr = ptr_bq_vectors->data(j);
                        
                        if (const auto fact = mvec_ptr[j]; std::fabs(fact) > 1.0e-15)
                        {
                            #pragma omp simd
                            for (size_t k = bstart; k < bend; k++)
                            {
                                jvec_ptr[k] += tint_ptr[k] * fact;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return jvec;
}
