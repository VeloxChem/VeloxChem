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

#include "RIFockGradDriver.hpp"

#include <ranges>

#include "ThreeCenterElectronRepulsionGeomX00Driver.hpp"
#include "ThreeCenterElectronRepulsionGeom0X0Driver.hpp"
#include "TwoCenterElectronRepulsionGeomX00Driver.hpp"

auto
CRIFockGradDriver::compute(const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const int                  iatom) const -> TPoint<double>
{
    return TPoint<double>(_comp_eri_grad(basis, aux_basis, molecule, gamma, density, iatom));
}

auto
CRIFockGradDriver::compute(const CT4CScreener&        screener,
                           const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const int                  iatom,
                           const int                  ithreshold) const -> TPoint<double>
{
    return TPoint<double>(_comp_eri_grad(screener, basis, aux_basis, molecule, gamma, density, iatom, ithreshold));
}

auto
CRIFockGradDriver::direct_compute(const CT4CScreener&        screener,
                                  const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& gamma,
                                  const CMatrix&             density,
                                  const int                  iatom,
                                  const int                  ithreshold) const -> TPoint<double>
{
    auto tg_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    const auto dvec = density.flat_values();
            
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto g1_xyz = t3cb_drv.compute(screener, basis, aux_basis, molecule, gamma, dvec, iatom, ithreshold);
    
    tg_xyz[0] += g1_xyz[0]; tg_xyz[1] += g1_xyz[1]; tg_xyz[2] += g1_xyz[2];
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto g2_xyz = t3ck_drv.compute(basis, aux_basis, molecule, gamma, dvec, iatom);
    
    tg_xyz[0] += g2_xyz[0]; tg_xyz[1] += g2_xyz[1]; tg_xyz[2] += g2_xyz[2];
    
    const auto t2c_drv = CTwoCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto g3_xyz = t2c_drv.compute(aux_basis, molecule, gamma, iatom);
    
    tg_xyz[0] -= g3_xyz[0]; tg_xyz[1] -= g3_xyz[1]; tg_xyz[2] -= g3_xyz[2];
    
    return TPoint<double>(tg_xyz);
}

auto
CRIFockGradDriver::direct_compute(const CT4CScreener&        screener,
                                  const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& bra_gamma,
                                  const std::vector<double>& ket_gamma,
                                  const CMatrix&             bra_density,
                                  const CMatrix&             ket_density,
                                  const int                  iatom,
                                  const int                  ithreshold) const -> TPoint<double>
{
    auto tg_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    const auto bra_dvec = bra_density.flat_values();
    
    const auto ket_dvec = ket_density.flat_values();
            
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto g1_xyz = t3cb_drv.compute(screener, basis, aux_basis, molecule,
                                         ket_gamma, bra_dvec, iatom, ithreshold);
    
    tg_xyz[0] += g1_xyz[0]; tg_xyz[1] += g1_xyz[1]; tg_xyz[2] += g1_xyz[2];
    
    const auto g2_xyz = t3cb_drv.compute(screener, basis, aux_basis, molecule,
                                         bra_gamma, ket_dvec, iatom, ithreshold);
    
    tg_xyz[0] += g2_xyz[0]; tg_xyz[1] += g2_xyz[1]; tg_xyz[2] += g2_xyz[2];
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto g3_xyz = t3ck_drv.compute(basis, aux_basis, molecule, ket_gamma, bra_dvec, iatom);
    
    tg_xyz[0] += g3_xyz[0]; tg_xyz[1] += g3_xyz[1]; tg_xyz[2] += g3_xyz[2];
    
    const auto g4_xyz = t3ck_drv.compute(basis, aux_basis, molecule, bra_gamma, ket_dvec, iatom);
    
    tg_xyz[0] += g4_xyz[0]; tg_xyz[1] += g4_xyz[1]; tg_xyz[2] += g4_xyz[2];
    
    tg_xyz[0] *= 0.25; tg_xyz[1] *= 0.25;  tg_xyz[2] *= 0.25; 
    
    const auto t2c_drv = CTwoCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto g5_xyz = t2c_drv.compute(aux_basis, molecule, bra_gamma, ket_gamma, iatom);
    
    tg_xyz[0] -= g5_xyz[0]; tg_xyz[1] -= g5_xyz[1]; tg_xyz[2] -= g5_xyz[2];
    
    return TPoint<double>(tg_xyz);
}

auto
CRIFockGradDriver::compute(const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const std::vector<int>     atoms) const -> std::vector<TPoint<double>>
{
    std::vector<TPoint<double>> grads(atoms.size(), TPoint<double>({0.0, 0.0, 0.0}));
    
    // prepare pointers for OMP parallel region

    auto ptr_gamma = &gamma;

    auto ptr_density = &density;

    auto ptr_molecule = &molecule;
    
    auto ptr_basis = &basis;
    
    auto ptr_aux_basis = &aux_basis;
    
    auto ptr_atoms = atoms.data();
    
    auto ptr_grads = grads.data();

    // execute OMP tasks with static scheduling
    
    const auto natoms = atoms.size();

#pragma omp parallel shared(ptr_gamma, ptr_density, ptr_molecule, ptr_basis, ptr_aux_basis, ptr_grads, ptr_atoms, natoms)
    {
#pragma omp single nowait
        {
            std::ranges::for_each(std::views::iota(size_t{0}, natoms), [&] (const auto index) {
#pragma omp task firstprivate(index)
                {
                    ptr_grads[index] = TPoint<double>(_comp_eri_grad(*ptr_basis, *ptr_aux_basis, *ptr_molecule, *ptr_gamma, *ptr_density, ptr_atoms[index])); 
                }
            });
        }
    }
    
    return grads;
}

auto
CRIFockGradDriver::compute(const CT4CScreener&        screener,
                           const CMolecularBasis&     basis,
                           const CMolecularBasis&     aux_basis,
                           const CMolecule&           molecule,
                           const std::vector<double>& gamma,
                           const CMatrix&             density,
                           const std::vector<int>     atoms,
                           const int                  ithreshold) const -> std::vector<TPoint<double>>
{
    std::vector<TPoint<double>> grads(atoms.size(), TPoint<double>({0.0, 0.0, 0.0}));
    
    // prepare pointers for OMP parallel region

    auto ptr_gamma = &gamma;

    auto ptr_density = &density;

    auto ptr_molecule = &molecule;
    
    auto ptr_basis = &basis;
    
    auto ptr_aux_basis = &aux_basis;
    
    auto ptr_atoms = atoms.data();
    
    auto ptr_grads = grads.data();
    
    auto ptr_screener = &screener;

    // execute OMP tasks with static scheduling
    
    const auto natoms = atoms.size();

#pragma omp parallel shared(ptr_screener, ptr_gamma, ptr_density, ptr_molecule, ptr_basis, ptr_aux_basis, ptr_grads, ptr_atoms, natoms, ithreshold)
    {
#pragma omp single nowait
        {
            std::ranges::for_each(std::views::iota(size_t{0}, natoms), [&] (const auto index) {
#pragma omp task firstprivate(index)
                {
                    ptr_grads[index] = TPoint<double>(_comp_eri_grad(*ptr_screener, *ptr_basis, *ptr_aux_basis, *ptr_molecule, *ptr_gamma, *ptr_density, ptr_atoms[index], ithreshold));
                }
            });
        }
    }
    
    return grads;
}

auto
CRIFockGradDriver::_comp_eri_grad(const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& gamma,
                                  const CMatrix&             density,
                                  const int                  iatom) const -> std::array<double, 3>
{
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    // set up flat density
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    // set up accumulation for aux. basis derivatives contribution
    
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gints = t3cb_drv.compute(basis, aux_basis, molecule, iatom);
    
    auto nelems = dvec.size();

    std::vector<double> gvec(nelems, 0.0);
    
    auto mask_indices = gints.mask_indices();

    for (size_t i = 0; i < gints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
     
        for (const auto [gidx, lidx] : mask_indices)
        {
            const double fact = gamma[gidx];

            auto eri_ptr = gints.data(i * mask_indices.size() + lidx);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < nelems; j++)
        {
            gsum += dvec_ptr[j] * gvec_ptr[j];
        }
        
        g_xyz[i] = 4.0 * gsum;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto tints = t3ck_drv.compute(basis, aux_basis, molecule, iatom);
    
    const auto indices = tints.indices();
    
    auto width = tints.width();
    
    mask_indices = tints.mask_indices();
    
    nelems = mask_indices.size() * width;

    gvec = std::vector<double>(nelems, 0.0);
    
    for (size_t i = 0; i < tints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
        
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[j];
            
            auto eri_ptr = tints.data(i * indices.size() + j);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (const auto [gidx, lidx] : mask_indices)
        {
            for (size_t j = 0; j < width; j++)
            {
                if (gidx < j)
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else if (gidx == j)
                {
                    gsum += 2.0 * gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(j, gidx, width)];
                }
            }
        }
        
        g_xyz[i] += 4.0 * gsum;
    }
    
    // set up accumulation of basis derivatives contribution
    
    const auto labels = std::vector<std::string>({"X", "Y", "Z"});
    
    const auto t2c_drv = CTwoCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gmats = t2c_drv.compute(aux_basis, molecule, iatom);
    
    width = indices.size();
    
    for (size_t i = 0; i < 3; i++)
    {
        const auto gmat = gmats.matrix(labels[i])->full_matrix();
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < width; j++)
        {
            for (size_t k = 0; k < width; k++)
            {
                gsum += gamma[j] * gamma[k] * (gmat.at({j, k}) + gmat.at({k, j}));
            }
        }
        
        g_xyz[i] -=  2.0 * gsum;
    }
    
    return g_xyz;
}

auto
CRIFockGradDriver::_comp_eri_grad(const CT4CScreener&        screener,
                                  const CMolecularBasis&     basis,
                                  const CMolecularBasis&     aux_basis,
                                  const CMolecule&           molecule,
                                  const std::vector<double>& gamma,
                                  const CMatrix&             density,
                                  const int                  iatom,
                                  const int                  ithreshold) const -> std::array<double, 3>
{
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    // set up flat density
    
    const auto dvec = density.flat_values();
    
    auto dvec_ptr = dvec.data();
    
    // set up accumulation for aux. basis derivatives contribution
    
    const auto t3cb_drv = CThreeCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gints = t3cb_drv.compute(screener, basis, aux_basis, molecule, iatom, ithreshold);
    
    auto nelems = dvec.size();

    std::vector<double> gvec(nelems, 0.0);
    
    auto mask_indices = gints.mask_indices();

    for (size_t i = 0; i < gints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
     
        for (const auto [gidx, lidx] : mask_indices)
        {
            const double fact = gamma[gidx];

            auto eri_ptr = gints.data(i * mask_indices.size() + lidx);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < nelems; j++)
        {
            gsum += dvec_ptr[j] * gvec_ptr[j];
        }
        
        g_xyz[i] = 4.0 * gsum;
    }
    
    //std::cout << " First : " << g_xyz[0] << " , " << g_xyz[1] << " , " << g_xyz[2] << std::endl;
    
    // set up accumulation of basis derivatives contribution
    
    const auto t3ck_drv = CThreeCenterElectronRepulsionGeom0X0Driver<1>();
    
    const auto tints = t3ck_drv.compute(basis, aux_basis, molecule, iatom);
    
    const auto indices = tints.indices();
    
    auto width = tints.width();
    
    mask_indices = tints.mask_indices();
    
    nelems = mask_indices.size() * width;

    gvec = std::vector<double>(nelems, 0.0);
    
    for (size_t i = 0; i < tints.aux_blocks(); i++)
    {
        std::ranges::fill(gvec, 0.0);
        
        auto gvec_ptr = gvec.data();
        
        for (size_t j = 0; j < indices.size(); j++)
        {
            const double fact = gamma[j];
            
            auto eri_ptr = tints.data(i * indices.size() + j);
        
            #pragma omp simd
            for (size_t k = 0; k < nelems; k++)
            {
                gvec_ptr[k] += fact * eri_ptr[k];
            }
        }
        
        double gsum = 0.0;
        
        for (const auto [gidx, lidx] : mask_indices)
        {
            for (size_t j = 0; j < width; j++)
            {
                if (gidx < j)
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else if (gidx == j)
                {
                    gsum += 2.0 * gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(gidx, j, width)];
                }
                else
                {
                    gsum += gvec_ptr[lidx * width + j] * dvec_ptr[mathfunc::uplo_rm_index(j, gidx, width)];
                }
            }
        }
        
        g_xyz[i] += 4.0 * gsum;
    }
    
    //std::cout << " Second : " << g_xyz[0] << " , " << g_xyz[1] << " , " << g_xyz[2] << std::endl;
    
    // set up accumulation of basis derivatives contribution
    
    const auto labels = std::vector<std::string>({"X", "Y", "Z"});
    
    const auto t2c_drv = CTwoCenterElectronRepulsionGeomX00Driver<1>();
    
    const auto gmats = t2c_drv.compute(aux_basis, molecule, iatom);
    
    width = indices.size();
    
    for (size_t i = 0; i < 3; i++)
    {
        const auto gmat = gmats.matrix(labels[i])->full_matrix();
        
        double gsum = 0.0;
        
        for (size_t j = 0; j < width; j++)
        {
            for (size_t k = 0; k < width; k++)
            {
                gsum += gamma[j] * gamma[k] * (gmat.at({j, k}) + gmat.at({k, j}));
            }
        }
        
        g_xyz[i] -=  2.0 * gsum;
    }
    
    //std::cout << " Total : " << g_xyz[0] << " , " << g_xyz[1] << " , " << g_xyz[2] << std::endl;
    
    return g_xyz;
}
