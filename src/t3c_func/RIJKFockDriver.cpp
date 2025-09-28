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
#include "ThreeCenterElectronRepulsionDriver.hpp"

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
    
    // compute atomic batches contributions to B^Q vectors
    
    CThreeCenterElectronRepulsionDriver eri_drv;
    
    for (int i = 0; i < nbatches; i++)
    {
        const auto atoms = omp::partition_atoms(natoms, i, nbatches);
        
        const auto tints = eri_drv.compute(basis, aux_basis, molecule, atoms);
     
        for (size_t j = 0; j < lindices.size(); j++)
        {
            auto ptr_bq_vec = _bq_vectors.data(j);
            
            for (const auto [gidx, lidx] : tints.mask_indices())
            {
                if (const auto fact = metric.at({lindices[j], gidx}); std::fabs(fact) > 1.0e-15)
                {
                    auto ptr_tints = tints.data(lidx);
                    
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
