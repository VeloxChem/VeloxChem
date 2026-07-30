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

#include "CorePotentialDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "MatrixFunc.hpp"
#include "ProjectedCorePotentialDriver.hpp"
#include "LocalCorePotentialDriver.hpp"

auto
CCorePotentialDriver::compute(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrix
{
    // set up ECP matrix

    auto ecp_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    ecp_mat.zero();

    // local core potential part
    
    CLocalCorePotentialDriver loc_drv;
    
    const auto loc_mat = loc_drv.compute(basis, molecule, atom_potential.get_local_potential(), 0);
    
    ecp_mat = ecp_mat + loc_mat;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialDriver proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mat = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i], 0);
            
            ecp_mat = ecp_mat + proj_mat;
        }
    }
    
    return ecp_mat;
}

auto
CCorePotentialDriver::compute(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int>& atoms) const -> CMatrix
{
    // set up ECP matrix

    auto ecp_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    ecp_mat.zero();
    
    // run over list of atoms
    
    for (const auto atom : atoms)
    {
        const auto loc_ecp = basis.get_ecp_potential(atom);
        
        auto loc_mat = compute(basis, molecule.shift_origin(atom), loc_ecp);
        
        ecp_mat = ecp_mat + loc_mat;
    }
    
    return ecp_mat;
}
