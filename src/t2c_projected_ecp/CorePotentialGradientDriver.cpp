//
//  CorePotentialGradientDriver.cpp
//  VeloxChem
//
//  Created by Zilvinas Rinkevicius on 2026-03-11.
//

#include "CorePotentialGradientDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "MatrixFunc.hpp"
#include "MatricesFunc.hpp"
#include "ProjectedCorePotentialGeomX00Driver.hpp"
#include "ProjectedCorePotentialGeom0X0Driver.hpp"
#include "LocalCorePotentialGeomX00Driver.hpp"
#include "LocalCorePotentialGeom0X0Driver.hpp"


auto
CCorePotentialGradientDriver::compute_bra_grad(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom) const -> CMatrices
{
    // set up ECP gradient matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{1}, basis, mat_t::general);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeomX00Driver<1> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential(), iatom);
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeomX00Driver<1> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i], iatom);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialGradientDriver::compute_bra_grad(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int>& atoms, const int iatom) const -> CMatrices
{
    // set up ECP gradient matrices

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{1}, basis, mat_t::general);

    ecp_mats.zero();
    
    // run over list of atoms
    
    for (const auto atom : atoms)
    {
        const auto loc_ecp = basis.get_ecp_potential(atom);
        
        auto loc_mats = compute_bra_grad(basis, molecule.shift_origin(atom), loc_ecp, iatom);
        
        ecp_mats = ecp_mats + loc_mats;
    }
    
    return ecp_mats;
}

auto
CCorePotentialGradientDriver::compute_pot_grad(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrices
{
    // set up ECP gradient matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{1}, basis, mat_t::symmetric);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeom0X0Driver<1> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential());
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeom0X0Driver<1> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i]);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialGradientDriver::compute_pot_grad(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices
{
    const auto loc_ecp = basis.get_ecp_potential(iatom);
    
    return compute_pot_grad(basis, molecule.shift_origin(iatom), loc_ecp);
}
