#include "CorePotentialHessianDriver.hpp"

#include "ProjectedCorePotentialGeomX00Driver.hpp"
#include "ProjectedCorePotentialGeom0X0Driver.hpp"
#include "ProjectedCorePotentialGeomXY0Driver.hpp"
#include "ProjectedCorePotentialGeomX0YDriver.hpp"
#include "LocalCorePotentialGeomX00Driver.hpp"
#include "LocalCorePotentialGeomX0YDriver.hpp"
#include "LocalCorePotentialGeom0X0Driver.hpp"
#include "LocalCorePotentialGeomXY0Driver.hpp"

auto
CCorePotentialHessianDriver::compute_geom_200(const CMolecularBasis    &basis,
                                              const CMolecule          &molecule,
                                              const CAtomCorePotential &atom_potential,
                                              const int                 iatom) const -> CMatrices
{
    // set up ECP hessian matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{2}, basis, mat_t::general);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeomX00Driver<2> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential(), iatom);
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeomX00Driver<2> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i], iatom);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_200(const CMolecularBasis   &basis,
                                              const CMolecule         &molecule,
                                              const std::vector<int>  &atoms,
                                              const int                iatom) const -> CMatrices
{
    // set up ECP gradient matrices

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{2}, basis, mat_t::general);

    ecp_mats.zero();
    
    // run over list of atoms
    
    for (const auto atom : atoms)
    {
        const auto loc_ecp = basis.get_ecp_potential(atom);
        
        auto loc_mats = compute_geom_200(basis, molecule.shift_origin(atom), loc_ecp, iatom);
        
        ecp_mats = ecp_mats + loc_mats;
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_101(const CMolecularBasis    &basis,
                                              const CMolecule          &molecule,
                                              const CAtomCorePotential &atom_potential,
                                              const int                 iatom,
                                              const int                 jatom) const -> CMatrices
{
    // set up ECP hessian matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 2>{1, 1}, basis, mat_t::general);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeomX0YDriver<1, 1> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential(), iatom, jatom);
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeomX0YDriver<1, 1> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i], iatom, jatom);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_101(const CMolecularBasis  &basis,
                                              const CMolecule        &molecule,
                                              const std::vector<int> &atoms,
                                              const int               iatom,
                                              const int               jatom) const -> CMatrices
{
    // set up ECP gradient matrices

    auto ecp_mats = matfunc::make_matrices(std::array<int, 2>{1, 1}, basis, mat_t::general);

    ecp_mats.zero();
    
    // run over list of atoms
    
    for (const auto atom : atoms)
    {
        const auto loc_ecp = basis.get_ecp_potential(atom);
        
        auto loc_mats = compute_geom_101(basis, molecule.shift_origin(atom), loc_ecp, iatom, jatom);
        
        ecp_mats = ecp_mats + loc_mats;
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_110(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential, const int iatom) const -> CMatrices
{
    // set up ECP hessian matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 2>{1, 1}, basis, mat_t::general);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeomXY0Driver<1, 1> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential(), iatom);
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeomXY0Driver<1, 1> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i], iatom);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_110(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom, const int jatom) const -> CMatrices
{
    const auto loc_ecp = basis.get_ecp_potential(jatom);
        
    return compute_geom_110(basis, molecule.shift_origin(jatom), loc_ecp, iatom);
}

auto
CCorePotentialHessianDriver::compute_geom_020(const CMolecularBasis &basis, const CMolecule &molecule, const CAtomCorePotential& atom_potential) const -> CMatrices
{
    // set up ECP gradient matrix

    auto ecp_mats = matfunc::make_matrices(std::array<int, 1>{2}, basis, mat_t::symmetric);

    ecp_mats.zero();

    // local core potential part
    
    CLocalCorePotentialGeom0X0Driver<2> loc_drv;
    
    const auto loc_mats = loc_drv.compute(basis, molecule, atom_potential.get_local_potential());
    
    ecp_mats = ecp_mats + loc_mats;
    
    // projected core potential part
    
    const auto ang_moms = atom_potential.get_angular_momentums();
    
    const auto proj_pots = atom_potential.get_projected_potentials();
    
    if (size_t nterms = ang_moms.size(); nterms > 0)
    {
        CProjectedCorePotentialGeom0X0Driver<2> proj_drv;
        
        for (size_t i = 0; i < nterms; i++)
        {
            const auto proj_mats = proj_drv.compute(basis, molecule, proj_pots[i], ang_moms[i]);
            
            ecp_mats = ecp_mats + proj_mats;
        }
    }
    
    return ecp_mats;
}

auto
CCorePotentialHessianDriver::compute_geom_020(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices
{
    const auto loc_ecp = basis.get_ecp_potential(iatom);
    
    return compute_geom_020(basis, molecule.shift_origin(iatom), loc_ecp);
}
