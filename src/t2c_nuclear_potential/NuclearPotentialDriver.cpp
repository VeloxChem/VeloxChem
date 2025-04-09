#include "NuclearPotentialDriver.hpp"

#include <algorithm>
#include <iostream>
#include <ranges>
#include <utility>

#include "GtoFunc.hpp"
#include "MatrixFunc.hpp"
#include "NuclearPotentialFunc.hpp"
#include "OpenMPFunc.hpp"
#include "DenseMatrixDistributor.hpp"
#include "T2CDistributor.hpp"
#include "TensorComponents.hpp"
#include "NuclearPotentialGridFunc.hpp"

auto
CNuclearPotentialDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix
{
    auto charges = molecule.charges();

    auto coordinates = molecule.coordinates("au");

    return compute(charges, coordinates, basis, molecule);
}

auto
CNuclearPotentialDriver::compute(const std::vector<double>&         charges,
                                 const std::vector<TPoint<double>>& coordinates,
                                 const CMolecularBasis&             basis,
                                 const CMolecule&                   molecule) const -> CMatrix
{
    // set up nuclear potential matrix

    auto npot_mat = matfunc::make_matrix(basis, mat_t::symmetric);

    npot_mat.zero();

    // prepare pointers for OMP parallel region

    auto ptr_charges = &charges;

    auto ptr_coordinates = &coordinates;

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_npot_mat = &npot_mat;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_charges, ptr_coordinates, ptr_basis, ptr_molecule, ptr_npot_mat)
    {
#pragma omp single nowait
        {
            const auto gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);

            const auto tasks = omp::make_work_tasks(gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto bra_gtos    = gto_blocks[task[0]];
                auto ket_gtos    = gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
                bool bkequal     = (task[0] == task[1]) && (task[2] == task[4]) && (task[3] == task[5]);
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal)
                {
                    CT2CDistributor<CMatrix> distributor(ptr_npot_mat, *ptr_coordinates, *ptr_charges);
                    npotfunc::compute(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                }
            });
        }
    }

    return npot_mat;
}

auto
CNuclearPotentialDriver::compute(CDenseMatrix&                   gmatrix,
                                 const std::vector<CGtoBlock>&   gto_blocks,
                                 const CDenseMatrix&             fmatrix,
                                 const std::map<size_t, size_t>& ao_mask,
                                 const size_t                    gindex,
                                 const double                    gpoint_x,
                                 const double                    gpoint_y,
                                 const double                    gpoint_z,
                                 const double                    gpoint_w) const -> void
{
    std::vector<double> charges({1.0, });
    
    std::vector<TPoint<double>> coordinates({TPoint<double>({gpoint_x, gpoint_y, gpoint_z}),});
    
    if (const auto nblocks = gto_blocks.size(); nblocks > 0)
    {
        for (size_t i = 0; i < nblocks; i++)
        {
            auto bra_indices = std::pair<size_t, size_t>{size_t{0}, gto_blocks[i].number_of_basis_functions()};
            
            CDenseMatrixDistributor distributor_ii(&gmatrix, coordinates, charges, &fmatrix, ao_mask, gpoint_w, gindex);
            
            npotfunc::compute(distributor_ii, gto_blocks[i], gto_blocks[i], bra_indices, bra_indices, true);
            
            for (size_t j = i + 1; j < nblocks; j++)
            {
                auto ket_indices = std::pair<size_t, size_t>{size_t{0}, gto_blocks[j].number_of_basis_functions()};
                
                CDenseMatrixDistributor distributor_ij(&gmatrix, coordinates, charges, &fmatrix, ao_mask, gpoint_w, gindex);
                
                npotfunc::compute(distributor_ij, gto_blocks[i], gto_blocks[j], bra_indices, ket_indices, false);
            }
        }
    }
}

auto
CNuclearPotentialDriver::compute(const std::vector<CGtoBlock>&   gto_blocks,
                                 const CDenseMatrix&             fmatrix,
                                 const std::map<size_t, size_t>& ao_mask,
                                 const std::vector<double>&      gcoords_x,
                                 const std::vector<double>&      gcoords_y,
                                 const std::vector<double>&      gcoords_z,
                                 const std::vector<double>&      gweights) const -> CDenseMatrix
{
    const auto npoints = static_cast<int>(gweights.size());
    
    const auto naos = static_cast<int>(ao_mask.size());
    
    CDenseMatrix mat_g(npoints, naos);
    
    mat_g.zero();
    
    if (const auto nblocks = gto_blocks.size(); nblocks > 0)
    {
        for (size_t i = 0; i < nblocks; i++)
        {
            const auto nbra_gtos = gto_blocks[i].number_of_basis_functions();
            
            const auto bra_mom = gto_blocks[i].angular_momentum();
            
            const auto bra_indices = gto_blocks[i].orbital_indices();
            
            // const auto nbra_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{bra_mom, });
            
            const auto nbra_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_mom, });
            
            for (size_t j = i; j < nblocks; j++)
            {
                const auto nket_gtos = gto_blocks[j].number_of_basis_functions();
                
                const auto ket_mom = gto_blocks[j].angular_momentum();
                
                const auto ket_indices = gto_blocks[j].orbital_indices();
                
                // const auto nket_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{ket_mom, });
                
                const auto nket_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_mom, });
                
                // allocate local buffers
                
                const std::array<size_t, 4> cdims{0, 0, _get_buffer_rows(bra_mom, ket_mom), static_cast<size_t>(npoints)};
                
                auto cbuffer = CSubMatrix(cdims);
                
                const std::array<size_t, 4> sdims{0, 0, static_cast<size_t>(nbra_spher_comps * nket_spher_comps), static_cast<size_t>(npoints)};
                
                auto sbuffer = CSubMatrix(sdims);
                
                // loop over GTO pairs
                
                for (int k = 0; k < nbra_gtos; k++)
                {
                    const int lstart = (i == j) ? k : 0;
                    
                    for (int l = lstart; l < nket_gtos; l++)
                    {
                        sbuffer.zero();
                        
                        cbuffer.zero();
                        
                        npotfunc::compute(sbuffer, cbuffer, gcoords_x, gcoords_y, gcoords_z, gweights, gto_blocks[i], gto_blocks[j], k, l);
                        
//                        std::cout << " *** Block *** (" << k << "," << l << ")" << std::endl;
//                        
//                        for (size_t m = 0; m < npoints; m++)
//                        {
//                            std::cout << " m = " << m << " val = " << cbuffer.at({6, m}) << std::endl; 
//                        }
                        
                        // FIX ME: Add distributor here....
                        
                        if ((bra_mom + ket_mom) == 0)
                        {
                            t2cfunc::distribute(mat_g, cbuffer, 6, fmatrix, gweights, ao_mask, bra_indices, ket_indices, bra_mom, ket_mom, k, l, i == j);
                        }
                        else
                        {
                            t2cfunc::distribute(mat_g, sbuffer, 0, fmatrix, gweights, ao_mask, bra_indices, ket_indices, bra_mom, ket_mom, k, l, i == j);
                        }
                    }
                }
            }
        }
    }
    
    return mat_g;
}

auto
CNuclearPotentialDriver::_get_buffer_rows(const int bra_angmom,
                                          const int ket_angmom) const -> size_t
{
    if ((bra_angmom == 0) && (ket_angmom == 0)) return 7;
    
    return 0;
}
