#include "FockDriver.hpp"

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "ElectronRepulsionFunc.hpp"
#include "T4CMatrixDistributor.hpp"

auto
CFockDriver::compute(const CMolecularBasis& basis,
                     const CMolecule&       molecule,
                     const CMatrix&         density,
                     const std::string&     label,
                     const double           exchange_factor,
                     const double           omega) const -> CMatrix
{
    // set up Fock matrix
    
    auto fock_mat = CMatrix(density);
    
    fock_mat.zero();
    
    // set up basis function pairs blocks
    
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    // prepare pointers for OMP parallel region
    
    auto ptr_gto_pair_blocks = &gto_pair_blocks;

    auto ptr_density = &density;
    
    auto ptr_fock = &fock_mat;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_density, ptr_fock, label, exchange_factor, omega)
    {
#pragma omp single nowait
        {
            const auto nblocks = ptr_gto_pair_blocks->size();
            
            auto ptr_gto_pairs_data = ptr_gto_pair_blocks->data();
            
            for (size_t i = 0; i < nblocks; i++)
            {
                for (size_t j = i; j < nblocks; j++)
                {
                    #pragma omp task firstprivate(i, j)
                    {
                        auto bra_gpairs = ptr_gto_pairs_data[i];
                        
                        auto ket_gpairs = ptr_gto_pairs_data[j];
                        
                        CT4CMatrixDistributor distributor(ptr_fock, ptr_density, label, exchange_factor, omega);
                        
                        distributor.set_indices(bra_gpairs, ket_gpairs); 
                        
                        auto bra_range = std::pair<size_t, size_t>(0, bra_gpairs.number_of_contracted_pairs());
                        
                        auto ket_range = std::pair<size_t, size_t>(0, ket_gpairs.number_of_contracted_pairs()); 
                        
                        erifunc::compute<CT4CMatrixDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range, i == j);
                        
                        //distributor.accumulate(bra_gpairs, ket_gpairs); 
                    }
                }
            }
        }
    }
    
    fock_mat.symmetrize(); 
    
    return fock_mat;
}
