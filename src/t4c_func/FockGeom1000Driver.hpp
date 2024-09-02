#ifndef FockGeom1000Driver_hpp
#define FockGeom1000Driver_hpp

#include <cstdint>
#include <vector>

/**
 Class CFockGeom1000Driver provides methods for computing Fock matrices
 using four center electron repulsion integrals.

 @author Z. Rinkevicius
 */
class CFockGeom1000Driver
{
   public:
    /**
     Creates a Fock matrices  driver.
     */
    CFockGeom1000Driver() = default;

    /**
     Computes Fock matrix for given density, basis and molecule.
     */
    template <class T, class U>
    auto
    compute(T* distributor, const CMolecularBasis& basis, const CMolecule& molecule) const -> void
    {
        for (auto& [atoms, mats] : distributor->get_map(0))
        {
            for (auto& [tints, ctasks] : distributor->integral_task(atoms))
            {
                const auto [a, b, c, d] = tints;

                const auto bra_gto_pair_blocks = gtofunc::make_gto_pair_blocks(gtofunc::makeGtoBlocks(basis,
                                                                                                      molecule,
                                                                                                      {
                                                                                                          a,
                                                                                                      }),
                                                                               gtofunc::makeGtoBlocks(basis, molecule));

                const auto ket_gto_pair_blocks = gtofunc::make_gto_pair_blocks(gtofunc::makeGtoBlocks(basis,
                                                                                                      molecule,
                                                                                                      {
                                                                                                          c,
                                                                                                      }),
                                                                               gtofunc::makeGtoBlocks(basis, molecule));

                const auto work_groups = omp::makeWorkGroup(bra_gto_pair_blocks, ket_gto_pair_blocks);

                // prepare pointers for OMP parallel region

                auto ptr_bra_gto_pair_blocks = bra_gto_pair_blocks.data();

                auto ptr_ket_gto_pair_blocks = ket_gto_pair_blocks.data();

                auto ptr_work_groups = work_groups.data();

                // execute OMP tasks with static scheduling

                omp::setStaticScheduler();
                const auto ntasks = work_groups.size();

#pragma omp parallel shared(ntasks, ptr_bra_gto_pair_blocks, ptr_ket_gto_pair_blocks, ptr_work_groups, )
                {
#pragma omp single nowait
                    {
                        for (size_t i = 0; i < ntasks; i++)
                        {
#pragma omp task firstprivate(i)
                            {
                                for (const auto& task : ptr_work_groups[i])
                                {
                                    const auto bra_gto_pair_block = ptr_bra_gto_pair_blocks[task[0]];  // task [2] not used

                                    const auto ket_gto_pair_block = ptr_ket_gto_pair_blocks[task[1]];

                                    eri_exp_g1010::compute(distributor, bra_gto_pair_block, ket_gto_pair_block, task[4], task[5]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif /* FockGeom1000Driver_hpp */
