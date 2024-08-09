#ifndef ExpGeom1000Driver_hpp
#define ExpGeom1000Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "OpenMPFunc.hpp"
#include "FockType.hpp"
#include "Matrix.hpp"
#include "Matrices.hpp"
#include "FockMatrix.hpp"
#include "FockMatrices.hpp"
#include "T4CScreener.hpp"
#include "TwoElExp.hpp"
#include "GeoIntArchetypeTask.hpp"
#include "EriGeom1000Func.hpp"

/**
 Class CExpGeom1000Driver provides methods for computing Fock matrices
 using four center electron repulsion integrals..

 @author Z. Rinkevicius
 */
class CExpGeom1000Driver
{
   public:
    /**
     Creates a Fock matrices  driver.
     */
    CExpGeom1000Driver() = default;

    /**
     Computes Fock matrix for given density, basis and molecule.
     */
    template <class T, class U>
    auto inline
    compute(      T* distributor
            const CMolecularBasis& basis,
            const CMolecule&       molecule,
            const U*               screener) const -> void
    {
        // construct ket side data

        const auto ket_gto_pair_blocks = screener->get_gto_pair_blocks();

        //
        //  Loop over atom indices in arch_task
        //

        for (auto& [atoms, mats] : results[0].getMap())
        {
            for (auto& [tints, ctasks] : task.getIntegralTask(atoms))
            {
                const auto [a, b, c, d ] = tints;

                const auto bra_gto_pair_blocks = gtofunc::make_gto_pair_blocks(gtofunc::makeGtoBlocks(basis, molecule, {a, }),
                                                                               gtofunc::makeGtoBlocks(basis, molecule));

                // MR: I think the "15" is screening
                const auto work_groups = omp::make_work_group(bra_gto_pair_blocks, ket_gto_pair_blocks, 15);

                // prepare pointers for OMP parallel region

                auto ptr_bra_gto_pair_blocks = bra_gto_pair_blocks.data();
                auto ptr_ket_gto_pair_blocks = ket_gto_pair_blocks.data();
                auto ptr_work_groups = work_groups.data();

                // execute OMP tasks with static scheduling

                omp::setStaticScheduler();

                const auto ntasks = work_groups.size();

#pragma omp parallel shared(ntasks, ptr_bra_gto_pair_blocks, ptr_ket_gto_pair_blocks, ptr_work_groups)
                {
#pragma omp single nowait
                    {
                        for (size_t i = 0; i < ntasks; i++)
                        {
#pragma omp task firstprivate(i)
                            {
                                for (const auto& task : ptr_work_groups[i])
                                {
                                    const auto bra_gto_pair_block = ptr_bra_gto_pair_blocks[task[0]];

                                    const auto ket_gto_pair_block = ptr_ket_gto_pair_blocks[task[1]].get_gto_pair_block(task[3]);

                                    eri_g1000::compute(distributor, bra_gto_pair_block, ket_gto_pair_block, {task[4], task[5]}, {task[6], task[7]});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif /* ExpGeom1000Driver_hpp */