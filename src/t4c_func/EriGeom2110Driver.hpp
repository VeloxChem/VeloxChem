#ifndef ExpGeom2110Driver_hpp
#define ExpGeom2110Driver_hpp

#include <vector>

#include "EriGeom2110Func.hpp"
#include "FockMatrices.hpp"
#include "FockMatrix.hpp"
#include "FockType.hpp"
#include "GeoIntArchetypeTask.hpp"
#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "OpenMPFunc.hpp"
#include "T4CScreener.hpp"
#include "TwoElExp.hpp"

/**
 Class CExpGeom2110Driver provides methods for computing Fock matrices
 using four center electron repulsion integrals..

 @author Z. Rinkevicius
 */
class CExpGeom2110Driver
{
   public:
    /**
     Creates a Fock matrices  driver.
     */
    CExpGeom2110Driver() = default;

    /**
     Computes Fock matrix for given density, basis and molecule.
     */
    template <class T, class U>
    auto inline compute(T* distributor const CMolecularBasis& basis, const CMolecule& molecule, const U* screener) const -> void
    {
        //
        //  Loop over atom indices in arch_task
        //

        for (auto& [atoms, mats] : results[0].getMap())
        {
            for (auto& [tints, ctasks] : task.getIntegralTask(atoms))
            {
                const auto [a, b, c, d] = tints;

                const auto bra_gto_pair_blocks = gtofunc::makeGtoPairBlocks(gtofunc::makeGtoBlocks(basis,
                                                                                                   molecule,
                                                                                                   {
                                                                                                       a,
                                                                                                   }),
                                                                            gtofunc::makeGtoBlocks(basis,
                                                                                                   molecule,
                                                                                                   {
                                                                                                       b,
                                                                                                   }));

                const auto ket_gto_pair_blocks = gtofunc::makeGtoPairBlocks(gtofunc::makeGtoBlocks(basis,
                                                                                                   molecule,
                                                                                                   {
                                                                                                       c,
                                                                                                   }),
                                                                            gtofunc::makeGtoBlocks(basis, molecule));

                // MR: I think the "15" is screening
                const auto work_groups = omp::make_work_group(bra_gto_pair_blocks, ket_gto_pair_blocks, 15);

                // prepare pointers for OMP parallel region

                auto ptr_bra_gto_pair_blocks = bra_gto_pair_blocks.data();
                auto ptr_ket_gto_pair_blocks = ket_gto_pair_blocks.data();
                auto ptr_work_groups         = work_groups.data();

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

                                    const auto ket_gto_pair_block = ptr_ket_gto_pair_blocks[task[1]];

                                    eri_g2110::compute(distributor, bra_gto_pair_block, ket_gto_pair_block, {task[4], task[5]}, {task[6], task[7]});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif /* ExpGeom2110Driver_hpp */