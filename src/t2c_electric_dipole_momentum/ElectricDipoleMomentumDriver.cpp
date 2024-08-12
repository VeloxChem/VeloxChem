#include "ElectricDipoleMomentumDriver.hpp"

#include "ElectricDipoleMomentumFunc.hpp"
#include "GtoFunc.hpp"
#include "MatricesFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T2CDistributor.hpp"

auto
CElectricDipoleMomentumDriver::compute(const CMolecularBasis& basis, const CMolecule& molecule, const TPoint<double>& origin) const -> CMatrices
{
    // set up electric dipole matrices

    auto dip_mats = matfunc::make_matrices(std::array<int, 1>{1}, basis, mat_t::symmetric);

    dip_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_dip_mats = &dip_mats;

    auto ptr_origin = &origin;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_dip_mats, ptr_origin)
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
                    const auto                 coords = std::vector<TPoint<double>>(1, *ptr_origin);
                    const auto                 data   = std::vector<double>();
                    CT2CDistributor<CMatrices> distributor(ptr_dip_mats, coords, data);
                    dipfunc::compute(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                }
            });
        }
    }

    return dip_mats;
}
