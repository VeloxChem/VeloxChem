#ifndef ThreeCenterElectronRepulsionGeomX00Driver_hpp
#define ThreeCenterElectronRepulsionGeomX00Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

/// @brief Class  CThreeCenterElectronRepulsionGeomX00Driver provides methods for computing arbitrary order three-center
/// electron repulsion integral derivatives with respect bra side.
template <int N>
class CThreeCenterElectronRepulsionGeomX00Driver
{
   public:
    /// @brief Creates an electron repulsion derivative integrals driver.
    CThreeCenterElectronRepulsionGeomX00Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap derivative integrals driver to be copied.
    CThreeCenterElectronRepulsionGeomX00Driver(const CThreeCenterElectronRepulsionGeomX00Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap derivative integrals driver  to be moved.
    CThreeCenterElectronRepulsionGeomX00Driver(CThreeCenterElectronRepulsionGeomX00Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterElectronRepulsionGeomX00Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap derivative integrals driver to be copy assigned.
    /// @return The assigned overlap derivative integrals driver.
    auto operator=(const CThreeCenterElectronRepulsionGeomX00Driver &other) -> CThreeCenterElectronRepulsionGeomX00Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap derivative integrals driver to be move assigned.
    /// @return The assigned overlap derivative integrals driver .
    auto operator=(CThreeCenterElectronRepulsionGeomX00Driver &&other) noexcept -> CThreeCenterElectronRepulsionGeomX00Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap derivative integrals driver  to be compared.
    /// @return True if overlap derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterElectronRepulsionGeomX00Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap derivative integrals driver to be compared.
    /// @return True if overlap derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterElectronRepulsionGeomX00Driver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom) const -> TPoint<double>;
};

template <int N>
auto
CThreeCenterElectronRepulsionGeomX00Driver<N>::compute(const CMolecularBasis &basis,
                                                       const CMolecularBasis &aux_basis,
                                                       const CMolecule       &molecule,
                                                       const int             iatom) const -> TPoint<double>
{
    // set up gradient values

    auto atom_grad = TPoint<double>({0.0, 0.0, 0.0});
    
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule, {iatom, });
    
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    size_t block_start = 0;
    
    for (const auto& bra_gto_block : bra_gto_blocks)
    {
        for (const auto& gto_pair_block : gto_pair_blocks)
        {
            
        }
        
        // std::cout << "Task : " << task[0] << " , " << task[1] << " , " << task[2] << std::endl;
    }

//    // prepare pointers for OMP parallel region
//
//    auto ptr_basis = &basis;
//
//    auto ptr_molecule = &molecule;
//
//    auto ptr_eri_mats = &eri_mats;
//
//    // execute OMP tasks with static scheduling
//
//    omp::set_static_scheduler();
//
//#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_eri_mats, iatom)
//    {
//#pragma omp single nowait
//        {
//            const auto bra_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule, {iatom, });
//            
//            const auto ket_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);
//
//            const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);
//
//            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
//                auto bra_gtos    = bra_gto_blocks[task[0]];
//                auto ket_gtos    = ket_gto_blocks[task[1]];
//                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
//                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
//#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices)
//                {
//                    CT2CDistributor<CMatrices> distributor(ptr_eri_mats);
//                    if constexpr (N == 1)
//                    {
//                        t2cerifunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
//                    }
//                }
//            });
//        }
//    }

    return atom_grad;
}

#endif /* ThreeCenterElectronRepulsionGeomX00Driver_hpp */
