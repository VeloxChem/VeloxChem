#ifndef ThreeCenterElectronRepulsionGeomX00Driver_hpp
#define ThreeCenterElectronRepulsionGeomX00Driver_hpp

#include <vector>
#include <ranges>
#include <array>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "TensorComponents.hpp"
#include "T3CUtils.hpp"
#include "T3CGeomX00Distributor.hpp"
#include "T3CGeom100SumDistributor.hpp"
#include "ThreeCenterElectronRepulsionGeom100Func.hpp"
#include "T4CScreener.hpp"
#include "OpenMPFunc.hpp"

#include <iostream>

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
    /// @return The electron repulsion tensor.
    auto compute(const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom) const -> CT3FlatBuffer<double>;
    
    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param screener The screener with basis function pairs data.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The electron repulsion tensor.
    auto compute(const CT4CScreener    &screener,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom,
                 const int             ithreshold) const -> CT3FlatBuffer<double>;
    
    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param screener The screener with basis function pairs data.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param density The flat density matrix values.
    /// @param iatom The index of atom.
    /// @return The electron repulsion gradient.
    auto compute(const CT4CScreener&        screener,
                 const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const std::vector<double>& density,
                 const int                  iatom,
                 const int                  ithreshold) const -> std::array<double, 3>;
};

template <int N>
auto
CThreeCenterElectronRepulsionGeomX00Driver<N>::compute(const CMolecularBasis &basis,
                                                       const CMolecularBasis &aux_basis,
                                                       const CMolecule       &molecule,
                                                       const int             iatom) const -> CT3FlatBuffer<double>
{
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule, {iatom, });
    
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    // set up composite flat tensor for integrals
    
    const auto mask_indices = t3cfunc::mask_indices(bra_gto_blocks);
    
    CT3FlatBuffer<double> buffer(mask_indices, basis.dimensions_of_basis(), 3);
    
    // set up distributor
    
    CT3CGeomX00Distributor distributor(&buffer);
    
    // main compute loop
    
    std::ranges::for_each(bra_gto_blocks, [&](const auto& gblock) {
        auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
        std::ranges::for_each(gto_pair_blocks, [&](const auto& gp_pairs) {
            if constexpr (N == 1)
            {
                t3cerifunc::compute_geom_100(distributor, gblock, gp_pairs, bra_range);
            }
        });
    });
    
    return buffer;
}

template <int N>
auto
CThreeCenterElectronRepulsionGeomX00Driver<N>::compute(const CT4CScreener    &screener,
                                                       const CMolecularBasis &basis,
                                                       const CMolecularBasis &aux_basis,
                                                       const CMolecule       &molecule,
                                                       const int             iatom,
                                                       const int             ithreshold) const -> CT3FlatBuffer<double>
{
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule, {iatom, });
    
    // set up composite flat tensor for integrals
    
    const auto mask_indices = t3cfunc::mask_indices(bra_gto_blocks);
    
    CT3FlatBuffer<double> buffer(mask_indices, basis.dimensions_of_basis(), 3);
    
    // set up distributor
    
    CT3CGeomX00Distributor distributor(&buffer);
    
    // main compute loop
    
    std::ranges::for_each(bra_gto_blocks, [&](const auto& gblock) {
        auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
        for (const auto& gto_pair_block : screener.gto_pair_blocks())
        {
            for (int i = 0; i <= ithreshold; i++)
            {
                if (gto_pair_block.is_empty_gto_pair_block(i)) continue;
                
                const auto gp_pairs = gto_pair_block.gto_pair_block(i);
                
                if constexpr (N == 1)
                {
                    t3cerifunc::compute_geom_100(distributor, gblock, gp_pairs, bra_range);
                }
            }
        }
    });
    
    return buffer;
}

template <int N>
auto
CThreeCenterElectronRepulsionGeomX00Driver<N>::compute(const CT4CScreener&        screener,
                                                       const CMolecularBasis&     basis,
                                                       const CMolecularBasis&     aux_basis,
                                                       const CMolecule&           molecule,
                                                       const std::vector<double>& gamma,
                                                       const std::vector<double>& density,
                                                       const int                  iatom,
                                                       const int                  ithreshold) const -> std::array<double, 3>
{
    // set up number of rows
    
    const auto nrows = basis.dimensions_of_basis();
    
    // set up GTOs data
    
    const auto aux_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule, {iatom, });
    
    const auto gto_pair_blocks = screener.gto_pair_blocks();
    
    const auto tasks = omp::make_aux_work_tasks(aux_gto_blocks, gto_pair_blocks, ithreshold);
    
    // set up accumulation vector
    
    auto gvec_xyz = std::vector<std::array<double, 3>>(tasks.size(), {0.0, 0.0, 0.0});
    
    // prepare pointers for OMP parallel region
    
    auto ptr_gto_pair_blocks = gto_pair_blocks.data();
    
    auto ptr_aux_gto_blocks = aux_gto_blocks.data();
    
    auto ptr_gamma = gamma.data();
    
    auto ptr_density = density.data();
    
    auto ptr_tasks = &tasks;
    
    auto ptr_gvec_xyz = gvec_xyz.data();
    
    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_aux_gto_blocks, ptr_gamma, ptr_density, ptr_tasks, ptr_gvec_xyz, nrows)
    {
#pragma omp single nowait
        {
            size_t index = 0;
            std::ranges::for_each(std::ranges::reverse_view(*ptr_tasks), [&](const auto& task) {
                auto aux_idx = task[0];
                auto ket_idx = task[1];
                auto gp_idx  = task[2];
#pragma omp task firstprivate(index, aux_idx, ket_idx, gp_idx)
                {
                    auto gblock = ptr_aux_gto_blocks[aux_idx];
                    
                    auto gp_pairs = ptr_gto_pair_blocks[ket_idx].gto_pair_block(static_cast<int>(gp_idx));
                              
                    CT3CGeom100SumDistributor distributor(ptr_gvec_xyz[index].data(), ptr_gamma, ptr_density, nrows);
                    
                    auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
                    
                    if constexpr (N == 1)
                    {
                        t3cerifunc::compute_geom_100(distributor, gblock, gp_pairs, bra_range);
                    }
                }
                index++;
            });
        }
    }
    
    // sum up contributions from accumulation vector
    
    auto g_xyz = std::array<double, 3>({0.0, 0.0, 0.0});
    
    for (const auto& t_xyz : gvec_xyz)
    {
        g_xyz[0] += t_xyz[0];
        
        g_xyz[1] += t_xyz[1];
        
        g_xyz[2] += t_xyz[2];
    }
    
    return g_xyz;
}




#endif /* ThreeCenterElectronRepulsionGeomX00Driver_hpp */
