#ifndef ThreeCenterElectronRepulsionGeom0X0Driver_hpp
#define ThreeCenterElectronRepulsionGeom0X0Driver_hpp

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
#include "T3RectFlatBuffer.hpp"
#include "T3CGeom0X0Distributor.hpp"
#include "T3CGeom010SumDistributor.hpp"
#include "ThreeCenterElectronRepulsionGeom010Func.hpp"
#include "OpenMPFunc.hpp"

#include <iostream>

/// @brief Class  CThreeCenterElectronRepulsionGeom0X0Driver provides methods for computing arbitrary order three-center
/// electron repulsion integral derivatives with respect bra side.
template <int N>
class CThreeCenterElectronRepulsionGeom0X0Driver
{
   public:
    /// @brief Creates an electron repulsion derivative integrals driver.
    CThreeCenterElectronRepulsionGeom0X0Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap derivative integrals driver to be copied.
    CThreeCenterElectronRepulsionGeom0X0Driver(const CThreeCenterElectronRepulsionGeom0X0Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap derivative integrals driver  to be moved.
    CThreeCenterElectronRepulsionGeom0X0Driver(CThreeCenterElectronRepulsionGeom0X0Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterElectronRepulsionGeom0X0Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap derivative integrals driver to be copy assigned.
    /// @return The assigned overlap derivative integrals driver.
    auto operator=(const CThreeCenterElectronRepulsionGeom0X0Driver &other) -> CThreeCenterElectronRepulsionGeom0X0Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap derivative integrals driver to be move assigned.
    /// @return The assigned overlap derivative integrals driver .
    auto operator=(CThreeCenterElectronRepulsionGeom0X0Driver &&other) noexcept -> CThreeCenterElectronRepulsionGeom0X0Driver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap derivative integrals driver  to be compared.
    /// @return True if overlap derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterElectronRepulsionGeom0X0Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap derivative integrals driver to be compared.
    /// @return True if overlap derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterElectronRepulsionGeom0X0Driver &other) const -> bool = delete;

    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The electron repulsion integrals tensor.
    auto compute(const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom) const -> CT3RectFlatBuffer<double>;
    
    /// @brief Computes electron repulsion matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param aux_basis The auxilary molecular basis for fiting of four-center repulsion integrals.
    /// @param molecule The molecule.
    /// @param gamma The transformed Gamma vector.
    /// @param density The flat density matrix values.
    /// @param iatom The index of atom.
    /// @return The electron repulsion integrals tensor.
    auto compute(const CMolecularBasis&     basis,
                 const CMolecularBasis&     aux_basis,
                 const CMolecule&           molecule,
                 const std::vector<double>& gamma,
                 const std::vector<double>& density,
                 const int                  iatom) const -> std::array<double, 3>;
};

template <int N>
auto
CThreeCenterElectronRepulsionGeom0X0Driver<N>::compute(const CMolecularBasis &basis,
                                                       const CMolecularBasis &aux_basis,
                                                       const CMolecule       &molecule,
                                                       const int             iatom) const -> CT3RectFlatBuffer<double>
{
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule);
    
    const auto ketr_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, {iatom, });
    
    const auto ketf_gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(ketr_gto_blocks, ketf_gto_blocks);
    
    // set up composite flat tensor for integrals
    
    std::vector<size_t> aux_indices(aux_basis.dimensions_of_basis());
    
    std::iota(aux_indices.begin(), aux_indices.end(), size_t{0});
    
    const auto red_indices = t3cfunc::unique_indices(ketr_gto_blocks);
    
    std::vector<size_t> full_indices(basis.dimensions_of_basis());
    
    std::iota(full_indices.begin(), full_indices.end(), size_t{0});
    
    const auto mask_indices = t3cfunc::mask_indices(ketr_gto_blocks);
    
    CT3RectFlatBuffer<double> buffer(aux_indices, mask_indices, basis.dimensions_of_basis(), 3);
    
    // set distributor
    
    CT3CGeom0X0Distributor distributor(&buffer);
    
    // main compute loop
    
    std::ranges::for_each(bra_gto_blocks, [&](const auto& gblock) {
        auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
        std::ranges::for_each(gto_pair_blocks, [&](const auto& gp_pairs) {
            if constexpr (N == 1)
            {
                t3cerifunc::compute_geom_010(distributor, gblock, gp_pairs, bra_range);
            }
        });
    });
    
    return buffer;
}

template <int N>
auto
CThreeCenterElectronRepulsionGeom0X0Driver<N>::compute(const CMolecularBasis&     basis,
                                                       const CMolecularBasis&     aux_basis,
                                                       const CMolecule&           molecule,
                                                       const std::vector<double>& gamma,
                                                       const std::vector<double>& density,
                                                       const int                  iatom) const -> std::array<double, 3>
{
    // set up number of rows
    
    const auto nrows = basis.dimensions_of_basis();
    
    // set up GTOs data
    
    const auto bra_gto_blocks = gtofunc::make_gto_blocks(aux_basis, molecule);
    
    const auto ketr_gto_blocks = gtofunc::make_gto_blocks(basis, molecule, {iatom, });
    
    const auto ketf_gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(ketr_gto_blocks, ketf_gto_blocks);
    
    // set up accumalation vector
    
    const auto nblocks = bra_gto_blocks.size();
    
    const auto ngpblocks = gto_pair_blocks.size();
    
    auto gvec_xyz = std::vector<std::array<double, 3>>(nblocks * ngpblocks, {0.0, 0.0, 0.0});
    
    // prepare pointers for OMP parallel region
    
    auto ptr_bra_gto_blocks = bra_gto_blocks.data();
    
    auto ptr_gto_pair_blocks = gto_pair_blocks.data();
    
    auto ptr_gamma = gamma.data();
    
    auto ptr_density = density.data();
    
    auto ptr_gvec_xyz = gvec_xyz.data();
    
    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_gto_pair_blocks, ptr_bra_gto_blocks, ptr_gamma, ptr_density, ptr_gvec_xyz, nrows, nblocks, ngpblocks)
    {
#pragma omp single nowait
        {
            size_t index = 0;
            
            for (size_t i = 0; i < nblocks; i++)
            {
                for (size_t j = 0; j < ngpblocks; j++)
                {
#pragma omp task firstprivate(index, i, j)
                    {
                        auto gblock = ptr_bra_gto_blocks[i];
                        
                        auto gp_pairs = ptr_gto_pair_blocks[j];
                        
                        CT3CGeom010SumDistributor distributor(ptr_gvec_xyz[index].data(), ptr_gamma, ptr_density, nrows);
                        
                        auto bra_range = std::pair<size_t, size_t>{size_t{0}, gblock.number_of_basis_functions()};
                        
                        if constexpr (N == 1)
                        {
                            t3cerifunc::compute_geom_010(distributor, gblock, gp_pairs, bra_range);
                        }
                    }
                    
                    index++;
                }
            }
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

#endif /* ThreeCenterElectronRepulsionGeom0X0Driver_hpp */
