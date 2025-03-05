#ifndef ThreeCenterElectronRepulsionGeomX00Driver_hpp
#define ThreeCenterElectronRepulsionGeomX00Driver_hpp

#include <vector>
#include <ranges>

#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "Point.hpp"
#include "TensorComponents.hpp"
#include "T3CUtils.hpp"
#include "T3CGeomX00Distributor.hpp"
#include "ThreeCenterElectronRepulsionGeom100Func.hpp"
#include "T4CScreener.hpp"

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
    /// @return The nuclear potential matrix.
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
    /// @return The nuclear potential matrix.
    auto compute(const CT4CScreener    &screener,
                 const CMolecularBasis &basis,
                 const CMolecularBasis &aux_basis,
                 const CMolecule       &molecule,
                 const int             iatom,
                 const int             ithreshold) const -> CT3FlatBuffer<double>;
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


#endif /* ThreeCenterElectronRepulsionGeomX00Driver_hpp */
