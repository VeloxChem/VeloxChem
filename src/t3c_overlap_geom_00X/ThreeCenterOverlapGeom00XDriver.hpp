#ifndef ThreeCenterOverlapGeom00XDriver_hpp
#define ThreeCenterOverlapGeom00XDriver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "ThreeCenterOverlapGeom001Func.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class  COverlapGeomX00Driver provides methods for computing arbitrary order two-center
/// overlap integral derivatives with respect bra side.
template <int N>
class CThreeCenterOverlapGeom00XDriver
{
   public:
    /// @brief Creates an overlap derivative integrals driver.
    CThreeCenterOverlapGeom00XDriver() = default;

    /// @brief The default copy constructor.
    /// @param other The overlap derivative integrals driver to be copied.
    CThreeCenterOverlapGeom00XDriver(const CThreeCenterOverlapGeom00XDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The overlap derivative integrals driver  to be moved.
    CThreeCenterOverlapGeom00XDriver(CThreeCenterOverlapGeom00XDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CThreeCenterOverlapGeom00XDriver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The overlap derivative integrals driver to be copy assigned.
    /// @return The assigned overlap derivative integrals driver.
    auto operator=(const CThreeCenterOverlapGeom00XDriver &other) -> CThreeCenterOverlapGeom00XDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The overlap derivative integrals driver to be move assigned.
    /// @return The assigned overlap derivative integrals driver .
    auto operator=(CThreeCenterOverlapGeom00XDriver &&other) noexcept -> CThreeCenterOverlapGeom00XDriver & = delete;

    /// @brief The equality operator.
    /// @param other The overlap derivative integrals driver  to be compared.
    /// @return True if overlap derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CThreeCenterOverlapGeom00XDriver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The overlap derivative integrals driver to be compared.
    /// @return True if overlap derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CThreeCenterOverlapGeom00XDriver &other) const -> bool = delete;

    /// @brief Computes overlapl matrix derivative fro molecule and molecular basis.
    /// @param exponents The vector of Gausian exponents of external centers.
    /// @param factors The vector of scaling factors of Gaussian exponents of external centers.
    /// @param coordinates The vector of external center coordinates.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const std::vector<double>         &exponents,
                 const std::vector<double>         &factors,
                 const std::vector<TPoint<double>> &coordinates,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule,
                 const int                         iatom) const -> CMatrices;
};

template <int N>
auto
CThreeCenterOverlapGeom00XDriver<N>::compute(const std::vector<double>         &exponents,
                                                     const std::vector<double>         &factors,
                                                     const std::vector<TPoint<double>> &coordinates,
                                                     const CMolecularBasis             &basis,
                                                     const CMolecule                   &molecule,
                                                     const int                         iatom) const -> CMatrices
{
    // set up operator derivatives matrices

    auto ovl_mats = matfunc::make_matrices(std::array<int, 1>{N,}, basis, mat_t::general);

    ovl_mats.zero();
    
    // prepare external Gaussians data
    
    auto gnorms = exponents;
    
    gnorms.insert(gnorms.end(), factors.cbegin(), factors.cend());

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_ovl_mats = &ovl_mats;
    
    auto ptr_coordinates = &coordinates;

    auto ptr_gnorms = &gnorms;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_ovl_mats, ptr_coordinates, ptr_gnorms,  iatom)
    {
#pragma omp single nowait
        {
            const auto bra_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);
            
            const auto ket_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule, {iatom, });
            
            const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto& task) {
                auto bra_gtos    = bra_gto_blocks[task[0]];
                auto ket_gtos    = ket_gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices)
                {
                    CT2CDistributor<CMatrices> distributor(ptr_ovl_mats, *ptr_coordinates, *ptr_gnorms);
                    if constexpr (N == 1)
                    {
                        t3ovlfunc::compute_geom_001(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return ovl_mats;
}

#endif /* ThreeCenterOverlapGeom00XDriver_hpp */
