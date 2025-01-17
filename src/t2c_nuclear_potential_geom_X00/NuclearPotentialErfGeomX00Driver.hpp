#ifndef NuclearPotentialErfGeomX00Driver_hpp
#define NuclearPotentialErfGeomX00Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "NuclearPotentialErfGeom100Func.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class  CNuclearPotentialErfGeomX00Driver provides methods for computing arbitrary order two-center
/// overlap integral derivatives with respect bra side.
template <int N>
class CNuclearPotentialErfGeomX00Driver
{
   public:
    /// @brief Creates an nuclear potential derivative integrals driver.
    CNuclearPotentialErfGeomX00Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The nuclear potential derivative integrals driver to be copied.
    CNuclearPotentialErfGeomX00Driver(const CNuclearPotentialErfGeomX00Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The nuclear potential derivative integrals driver  to be moved.
    CNuclearPotentialErfGeomX00Driver(CNuclearPotentialErfGeomX00Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialErfGeomX00Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be copy assigned.
    /// @return The assigned nuclear potential derivative integrals driver.
    auto operator=(const CNuclearPotentialErfGeomX00Driver &other) -> CNuclearPotentialErfGeomX00Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be move assigned.
    /// @return The assigned nuclear potential derivative integrals driver .
    auto operator=(CNuclearPotentialErfGeomX00Driver &&other) noexcept -> CNuclearPotentialErfGeomX00Driver & = delete;

    /// @brief The equality operator.
    /// @param other The nuclear potential derivative integrals driver  to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialErfGeomX00Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The nuclear potential derivative integrals driver to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialErfGeomX00Driver &other) const -> bool = delete;

    /// @brief Computes overlapl matrix derivative fro molecule and molecular basis.
    /// @param omegas The vector of range-separation factors.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const std::vector<double> &omegas, const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices;
};

template <int N>
auto
CNuclearPotentialErfGeomX00Driver<N>::compute(const std::vector<double> &omegas, const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices
{
    // set up operator derivatives matrices

    auto npot_mats = matfunc::make_matrices(std::array<int, 1>{N}, basis, mat_t::general);

    npot_mats.zero();

    // atoms charges and coordinates

    auto charges = molecule.charges();

    auto coordinates = molecule.coordinates("au");

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;
    
    auto ptr_omegas = &omegas;

    auto ptr_npot_mats = &npot_mats;

    auto ptr_charges = &charges;

    auto ptr_coordinates = &coordinates;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_omegas, ptr_npot_mats, ptr_charges, ptr_coordinates, iatom)
    {
#pragma omp single nowait
        {
            const auto bra_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis,
                                                                 *ptr_molecule,
                                                                 {
                                                                     iatom,
                                                                 });

            const auto ket_gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);

            const auto tasks = omp::make_work_tasks(bra_gto_blocks, ket_gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto &task) {
                auto bra_gtos    = bra_gto_blocks[task[0]];
                auto ket_gtos    = ket_gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices)
                {
                    CT2CDistributor<CMatrices> distributor(ptr_npot_mats, *ptr_coordinates, *ptr_charges);
                    if constexpr (N == 1)
                    {
                        npotfunc::compute_geom_100(distributor, *ptr_omegas, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return npot_mats;
}

#endif /* NuclearPotentialErfGeomX00Driver_hpp */
