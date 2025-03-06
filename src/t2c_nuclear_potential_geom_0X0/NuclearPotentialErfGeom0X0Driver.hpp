#ifndef NuclearPotentialErfGeom0X0Driver_hpp
#define NuclearPotentialErfGeom0X0Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "NuclearPotentialErfGeom010Func.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class CNuclearPotentialGeom0X0Driver provides methods for computing arbitrary order two-center
/// nuclear potential derivative integrals derivatives with respect to operator.
template <int N>
class CNuclearPotentialErfGeom0X0Driver
{
   public:
    /// @brief Creates an nuclear potential derivative integrals driver.
    CNuclearPotentialErfGeom0X0Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The nuclear potential derivative integrals driver to be copied.
    CNuclearPotentialErfGeom0X0Driver(const CNuclearPotentialErfGeom0X0Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The nuclear potential derivative integrals driver  to be moved.
    CNuclearPotentialErfGeom0X0Driver(CNuclearPotentialErfGeom0X0Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialErfGeom0X0Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be copy assigned.
    /// @return The assigned nuclear potential derivative integrals driver.
    auto operator=(const CNuclearPotentialErfGeom0X0Driver &other) -> CNuclearPotentialErfGeom0X0Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be move assigned.
    /// @return The assigned nuclear potential derivative integrals driver .
    auto operator=(CNuclearPotentialErfGeom0X0Driver &&other) noexcept -> CNuclearPotentialErfGeom0X0Driver & = delete;

    /// @brief The equality operator.
    /// @param other The nuclear potential derivative integrals driver  to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialErfGeom0X0Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The nuclear potential derivative integrals driver to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialErfGeom0X0Driver &other) const -> bool = delete;

    /// @brief Computes nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param multipoles The vector of external multipoles.
    /// @param coordinates The vector of coordinates of external charges.
    /// @param omegas The vector of range-separation factors.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @return The nuclear potential matrix.
    auto compute(const std::vector<double>         &multipoles,
                 const std::vector<TPoint<double>> &coordinates,
                 const std::vector<double>         &omegas,
                 const CMolecularBasis             &basis,
                 const CMolecule                   &molecule) const -> CMatrices;

    /// @brief Computes nuclear potential matrix for given set of external charges,  molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param omega The range-separation factor.
    /// @param iatom The index of selected atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const double omega, const int iatom) const -> CMatrices;
};

template <int N>
auto
CNuclearPotentialErfGeom0X0Driver<N>::compute(const std::vector<double>         &multipoles,
                                              const std::vector<TPoint<double>> &coordinates,
                                              const std::vector<double>         &omegas,
                                              const CMolecularBasis             &basis,
                                              const CMolecule                   &molecule) const -> CMatrices
{
    // set up operator derivatives matrices

    auto op_mats = matfunc::make_matrices(std::array<int, 1>{N}, basis, mat_t::symmetric);

    op_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_op_mats = &op_mats;

    auto ptr_coordinates = &coordinates;

    auto ptr_multipoles = &multipoles;
    
    auto ptr_omegas = &omegas;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_op_mats, ptr_omegas, ptr_multipoles, ptr_coordinates)
    {
#pragma omp single nowait
        {
            const auto gto_blocks = gtofunc::make_gto_blocks(*ptr_basis, *ptr_molecule);

            const auto tasks = omp::make_work_tasks(gto_blocks);

            std::ranges::for_each(std::ranges::reverse_view(tasks), [&](const auto &task) {
                auto bra_gtos    = gto_blocks[task[0]];
                auto ket_gtos    = gto_blocks[task[1]];
                auto bra_indices = std::pair<size_t, size_t>{task[2], task[3]};
                auto ket_indices = std::pair<size_t, size_t>{task[4], task[5]};
                bool bkequal     = (task[0] == task[1]) && (task[2] == task[4]) && (task[3] == task[5]);
#pragma omp task firstprivate(bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal)
                {
                    CT2CDistributor<CMatrices> distributor(ptr_op_mats, *ptr_coordinates, *ptr_multipoles);
                    if constexpr (N == 1)
                    {
                        npotfunc::compute_geom_010(distributor, *ptr_omegas, bra_gtos, ket_gtos, bra_indices, ket_indices, bkequal);
                    }
                    // TODO: Add other order here...
                }
            });
        }
    }

    return op_mats;
}

template <int N>
auto
CNuclearPotentialErfGeom0X0Driver<N>::compute(const CMolecularBasis &basis, const CMolecule &molecule, const double omega, const int iatom) const -> CMatrices
{
    const auto charge = molecule.charges()[iatom];

    std::vector<double> multipoles;

    if constexpr (N == 1)
    {
        multipoles = std::vector<double>({charge, charge, charge});
    }

    auto coords = molecule.atom_coordinates(iatom, "au");

    return compute(multipoles,
                   {
                       coords,
                   },
                   {
                        omega,
                   },
                   basis,
                   molecule);
}

#endif /* NuclearPotentialErfGeom0X0Driver_hpp */
