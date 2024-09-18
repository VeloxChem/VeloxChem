#ifndef NuclearPotentialGeomXY0Driver_hpp
#define NuclearPotentialGeomXY0Driver_hpp

#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "NuclearPotentialGeom110Func.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class  CNuclearPotentialGeomXY0Driver provides methods for computing arbitrary order two-center
/// nuclear potential derivative integrals derivatives with respect to operator.
template <int N, int M>
class CNuclearPotentialGeomXY0Driver
{
   public:
    /// @brief Creates an nuclear potential derivative integrals driver.
    CNuclearPotentialGeomXY0Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The nuclear potential derivative integrals driver to be copied.
    CNuclearPotentialGeomXY0Driver(const CNuclearPotentialGeomXY0Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The nuclear potential derivative integrals driver  to be moved.
    CNuclearPotentialGeomXY0Driver(CNuclearPotentialGeomXY0Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CNuclearPotentialGeomXY0Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be copy assigned.
    /// @return The assigned nuclear potential derivative integrals driver.
    auto operator=(const CNuclearPotentialGeomXY0Driver &other) -> CNuclearPotentialGeomXY0Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The nuclear potential derivative integrals driver to be move assigned.
    /// @return The assigned nuclear potential derivative integrals driver .
    auto operator=(CNuclearPotentialGeomXY0Driver &&other) noexcept -> CNuclearPotentialGeomXY0Driver & = delete;

    /// @brief The equality operator.
    /// @param other The nuclear potential derivative integrals driver  to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CNuclearPotentialGeomXY0Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The nuclear potential derivative integrals driver to be compared.
    /// @return True if nuclear potential derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CNuclearPotentialGeomXY0Driver &other) const -> bool = delete;

    /// @brief Computes overlapl matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @param jatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom, const int jatom) const -> CMatrices;
};

template <int N, int M>
auto
CNuclearPotentialGeomXY0Driver<N, M>::compute(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom, const int jatom) const
    -> CMatrices
{
    // set up operator derivatives matrices

    auto npot_mats = matfunc::make_matrices(std::array<int, 2>{N, M}, basis, mat_t::general);

    npot_mats.zero();

    // atoms charges and coordinates

    const auto charge = molecule.charges()[jatom];

    std::vector<double> multipoles;

    if constexpr (M == 1)
    {
        multipoles = std::vector<double>({charge, charge, charge});
    }

    auto coordinates = std::vector<TPoint<double>>(1, molecule.atom_coordinates(jatom, "au"));

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_npot_mats = &npot_mats;

    auto ptr_multipoles = &multipoles;

    auto ptr_coordinates = &coordinates;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_npot_mats, ptr_multipoles, ptr_coordinates, iatom)
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
                    CT2CDistributor<CMatrices> distributor(ptr_npot_mats, *ptr_coordinates, *ptr_multipoles);
                    if constexpr ((N == 1) && (M == 1))
                    {
                        npotfunc::compute_geom_110(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return npot_mats;
}

#endif /* NuclearPotentialGeomXY0Driver_hpp */