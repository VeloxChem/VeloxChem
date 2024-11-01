#ifndef FockGeomX000Driver_hpp
#define FockGeomX000Driver_hpp

#include <string>

#include "ElectronRepulsionGeom1000Func.hpp"
#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OpenMPFunc.hpp"
#include "T4CGeomX0MatricesDistributor.hpp"

/// Class CFockGeomX000Driver provides methods for computing Fock matrices
/// using derivatives of four center electron repulsion integrals.
template <int N>
class CFockGeomX000Driver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockGeomX000Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CFockGeomX000Driver(const CFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CFockGeomX000Driver(CFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CFockGeomX000Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CFockGeomX000Driver &other) -> CFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CFockGeomX000Driver &&other) noexcept -> CFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CFockGeomX000Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CFockGeomX000Driver &other) const -> bool = delete;

    /// @brief Computes Fock matrix derivatives for given density, basis and molecule (N^4 scaling).
    /// @param iatom The index of atom to compute derivatives of Fock matrix.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param density The density matrix to construct Fock matrix.
    /// @param label The label of Fock matrix type.
    /// @param exchange_factor The exchange-correlation factors.
    /// @param omega The range separation factor.
    /// @return The Fock matrices.
    auto compute(const CMolecularBasis &basis,
                 const CMolecule       &molecule,
                 const CMatrix         &density,
                 const int              iatom,
                 const std::string     &label,
                 const double           exchange_factor,
                 const double           omega) const -> CMatrices;

    auto compute_with_screening(const CMolecularBasis &basis,
                                const CMolecule       &molecule,
                                const CT4CScreener&   screener,
                                const CT4CScreener&   screener_atom,
                                const CMatrix         &density,
                                const int              iatom,
                                const std::string     &label,
                                const double           exchange_factor,
                                const double           omega,
                                const int              ithreshold) const -> CMatrices;

   private:
    int _block_size_factor = 16;
};

template <int N>
auto
CFockGeomX000Driver<N>::compute(const CMolecularBasis &basis,
                                const CMolecule       &molecule,
                                const CMatrix         &density,
                                const int              iatom,
                                const std::string     &label,
                                const double           exchange_factor,
                                const double           omega) const -> CMatrices
{
    // set up Fock matrices

    auto fock_mats = matfunc::make_matrices(
        std::array<int, 1>{
            1,
        },
        basis,
        mat_t::general);

    fock_mats.zero();

    // set basis function pair blocks

    const auto bra_pairs = gtofunc::make_gto_blocks(basis,
                                                    molecule,
                                                    {
                                                        iatom,
                                                    });

    const auto ket_pairs = gtofunc::make_gto_blocks(basis, molecule);

    const auto bra_gto_pair_blocks = gtofunc::make_gto_pair_blocks(bra_pairs, ket_pairs);

    const auto ket_gto_pair_blocks = gtofunc::make_gto_pair_blocks(ket_pairs);

    // prepare pointers for OMP parallel region

    auto ptr_bra_gto_pair_blocks = &bra_gto_pair_blocks;

    auto ptr_ket_gto_pair_blocks = &ket_gto_pair_blocks;

    auto ptr_density = &density;

    auto ptr_focks = &fock_mats;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_bra_gto_pair_blocks, ptr_ket_gto_pair_blocks, ptr_density, ptr_focks, label, exchange_factor, omega)
    {
#pragma omp single nowait
        {
            const auto n_bra_blocks = ptr_bra_gto_pair_blocks->size();

            const auto n_ket_blocks = ptr_ket_gto_pair_blocks->size();

            auto ptr_bra_gto_pairs_data = ptr_bra_gto_pair_blocks->data();

            auto ptr_ket_gto_pairs_data = ptr_ket_gto_pair_blocks->data();

            std::ranges::for_each(views::rectangular(n_bra_blocks, n_ket_blocks) | std::views::reverse, [&](const auto &index) {
                const size_t i = index.first;
                const size_t j = index.second;
#pragma omp task firstprivate(i, j)
                {
                    auto                          bra_gpairs = ptr_bra_gto_pairs_data[i];
                    auto                          ket_gpairs = ptr_ket_gto_pairs_data[j];
                    CT4CGeomX0MatricesDistributor distributor(ptr_focks, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    auto bra_range = std::pair<size_t, size_t>(0, bra_gpairs.number_of_contracted_pairs());
                    auto ket_range = std::pair<size_t, size_t>(0, ket_gpairs.number_of_contracted_pairs());
                    erifunc::compute_geom_1000<CT4CGeomX0MatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    return fock_mats;
}

template <int N>
auto
CFockGeomX000Driver<N>::compute_with_screening(const CMolecularBasis &basis,
                                               const CMolecule       &molecule,
                                               const CT4CScreener&   screener,
                                               const CT4CScreener&   screener_atom,
                                               const CMatrix         &density,
                                               const int              iatom,
                                               const std::string     &label,
                                               const double           exchange_factor,
                                               const double           omega,
                                               const int             ithreshold) const -> CMatrices
{
    const auto nao = density.number_of_rows();

    int bsfac = _block_size_factor;

    if (nao < 450) { bsfac = 16 * _block_size_factor; }
    else if (nao < 900) { bsfac = 8 * _block_size_factor; }
    else if (nao < 1800) { bsfac = 4 * _block_size_factor; }
    else if (nao < 3600) { bsfac = 2 * _block_size_factor; }
    else { bsfac = _block_size_factor; }

    // set up Fock matrices

    auto fock_mats = matfunc::make_matrices(
        std::array<int, 1>{
            1,
        },
        basis,
        mat_t::general);

    fock_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_screener_atom = &screener_atom;

    auto ptr_density = &density;

    auto ptr_focks = &fock_mats;

    /*
    // set basis function pair blocks

    const auto bra_pairs = gtofunc::make_gto_blocks(basis,
                                                    molecule,
                                                    {
                                                        iatom,
                                                    });

    const auto ket_pairs = gtofunc::make_gto_blocks(basis, molecule);

    const auto bra_gto_pair_blocks = gtofunc::make_gto_pair_blocks(bra_pairs, ket_pairs);

    const auto ket_gto_pair_blocks = gtofunc::make_gto_pair_blocks(ket_pairs);

    // prepare pointers for OMP parallel region

    auto ptr_bra_gto_pair_blocks = &bra_gto_pair_blocks;

    auto ptr_ket_gto_pair_blocks = &ket_gto_pair_blocks;

    auto ptr_density = &density;

    auto ptr_focks = &fock_mats;
    */

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_screener, ptr_screener_atom, ptr_density, ptr_focks, label, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto bra_gto_pair_blocks = ptr_screener_atom->gto_pair_blocks();

            auto ket_gto_pair_blocks = ptr_screener->gto_pair_blocks();

            const auto work_tasks = omp::make_bra_ket_work_group(bra_gto_pair_blocks, ket_gto_pair_blocks, ithreshold, bsfac);

            std::ranges::for_each(std::views::reverse(work_tasks), [&](const auto& task) {
                const auto bra_gpairs = bra_gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = ket_gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
                //const bool diagonal   = (task[0] == task[1]) && (task[2] == task[3]) && (bra_range == ket_range);
                const bool diagonal   = false;
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range, diagonal)
                {
                    CT4CGeomX0MatricesDistributor distributor(ptr_focks, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute_geom_1000<CT4CGeomX0MatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range);
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    return fock_mats;
}

#endif /* FockGeomX000Driver_hpp */
