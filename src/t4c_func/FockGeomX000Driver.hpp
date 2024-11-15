#ifndef FockGeomX000Driver_hpp
#define FockGeomX000Driver_hpp

#include <string>
#include <vector>

#include "ElectronRepulsionGeom1000Func.hpp"
#include "DenseMatrix.hpp"
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
    auto compute(const CMolecularBasis& basis,
                 const CMolecule&       molecule,
                 const CMatrix&         density,
                 const int              iatom,
                 const std::string&     label,
                 const double           exchange_factor,
                 const double           omega) const -> CMatrices;

    auto compute(const CMolecularBasis& basis,
                 const CT4CScreener&    screener_atom,
                 const CT4CScreener&    screener,
                 const CMatrix&         density,
                 const int              iatom,
                 const std::string&     label,
                 const double           exchange_factor,
                 const double           omega,
                 const int              ithreshold) const -> CMatrices;

    auto compute(const CMolecularBasis& basis,
                 const CT4CScreener&    screener_atom,
                 const CT4CScreener&    screener,
                 const CMatrix&         density,
                 const CMatrix&         density2,
                 const int              iatom,
                 const std::string&     label,
                 const double           exchange_factor,
                 const double           omega,
                 const int              ithreshold) const -> std::vector<double>;

    auto set_block_size_factor(const int factor) -> void;

   private:
    int _block_size_factor = 4;

    auto _get_nao(const CMatrix& mat) const -> int;

    auto _determine_block_size_factor(const int nao) const -> int;
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
CFockGeomX000Driver<N>::compute(const CMolecularBasis& basis,
                                const CT4CScreener&    screener_atom,
                                const CT4CScreener&    screener,
                                const CMatrix&         density,
                                const int              iatom,
                                const std::string&     label,
                                const double           exchange_factor,
                                const double           omega,
                                const int              ithreshold) const -> CMatrices
{
    auto bsfac = _determine_block_size_factor(_get_nao(density));

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
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range)
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

template <int N>
auto
CFockGeomX000Driver<N>::compute(const CMolecularBasis& basis,
                                const CT4CScreener&    screener_atom,
                                const CT4CScreener&    screener,
                                const CMatrix&         density,
                                const CMatrix&         density2,
                                const int              iatom,
                                const std::string&     label,
                                const double           exchange_factor,
                                const double           omega,
                                const int              ithreshold) const -> std::vector<double>
{
    auto bsfac = _determine_block_size_factor(_get_nao(density));

    // prepare pointers for OMP parallel region

    auto ptr_screener = &screener;

    auto ptr_screener_atom = &screener_atom;

    auto ptr_density = &density;

    auto ptr_density_2 = &density2;

    auto nthreads = omp_get_max_threads();

    CDenseMatrix omp_values(nthreads, 3);

    omp_values.zero();

    auto ptr_omp_values = &omp_values;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_omp_values, ptr_screener, ptr_screener_atom, ptr_density, ptr_density_2, label, exchange_factor, omega, ithreshold)
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
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range)
                {
                    CT4CGeomX0MatricesDistributor distributor(ptr_density, ptr_density_2, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    erifunc::compute_geom_1000<CT4CGeomX0MatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range);
                    auto values = distributor.get_values();
                    auto thread_id = omp_get_thread_num();
                    ptr_omp_values->row(thread_id)[0] += values[0];
                    ptr_omp_values->row(thread_id)[1] += values[1];
                    ptr_omp_values->row(thread_id)[2] += values[2];
                }
            });
        }
    }

    std::vector<double> values(3, 0.0);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int d = 0; d < 3; d++)
        {
            values[d] += omp_values.row(thread_id)[d];
        }
    }

    return values;
}

template <int N>
auto
CFockGeomX000Driver<N>::set_block_size_factor(const int factor) -> void
{
    _block_size_factor = factor;
}

template <int N>
auto
CFockGeomX000Driver<N>::_get_nao(const CMatrix& mat) const -> int
{
    return mat.number_of_rows();
}

template <int N>
auto
CFockGeomX000Driver<N>::_determine_block_size_factor(const int nao) const -> int
{
    if (nao < 300)
    {
        return 32 * _block_size_factor;
    }
    else if (nao < 1500)
    {
        return 16 * _block_size_factor;
    }
    else if (nao < 4500)
    {
        return 8 * _block_size_factor;
    }
    else if (nao < 6000)
    {
        return 4 * _block_size_factor;
    }
    else if (nao < 10000)
    {
        return 2 * _block_size_factor;
    }
    else
    {
        return _block_size_factor;
    }
}

#endif /* FockGeomX000Driver_hpp */
