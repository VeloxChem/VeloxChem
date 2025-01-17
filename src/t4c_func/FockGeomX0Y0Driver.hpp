#ifndef FockGeomX0Y0Driver_hpp
#define FockGeomX0Y0Driver_hpp

#include <string>

#include "ElectronRepulsionGeom1010Func.hpp"
#include "GtoFunc.hpp"
#include "GtoPairBlockFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "T4CGeomXYMatricesDistributor.hpp"

/// Class CFockGeomX0Y0Driver provides methods for computing Fock matrices
/// using derivatives of four center electron repulsion integrals.
template <int N, int M>
class CFockGeomX0Y0Driver
{
   public:
    /// Creates a Fock matrices  driver.
    CFockGeomX0Y0Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The Fock matrices driver to be copied.
    CFockGeomX0Y0Driver(const CFockDriver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The Fock matrices driver  to be moved.
    CFockGeomX0Y0Driver(CFockDriver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CFockGeomX0Y0Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The Fock matrices driver to be copy assigned.
    /// @return The assigned Fock matrices driver.
    auto operator=(const CFockGeomX0Y0Driver &other) -> CFockDriver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The Fock matrices driver to be move assigned.
    /// @return The assigned Fock matrices driver .
    auto operator=(CFockGeomX0Y0Driver &&other) noexcept -> CFockDriver & = delete;

    /// @brief The equality operator.
    /// @param other The Fock matrices driver  to be compared.
    /// @return True if Fock matrices drivers  are equal, False otherwise.
    auto operator==(const CFockGeomX0Y0Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The Fock matrices driver to be compared.
    /// @return True if Fock matrices drivers  are not equal, False otherwise.
    auto operator!=(const CFockGeomX0Y0Driver &other) const -> bool = delete;

    /// @brief Computes Fock matrix derivatives for given density, basis and molecule (N^4 scaling).
    /// @param iatom The index of atom to compute derivatives of Fock matrix.
    /// @param jatom The index of atom to compute derivatives of Fock matrix.
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
                 const int              jatom,
                 const std::string     &label,
                 const double           exchange_factor,
                 const double           omega) const -> CMatrices;

    auto compute(const CMolecularBasis& basis,
                 const CT4CScreener&    screener_atom_i,
                 const CT4CScreener&    screener_atom_j,
                 const CMatrix&         density,
                 const CMatrix&         density2,
                 const int              iatom,
                 const int              jatom,
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

template <int N, int M>
auto
CFockGeomX0Y0Driver<N, M>::compute(const CMolecularBasis &basis,
                                   const CMolecule       &molecule,
                                   const CMatrix         &density,
                                   const int              iatom,
                                   const int              jatom,
                                   const std::string     &label,
                                   const double           exchange_factor,
                                   const double           omega) const -> CMatrices
{
    // set up Fock matrices

    auto fock_mats = matfunc::make_matrices(
        std::array<int, 2>{
            N, M
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
    
    const auto ket_pairs = gtofunc::make_gto_blocks(basis,
                                                    molecule,
                                                    {
                                                        jatom,
                                                    });
    
    const auto aux_pairs = gtofunc::make_gto_blocks(basis, molecule);

    const auto bra_gto_pair_blocks = gtofunc::make_gto_pair_blocks(bra_pairs, aux_pairs);

    const auto ket_gto_pair_blocks = gtofunc::make_gto_pair_blocks(ket_pairs, aux_pairs);
    
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
                    CT4CGeomXYMatricesDistributor distributor(ptr_focks, ptr_density, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    auto bra_range = std::pair<size_t, size_t>(0, bra_gpairs.number_of_contracted_pairs());
                    auto ket_range = std::pair<size_t, size_t>(0, ket_gpairs.number_of_contracted_pairs());
                    if constexpr ((N == 1) && (M == 1))
                    {
                        erifunc::compute_geom_1010<CT4CGeomXYMatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range);
                    }
                    distributor.accumulate(bra_gpairs, ket_gpairs);
                }
            });
        }
    }

    return fock_mats;
}

template <int N, int M>
auto
CFockGeomX0Y0Driver<N, M>::compute(const CMolecularBasis& basis,
                                   const CT4CScreener&    screener_atom_i,
                                   const CT4CScreener&    screener_atom_j,
                                   const CMatrix&         density,
                                   const CMatrix&         density2,
                                   const int              iatom,
                                   const int              jatom,
                                   const std::string&     label,
                                   const double           exchange_factor,
                                   const double           omega,
                                   const int              ithreshold) const -> std::vector<double>
{
    auto bsfac = _determine_block_size_factor(_get_nao(density));

    // prepare pointers for OMP parallel region

    auto ptr_screener_atom_i = &screener_atom_i;

    auto ptr_screener_atom_j = &screener_atom_j;

    auto ptr_density = &density;

    auto ptr_density_2 = &density2;

    auto nthreads = omp_get_max_threads();

    CDenseMatrix omp_values;

    if constexpr ((N == 1) && (M == 1))
    {
        omp_values = CDenseMatrix (nthreads, 9);
    }

    omp_values.zero();

    auto ptr_omp_values = &omp_values;

    // execute OMP tasks with static scheduling

#pragma omp parallel shared(ptr_omp_values, ptr_screener_atom_i, ptr_screener_atom_j, ptr_density, ptr_density_2, label, exchange_factor, omega, ithreshold)
    {
#pragma omp single nowait
        {
            auto bra_gto_pair_blocks = ptr_screener_atom_i->gto_pair_blocks();

            auto ket_gto_pair_blocks = ptr_screener_atom_j->gto_pair_blocks();

            const auto work_tasks = omp::make_bra_ket_work_group(bra_gto_pair_blocks, ket_gto_pair_blocks, ithreshold, bsfac);

            std::ranges::for_each(std::views::reverse(work_tasks), [&](const auto& task) {
                const auto bra_gpairs = bra_gto_pair_blocks[task[0]].gto_pair_block(static_cast<int>(task[2]));
                const auto ket_gpairs = ket_gto_pair_blocks[task[1]].gto_pair_block(static_cast<int>(task[3]));
                const auto bra_range  = std::pair<size_t, size_t>{task[4], task[5]};
                const auto ket_range  = std::pair<size_t, size_t>{task[6], task[7]};
#pragma omp task firstprivate(bra_gpairs, ket_gpairs, bra_range, ket_range)
                {
                    CT4CGeomXYMatricesDistributor distributor(ptr_density, ptr_density_2, label, exchange_factor, omega);
                    distributor.set_indices(bra_gpairs, ket_gpairs);
                    if constexpr ((N == 1) && (M == 1))
                    {
                        distributor.set_num_values(9);
                        erifunc::compute_geom_1010<CT4CGeomXYMatricesDistributor>(distributor, bra_gpairs, ket_gpairs, bra_range, ket_range);
                    }
                    auto values = distributor.get_values();
                    auto thread_id = omp_get_thread_num();
                    for (int idx = 0; idx < static_cast<int>(values.size()); idx++)
                    {
                        ptr_omp_values->row(thread_id)[idx] += values[idx];
                    }
                }
            });
        }
    }

    auto n_components = omp_values.getNumberOfColumns();

    std::vector<double> values(n_components, 0.0);

    for (int thread_id = 0; thread_id < nthreads; thread_id++)
    {
        for (int d = 0; d < n_components; d++)
        {
            values[d] += omp_values.row(thread_id)[d];
        }
    }

    return values;
}

template <int N, int M>
auto
CFockGeomX0Y0Driver<N, M>::set_block_size_factor(const int factor) -> void
{
    _block_size_factor = factor;
}

template <int N, int M>
auto
CFockGeomX0Y0Driver<N, M>::_get_nao(const CMatrix& mat) const -> int
{
    return mat.number_of_rows();
}

template <int N, int M>
auto
CFockGeomX0Y0Driver<N, M>::_determine_block_size_factor(const int nao) const -> int
{
    if (nao < 300)
    {
        return 32 * _block_size_factor;
    }
    else if (nao < 1500)
    {
        return 16 * _block_size_factor;
    }
    else if (nao < 3000)
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

#endif /* FockGeomX0Y0Driver_hpp */
