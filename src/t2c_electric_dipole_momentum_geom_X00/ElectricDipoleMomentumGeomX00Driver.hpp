#ifndef ElectricDipoleMomentumGeomX00Driver_hpp
#define ElectricDipoleMomentumGeomX00Driver_hpp

#include <vector>

#include "ElectricDipoleMomentumGeom100Func.hpp"
#include "GtoFunc.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class  CElectricDipoleMomentumGeomX00Driver provides methods for computing arbitrary order two-center
/// overlap integral derivatives with respect bra side.
template <int N>
class CElectricDipoleMomentumGeomX00Driver
{
   public:
    /// @brief Creates an electric dipole momentum derivative integrals driver.
    CElectricDipoleMomentumGeomX00Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The electric dipole momentum derivative integrals driver to be copied.
    CElectricDipoleMomentumGeomX00Driver(const CElectricDipoleMomentumGeomX00Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The electric dipole momentum derivative integrals driver  to be moved.
    CElectricDipoleMomentumGeomX00Driver(CElectricDipoleMomentumGeomX00Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CElectricDipoleMomentumGeomX00Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The electric dipole momentum derivative integrals driver to be copy assigned.
    /// @return The assigned electric dipole momentum derivative integrals driver.
    auto operator=(const CElectricDipoleMomentumGeomX00Driver &other) -> CElectricDipoleMomentumGeomX00Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The electric dipole momentum derivative integrals driver to be move assigned.
    /// @return The assigned electric dipole momentum derivative integrals driver .
    auto operator=(CElectricDipoleMomentumGeomX00Driver &&other) noexcept -> CElectricDipoleMomentumGeomX00Driver & = delete;

    /// @brief The equality operator.
    /// @param other The electric dipole momentum derivative integrals driver  to be compared.
    /// @return True if electric dipole momentum derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CElectricDipoleMomentumGeomX00Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The electric dipole momentum derivative integrals driver to be compared.
    /// @return True if electric dipole momentum derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CElectricDipoleMomentumGeomX00Driver &other) const -> bool = delete;

    /// @brief Computes overlapl matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param origin The origin of electric dipole momentum.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const TPoint<double> &origin, const int iatom) const -> CMatrices;
};

template <int N>
auto
CElectricDipoleMomentumGeomX00Driver<N>::compute(const CMolecularBasis &basis,
                                                 const CMolecule       &molecule,
                                                 const TPoint<double>  &origin,
                                                 const int              iatom) const -> CMatrices
{
    // set up operator derivatives matrices

    auto dip_mats = matfunc::make_matrices(std::array<int, 2>{N, 1}, basis, mat_t::general);

    dip_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_origin = &origin;

    auto ptr_dip_mats = &dip_mats;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_origin, ptr_dip_mats, iatom)
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
                    const auto                 coords = std::vector<TPoint<double>>(1, *ptr_origin);
                    const auto                 data   = std::vector<double>();
                    CT2CDistributor<CMatrices> distributor(ptr_dip_mats, coords, data);
                    if constexpr (N == 1)
                    {
                        dipfunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return dip_mats;
}

#endif /* ElectricDipoleMomentumGeomX00Driver_hpp */
