//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef KineticEnergyGeomX00Driver_hpp
#define KineticEnergyGeomX00Driver_hpp

#include <vector>

#include "GtoFunc.hpp"
#include "KineticEnergyGeom100Func.hpp"
#include "KineticEnergyGeom200Func.hpp"
#include "Matrices.hpp"
#include "MatricesFunc.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OpenMPFunc.hpp"
#include "Point.hpp"
#include "T2CDistributor.hpp"

/// @brief Class  CKineticEnergyGeomX00Driver provides methods for computing arbitrary order two-center
/// overlap integral derivatives with respect bra side.
template <int N>
class CKineticEnergyGeomX00Driver
{
   public:
    /// @brief Creates an kinetic energy derivative integrals driver.
    CKineticEnergyGeomX00Driver() = default;

    /// @brief The default copy constructor.
    /// @param other The kinetic energy derivative integrals driver to be copied.
    CKineticEnergyGeomX00Driver(const CKineticEnergyGeomX00Driver &other) = delete;

    /// @brief The default move constructor.
    /// @param other The kinetic energy derivative integrals driver  to be moved.
    CKineticEnergyGeomX00Driver(CKineticEnergyGeomX00Driver &&other) noexcept = delete;

    /// @brief The default destructor.
    ~CKineticEnergyGeomX00Driver() = default;

    /// @brief The default copy assignment operator.
    /// @param other The kinetic energy derivative integrals driver to be copy assigned.
    /// @return The assigned kinetic energy derivative integrals driver.
    auto operator=(const CKineticEnergyGeomX00Driver &other) -> CKineticEnergyGeomX00Driver & = delete;

    /// @brief The default move assignment operator.
    /// @param other The kinetic energy derivative integrals driver to be move assigned.
    /// @return The assigned kinetic energy derivative integrals driver .
    auto operator=(CKineticEnergyGeomX00Driver &&other) noexcept -> CKineticEnergyGeomX00Driver & = delete;

    /// @brief The equality operator.
    /// @param other The kinetic energy derivative integrals driver  to be compared.
    /// @return True if kinetic energy derivative integrals drivers  are equal, False otherwise.
    auto operator==(const CKineticEnergyGeomX00Driver &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The kinetic energy derivative integrals driver to be compared.
    /// @return True if kinetic energy derivative integrals drivers  are not equal, False otherwise.
    auto operator!=(const CKineticEnergyGeomX00Driver &other) const -> bool = delete;

    /// @brief Computes overlapl matrix derivative fro molecule and molecular basis.
    /// @param basis The molecular basis.
    /// @param molecule The molecule.
    /// @param iatom The index of atom.
    /// @return The nuclear potential matrix.
    auto compute(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices;
};

template <int N>
auto
CKineticEnergyGeomX00Driver<N>::compute(const CMolecularBasis &basis, const CMolecule &molecule, const int iatom) const -> CMatrices
{
    // set up operator derivatives matrices

    auto kin_mats = matfunc::make_matrices(std::array<int, 1>{N}, basis, mat_t::general);

    kin_mats.zero();

    // prepare pointers for OMP parallel region

    auto ptr_basis = &basis;

    auto ptr_molecule = &molecule;

    auto ptr_kin_mats = &kin_mats;

    // execute OMP tasks with static scheduling

    omp::set_static_scheduler();

#pragma omp parallel shared(ptr_basis, ptr_molecule, ptr_kin_mats, iatom)
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
                    CT2CDistributor<CMatrices> distributor(ptr_kin_mats);
                    if constexpr (N == 1)
                    {
                        kinfunc::compute_geom_100(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                    if constexpr (N == 2)
                    {
                        kinfunc::compute_geom_200(distributor, bra_gtos, ket_gtos, bra_indices, ket_indices, false);
                    }
                }
            });
        }
    }

    return kin_mats;
}

#endif /* KineticEnergyGeomX00Driver_hpp */
