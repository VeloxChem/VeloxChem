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

#ifndef AtomCorePotential_hpp
#define AtomCorePotential_hpp

#include <cstddef>
#include <vector>

#include "BaseCorePotential.hpp"

/// @brief Class CAtomCorePotential stores data about atom core potential and
/// provides set of methods for handling of atom core potential data.
class CAtomCorePotential
{
   public:
    /// @brief The default constructor.
    CAtomCorePotential();

    /// @brief The constructor with local and projected potentials.
    /// @param local_potential The local potential of atom's core potential.
    /// @param projected_potentials The vector of projected potentials of atom's core potential.
    /// @param angular_momentums The vector of angular momentums associated with projected potentials.
    /// @param core_electrons The number of core electrons.
    CAtomCorePotential(const CBaseCorePotential              &local_potential,
                       const std::vector<CBaseCorePotential> &projected_potentials,
                       const std::vector<int>                &angular_momentums,
                       const int                             core_electrons);

    /// @brief The default copy constructor.
    /// @param other The atom core potential to be copied.
    CAtomCorePotential(const CAtomCorePotential &other);

    /// @brief The default move constructor.
    /// @param other The atom core potential to be moved.
    CAtomCorePotential(CAtomCorePotential &&other) noexcept;

    /// @brief The default destructor.
    ~CAtomCorePotential() = default;

    /// @brief The default copy assignment operator.
    /// @param other The atom core potential to be copy assigned.
    /// @return The assigned atom core potential.
    auto operator=(const CAtomCorePotential &other) -> CAtomCorePotential&;

    /// @brief The default move assignment operator.
    /// @param other The atom core potential to be move assigned.
    /// @return The assigned atom core potential.
    auto operator=(CAtomCorePotential &&other) noexcept -> CAtomCorePotential &;

    /// @brief The equality operator.
    /// @param other The atom core potential to be compared.
    /// @return True if atom core potentials are equal, False otherwise.
    auto operator==(const CAtomCorePotential &other) const -> bool;

    /// @brief The equality operator.
    /// @param other The atom core potential to be compared.
    /// @return True if atom core potentials are not equal, False otherwise.
    auto operator!=(const CAtomCorePotential &other) const -> bool;

    /// @brief Sets local potential in atom core potential.
    /// @param local_potential The local potential.
    auto set_local_potential(const CBaseCorePotential &local_potential) -> void;

    /// @brief Sets projected core potentials and associated angular momentums in atom core potential.
    /// @param projected_potentials The vector of projected core potentials.
    /// @param angular_momentums The vector of angular momentums.
    auto set_projected_potentials(const std::vector<CBaseCorePotential> &projected_potentials,
                                  const std::vector<int>                &angular_momentums) -> void;
   
    /// @brief Adds projected potential to atom core potential.
    /// @param projected_potential The projected potential.
    /// @param angular_momentum The angular momentum of projected potential.
    auto add_projected_potential(const CBaseCorePotential &projected_potential,
                                 const int                 angular_momentum) -> void;

    /// @brief Sets number of core electrons..
    /// @param core_electrons The number of core electrons..
    auto set_number_core_electrons(const int core_electrons) -> void;
    
    /// @brief Gets local potential from atom core potential.
    /// @return The local potential.
    auto get_local_potential() const -> CBaseCorePotential;

    /// @brief Gets vector of projected potentials from atom core potential.
    /// @return The vector of projected potentials.
    auto get_projected_potentials() const -> std::vector<CBaseCorePotential>;

    /// @brief Gets vector of angular momentums associated with projected potentials in atom core potential.
    /// @return The vector of angular momentums.
    auto get_angular_momentums() const -> std::vector<int>;

    /// @brief Gets number of core electrons in atom core potential.
    /// @return The number of core electrons.
    auto number_of_core_electrons() const -> int;

   private:
    
    /// @brief The local core potential.
    CBaseCorePotential _local_potential;

    /// @brief The vector of projected core potentials.
    std::vector<CBaseCorePotential> _projected_potentials;

    /// @brief The vector of angular momentums associated with projected core potentials.
    std::vector<int> _angular_momentums;
    
    /// @brief The number of core electrons.
    int _core_electrons;
};

#endif /* AtomCorePotential_hpp */
