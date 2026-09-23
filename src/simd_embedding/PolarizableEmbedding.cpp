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


#include "PolarizableEmbedding.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "EmbeddingError.hpp"

CPolarizableEmbedding::CPolarizableEmbedding()

    : _polarizable(true)

    , _nonpolarizable(false)
{
}

auto
CPolarizableEmbedding::get_polarizable_region() const -> const CEmbeddingRegion &
{
    return _polarizable;
}

auto
CPolarizableEmbedding::polarizable_region() -> CEmbeddingRegion &
{
    return _polarizable;
}

auto
CPolarizableEmbedding::get_nonpolarizable_region() const -> const CEmbeddingRegion &
{
    return _nonpolarizable;
}

auto
CPolarizableEmbedding::nonpolarizable_region() -> CEmbeddingRegion &
{
    return _nonpolarizable;
}

auto
CPolarizableEmbedding::number_of_molecules() const -> int
{
    return _polarizable.number_of_molecules() + _nonpolarizable.number_of_molecules();
}

auto
CPolarizableEmbedding::number_of_sites() const -> int
{
    return _polarizable.number_of_sites() + _nonpolarizable.number_of_sites();
}

auto
CPolarizableEmbedding::number_of_polarizable_sites() const -> int
{
    return _polarizable.number_of_polarizable_sites();
}

auto
CPolarizableEmbedding::is_polarizable() const -> bool
{
    return _polarizable.is_polarizable();
}

auto
CPolarizableEmbedding::permanent_multipoles(const int order) const -> const TPermanentMultipoles &
{
    const auto &polarizable = _polarizable.permanent_multipoles(order);

    const auto &nonpolarizable = _nonpolarizable.permanent_multipoles(order);

    const auto versions = std::make_pair(_polarizable.version(), _nonpolarizable.version());

    auto &gathered = _multipoles[static_cast<size_t>(order)];

    // NOTE: the two regions hold their own multipoles and rebuild them when
    // they are added to. This only puts the two together, and needs to do that
    // again when either of them has changed. Their versions say that and the
    // number of sites does not: a molecule whose sites carry nothing of this
    // order changes the region without changing how many of them there are.

    if (_has_multipoles[static_cast<size_t>(order)] && (_multipoles_version[static_cast<size_t>(order)] == versions))
    {
        return gathered;
    }

    gathered.order = order;

    gathered.coordinates.clear();

    gathered.values.clear();

    gathered.coordinates.reserve(polarizable.coordinates.size() + nonpolarizable.coordinates.size());

    gathered.values.reserve(polarizable.values.size() + nonpolarizable.values.size());

    for (const auto *region : {&polarizable, &nonpolarizable})
    {
        gathered.coordinates.insert(gathered.coordinates.end(), region->coordinates.begin(), region->coordinates.end());

        gathered.values.insert(gathered.values.end(), region->values.begin(), region->values.end());
    }

    _multipoles_version[static_cast<size_t>(order)] = versions;

    _has_multipoles[static_cast<size_t>(order)] = true;

    return gathered;
}

auto
CPolarizableEmbedding::permanent_nuclear_energy(const CMolecule &molecule, const CMolecularBasis &basis) const -> double
{
    // NOTE: checked here, where this library throws, rather than left to the
    // molecule, where the convention is to assert and a basis built for some
    // other molecule would take the interpreter down with it.

    embedding::require(static_cast<int>(basis.basis_sets_indices().size()) == molecule.number_of_atoms(),
                       std::string("PolarizableEmbedding: The basis describes ") +
                           std::to_string(basis.basis_sets_indices().size()) + std::string(" atoms and the molecule has ") +
                           std::to_string(molecule.number_of_atoms()));

    const auto charges = molecule.effective_charges(basis);

    const auto &points = molecule.coordinates();

    const auto &gathered = permanent_multipoles(0);

    const auto nsites = static_cast<int>(gathered.number_of_sites());

    const auto nnuclei = static_cast<int>(charges.size());

    const auto *sites = gathered.coordinates.data();

    const auto *values = gathered.values.data();

    const auto *nuclei = points.data();

    const auto *nuclear_charges = charges.data();

    double energy = 0.0;

    double closest = std::numeric_limits<double>::max();

    // NOTE: the sites are the millions and the nuclei are the tens, so the
    // sites are streamed once with the nuclei sitting in cache rather than the
    // other way round. Nothing is allocated: the distances are used where they
    // are made and never gathered into an array of their own.

#pragma omp parallel for reduction(+ : energy) reduction(min : closest) if (nsites > 1024)
    for (int isite = 0; isite < nsites; isite++)
    {
        const double x = sites[3 * isite + 0];

        const double y = sites[3 * isite + 1];

        const double z = sites[3 * isite + 2];

        double potential = 0.0;

        for (int inuc = 0; inuc < nnuclei; inuc++)
        {
            const auto xyz = nuclei[inuc].coordinates();

            const double dx = x - xyz[0];

            const double dy = y - xyz[1];

            const double dz = z - xyz[2];

            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (r < closest) closest = r;

            potential += nuclear_charges[inuc] / r;
        }

        energy += values[isite] * potential;
    }

    // NOTE: refused here and not inside the loop, which is a parallel region
    // and no place to throw from.

    embedding::require(closest > 1.0e-6,
                       std::string("PolarizableEmbedding: A site of the environment sits ") + std::to_string(closest) +
                           std::string(" bohr from a nucleus of the quantum region"));

    return energy;
}
