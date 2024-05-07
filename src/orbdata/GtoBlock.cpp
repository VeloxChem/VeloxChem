//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <cmath>

#include "ErrorHandler.hpp"
#include "GtoBlock.hpp"

CGtoBlock::CGtoBlock(const std::vector<TPoint3D>& coordinates,
                     const std::vector<double>&   exponents,
                     const std::vector<double>&   norms,
                     const std::vector<int64_t>&  orb_indexes,
                     const std::vector<int64_t>&  atm_indexes,
                     const int64_t                angmom,
                     const int64_t                npgtos)

    : _coordinates(coordinates)

    , _exponents(exponents)

    , _norms(norms)

    , _orb_indexes(orb_indexes)

    , _atm_indexes(atm_indexes)

    , _angmom(angmom)

    , _npgtos(npgtos)
{
}

CGtoBlock::CGtoBlock(const CMolecularBasis& basis, const CMolecule& molecule, const int64_t angmom, const int64_t npgtos)
{
    if (const auto gtos = basis.getBasisFunctions(angmom, npgtos); !gtos.empty())
    {
        _angmom = angmom;

        _npgtos = npgtos;

        _orb_indexes = basis.getIndexMap(angmom, npgtos);

        _atm_indexes = basis.getAtomicIndexes(angmom, npgtos);

        if (const auto natoms = _atm_indexes.size(); natoms > 0)
        {
            for (size_t i = 0; i < natoms; i++)
            {
                _coordinates.push_back(molecule.getAtomCoordinates(_atm_indexes[i]));
            }
        }

        const auto ncgtos = static_cast<int64_t>(gtos.size());

        _exponents = std::vector<double>(ncgtos * npgtos, 0.0);

        _norms = std::vector<double>(ncgtos * npgtos, 0.0);

        for (int64_t i = 0; i < ncgtos; i++)
        {
            const auto fexps = gtos[i].getExponents();

            const auto fnorms = gtos[i].getNormalizationFactors();

            for (int64_t j = 0; j < npgtos; j++)
            {
                _exponents[j * ncgtos + i] = fexps[j];

                _norms[j * ncgtos + i] = fnorms[j];
            }
        }
    }
}

CGtoBlock::CGtoBlock(const CMolecularBasis&      basis,
                     const CMolecule&            molecule,
                     const std::vector<int64_t>& atoms,
                     const int64_t               angmom,
                     const int64_t               npgtos)
{
    if (const auto gtos = basis.getBasisFunctions(atoms, angmom, npgtos); !gtos.empty())
    {
        _angmom = angmom;

        _npgtos = npgtos;

        _orb_indexes = basis.getIndexMap(atoms, angmom, npgtos);

        _atm_indexes = basis.getAtomicIndexes(atoms, angmom, npgtos);

        if (const auto natoms = _atm_indexes.size(); natoms > 0)
        {
            for (size_t i = 0; i < natoms; i++)
            {
                _coordinates.push_back(molecule.getAtomCoordinates(_atm_indexes[i]));
            }
        }

        const auto ncgtos = static_cast<int64_t>(gtos.size());

        _exponents = std::vector<double>(ncgtos * npgtos, 0.0);

        _norms = std::vector<double>(ncgtos * npgtos, 0.0);

        for (int64_t i = 0; i < ncgtos; i++)
        {
            const auto fexps = gtos[i].getExponents();

            const auto fnorms = gtos[i].getNormalizationFactors();

            for (int64_t j = 0; j < npgtos; j++)
            {
                const auto joff = j * ncgtos + i;
                
                _exponents[joff] = fexps[j];

                _norms[joff] = fnorms[j];
            }
        }
    }
}

auto
CGtoBlock::getCoordinates() const -> std::vector<TPoint3D>
{
    return _coordinates;
}

auto
CGtoBlock::getExponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CGtoBlock::getNormalizationFactors() const -> std::vector<double>
{
    return _norms;
}

auto
CGtoBlock::getOrbitalIndexes() const -> std::vector<int64_t>
{
    return _orb_indexes;
}

auto
CGtoBlock::getAtomicIndexes() const -> std::vector<int64_t>
{
    return _atm_indexes;
}

auto
CGtoBlock::getAtomicOrbitalsIndexes() const -> std::vector<int64_t>
{
    std::vector<int64_t> ao_inds;

    // go through spherical harmonics components

    for (int64_t comp = 0; comp < _angmom * 2 + 1; comp++)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            ao_inds.push_back(comp * _orb_indexes[0] + _orb_indexes[ind]);
        }
    }

    return ao_inds;
}

auto
CGtoBlock::getAtomicOrbitalsIndexesForCartesian() const -> std::vector<int64_t>
{
    errors::assertMsgCritical(_angmom <= 2, std::string("GtoBlock: getAtomicOrbitalsIndexesForCartesian only supports up to d-orbitals"));

    std::vector<int64_t> ao_inds;

    // go through Cartesian components

    for (int64_t comp = 0; comp < (_angmom + 1) * (_angmom + 2) / 2; comp++)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            ao_inds.push_back(comp * _orb_indexes[0] + _orb_indexes[ind]);
        }
    }

    return ao_inds;
}

auto
CGtoBlock::getSphericalToCartesianMapping() const -> std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>>
{
    errors::assertMsgCritical((0 <= _angmom) && (_angmom <= 2), std::string("GtoBlock: getSphericalToCartesianMapping only works for s/p/d orbitals"));

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_p;

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_comp_map;

    if (_angmom == 0)
    {
        // s (0) -> s (0)

        sph_cart_comp_map[0] = std::vector<std::pair<int64_t, double>>({{0, 1.0}});
    }
    else if (_angmom == 1)
    {
        // p-1 (0) -> py (1)
        // p_0 (1) -> pz (2)
        // p+1 (2) -> px (0)

        sph_cart_comp_map[0] = std::vector<std::pair<int64_t, double>>({{1, 1.0}});
        sph_cart_comp_map[1] = std::vector<std::pair<int64_t, double>>({{2, 1.0}});
        sph_cart_comp_map[2] = std::vector<std::pair<int64_t, double>>({{0, 1.0}});
    }
    else if (_angmom == 2)
    {
        const double f2_3 = 2.0 * std::sqrt(3.0);

        // d-2 (0) <- dxy(1) * f2_3
        // d-1 (1) <- dyz(4) * f2_3
        // d_0 (2) <- dzz(5) * 2.0 - dxx(0) - dyy(3)
        // d+1 (3) <- dxz(2) * f2_3
        // d+2 (4) <- (dxx(0) - dyy(3)) * 0.5 * f2_3

        sph_cart_comp_map[0] = std::vector<std::pair<int64_t, double>>({{1, f2_3}});
        sph_cart_comp_map[1] = std::vector<std::pair<int64_t, double>>({{4, f2_3}});
        sph_cart_comp_map[2] = std::vector<std::pair<int64_t, double>>({{0, -1.0}, {3, -1.0}, {5, 2.0}});
        sph_cart_comp_map[3] = std::vector<std::pair<int64_t, double>>({{2, f2_3}});
        sph_cart_comp_map[4] = std::vector<std::pair<int64_t, double>>({{0, 0.5 * f2_3}, {3, -0.5 * f2_3}});
    }

    for (const auto& [sph_comp, cart_comp_coef_vec] : sph_cart_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            auto sph_ind = sph_comp * _orb_indexes[0] + _orb_indexes[ind];

            sph_cart_p[sph_ind] = std::vector<std::pair<int64_t, double>>();

            for (const auto& cart_comp_coef : cart_comp_coef_vec)
            {
                auto cart_comp = cart_comp_coef.first;
                auto cart_coef = cart_comp_coef.second;

                auto cart_ind = cart_comp * _orb_indexes[0] + _orb_indexes[ind];

                sph_cart_p[sph_ind].push_back(std::pair<int64_t, double>({cart_ind, cart_coef}));
            }
        }
    }

    return sph_cart_p;
}

auto
CGtoBlock::getCartesianToSphericalMappingForP() const -> std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>>
{
    errors::assertMsgCritical(_angmom == 1, std::string("GtoBlock: getCartesianToSphericalMappingForP only works for p-orbitals"));

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_p;

    // p-1 (0) <- py (1)
    // p_0 (1) <- pz (2)
    // p+1 (2) <- px (0)

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_comp_map;

    cart_sph_comp_map[0] = std::vector<std::pair<int64_t, double>>({{2, 1.0}});
    cart_sph_comp_map[1] = std::vector<std::pair<int64_t, double>>({{0, 1.0}});
    cart_sph_comp_map[2] = std::vector<std::pair<int64_t, double>>({{1, 1.0}});

    for (const auto& [cart_comp, sph_comp_coef_vec] : cart_sph_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            auto cart_ind = cart_comp * _orb_indexes[0] + _orb_indexes[ind];

            cart_sph_p[cart_ind] = std::vector<std::pair<int64_t, double>>();

            for (const auto& sph_comp_coef : sph_comp_coef_vec)
            {
                auto sph_comp = sph_comp_coef.first;
                auto sph_coef = sph_comp_coef.second;

                auto sph_ind = sph_comp * _orb_indexes[0] + _orb_indexes[ind];

                cart_sph_p[cart_ind].push_back(std::pair<int64_t, double>({sph_ind, sph_coef}));
            }
        }
    }

    return cart_sph_p;
}

auto
CGtoBlock::getSphericalToCartesianMappingForP() const -> std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>>
{
    errors::assertMsgCritical(_angmom == 1, std::string("GtoBlock: getSphericalToCartesianMappingForP only works for p-orbitals"));

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_p;

    // p-1 (0) -> py (1)
    // p_0 (1) -> pz (2)
    // p+1 (2) -> px (0)

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> sph_cart_comp_map;

    sph_cart_comp_map[0] = std::vector<std::pair<int64_t, double>>({{1, 1.0}});
    sph_cart_comp_map[1] = std::vector<std::pair<int64_t, double>>({{2, 1.0}});
    sph_cart_comp_map[2] = std::vector<std::pair<int64_t, double>>({{0, 1.0}});

    for (const auto& [sph_comp, cart_comp_coef_vec] : sph_cart_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            auto sph_ind = sph_comp * _orb_indexes[0] + _orb_indexes[ind];

            sph_cart_p[sph_ind] = std::vector<std::pair<int64_t, double>>();

            for (const auto& cart_comp_coef : cart_comp_coef_vec)
            {
                auto cart_comp = cart_comp_coef.first;
                auto cart_coef = cart_comp_coef.second;

                auto cart_ind = cart_comp * _orb_indexes[0] + _orb_indexes[ind];

                sph_cart_p[sph_ind].push_back(std::pair<int64_t, double>({cart_ind, cart_coef}));
            }
        }
    }

    return sph_cart_p;
}

auto
CGtoBlock::getCartesianToSphericalMappingForD() const -> std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>>
{
    errors::assertMsgCritical(_angmom == 2, std::string("GtoBlock: getCartesianToSphericalMappingForD only works for d-orbitals"));

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_d;

    const double f2_3 = 2.0 * std::sqrt(3.0);

    // d-2 (0) <- dxy(1) * f2_3
    // d-1 (1) <- dyz(4) * f2_3
    // d_0 (2) <- dzz(5) * 2.0 - dxx(0) - dyy(3)
    // d+1 (3) <- dxz(2) * f2_3
    // d+2 (4) <- (dxx(0) - dyy(3)) * 0.5 * f2_3

    std::unordered_map<int64_t, std::vector<std::pair<int64_t, double>>> cart_sph_comp_map;

    cart_sph_comp_map[0] = std::vector<std::pair<int64_t, double>>({{2, -1.0}, {4, 0.5 * f2_3}});
    cart_sph_comp_map[1] = std::vector<std::pair<int64_t, double>>({{0, f2_3}});
    cart_sph_comp_map[2] = std::vector<std::pair<int64_t, double>>({{3, f2_3}});
    cart_sph_comp_map[3] = std::vector<std::pair<int64_t, double>>({{2, -1.0}, {4, -0.5 * f2_3}});
    cart_sph_comp_map[4] = std::vector<std::pair<int64_t, double>>({{1, f2_3}});
    cart_sph_comp_map[5] = std::vector<std::pair<int64_t, double>>({{2, 2.0}});

    for (const auto& [cart_comp, sph_comp_coef_vec] : cart_sph_comp_map)
    {
        // go through CGTOs in this block
        // note that ind starts from 1
        // because _orb_indexes[0] is the total number of CGTOs of _angmom
        // which could be larger than the number of CGTOs in this block

        for (size_t ind = 1; ind < _orb_indexes.size(); ind++)
        {
            auto cart_ind = cart_comp * _orb_indexes[0] + _orb_indexes[ind];

            cart_sph_d[cart_ind] = std::vector<std::pair<int64_t, double>>();

            for (const auto& sph_comp_coef : sph_comp_coef_vec)
            {
                auto sph_comp = sph_comp_coef.first;
                auto sph_coef = sph_comp_coef.second;

                auto sph_ind = sph_comp * _orb_indexes[0] + _orb_indexes[ind];

                cart_sph_d[cart_ind].push_back(std::pair<int64_t, double>({sph_ind, sph_coef}));
            }
        }
    }

    return cart_sph_d;
}

auto
CGtoBlock::getAngularMomentum() const -> int64_t
{
    return _angmom;
}

auto
CGtoBlock::getNumberOfPrimitives() const -> int64_t
{
    return _npgtos;
}

auto
CGtoBlock::getNumberOfBasisFunctions() const -> int64_t
{
    return static_cast<int64_t>(_coordinates.size());
}
