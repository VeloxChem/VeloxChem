//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#include "CommonNeighbors.hpp"

#include <cmath>
#include <set>
#include <algorithm>

#include <iostream>

CCommonNeighbors::CCommonNeighbors()
    
    :  _cutRadius(0.0)

    , _adjacencies(CMemBlock<int32_t>())

    , _bonds(CDenseMatrix())

    , _signatures(std::vector<CThreeIndexes>())

    , _repetitions(std::vector<int32_t>())
{
    
}

CCommonNeighbors::CCommonNeighbors(const double                      cutRadius,
                                   const CMemBlock<int32_t>&         adjacencies,
                                   const CDenseMatrix&               bonds,
                                   const std::vector<CThreeIndexes>& signatures,
                                   const std::vector<int32_t>&       repetitions)
    :  _cutRadius(cutRadius)

    , _adjacencies(adjacencies)

    , _bonds(bonds)

    , _signatures(signatures)

    , _repetitions(repetitions)
{
    
}

CCommonNeighbors::CCommonNeighbors(const CMolecule& molecule,
                                   const double     cutRadius)
{
    _cutRadius = cutRadius;
    
    _computeBonds(molecule);
    
    _computeAdjacencies();
}

CCommonNeighbors::CCommonNeighbors(const CCommonNeighbors& source)

    : _cutRadius(source._cutRadius)

    , _adjacencies(source._adjacencies)
    
    , _bonds(source._bonds)

    , _signatures(source._signatures)

    , _repetitions(source._repetitions)
{
    
}

CCommonNeighbors::CCommonNeighbors(CCommonNeighbors&& source) noexcept

    : _cutRadius(std::move(source._cutRadius))

    , _adjacencies(std::move(source._adjacencies))
    
    , _bonds(std::move(source._bonds))

    , _signatures(std::move(source._signatures))

    , _repetitions(std::move(source._repetitions))
{
    
}

CCommonNeighbors::~CCommonNeighbors()
{
    
}

CCommonNeighbors&
CCommonNeighbors::operator=(const CCommonNeighbors& source)
{
    if (this == &source) return *this;

    _cutRadius = source._cutRadius;

    _adjacencies = source._adjacencies;
       
    _bonds = source._bonds;

    _signatures = source._signatures;

    _repetitions = source._repetitions;

    return *this;
}

CCommonNeighbors&
CCommonNeighbors::operator=(CCommonNeighbors&& source) noexcept
{
    if (this == &source) return *this;

    _cutRadius = std::move(source._cutRadius);

    _adjacencies = std::move(source._adjacencies);
    
    _bonds = std::move(source._bonds);

    _signatures = std::move(source._signatures);

    _repetitions = std::move(source._repetitions);

    return *this;
}

bool
CCommonNeighbors::operator==(const CCommonNeighbors& other) const
{
    if (std::fabs(_cutRadius - other._cutRadius) > 1.0e-13) return false;
    
    if (_adjacencies != other._adjacencies) return false;
    
    if (_bonds != other._bonds) return false;

    if (_signatures != other._signatures) return false;
   
    if (_repetitions != other._repetitions) return false;
    
    return true;
}

bool
CCommonNeighbors::operator!=(const CCommonNeighbors& other) const
{
    return !(*this == other);
}

void
CCommonNeighbors::generate(const double radius)
{
    for (const auto& tidx : getBondPairs())
    {
        const auto atoms = _getCommonAtoms(tidx, radius);
        
        std::cout << "index: " << tidx.first() << ":"  << tidx.second() << ": " << std::endl;
        
        for (size_t i = 0; i < atoms.size(); i++)
        {
            std::cout << atoms[i] << " ";
        }
        
        std::cout << std::endl;
    }
}

CDenseMatrix
CCommonNeighbors::getBonds() const
{
    return _bonds;
}

CMemBlock<int32_t>
CCommonNeighbors::getAdjacencies() const
{
    return _adjacencies;
}

std::vector<CTwoIndexes>
CCommonNeighbors::getBondPairs() const
{
    std::vector<CTwoIndexes> bpairs;
    
    const auto natoms = _bonds.getNumberOfRows();
    
    auto edges = _adjacencies.data();
    
    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            if (edges[i * natoms +  j] == 1)
            {
                bpairs.push_back(CTwoIndexes(i, j));
            }
        }
    }
    
    return bpairs;
}

void
CCommonNeighbors::_computeBonds(const CMolecule& molecule)
{
    const auto natoms = molecule.getNumberOfAtoms();
    
    _bonds = CDenseMatrix(natoms);
    
    // set up pointers to molecular coordinates
    
    auto rx = molecule.getCoordinatesX();
    
    auto ry = molecule.getCoordinatesY();
    
    auto rz = molecule.getCoordinatesZ();
    
    // compute pairwise atomic distances
    
    auto rab = _bonds.values();
    
    for (int32_t i = 0; i < natoms; i++)
    {
        const auto ax = rx[i];
        
        const auto ay = ry[i];
        
        const auto az = rz[i];
        
        rab[i * natoms + i] = 0.0;
       
        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto frad = std::sqrt((ax - rx[j]) * (ax - rx[j]) +
                               
                                  (ay - ry[j]) * (ay - ry[j]) +
                               
                                  (az - rz[j]) * (az - rz[j]));
           
            rab[i * natoms + j] = frad;
            
            rab[j * natoms + i] = frad;
        }
    }
}

void
CCommonNeighbors::_computeAdjacencies()
{
    const auto natoms = _bonds.getNumberOfRows();
    
    auto rab = _bonds.values();
    
    _adjacencies = CMemBlock<int32_t>(natoms * natoms);
    
    auto edges = _adjacencies.data();
    
    for (int32_t i = 0; i < natoms; i++)
    {
        edges[i * natoms + i]  = 0;
        
        for (int32_t j = i + 1; j < natoms; j++)
        {
            int32_t fact = (rab[i * natoms + j] <= _cutRadius) ? 1 : 0;
            
            edges[i * natoms +  j] = fact;
            
            edges[j * natoms +  i] = fact;
        }
    }
}

std::vector<int32_t>
CCommonNeighbors::_getCommonAtoms(const CTwoIndexes& atomsPair,
                                  const double       radius)
{
    const auto natoms = _bonds.getNumberOfRows();
    
    // i-th atom environment
    
    const auto iatom = atomsPair.first();
    
    auto rab = _bonds.row(iatom);
    
    std::set<int32_t> idxa;
    
    for (int32_t i = 0; i < natoms; i++)
    {
        if (i != iatom)
        {
            if (rab[i] < radius) idxa.insert(i);
        }
    }
    
    // j-th atom environment
    
    const auto jatom = atomsPair.second();
    
    rab = _bonds.row(jatom);
    
    std::set<int32_t> idxb;
    
    for (int32_t i = 0; i < natoms; i++)
    {
        if (i != jatom)
        {
            if (rab[i] < radius) idxb.insert(i);
        }
    }
    
    std::vector<int32_t> atoms(natoms);
    
    auto it=std::set_intersection(idxa.begin(), idxa.end(),
                                  idxb.begin(), idxb.end(),
                                  atoms.begin());
    
    atoms.resize(it - atoms.begin());
    
    return atoms;
}

std::ostream&
operator<<(      std::ostream&     output,
           const CCommonNeighbors& source)
{
    output << std::endl;

    output << "[CCommonNeighbors (Object):" << &source << "]" << std::endl;
    
    output << "_cutRadius: " << source._cutRadius << std::endl;

    output << "_adjacencies: " << source._adjacencies << std::endl;
    
    output << "_bonds: " << source._bonds << std::endl;
    
    output << "_signatures: " << source._signatures.size() << std::endl;
    
    for (size_t i = 0; i < source._signatures.size(); i++)
    {
        output << source._signatures[i] << std::endl;
    }
    
    output << "_repetitions: " << source._repetitions.size() << std::endl;
    
    for (size_t i = 0; i < source._repetitions.size(); i++)
    {
        output << source._repetitions[i] << std::endl;
    }

    return output;
}

