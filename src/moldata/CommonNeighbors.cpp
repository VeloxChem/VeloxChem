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
#include <iterator>
#include <sstream>

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"

CCommonNeighbors::CCommonNeighbors()
    
    :  _cutRadius(0.0)

    , _adjacencies(CMemBlock<int32_t>())

    , _bonds(CDenseMatrix())

    , _signatures(std::vector<CFourIndexes>())

    , _repetitions(std::vector<int32_t>())

    , _idAtomic(std::vector<int32_t>())

    , _composition(std::set<int32_t>())

    , _bondLabels(std::vector<std::string>())
{
    
}

CCommonNeighbors::CCommonNeighbors(const double                     cutRadius,
                                   const CMemBlock<int32_t>&        adjacencies,
                                   const CDenseMatrix&              bonds,
                                   const std::vector<CFourIndexes>& signatures,
                                   const std::vector<int32_t>&      repetitions,
                                   const std::vector<int32_t>&      idAtomic,
                                   const std::set<int32_t>&         composition,
                                   const std::vector<std::string>&  bondLabels)
    :  _cutRadius(cutRadius)

    , _adjacencies(adjacencies)

    , _bonds(bonds)

    , _signatures(signatures)

    , _repetitions(repetitions)

    , _idAtomic(idAtomic)

    , _composition(composition)

    , _bondLabels(bondLabels)
{
    
}

CCommonNeighbors::CCommonNeighbors(const CMolecule& molecule,
                                   const double     cutRadius)
{
    _cutRadius = cutRadius;
    
    _computeBonds(molecule);
    
    _computeAdjacencies();
    
    _setBondLabels();
}

CCommonNeighbors::CCommonNeighbors(const CCommonNeighbors& source)

    : _cutRadius(source._cutRadius)

    , _adjacencies(source._adjacencies)
    
    , _bonds(source._bonds)

    , _signatures(source._signatures)

    , _repetitions(source._repetitions)

    , _idAtomic(source._idAtomic)

    , _composition(source._composition)

    , _bondLabels(source._bondLabels)
{
    
}

CCommonNeighbors::CCommonNeighbors(CCommonNeighbors&& source) noexcept

    : _cutRadius(std::move(source._cutRadius))

    , _adjacencies(std::move(source._adjacencies))
    
    , _bonds(std::move(source._bonds))

    , _signatures(std::move(source._signatures))

    , _repetitions(std::move(source._repetitions))

    , _idAtomic(std::move(source._idAtomic))

    , _composition(std::move(source._composition))

    , _bondLabels(std::move(source._bondLabels))
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
    
    _idAtomic = source._idAtomic;
    
    _composition = source._composition;
    
    _bondLabels = source._bondLabels;

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
    
    _idAtomic = std::move(source._idAtomic);
    
    _composition = std::move(source._composition);
    
    _bondLabels = std::move(source._bondLabels);

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
    
    if (_idAtomic != other._idAtomic) return false;
    
    if (_composition != other._composition) return false;
    
    if (_bondLabels != other._bondLabels) return false;
    
    return true;
}

bool
CCommonNeighbors::operator!=(const CCommonNeighbors& other) const
{
    return !(*this == other);
}

#include <iostream>

void
CCommonNeighbors::generate(const double radius)
{
    _signatures.clear();
    
    _repetitions.clear();
    
    for (const auto& tidx : getBondPairs())
    {
        const auto idbond = _getBondIdentifier(tidx);
        
        const auto catoms = _getCommonAtoms(tidx, radius);
        
        const auto cbonds = _getNumberOfCommonBonds(catoms);
        
        const auto mcbond = _getLongestCommonBond(catoms);
        
        _add(CFourIndexes(idbond, static_cast<int32_t>(catoms.size()), cbonds, mcbond));
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

std::vector<CFourIndexes>
CCommonNeighbors::getSignatures() const
{
    return _signatures;
}

std::vector<int32_t>
CCommonNeighbors::getRepetitions() const
{
    return _repetitions;
}

double
CCommonNeighbors::compJaccardIndex(const CCommonNeighbors& other)
{
    // compute intersection
    
    std::vector<CFourIndexes> ivecab;
    
    for (const auto& tval : _signatures)
    {
        if (other.find(tval) != -1)
        {
            ivecab.push_back(tval);
        }
    }
    
    int32_t ifact = 0;
    
    for (const auto& tval : ivecab)
    {
        const auto repa = _repetitions[find(tval)];
        
        const auto repb = other._repetitions[other.find(tval)];
        
        ifact += (repa > repb) ? repb : repa;
    }
    
    // compute union
    
    int32_t ufact = 0;
    
    for (size_t i = 0; i < _signatures.size(); i++)
    {
        const auto idx = other.find(_signatures[i]);
        
        const auto repa = _repetitions[i];
        
        if (idx == -1)
        {
            ufact += repa;
        }
        else
        {
            const auto repb = other._repetitions[idx];
            
            ufact += (repa > repb) ? repa : repb;
        }
    }
    
    for (int32_t i = 0; i < other._signatures.size(); i++)
    {
        const auto idx = find(other._signatures[i]);
        
        if (idx == - 1)
        {
            ufact += other._repetitions[i];
        }
    }
    
    return  static_cast<double>(ifact) / static_cast<double>(ufact);
}

int32_t
CCommonNeighbors::find(const CFourIndexes& signature) const
{
    for (size_t i = 0; i < _signatures.size(); i++)
    {
        if (_signatures[i] == signature)
        {
            return static_cast<int32_t>(i);
        }
    }
    
    return -1;
}

std::string
CCommonNeighbors::getSignaturesRepr() const
{
    std::stringstream ss;
    
    for (size_t i = 0; i < _signatures.size(); i++)
    {
        const auto tsig = _signatures[i];
        
        ss << fstr::format(_bondLabels[tsig.first()], 16, fmt::left);
        
        auto str = "(" + std::to_string(tsig.second())
                 
                 + "," + std::to_string(tsig.third())
                
                 + "," + std::to_string(tsig.fourth()) + ")";
        
        ss << fstr::format(str, 18, fmt::left);
        
        ss << fstr::format(std::to_string(_repetitions[i]), 8, fmt::left);
    
        ss << "\n";
    }
    
    return ss.str();
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
    
    // set up pointer to atomic identifiers
    
    auto eleids = molecule.getIdsElemental();
    
    for (int32_t i = 0; i < natoms; i++)
    {
        // atomic identifiers
        
        _idAtomic.push_back(eleids[i]);
        
        // compute bond distances
        
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
    
    // set up elemental composition
    
    _composition.insert(_idAtomic.begin(), _idAtomic.end());
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

int32_t
CCommonNeighbors::_getNumberOfCommonBonds(const std::vector<int32_t>& atoms)
{
    const auto natoms = _bonds.getNumberOfRows();
    
    auto edges = _adjacencies.data();
    
    const auto ndim = atoms.size();
    
    int32_t nbonds = 0;
    
    for (size_t i = 0; i < ndim; i++)
    {
        for (size_t j = i + 1; j < ndim; j++)
        {
            if (edges[atoms[i] * natoms + atoms[j]] == 1) nbonds++;
        }
    }
    
    return nbonds; 
}

int32_t
CCommonNeighbors::_getLongestCommonBond(const std::vector<int32_t>& atoms)
{
    const auto natoms = _bonds.getNumberOfRows();
    
    auto edges = _adjacencies.data();
    
    const auto ndim = atoms.size();
    
    // select all bonded atoms
    
    std::set<int32_t> bpairs;
    
    for (size_t i = 0; i < ndim; i++)
    {
        for (size_t j = i + 1; j < ndim; j++)
        {
            if (edges[atoms[i] * natoms + atoms[j]] == 1)
            {
                bpairs.insert(atoms[i]);
                
                bpairs.insert(atoms[j]);
            }
        }
    }
    
    // compute longest bond
    
    int32_t nbonds = 0;
    
    if (bpairs.size() > 1)
    {
        std::vector<int32_t> bvec(bpairs.begin(), bpairs.end());
        
        int32_t mbonds = static_cast<int32_t>(bpairs.size() - 1);
        
        do {
            auto cbonds = _getLongestBondForPath(bvec);
            
            if (cbonds > nbonds)  nbonds = cbonds;
            
            if (nbonds == mbonds) break;
            
        } while (std::next_permutation(bvec.begin(), bvec.end()));
    }
    
    return nbonds;
}

int32_t
CCommonNeighbors::_getLongestBondForPath(const std::vector<int32_t>& path)
{
    const auto natoms = _bonds.getNumberOfRows();
    
    auto edges = _adjacencies.data();
    
    const auto ndim = path.size() - 1;
    
    int32_t mbonds = 0;
    
    int32_t cbonds = 0;
    
    for (size_t i = 0; i < ndim; i++)
    {
        if (edges[path[i] * natoms + path[i + 1]] == 1)
        {
            cbonds++;
        }
        else
        {
            if (cbonds > mbonds) mbonds = cbonds;
            
            cbonds = 0;
        }
    }
    
    return mbonds;
}

int32_t
CCommonNeighbors::_getBondIdentifier(const CTwoIndexes& atomsPair)
{
    const auto ndim = _composition.size();
        
    auto it = _composition.find(_idAtomic[atomsPair.first()]);
    
    const auto idxa = std::distance(_composition.begin(), it);
    
    it = _composition.find(_idAtomic[atomsPair.second()]);
    
    const auto idxb = std::distance(_composition.begin(), it);
    
    if (idxa > idxb)
    {
        return static_cast<int32_t>(idxb * ndim + idxa);
    }
    else
    {
        return static_cast<int32_t>(idxa * ndim + idxb);
    }
}

void
CCommonNeighbors::_add(const CFourIndexes& signature)
{
    const auto idx = find(signature);
    
    if (idx == -1)
    {
        _signatures.push_back(signature);
        
        _repetitions.push_back(1);
    }
    else
    {
        _repetitions[idx] = _repetitions[idx] + 1;
    }
}

void
CCommonNeighbors::_setBondLabels()
{
    auto bele = CChemicalElement();
    
    auto kele = CChemicalElement();
    
    for (const auto& bval : _composition)
    {
        bele.setAtomType(bval);
        
        for (const auto& kval : _composition)
        {
            kele.setAtomType(kval);
            
            const auto str = bele.getName() + "-"
                            
                           + kele.getName();
            
            _bondLabels.push_back(str);
        }
    }
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

    output << "_idAtomic: " << source._idAtomic.size() << std::endl;
    
    for (size_t i = 0; i < source._idAtomic.size(); i++)
    {
        output << source._idAtomic[i] << std::endl;
    }
    
    output << "_composition: " << source._composition.size() << std::endl;
       
    for (const auto& tval : source._composition)
    {
           output << tval << std::endl;
    }
    
    output << "_bondLabels: " << source._bondLabels.size() << std::endl;
       
    for (const auto& tval : source._bondLabels)
    {
           output << tval << std::endl;
    }
    
    return output;
}
