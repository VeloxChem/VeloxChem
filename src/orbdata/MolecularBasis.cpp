//
//                           VELOXCHEM 1.0-RC
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

#include "MolecularBasis.hpp"

#include <array>
#include <sstream>

#include <mpi.h>

#include "AngularMomentum.hpp"
#include "ChemicalElement.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

CMolecularBasis::CMolecularBasis()

    : _maxAngularMomentum(-1)
{
}

CMolecularBasis::CMolecularBasis(const CMolecularBasis& source)

    : _atomicBasisSets(source._atomicBasisSets)

    , _maxAngularMomentum(source._maxAngularMomentum)

    , _label(source._label)
{
}

CMolecularBasis::CMolecularBasis(CMolecularBasis&& source) noexcept

    : _atomicBasisSets(std::move(source._atomicBasisSets))

    , _maxAngularMomentum(std::move(source._maxAngularMomentum))

    , _label(std::move(source._label))
{
}

CMolecularBasis::~CMolecularBasis()
{
}

CMolecularBasis&
CMolecularBasis::operator=(const CMolecularBasis& source)
{
    if (this == &source) return *this;

    _atomicBasisSets = source._atomicBasisSets;

    _maxAngularMomentum = source._maxAngularMomentum;

    _label = source._label;

    return *this;
}

CMolecularBasis&
CMolecularBasis::operator=(CMolecularBasis&& source) noexcept
{
    if (this == &source) return *this;

    _atomicBasisSets = std::move(source._atomicBasisSets);

    _maxAngularMomentum = std::move(source._maxAngularMomentum);

    _label = std::move(source._label);

    return *this;
}

bool
CMolecularBasis::operator==(const CMolecularBasis& other) const
{
    if (_atomicBasisSets.size() != other._atomicBasisSets.size()) return false;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        if (_atomicBasisSets[i] != other._atomicBasisSets[i]) return false;
    }

    if (_maxAngularMomentum != other._maxAngularMomentum) return false;

    if (_label != other._label) return false;

    return true;
}

bool
CMolecularBasis::operator!=(const CMolecularBasis& other) const
{
    return !(*this == other);
}

void
CMolecularBasis::setMaxAngularMomentum(const int32_t maxAngularMomentum)
{
    _maxAngularMomentum = maxAngularMomentum;
}

void
CMolecularBasis::setLabel(const std::string& label)
{
    _label = label;
}

void
CMolecularBasis::addAtomBasis(const CAtomBasis& atomBasis)
{
    _atomicBasisSets.push_back(atomBasis);

    auto mAngularMomentum = atomBasis.getMaxAngularMomentum();

    if (mAngularMomentum > _maxAngularMomentum)
    {
        _maxAngularMomentum = mAngularMomentum;
    }
}

CMolecularBasis
CMolecularBasis::reduceToValenceBasis() const
{
    CMolecularBasis molbas;

    auto strlbl = _label;

    strlbl.append("(VAL)");

    molbas.setLabel(strlbl);

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        molbas.addAtomBasis(_atomicBasisSets[i].reduceToValenceBasis());
    }

    return molbas;
}

int32_t
CMolecularBasis::getMaxAngularMomentum() const
{
    return _maxAngularMomentum;
}

int32_t
CMolecularBasis::getMaxAngularMomentum(const int32_t idElemental) const
{
    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        if (_atomicBasisSets[i].getIdElemental() == idElemental)
        {
            return _atomicBasisSets[i].getMaxAngularMomentum();
        }
    }

    return -1;
}

int32_t
CMolecularBasis::getMolecularMaxAngularMomentum(const CMolecule& molecule) const
{
    int32_t max_angl = 0;

    auto elmlst = molecule.getElementalComposition();

    for (auto i = elmlst.cbegin(); i != elmlst.cend(); ++i)
    {
        int32_t elem_angl = getMaxAngularMomentum(*i);

        if (max_angl < elem_angl)
        {
            max_angl = elem_angl;
        }
    }

    return max_angl;
}

std::string
CMolecularBasis::getLabel() const
{
    return _label;
}

int32_t
CMolecularBasis::getNumberOfBasisFunctions(const int32_t idElemental, const int32_t angularMomentum) const
{
    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        if (_atomicBasisSets[i].getIdElemental() == idElemental)
        {
            return _atomicBasisSets[i].getNumberOfBasisFunctions(angularMomentum);
        }
    }

    return 0;
}

int32_t
CMolecularBasis::getNumberOfBasisFunctions(const CMolecule& molecule, const int32_t angularMomentum) const
{
    return getNumberOfBasisFunctions(molecule, 0, molecule.getNumberOfAtoms(), angularMomentum);
}

int32_t
CMolecularBasis::getNumberOfBasisFunctions(const CMolecule& molecule, const int32_t iAtom, const int32_t nAtoms, const int32_t angularMomentum) const
{
    int32_t nbfuncs = 0;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        nbfuncs += molecule.getNumberOfAtoms(iAtom, nAtoms, _atomicBasisSets[i].getIdElemental())

                   * _atomicBasisSets[i].getNumberOfBasisFunctions(angularMomentum);
    }

    return nbfuncs;
}

int32_t
CMolecularBasis::getNumberOfPrimitiveBasisFunctions(const CMolecule& molecule, const int32_t angularMomentum) const
{
    return getNumberOfPrimitiveBasisFunctions(molecule, 0, molecule.getNumberOfAtoms(), angularMomentum);
}

int32_t
CMolecularBasis::getNumberOfPrimitiveBasisFunctions(const CMolecule& molecule,
                                                    const int32_t    iAtom,
                                                    const int32_t    nAtoms,
                                                    const int32_t    angularMomentum) const
{
    int32_t npfuncs = 0;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        npfuncs += molecule.getNumberOfAtoms(iAtom, nAtoms, _atomicBasisSets[i].getIdElemental())

                   * _atomicBasisSets[i].getNumberOfPrimitiveFunctions(angularMomentum);
    }

    return npfuncs;
}

int32_t
CMolecularBasis::getDimensionsOfBasis(const CMolecule& molecule) const
{
    int32_t ndim = 0;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        auto idelem = _atomicBasisSets[i].getIdElemental();

        auto natoms = molecule.getNumberOfAtoms(idelem);

        auto mang = _atomicBasisSets[i].getMaxAngularMomentum();

        for (int32_t j = 0; j <= mang; j++)
        {
            ndim += natoms * angmom::to_SphericalComponents(j)

                    * _atomicBasisSets[i].getNumberOfBasisFunctions(j);
        }
    }

    return ndim;
}

int32_t
CMolecularBasis::getPartialDimensionsOfBasis(const CMolecule& molecule, const int32_t angularMomentum) const
{
    int32_t ndim = 0;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        auto idelem = _atomicBasisSets[i].getIdElemental();

        auto natoms = molecule.getNumberOfAtoms(idelem);

        for (int32_t j = 0; j < angularMomentum; j++)
        {
            ndim += natoms * angmom::to_SphericalComponents(j)

                    * _atomicBasisSets[i].getNumberOfBasisFunctions(j);
        }
    }

    return ndim;
}

int32_t
CMolecularBasis::getDimensionsOfPrimitiveBasis(const CMolecule& molecule) const
{
    int32_t ndim = 0;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        auto idelem = _atomicBasisSets[i].getIdElemental();

        auto natoms = molecule.getNumberOfAtoms(idelem);

        auto mang = _atomicBasisSets[i].getMaxAngularMomentum();

        for (int32_t j = 0; j <= mang; j++)
        {
            ndim += natoms * angmom::to_SphericalComponents(j)

                    * _atomicBasisSets[i].getNumberOfPrimitiveFunctions(j);
        }
    }

    return ndim;
}

CAtomBasis
CMolecularBasis::getAtomBasis(const int32_t idElemental) const
{
    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        if (_atomicBasisSets[i].getIdElemental() == idElemental)
        {
            return _atomicBasisSets[i];
        }
    }

    return CAtomBasis();
}

std::vector<CBasisFunction>
CMolecularBasis::getBasisFunctions(const int32_t idElemental, const int32_t angularMomentum) const
{
    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        if (_atomicBasisSets[i].getIdElemental() == idElemental)
        {
            return _atomicBasisSets[i].getBasisFunctions(angularMomentum);
        }
    }

    return std::vector<CBasisFunction>();
}

std::vector<std::string>
CMolecularBasis::getAOBasisMap(const CMolecule& molecule) const
{
    std::vector<std::string> strmap;

    auto natoms = molecule.getNumberOfAtoms();

    auto idselm = molecule.getIdsElemental();

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        for (int32_t j = 0; j < angmom::to_SphericalComponents(i); j++)
        {
            for (int32_t k = 0; k < natoms; k++)
            {
                auto gtos = getBasisFunctions(idselm[k], i);

                auto ngtos = static_cast<int32_t>(gtos.size());

                for (int32_t l = 0; l < ngtos; l++)
                {
                    std::stringstream st;

                    st.setf(std::ios::fixed);

                    st.width(4);

                    st << k + 1;

                    st << " ";

                    auto lbl = molecule.getLabel(k);

                    st << lbl;

                    if (lbl.size() == 1) st << " ";

                    st << " ";

                    st.setf(std::ios::fixed);

                    st.width(2);

                    st << l + 1;

                    st << angmom::getStringOfAngularMomentum(i, j);

                    strmap.push_back(st.str());
                }
            }
        }
    }

    return strmap;
}

CMemBlock2D<int32_t>
CMolecularBasis::getIndexMapForDalton(const CMolecule &molecule) const
{
    // allocates index map
    
    auto natoms = molecule.getNumberOfAtoms();
    
    CMemBlock2D<int32_t> idsmap(getDimensionsOfBasis(molecule), 2);
    
    // set up pointers to index map
    
    auto vlxidx = idsmap.data(0);
    
    auto dalidx = idsmap.data(1);
    
    // set up pointer to chemical element data
    
    auto idselm = molecule.getIdsElemental();
    
    // indexing map for p-functions
    
    std::array<int32_t, 3> pord({2, 0, 1});
    
    // loop over atoms
    
    int32_t curgto = 0;
    
    for (int32_t i = 0; i < natoms; i++)
    {
        auto atmbas = getAtomBasis(idselm[i]);
        
        for (int32_t j = 0; j <= atmbas.getMaxAngularMomentum(); j++)
        {
            auto gtopos = getPartialDimensionsOfBasis(molecule, j);
            
            auto blkpos = getPositionInAngularBlock(molecule, i, j);
            
            auto ngtos = getNumberOfBasisFunctions(molecule, j); 
            
            for (int32_t k = 0; k < atmbas.getNumberOfBasisFunctions(j); k++)
            {
                for (int32_t l = 0; l < angmom::to_SphericalComponents(j); l++)
                {
                    dalidx[curgto] = curgto;
                    
                    if (j == 1)
                    {
                        vlxidx[curgto] = gtopos + ngtos * pord[l] + blkpos + k;
                    }
                    else
                    {
                        vlxidx[curgto] = gtopos + ngtos * l + blkpos + k;
                    }
                    
                    curgto++;
                }
            }
        }
    }
    
    return idsmap;
}

int32_t
CMolecularBasis::getPositionInAngularBlock(const CMolecule& molecule,
                                           const int32_t    iAtom,
                                           const int32_t    angularMomentum) const
{
    if ((iAtom < molecule.getNumberOfAtoms()) && (angularMomentum <= getMaxAngularMomentum()))
    {
        auto idselm = molecule.getIdsElemental();
        
        int32_t bfpos = 0;
        
        for (int32_t i = 0; i < iAtom; i++)
        {
            bfpos += getNumberOfBasisFunctions(idselm[i], angularMomentum);
        }
        
        return bfpos; 
    }

    return -1;
}

std::string
CMolecularBasis::printBasis(const std::string& title, const CMolecule& molecule) const
{
    std::string str = "Molecular Basis (" + title + ")";

    std::stringstream ss;

    ss << str << "\n";

    ss << std::string(str.size() + 2, '=') << "\n\n";

    str.assign("Basis: ");

    str.append(_label);

    ss << fstr::format(str, 54, fmt::left) << "\n\n";

    ss << "  Atom ";

    ss << fstr::format(std::string("Contracted GTOs"), 25, fmt::left);

    ss << fstr::format(std::string("Primitive GTOs"), 25, fmt::left);

    ss << "\n\n";

    for (auto i = _atomicBasisSets.cbegin(); i != _atomicBasisSets.cend(); ++i)
    {
        std::string lbl("  ");

        CChemicalElement ce;

        ce.setAtomType(i->getIdElemental());

        lbl.append(ce.getName());

        ss << fstr::format(lbl, 6, fmt::left);

        ss << fstr::format(i->getContractionString(), 25, fmt::left);

        ss << fstr::format(i->getPrimitivesString(), 25, fmt::left);

        ss << "\n";
    }

    ss << "\n";

    str.assign("Contracted Basis Functions : ");

    str.append(std::to_string(getDimensionsOfBasis(molecule)));

    ss << fstr::format(str, 54, fmt::left) << "\n";

    str.assign("Primitive Basis Functions  : ");

    str.append(std::to_string(getDimensionsOfPrimitiveBasis(molecule)));

    ss << fstr::format(str, 54, fmt::left) << "\n";

    ss << "\n";

    return ss.str();
}

void
CMolecularBasis::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_maxAngularMomentum, comm);

        mpi::bcast(_label, rank, comm);

        int32_t natombases = static_cast<int32_t>(_atomicBasisSets.size());

        mpi::bcast(natombases, comm);

        for (int32_t i = 0; i < natombases; i++)
        {
            CAtomBasis atmbasis;

            if (rank == mpi::master()) atmbasis = _atomicBasisSets[i];

            atmbasis.broadcast(rank, comm);

            if (rank != mpi::master()) addAtomBasis(atmbasis);

            MPI_Barrier(comm);
        }
    }
}

std::string
CMolecularBasis::repr() const
{
    std::ostringstream os;

    os << std::endl;

    os << "[CMolecularBasis (Object):" << this << "]" << std::endl;

    os << "_label: " << _label << std::endl;

    os << "_maxAngularMomentum: " << _maxAngularMomentum;

    os << std::endl;

    os << "_atomicBasisSets: " << std::endl;

    for (size_t i = 0; i < _atomicBasisSets.size(); i++)
    {
        os << "_atomicBasisSets[" << i << "]: " << std::endl;

        os << _atomicBasisSets[i] << std::endl;
    }

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CMolecularBasis& source)
{
    return (output << source.repr());
}
