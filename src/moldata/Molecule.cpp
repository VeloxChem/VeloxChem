//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "Molecule.hpp"

#include <cmath>

#include "StringFormat.hpp"
#include "Codata.hpp"
#include "MathFunc.hpp"

CMolecule::CMolecule()

    : _charge(0.0)

    , _multiplicity(1)
{

}

CMolecule::CMolecule(const std::vector<double>&      atomCoordinates,
                     const std::vector<double>&      atomCharges,
                     const std::vector<double>&      atomMasses,
                     const std::vector<std::string>& atomLabels,
                     const std::vector<int32_t>&     idsElemental)

    : _charge(0.0)

    , _multiplicity(1)
{
    auto natoms = static_cast<int32_t>(idsElemental.size());
    
    // set up atom's properties
    
    _atomCoordinates = CMemBlock2D<double>(atomCoordinates, natoms, 3);
    
    _atomCharges = CMemBlock<double>(atomCharges);
    
    _atomMasses = CMemBlock<double>(atomMasses);
    
    _atomLabels = atomLabels;
    
    _idsElemental = CMemBlock<int32_t>(idsElemental);
    
    // set up default indexing of atoms in molecule
    
    setAtomicIndexes(0);
}

CMolecule::CMolecule(const CMolecule& source)

    : _charge(source._charge)

    , _multiplicity(source._multiplicity)

    , _atomCoordinates(source._atomCoordinates)

    , _atomCharges(source._atomCharges)

    , _atomMasses(source._atomMasses)

    , _atomLabels(source._atomLabels)

    , _idsAtomic(source._idsAtomic)

    , _idsElemental(source._idsElemental)
{

}

CMolecule::CMolecule(CMolecule&& source) noexcept

    : _charge(std::move(source._charge))

    , _multiplicity(std::move(source._multiplicity))

    , _atomCoordinates(std::move(source._atomCoordinates))

    , _atomCharges(std::move(source._atomCharges))

    , _atomMasses(std::move(source._atomMasses))

    , _atomLabels(std::move(source._atomLabels))

    , _idsAtomic(std::move(source._idsAtomic))

    , _idsElemental(std::move(source._idsElemental))
{

}

CMolecule::~CMolecule()
{

}

CMolecule
CMolecule::getSubMolecule(int32_t start_index, int32_t num_atoms)
{
    std::vector<double>      atomCoordinates;
    std::vector<double>      atomCharges;
    std::vector<double>      atomMasses;
    std::vector<std::string> atomLabels;
    std::vector<int32_t>     idsElemental;

    for (int i = start_index; i < start_index + num_atoms; i++) {
        atomCoordinates.push_back(_atomCoordinates.data(0)[i]);
    }
    for (int i = start_index; i < start_index + num_atoms; i++) {
        atomCoordinates.push_back(_atomCoordinates.data(1)[i]);
    }
    for (int i = start_index; i < start_index + num_atoms; i++) {
        atomCoordinates.push_back(_atomCoordinates.data(2)[i]);
    }

    for (int i = start_index; i < start_index + num_atoms; i++) {
        atomCharges.push_back(_atomCharges.data()[i]);
        atomMasses.push_back(_atomMasses.data()[i]);
        atomLabels.push_back(_atomLabels.data()[i]);
        idsElemental.push_back(_idsElemental.data()[i]);
    }

    return CMolecule(atomCoordinates,
                     atomCharges,
                     atomMasses,
                     atomLabels,
                     idsElemental);
}

CMolecule&
CMolecule::operator=(const CMolecule& source)
{
    if (this == &source) return *this;

    _charge = source._charge;

    _multiplicity = source._multiplicity;

    _atomCoordinates = source._atomCoordinates;
    
    _atomCharges = source._atomCharges;
    
    _atomMasses = source._atomMasses;
    
    _atomLabels = source._atomLabels;
    
    _idsAtomic = source._idsAtomic;
    
    _idsElemental = source._idsElemental;

    return *this;
}

CMolecule&
CMolecule::operator=(CMolecule&& source) noexcept
{
    if (this == &source) return *this;

    _charge = std::move(source._charge);

    _multiplicity = std::move(source._multiplicity);

    _atomCoordinates = std::move(source._atomCoordinates);
    
    _atomCharges = std::move(source._atomCharges);
    
    _atomMasses = std::move(source._atomMasses);
    
    _atomLabels = std::move(source._atomLabels);
    
    _idsAtomic = std::move(source._idsAtomic);
    
    _idsElemental = std::move(source._idsElemental);

    return *this;
}

bool
CMolecule::operator==(const CMolecule& other) const
{
    if (std::fabs(_charge - other._charge) > 1.0e-13) return false;

    if (_multiplicity != other._multiplicity) return false;

    if (_atomCoordinates != other._atomCoordinates) return false;
    
    if (_atomCharges != other._atomCharges) return false;
    
    if (_atomMasses != other._atomMasses) return false;
    
    if (_atomLabels.size() != other._atomLabels.size()) return false;
    
    for (size_t i = 0; i < _atomLabels.size(); i++)
    {
        if (_atomLabels[i] != other._atomLabels[i]) return false;
    }
    
    if (_idsAtomic != other._idsAtomic) return false;
    
    if (_idsElemental != other._idsElemental) return false;
   
    return true;
}

bool
CMolecule::operator!=(const CMolecule& other) const
{
    return !(*this == other);
}

void
CMolecule::setAtomicIndexes(const int32_t startIndex)
{
    for (int32_t i = 0; i < _idsAtomic.size(); i++)
    {
        _idsAtomic.at(i) = startIndex + i;
    }
}

void
CMolecule::setCharge(const double charge)
{
    _charge = charge;
}

void
CMolecule::setMultiplicity(const int32_t multiplicity)
{
    _multiplicity = multiplicity; 
}

double
CMolecule::getCharge() const
{
    return _charge;
}

int32_t
CMolecule::getMultiplicity() const
{
    return _multiplicity;
}

int32_t
CMolecule::getNumberOfAtoms() const
{
    return _idsElemental.size();
}

int32_t
CMolecule::getNumberOfAtoms(const int32_t idElemental) const
{
    return getNumberOfAtoms(0, _idsElemental.size(), idElemental);
}

int32_t
CMolecule::getNumberOfAtoms(const int32_t iAtom,
                            const int32_t nAtoms,
                            const int32_t idElemental) const
{
    int32_t natoms = 0;
    
    for (int32_t i = iAtom; i < (iAtom + nAtoms); i++)
    {
        if (_idsElemental.at(i) == idElemental) natoms++;
    }
    
    return natoms;
}

std::set<int32_t>
CMolecule::getElementalComposition() const
{
    std::set<int32_t> elemset;
    
    for (int32_t i = 0; i < _idsElemental.size(); i++)
    {
        elemset.insert(_idsElemental.at(i));
    }
    
    return elemset;
}

int32_t
CMolecule::getNumberOfElectrons() const
{
    double nelectrons = -_charge;
    
    for (int32_t i = 0; i < _atomCharges.size(); i++)
    {
        nelectrons += _atomCharges.at(i);
    }
    
    return static_cast<int32_t>(nelectrons);
}

const int32_t*
CMolecule::getIdsElemental() const
{
    return _idsElemental.data();
}

const double*
CMolecule::getCoordinatesX() const
{
    return _atomCoordinates.data(0);
}

const double*
CMolecule::getCoordinatesY() const
{
    return _atomCoordinates.data(1);
}

const double*
CMolecule::getCoordinatesZ() const
{
    return _atomCoordinates.data(2);
}

CMemBlock2D<double>
CMolecule::getCoordinates() const
{
    return _atomCoordinates;
}

CMemBlock<double>
CMolecule::getCharges() const
{
    return _atomCharges; 
}

CMemBlock<double>
CMolecule::getMinDistances() const
{
    // allocate and initialize distances
    
    auto natoms = getNumberOfAtoms();
    
    CMemBlock<double> mdists(natoms);
    
    auto rmin = mdists.data();
    
    mathfunc::set_to(rmin, 1.0e24, natoms);
    
    // set pointers to coordinates
    
    auto rx = getCoordinatesX();
    
    auto ry = getCoordinatesY();
    
    auto rz = getCoordinatesZ();
    
    // determine distance to closest neighbouring atom
    
    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto rab = mathfunc::distance(rx[i], ry[i], rz[i],
                                          rx[j], ry[j], rz[j]);
            
            if (rab < rmin[i]) rmin[i] = rab;
            
            if (rab < rmin[j]) rmin[j] = rab;
        }
    }
    
    return mdists;
}

void
CMolecule::printGeometry(COutputStream& oStream) const
{
    oStream << fmt::header;
    
    oStream << "Molecular Geometry (Angstroms)" << fmt::end;
    
    oStream << std::string(32, '=') << fmt::end << fmt::blank;
    
    oStream << "  Atom ";
    
    oStream << fstr::format(std::string("Coordinate X"), 20, fmt::right);
    
    oStream << "  ";
    
    oStream << fstr::format(std::string("Coordinate Y"), 20, fmt::right);
    
    oStream << "  ";
    
    oStream << fstr::format(std::string("Coordinate Z"), 20, fmt::right);
    
    oStream << "  " << fmt::end << fmt::blank;
    
    auto factor = units::getBohrValueInAngstroms();
    
    auto coordx = _atomCoordinates.data(0);
    
    auto coordy = _atomCoordinates.data(1);
    
    auto coordz = _atomCoordinates.data(2);
    
    for (size_t i = 0; i < _atomCoordinates.size(0); i++)
    {
        std::string label("  ");
        
        label.append(_atomLabels.at(i));
        
        oStream << fstr::format(label, 6, fmt::left);
        
        oStream << fstr::to_string(factor * coordx[i], 12, 22, fmt::right);
        
        oStream << fstr::to_string(factor * coordy[i], 12, 22, fmt::right);
        
        oStream << fstr::to_string(factor * coordz[i], 12, 22, fmt::right);
        
        oStream << fmt::end;
    }
    
    oStream << fmt::blank;
}

bool
CMolecule::checkProximity(const double        minDistance,
                               COutputStream& oStream) const
{
    auto natoms = _atomCoordinates.size(0);
    
    auto coordx = _atomCoordinates.data(0);
    
    auto coordy = _atomCoordinates.data(1);
    
    auto coordz = _atomCoordinates.data(2);
    
    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            auto rab = mathfunc::distance(coordx[i], coordy[i], coordz[i],
                                          coordx[j], coordy[j], coordz[j]);
            
            if (rab < minDistance)
            {
                _errorAtomsToClose(i, j, minDistance, oStream);
                
                return false;
            }
        }
    }
    
    return true;
}

void
CMolecule::broadcast(int32_t  rank,
                     MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_charge, comm);
        
        mpi::bcast(_multiplicity, comm);
        
        _atomCoordinates.broadcast(rank, comm);
        
        _atomCharges.broadcast(rank, comm);
        
        _atomMasses.broadcast(rank, comm);
        
        mpi::bcast(_atomLabels, rank, comm);
        
        _idsAtomic.broadcast(rank, comm);
        
        _idsElemental.broadcast(rank, comm);
    }
}

void
CMolecule::_errorAtomsToClose(const int32_t        iAtomA,
                              const int32_t        iAtomB,
                              const double         minDistance,
                                    COutputStream& oStream) const
{
    oStream << fmt::cerror << "Distance between atoms:" << fmt::end;
    
    oStream << _atomLabels[iAtomA] << "(" << std::to_string(iAtomA);
    
    oStream << ") and " << _atomLabels[iAtomB] << "(" << std::to_string(iAtomB);
    
    oStream << ") is smaller the minimal separation distance of ";
    
    oStream << std::to_string(minDistance) << " Bohr!" << fmt::end;
    
    oStream << "Please correct Your input!" << fmt::end << fmt::blank;
}

std::ostream&
operator<<(      std::ostream& output,
           const CMolecule&    source)
{
    output << std::endl;

    output << "[CMolecule (Object):" << &source << "]" << std::endl;

    output << "_charge: " << source._charge <<  std::endl;

    output << "_multiplicity: " << source._multiplicity <<  std::endl;
    
    output << "_atomCoordinates: " << source._atomCoordinates <<  std::endl;
    
    output << "_atomCharges: " << source._atomCharges << std::endl;
    
    output << "_atomMasses: " << source._atomMasses << std::endl;
    
    output << "_atomLabels: " << std::endl;

    for (size_t i = 0; i < source._atomLabels.size(); i++)
    {
        output << "atomsLabels_[" << i << "]" << std::endl;

        output << source._atomLabels[i] << std::endl;
    }

    output << "_idsAtomic: " << source._idsAtomic << std::endl;
    
    output << "_idsElemental: " << source._idsElemental << std::endl;
    
    return output;
}
