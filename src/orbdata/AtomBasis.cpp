//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AtomBasis.hpp"

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"
#include "MpiFunc.hpp"

CAtomBasis::CAtomBasis()

    : _idElemental(-1)

    , _maxAngularMomentum(-1)
{

}

CAtomBasis::CAtomBasis(const CAtomBasis& source)

    : _basisFunctions(source._basisFunctions)

    , _idElemental(source._idElemental)

    , _maxAngularMomentum(source._maxAngularMomentum)
{

}

CAtomBasis::CAtomBasis(CAtomBasis&& source) noexcept

    : _basisFunctions(std::move(source._basisFunctions))

    , _idElemental(std::move(source._idElemental))

    , _maxAngularMomentum(std::move(source._maxAngularMomentum))
{

}

CAtomBasis::~CAtomBasis()
{

}

CAtomBasis&
CAtomBasis::operator=(const CAtomBasis& source)
{
    if (this == &source) return *this;

    _basisFunctions = source._basisFunctions;

    _idElemental = source._idElemental;

    _maxAngularMomentum = source._maxAngularMomentum;

    return *this;
}

CAtomBasis&
CAtomBasis::operator=(CAtomBasis&& source) noexcept
{
    if (this == &source) return *this;

    _basisFunctions = std::move(source._basisFunctions);

    _idElemental = std::move(source._idElemental);

    _maxAngularMomentum = std::move(source._maxAngularMomentum);

    return *this;
}

bool
CAtomBasis::operator==(const CAtomBasis& other) const
{
    if (_basisFunctions.size() != other._basisFunctions.size()) return false;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i] != other._basisFunctions[i]) return false;
    }

    if (_idElemental != other._idElemental) return false;

    if (_maxAngularMomentum != other._maxAngularMomentum) return false;

    return true;
}

bool
CAtomBasis::operator!=(const CAtomBasis& other) const
{
    return !(*this == other);
}

void
CAtomBasis::setIdElemental(const int32_t idElemental)
{
    _idElemental = idElemental;
}

void
CAtomBasis::setMaxAngularMomentum(const int32_t maxAngularMomentum)
{
    _maxAngularMomentum = maxAngularMomentum;
}

void
CAtomBasis::addBasisFunction(const CBasisFunction& basisFunction)
{
    _basisFunctions.push_back(basisFunction);

    auto bAngularMomentum = basisFunction.getAngularMomentum();

    if (bAngularMomentum > _maxAngularMomentum)
    {
         _maxAngularMomentum = bAngularMomentum;
     }
}

int32_t
CAtomBasis::getIdElemental() const
{
    return _idElemental;
}

int32_t
CAtomBasis::getMaxAngularMomentum() const
{
    return _maxAngularMomentum;
}

int32_t
CAtomBasis::getNumberOfBasisFunctions(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;

    int32_t nbfuncs = 0;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            nbfuncs += _basisFunctions[i].getNumberOfContractedFunctions();
        }
    }

    return nbfuncs;
}

int32_t
CAtomBasis::getNumberOfReducedBasisFunctions(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;
    
    int32_t nbfuncs = 0;
    
    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            nbfuncs++;
        }
    }
    
    return nbfuncs;
}

int32_t
CAtomBasis::getNumberOfPrimitiveFunctions(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;

    int32_t npfuncs = 0;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
         if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
         {
             npfuncs += _basisFunctions[i].getNumberOfPrimitiveFunctions();
         }
    }

    return npfuncs;
}

int32_t
CAtomBasis::getNumberOfNormalizationFactors(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;
    
    int32_t nfacts = 0;
    
    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            nfacts += _basisFunctions[i].getNumberOfNormalizationFactors();
        }
    }
    
    return nfacts;
}

std::string
CAtomBasis::getContractionString() const
{
    std::string str("(");

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        str.append(std::to_string(getNumberOfBasisFunctions(i)));

        str.append(fstr::to_AngularMomentum(i));

        if (i != _maxAngularMomentum) str.append(",");
    }

    str.append(")");

    return str;
}

std::string
CAtomBasis::getPrimitivesString() const
{
    std::string str("(");

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        str.append(std::to_string(getNumberOfPrimitiveFunctions(i)));

        str.append(fstr::to_AngularMomentum(i));

        if (i != _maxAngularMomentum) str.append(",");
    }

    str.append(")");

    return str;
}

std::vector<CBasisFunction>
CAtomBasis::getBasisFunctions(const int32_t angularMomentum) const
{
    std::vector<CBasisFunction> basvector;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            basvector.push_back(_basisFunctions[i]);
        }
    }

    return basvector;
}

CAtomBasis
CAtomBasis::reduceToValenceBasis() const
{
    // set atomic shell max. angular momentum
    
    CChemicalElement chemele;
    
    chemele.setAtomType(_idElemental);
    
    auto mang = chemele.getMaxAngularMomentum();
    
    // generate valence basis
    
    CAtomBasis valbas;
    
    valbas.setIdElemental(_idElemental);
    
    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() <= mang)
        {
            valbas.addBasisFunction(_basisFunctions[i]);
        }
    }
    
    return valbas;
}

void
CAtomBasis::broadcast(int32_t  rank,
                      MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_idElemental, comm);

        mpi::bcast(_maxAngularMomentum, comm);

        int32_t nbasfuncs = static_cast<int32_t>(_basisFunctions.size());

        mpi::bcast(nbasfuncs, comm);

        for (int32_t i = 0; i < nbasfuncs; i++)
        {
            CBasisFunction bfunc;

            if (rank == mpi::master()) bfunc = _basisFunctions[i];

            bfunc.broadcast(rank, comm);

            if (rank != mpi::master()) addBasisFunction(bfunc);
            
            MPI_Barrier(comm); 
        }
    }
}

std::ostream&
operator<<(      std::ostream& output,
           const CAtomBasis&   source)
{
    output << std::endl;

    output << "[CAtomBasis (Object):" << &source << "]" << std::endl;

    output << "_idElemental: " << source._idElemental << std::endl;

    output << "_maxAngularMomentum: " << source._maxAngularMomentum;

    output << std::endl;

    output << "_basisFunctions: " << std::endl;

    for (size_t i = 0; i < source._basisFunctions.size(); i++)
    {
        output << "_basisFunctions[" << i << "]: "<< std::endl;

        output << source._basisFunctions[i] << std::endl;
    }

    return output;
}
