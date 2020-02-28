//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "MolecularOrbitals.hpp"

#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "ErrorHandler.hpp"

CMolecularOrbitals::CMolecularOrbitals()

    : _orbitalsType(molorb::rest)

    , _orbitals(std::vector<CDenseMatrix>())

    , _energies(std::vector<CMemBlock<double>>())
{
}

CMolecularOrbitals::CMolecularOrbitals(const std::vector<CDenseMatrix>&      orbitals,
                                       const std::vector<CMemBlock<double>>& energies,
                                       const molorb                          orbitalsType)

    : _orbitalsType(orbitalsType)

    , _orbitals(orbitals)

    , _energies(energies)
{
}

CMolecularOrbitals::CMolecularOrbitals(const CMolecularOrbitals& source)

    : _orbitalsType(source._orbitalsType)

    , _orbitals(source._orbitals)

    , _energies(source._energies)
{
}

CMolecularOrbitals::CMolecularOrbitals(CMolecularOrbitals&& source) noexcept

    : _orbitalsType(std::move(source._orbitalsType))

    , _orbitals(std::move(source._orbitals))

    , _energies(std::move(source._energies))
{
}

CMolecularOrbitals::~CMolecularOrbitals()
{
}

CMolecularOrbitals&
CMolecularOrbitals::operator=(const CMolecularOrbitals& source)
{
    if (this == &source) return *this;

    _orbitalsType = source._orbitalsType;

    _orbitals = source._orbitals;

    _energies = source._energies;

    return *this;
}

CMolecularOrbitals&
CMolecularOrbitals::operator=(CMolecularOrbitals&& source) noexcept
{
    if (this == &source) return *this;

    _orbitalsType = std::move(source._orbitalsType);

    _orbitals = std::move(source._orbitals);

    _energies = std::move(source._energies);

    return *this;
}

bool
CMolecularOrbitals::operator==(const CMolecularOrbitals& other) const
{
    if (_orbitalsType != other._orbitalsType) return false;

    if (_orbitals.size() != other._orbitals.size()) return false;

    for (size_t i = 0; i < _orbitals.size(); i++)
    {
        if (_orbitals[i] != other._orbitals[i]) return false;
    }

    if (_energies.size() != other._energies.size()) return false;

    for (size_t i = 0; i < _energies.size(); i++)
    {
        if (_energies[i] != other._energies[i]) return false;
    }

    return true;
}

bool
CMolecularOrbitals::operator!=(const CMolecularOrbitals& other) const
{
    return !(*this == other);
}

CMolecularOrbitals
CMolecularOrbitals::insert(const CMolecule& molecule, const CMolecularBasis& aoBasis, const CMolecularBasis& minBasis) const
{
    std::vector<CDenseMatrix> orbvec;

    // create orbital coeficients vector

    auto naos = aoBasis.getDimensionsOfBasis(molecule);

    for (size_t i = 0; i < _orbitals.size(); i++)
    {
        CDenseMatrix cmat(naos, _orbitals[i].getNumberOfColumns());

        cmat.zero();

        orbvec.push_back(cmat);
    }

    // set up pointer to chemical elements data

    auto idselm = molecule.getIdsElemental();

    // insert molecular orbitals coeficients

    int32_t midx = 0;

    int32_t cidx = 0;

    for (int32_t i = 0; i <= minBasis.getMaxAngularMomentum(); i++)
    {
        for (int32_t j = 0; j < angmom::to_SphericalComponents(i); j++)
        {
            for (int32_t k = 0; k < molecule.getNumberOfAtoms(); k++)
            {
                auto mbfs = minBasis.getNumberOfBasisFunctions(idselm[k], i);

                auto cbfs = aoBasis.getNumberOfBasisFunctions(idselm[k], i);

                // copy orbital coeficients

                if (mbfs > 0)
                {
                    for (size_t l = 0; l < _orbitals.size(); l++)
                    {
                        auto morbs = _orbitals[l].values();

                        auto corbs = orbvec[l].values();

                        auto mdim = _orbitals[l].getNumberOfColumns();

                        for (int32_t m = 0; m < mbfs; m++)
                        {
                            auto moff = (midx + m) * mdim;

                            auto coff = (cidx + m) * mdim;

                            for (int32_t n = 0; n < mdim; n++)
                            {
                                corbs[coff + n] = morbs[moff + n];
                            }
                        }
                    }

                    midx += mbfs;
                }

                cidx += cbfs;
            }
        }
    }

    return CMolecularOrbitals(orbvec, _energies, _orbitalsType);
}

molorb
CMolecularOrbitals::getOrbitalsType() const
{
    return _orbitalsType;
}

int32_t
CMolecularOrbitals::getNumberOfRows() const
{
    return _orbitals[0].getNumberOfRows();
}

int32_t
CMolecularOrbitals::getNumberOfColumns() const
{
    return _orbitals[0].getNumberOfColumns();
}

const double*
CMolecularOrbitals::alphaOrbitals() const
{
    return _orbitals[0].values();
}

const double*
CMolecularOrbitals::betaOrbitals() const
{
    if (_orbitalsType == molorb::rest)
    {
        return alphaOrbitals();
    }

    if (_orbitalsType == molorb::unrest)
    {
        return _orbitals[1].values();
    }

    return nullptr;
}

CDenseMatrix
CMolecularOrbitals::alphaOrbitals(const int32_t iMolecularOrbital, const int32_t nMolecularOrbitals) const
{
    return _orbitals[0].slice(0, iMolecularOrbital, getNumberOfRows(), nMolecularOrbitals);
}

CDenseMatrix
CMolecularOrbitals::betaOrbitals(const int32_t iMolecularOrbital, const int32_t nMolecularOrbitals) const
{
    if (_orbitalsType == molorb::rest)
    {
        return alphaOrbitals(iMolecularOrbital, nMolecularOrbitals);
    }

    if (_orbitalsType == molorb::unrest)
    {
        return _orbitals[1].slice(0, iMolecularOrbital, getNumberOfRows(), nMolecularOrbitals);
    }

    return CDenseMatrix();
}

CDenseMatrix
CMolecularOrbitals::alphaOrbitals(const std::vector<int32_t>& iMolecularOrbitals) const
{
    return _orbitals[0].selectByColumn(iMolecularOrbitals);
}

CDenseMatrix
CMolecularOrbitals::betaOrbitals(const std::vector<int32_t>& iMolecularOrbitals) const
{
    if (_orbitalsType == molorb::rest)
    {
        return alphaOrbitals(iMolecularOrbitals);
    }

    if (_orbitalsType == molorb::unrest)
    {
        return _orbitals[1].selectByColumn(iMolecularOrbitals);
    }

    return CDenseMatrix();
}

const double*
CMolecularOrbitals::alphaEnergies() const
{
    return _energies[0].data();
}

const double*
CMolecularOrbitals::betaEnergies() const
{
    if (_orbitalsType == molorb::rest)
    {
        return alphaEnergies();
    }

    if (_orbitalsType == molorb::unrest)
    {
        return _energies[1].data();
    }

    return nullptr;
}

std::string
CMolecularOrbitals::getString() const
{
    std::string orb_str;

    orb_str += "Orbitals Type: " + to_string(_orbitalsType) + "\n";

    for (size_t i = 0; i < _orbitals.size(); i++)
    {
        orb_str += _orbitals[i].getString();
    }

    return orb_str;
}

CAODensityMatrix
CMolecularOrbitals::getAODensity(const int32_t nElectrons) const
{
    if ((nElectrons % 2) == 0)
    {
        auto ndim = nElectrons / 2;

        auto nrow = _orbitals[0].getNumberOfRows();

        auto ncol = _orbitals[0].getNumberOfColumns();

        if (ndim <= ncol)
        {
            auto cmo = _orbitals[0].slice(0, 0, nrow, ndim);

            auto den = denblas::multABt(cmo, cmo);

            return CAODensityMatrix({den}, denmat::rest);
        }
    }

    return CAODensityMatrix();
}

CAODensityMatrix
CMolecularOrbitals::getAODensity(const int32_t nAlphaElectrons, const int32_t nBetaElectrons) const
{
    auto ndima = nAlphaElectrons;

    auto ndimb = nBetaElectrons;

    auto nrowa = _orbitals[0].getNumberOfRows();

    auto ncola = _orbitals[0].getNumberOfColumns();

    auto nrowb = _orbitals[1].getNumberOfRows();

    auto ncolb = _orbitals[1].getNumberOfColumns();

    if ((ndima <= ncola) && (ndimb <= ncolb))
    {
        auto cmoa = _orbitals[0].slice(0, 0, nrowa, ndima);

        auto dena = denblas::multABt(cmoa, cmoa);

        auto cmob = _orbitals[1].slice(0, 0, nrowb, ndimb);

        auto denb = denblas::multABt(cmob, cmob);

        return CAODensityMatrix({dena, denb}, denmat::unrest);
    }

    return CAODensityMatrix();
}

CAODensityMatrix
CMolecularOrbitals::getRestrictedPairDensity(const int32_t iMolecularOrbital, const int32_t jMolecularOrbital) const
{
    auto nrow = _orbitals[0].getNumberOfRows();

    auto ncol = _orbitals[0].getNumberOfColumns();

    if ((iMolecularOrbital < ncol) && (jMolecularOrbital < ncol))
    {
        auto cmi = _orbitals[0].slice(0, iMolecularOrbital, nrow, 1);

        auto cmj = _orbitals[0].slice(0, jMolecularOrbital, nrow, 1);

        auto den = denblas::multABt(cmi, cmj);

        return CAODensityMatrix({den}, denmat::rmoij);
    }

    return CAODensityMatrix();
}

CAODensityMatrix
CMolecularOrbitals::getRestrictedPairDensity(const std::vector<int32_t>& iMolecularOrbitals, const std::vector<int32_t>& jMolecularOrbitals) const
{
    auto nrow = _orbitals[0].getNumberOfRows();

    auto ncol = _orbitals[0].getNumberOfColumns();

    if (iMolecularOrbitals.size() == jMolecularOrbitals.size())
    {
        std::vector<CDenseMatrix> denvec;

        for (size_t i = 0; i < iMolecularOrbitals.size(); i++)
        {
            auto icol = iMolecularOrbitals[i];

            auto jcol = jMolecularOrbitals[i];

            if ((icol < ncol) && (jcol < ncol))
            {
                auto cmi = _orbitals[0].slice(0, icol, nrow, 1);

                auto cmj = _orbitals[0].slice(0, jcol, nrow, 1);

                denvec.push_back(denblas::multABt(cmi, cmj));
            }
            else
            {
                return CAODensityMatrix();
            }
        }

        return CAODensityMatrix(denvec, denmat::rmoij);
    }

    return CAODensityMatrix();
}

CDenseMatrix
CMolecularOrbitals::transform(const CDenseMatrix& aoMatrix, const szblock spinPair) const
{
    // alpha - alpha

    if (spinPair == szblock::aa)
    {
        return denblas::multAtB(_orbitals[0], denblas::multAB(aoMatrix, _orbitals[0]));
    }

    // alpha - beta

    if (spinPair == szblock::ab)
    {
        return denblas::multAtB(_orbitals[0], denblas::multAB(aoMatrix, _orbitals[1]));
    }

    // beta - alpha

    if (spinPair == szblock::ba)
    {
        return denblas::multAtB(_orbitals[1], denblas::multAB(aoMatrix, _orbitals[0]));
    }

    // beta - beta

    if (spinPair == szblock::bb)
    {
        return denblas::multAtB(_orbitals[1], denblas::multAB(aoMatrix, _orbitals[1]));
    }

    return CDenseMatrix();
}

void
CMolecularOrbitals::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        // broadcast molecular orbital type

        int32_t motyp = 0;

        if (rank == mpi::master()) motyp = to_int(_orbitalsType);

        mpi::bcast(motyp, comm);

        if (rank != mpi::master()) _orbitalsType = to_molorb(motyp);

        // broadcast molecular orbitals

        int32_t norbs = static_cast<int32_t>(_orbitals.size());

        mpi::bcast(norbs, comm);

        for (int32_t i = 0; i < norbs; i++)
        {
            CDenseMatrix morb;

            if (rank == mpi::master()) morb = _orbitals[i];

            morb.broadcast(rank, comm);

            if (rank != mpi::master()) _orbitals.push_back(morb);
        }

        // broadcast orbital energies

        int32_t nvecs = static_cast<int32_t>(_energies.size());

        mpi::bcast(nvecs, comm);

        for (int32_t i = 0; i < nvecs; i++)
        {
            CMemBlock<double> vec;

            if (rank == mpi::master()) vec = _energies[i];

            vec.broadcast(rank, comm);

            if (rank != mpi::master()) _energies.push_back(vec);
        }
    }
}

std::ostream&
operator<<(std::ostream& output, const CMolecularOrbitals& source)
{
    output << std::endl;

    output << "[CMolecularOrbitals (Object):" << &source << "]" << std::endl;

    output << "_orbitalsType: " << to_string(source._orbitalsType) << std::endl;

    output << "_orbitals: " << std::endl;

    for (size_t i = 0; i < source._orbitals.size(); i++)
    {
        output << "_orbitals[" << i << "]: " << std::endl;

        output << source._orbitals[i] << std::endl;
    }

    output << "_energies: " << std::endl;

    for (size_t i = 0; i < source._energies.size(); i++)
    {
        output << "_energies[" << i << "]: " << std::endl;

        output << source._energies[i] << std::endl;
    }

    return output;
}
