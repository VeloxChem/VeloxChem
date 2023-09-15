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

        for (int64_t ind = 1; ind < static_cast<int64_t>(_orb_indexes.size()); ind++)
        {
            ao_inds.push_back(comp * _orb_indexes[0] + _orb_indexes[ind]);
        }
    }

    return ao_inds;
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
