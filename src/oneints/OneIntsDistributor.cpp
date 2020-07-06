//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "OneIntsDistributor.hpp"

#include "AngularMomentum.hpp"

COneIntsDistribution::COneIntsDistribution()

    : _distPattern(dist1e::batch)

    , _nRows(0)

    , _nColumns(0)

    , _intsData(nullptr)
{
}

COneIntsDistribution::COneIntsDistribution(double* intsData, const int32_t nRows, const int32_t nColumns, const dist1e distPattern)

    : _distPattern(distPattern)

    , _nRows(nRows)

    , _nColumns(nColumns)

    , _intsData(intsData)
{
}

COneIntsDistribution::COneIntsDistribution(const COneIntsDistribution& source)

    : _distPattern(source._distPattern)

    , _nRows(source._nRows)

    , _nColumns(source._nColumns)
{
    _intsData = source._intsData;
}

COneIntsDistribution::~COneIntsDistribution()
{
}

COneIntsDistribution&
COneIntsDistribution::operator=(const COneIntsDistribution& source)
{
    if (this == &source) return *this;

    _distPattern = source._distPattern;

    _nRows = source._nRows;

    _nColumns = source._nColumns;

    _intsData = source._intsData;

    return *this;
}

bool
COneIntsDistribution::operator==(const COneIntsDistribution& other) const
{
    if (_distPattern != other._distPattern) return false;

    if (_nRows != other._nRows) return false;

    if (_nColumns != other._nColumns) return false;

    if (_intsData != other._intsData) return false;

    return true;
}

bool
COneIntsDistribution::operator!=(const COneIntsDistribution& other) const
{
    return !(*this == other);
}

void
COneIntsDistribution::distribute(const CMemBlock2D<double>& spherInts,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const bool                 isBraEqualKet,
                                 const int32_t              iContrGto)
{
    // distribute one electron integrals into data batch

    if (_distPattern == dist1e::batch)
    {
        _distSpherIntsIntoBatch(spherInts, braGtoBlock, ketGtoBlock, iContrGto);

        return;
    }

    // distribute one electron integrals into symmetric square matrix

    if (_distPattern == dist1e::symsq)
    {
        _distSpherIntsIntoSymMatrix(spherInts, braGtoBlock, ketGtoBlock, isBraEqualKet, iContrGto);

        return;
    }

    // distribute one electron integrals into ant-symmetric square matrix

    if (_distPattern == dist1e::antisq)
    {
        _distSpherIntsIntoAntiSymMatrix(spherInts, braGtoBlock, ketGtoBlock, isBraEqualKet, iContrGto);

        return;
    }

    // distribute one electron integrals into data batch

    if (_distPattern == dist1e::rect)
    {
        _distSpherIntsIntoRectMatrix(spherInts, braGtoBlock, ketGtoBlock, iContrGto);

        return;
    }
}

void
COneIntsDistribution::_distSpherIntsIntoBatch(const CMemBlock2D<double>& spherInts,
                                              const CGtoBlock&           braGtoBlock,
                                              const CGtoBlock&           ketGtoBlock,
                                              const int32_t              iContrGto)
{
    // set up number of angular components on bra and ket sides

    auto bcomp = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());

    auto kcomp = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());

    // set up number contracted GTOs on bra and ket sides

    auto nbgtos = braGtoBlock.getNumberOfContrGtos();

    auto nkgtos = ketGtoBlock.getNumberOfContrGtos();

    for (int32_t i = 0; i < bcomp; i++)
    {
        // offset in rows indexing space

        auto ioff = (i * nbgtos + iContrGto) * _nColumns;

        for (int32_t j = 0; j < kcomp; j++)
        {
            // set up pointer to one electron integrals

            auto fvals = spherInts.data(i * kcomp + j);

            // offset in full indexing space

            auto ijoff = ioff + j * nkgtos;

            // distribute integrals

            for (int32_t k = 0; k < nkgtos; k++)
            {
                _intsData[ijoff + k] = fvals[k];
            }
        }
    }
}

void
COneIntsDistribution::_distSpherIntsIntoSymMatrix(const CMemBlock2D<double>& spherInts,
                                                  const CGtoBlock&           braGtoBlock,
                                                  const CGtoBlock&           ketGtoBlock,
                                                  const bool                 isBraEqualKet,
                                                  const int32_t              iContrGto)
{
    // set up angular momentum components on bra and ket sides

    auto bcomp = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());

    auto kcomp = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());

    // set up number of contracted GTOs on ket side

    auto kdim = ketGtoBlock.getNumberOfContrGtos();

    for (int32_t i = 0; i < bcomp; i++)
    {
        auto bidx = (braGtoBlock.getIdentifiers(i))[iContrGto];

        // loop over ket components

        for (int32_t j = 0; j < kcomp; j++)
        {
            // set up pointer to integrals

            auto fvals = spherInts.data(i * kcomp + j);

            auto kidx = ketGtoBlock.getIdentifiers(j);

            // loop over integrals

            if (isBraEqualKet)
            {
                for (int32_t k = 0; k < kdim; k++)
                {
                    _intsData[bidx * _nColumns + kidx[k]] = fvals[k];
                }
            }
            else
            {
                for (int32_t k = 0; k < kdim; k++)
                {
                    _intsData[bidx * _nColumns + kidx[k]] = fvals[k];

                    _intsData[kidx[k] * _nColumns + bidx] = fvals[k];
                }
            }
        }
    }
}

void
COneIntsDistribution::_distSpherIntsIntoAntiSymMatrix(const CMemBlock2D<double>& spherInts,
                                                      const CGtoBlock&           braGtoBlock,
                                                      const CGtoBlock&           ketGtoBlock,
                                                      const bool                 isBraEqualKet,
                                                      const int32_t              iContrGto)
{
    // set up angular momentum components on bra and ket sides

    auto bcomp = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());

    auto kcomp = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());

    // set up number of contracted GTOs on ket side

    auto kdim = ketGtoBlock.getNumberOfContrGtos();

    for (int32_t i = 0; i < bcomp; i++)
    {
        auto bidx = (braGtoBlock.getIdentifiers(i))[iContrGto];

        // loop over ket components

        for (int32_t j = 0; j < kcomp; j++)
        {
            // set up pointer to integrals

            auto fvals = spherInts.data(i * kcomp + j);

            auto kidx = ketGtoBlock.getIdentifiers(j);

            // loop over integrals

            if (isBraEqualKet)
            {
                for (int32_t k = 0; k < kdim; k++)
                {
                    _intsData[bidx * _nColumns + kidx[k]] = fvals[k];
                }
            }
            else
            {
                // NOTE: assummes upper triangle computation in one-electron
                //       integrals driver

                for (int32_t k = 0; k < kdim; k++)
                {
                    _intsData[bidx * _nColumns + kidx[k]] = fvals[k];

                    _intsData[kidx[k] * _nColumns + bidx] = -fvals[k];
                }
            }
        }
    }
}

void
COneIntsDistribution::_distSpherIntsIntoRectMatrix(const CMemBlock2D<double>& spherInts,
                                                   const CGtoBlock&           braGtoBlock,
                                                   const CGtoBlock&           ketGtoBlock,
                                                   const int32_t              iContrGto)
{
    // set up angular momentum components on bra and ket sides

    auto bcomp = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum());

    auto kcomp = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum());

    // set up number of contracted GTOs on ket side

    auto kdim = ketGtoBlock.getNumberOfContrGtos();

    for (int32_t i = 0; i < bcomp; i++)
    {
        auto bidx = (braGtoBlock.getIdentifiers(i))[iContrGto];

        // loop over ket components

        for (int32_t j = 0; j < kcomp; j++)
        {
            // set up pointer to integrals

            auto fvals = spherInts.data(i * kcomp + j);

            auto kidx = ketGtoBlock.getIdentifiers(j);

            // loop over integrals

            for (int32_t k = 0; k < kdim; k++)
            {
                _intsData[bidx * _nColumns + kidx[k]] = fvals[k];
            }
        }
    }
}

std::ostream&
operator<<(std::ostream& output, const COneIntsDistribution& source)
{
    output << std::endl;

    output << "[COneIntsDistribution (Object):" << &source << "]" << std::endl;

    output << "_distPattern: " << to_string(source._distPattern) << std::endl;

    output << "_nRows: " << source._nRows << std::endl;

    output << "_nColumns: " << source._nColumns << std::endl;

    output << "_intsData: " << source._intsData << std::endl;

    return output;
}
