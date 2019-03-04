//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoIntsDistributor.hpp"

#include <cmath>

#include "StringFormat.hpp"
#include "AngularMomentum.hpp"
#include "DistFock.hpp"
#include "DistMaxDensity.hpp"

CTwoIntsDistribution::CTwoIntsDistribution()

    : _distPattern(dist2e::batch)

    , _nRows(0)

    , _nColumns(0)

    , _idGtoPair(-1)

    , _intsData(nullptr)

    , _aoDensity(nullptr)

    , _aoFock(nullptr)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(      double* intsData,
                                           const int32_t nRows,
                                           const int32_t nColumns,
                                           const dist2e  distPattern)

    : CTwoIntsDistribution(intsData, nRows, nColumns, -1, distPattern)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(      double* intsData,
                                           const int32_t nRows,
                                           const int32_t nColumns,
                                           const int32_t idGtoPair,
                                           const dist2e  distPattern)
    : _distPattern(distPattern)

    , _nRows(nRows)

    , _nColumns(nColumns)

    , _idGtoPair(idGtoPair)

    , _intsData(intsData)

    , _aoDensity(nullptr)

    , _aoFock(nullptr)
{
    // override identifier of GTOs pair if not used in distributing integrals
    
    if (_distPattern != dist2e::qvalues) _idGtoPair = -1;
}

CTwoIntsDistribution::CTwoIntsDistribution(      CAOFockMatrix*    aoFock,
                                           const CAODensityMatrix* aoDensity)

    : _distPattern(dist2e::fock)

    , _nRows(0)

    , _nColumns(0)

    , _idGtoPair(-1)

    , _intsData(nullptr)

    , _aoDensity(aoDensity)

    , _aoFock(aoFock)
{
    
}

CTwoIntsDistribution::CTwoIntsDistribution(const CTwoIntsDistribution& source)

    : _distPattern(source._distPattern)

    , _nRows(source._nRows)

    , _nColumns(source._nColumns)

    , _idGtoPair(source._idGtoPair)

    , _intsData(source._intsData)

    , _aoDensity(source._aoDensity)

    , _aoFock(source._aoFock)
{
}

CTwoIntsDistribution::~CTwoIntsDistribution()
{
    
}

CTwoIntsDistribution&
CTwoIntsDistribution::operator=(const CTwoIntsDistribution& source)
{
    if (this == &source) return *this;
    
    _distPattern = source._distPattern;
    
    _nRows = source._nRows;
    
    _nColumns = source._nColumns;
    
    _idGtoPair = source._idGtoPair;
    
    _intsData = source._intsData;
    
    _aoDensity = source._aoDensity;
    
    _aoFock = source._aoFock;
    
    _fockContainer = source._fockContainer;
    
    return *this;
}

bool
CTwoIntsDistribution::operator==(const CTwoIntsDistribution& other) const
{
    if (_distPattern != other._distPattern) return false;
    
    if (_nRows != other._nRows) return false;
    
    if (_nColumns != other._nColumns) return false;
    
    if (_idGtoPair != other._idGtoPair) return false;
    
    if (_intsData != other._intsData) return false;
    
    if (_aoDensity != other._aoDensity) return false;
    
    if (_aoFock != other._aoFock) return false;
    
    if (_fockContainer != other._fockContainer) return false;
    
    return true;
}

bool
CTwoIntsDistribution::operator!=(const CTwoIntsDistribution& other) const
{
    return !(*this == other);
}

void
CTwoIntsDistribution::setFockContainer(const CGtoPairsBlock&      braGtoPairsBlock,
                                       const CGtoPairsBlock&      ketGtoPairsBlock)
{
    if (_distPattern == dist2e::fock)
    {
        _fockContainer = CFockContainer(_aoFock, braGtoPairsBlock,
                                        ketGtoPairsBlock);
    }
}

void
CTwoIntsDistribution::distribute(const CMemBlock2D<double>& spherInts,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const bool                 isBraEqualKet,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
{
    // distribute two electron integrals into data batch
    
    if (_distPattern == dist2e::batch)
    {
        _distSpherIntsIntoBatch(spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                isBraEqualKet, iContrPair);
        
        return; 
    }
    
    // distribute two electron integrals into qvalues vector
    
    if (_distPattern == dist2e::qvalues)
    {
        _distSpherIntsIntoQValues(spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                  isBraEqualKet, iContrPair);
        
        return;
    }
    
    if (_distPattern == dist2e::fock)
    {
        _distSpherIntsIntoFock(spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                               isBraEqualKet, nKetContrPairs, iContrPair);
        
        return;
    }
}

void
CTwoIntsDistribution::accumulate()
{
    _fockContainer.accumulate(_aoFock);
}

void
CTwoIntsDistribution::getMaxDensityElements(      CMemBlock<double>& maxDensityElements,
                                            const CGtoPairsBlock&    braGtoPairsBlock,
                                            const CGtoPairsBlock&    ketGtoPairsBlock,
                                            const bool               isBraEqualKet,
                                            const int32_t            nKetContrPairs,
                                            const int32_t            iContrPair) const
{
    // initialize vector of max. density elements
    
    maxDensityElements.zero();
    
    // set up number of AO Fock matrices
    
    auto nfock = _aoFock->getNumberOfFockMatrices();
    
    for (int32_t i = 0; i < nfock; i++)
    {
        // set up fock matrix type and origin
        
        auto fcktyp = _aoFock->getFockType(i);
        
        auto idden = _aoFock->getDensityIdentifier(i);
        
        // closed shell restricted Hatree-Fock: 2J + K
        
        if (fcktyp == fockmat::restjk)
        {
            distmaxden::getMaxRestDenJK(maxDensityElements,
                                        _aoDensity->totalDensity(idden),
                                        _aoDensity->getNumberOfColumns(idden),
                                        braGtoPairsBlock, ketGtoPairsBlock,
                                        nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted Hatree-Fock: J
        
        if (fcktyp == fockmat::restj)
        {
            distmaxden::getMaxRestDenJ(maxDensityElements,
                                       _aoDensity->totalDensity(idden),
                                       _aoDensity->getNumberOfColumns(idden),
                                       braGtoPairsBlock, ketGtoPairsBlock,
                                       nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted Hatree-Fock: K
        
        if (fcktyp == fockmat::restk)
        {
            distmaxden::getMaxRestDenK(maxDensityElements,
                                       _aoDensity->totalDensity(idden),
                                       _aoDensity->getNumberOfColumns(idden),
                                       braGtoPairsBlock, ketGtoPairsBlock,
                                       nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general Coulomb: J
        
        if (fcktyp == fockmat::rgenj)
        {
            distmaxden::getMaxRestGenDenJ(maxDensityElements,
                                          _aoDensity->getDensity(idden),
                                          _aoDensity->getNumberOfColumns(idden),
                                          braGtoPairsBlock, ketGtoPairsBlock,
                                          nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general exchange: K
        
        if (fcktyp == fockmat::rgenk)
        {
            distmaxden::getMaxRestGenDenK(maxDensityElements,
                                          _aoDensity->getDensity(idden),
                                          _aoDensity->getNumberOfColumns(idden),
                                          braGtoPairsBlock, ketGtoPairsBlock,
                                          nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general Fock: 2J - K
        
        if (fcktyp == fockmat::rgenjk)
        {
            distmaxden::getMaxRestGenDenJK(maxDensityElements,
                                           _aoDensity->getDensity(idden),
                                           _aoDensity->getNumberOfColumns(idden),
                                           braGtoPairsBlock, ketGtoPairsBlock,
                                           nKetContrPairs, iContrPair);
        }
        
    }
}

void
CTwoIntsDistribution::_distSpherIntsIntoBatch(const CMemBlock2D<double>& spherInts,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const bool                 isBraEqualKet,
                                              const int32_t              iContrPair)
{
    // determine number of angular components in shhell
    
    auto ncomp = angmom::to_SphericalComponents(braGtoPairsBlock.getBraAngularMomentum(),
                                                braGtoPairsBlock.getKetAngularMomentum())
    
               * angmom::to_SphericalComponents(ketGtoPairsBlock.getBraAngularMomentum(),
                                                ketGtoPairsBlock.getKetAngularMomentum());
    
    // determine starting poisition in integrals batch
    
    auto spos = _getStartIndexForBatch(ncomp, iContrPair, isBraEqualKet);
    
    // determine dimensions of GTOs pairs vector on ket side
    
    auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();

    if (isBraEqualKet) kdim = iContrPair + 1;
    
    // store integrals into batch of integrals
    
    for (int32_t i = 0; i < kdim; i++)
    {
        for (int32_t j = 0; j < ncomp; j++)
        {
            _intsData[spos + j] = (spherInts.data(j))[i];
        }
        
        spos += ncomp;
    }
}

void
CTwoIntsDistribution::_distSpherIntsIntoQValues(const CMemBlock2D<double>& spherInts,
                                                const CGtoPairsBlock&      braGtoPairsBlock,
                                                const CGtoPairsBlock&      ketGtoPairsBlock,
                                                const bool                 isBraEqualKet,
                                                const int32_t              iContrPair)
{
    // check if dimensions on bra and ket sides are correct
    
    if ((braGtoPairsBlock.getNumberOfScreenedContrPairs() != 1) ||
        (ketGtoPairsBlock.getNumberOfScreenedContrPairs() != 1))
    {
        // error handling code...
        
        return;
    }
    
    // determine number of angular components in shhell
    
    auto ncomp = angmom::to_SphericalComponents(braGtoPairsBlock.getBraAngularMomentum(),
                                                braGtoPairsBlock.getKetAngularMomentum())
    
               * angmom::to_SphericalComponents(ketGtoPairsBlock.getBraAngularMomentum(),
                                                ketGtoPairsBlock.getKetAngularMomentum());
    
    //  determine maximum integral value in shell
    
    auto mval = (spherInts.data(0))[0];
    
    for (int32_t j = 1; j < ncomp; j++)
    {
        auto fval = (spherInts.data(j))[0];
        
        if (fval > mval) mval = fval;
    }
    
    // compute Q value for contracted pair
    
    _intsData[_idGtoPair] = std::sqrt(mval);
}

void
CTwoIntsDistribution::_distSpherIntsIntoFock(const CMemBlock2D<double>& spherInts,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const bool                 isBraEqualKet,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair)
{
    // set up number of AO Fock matrices
    
    auto nfock = _aoFock->getNumberOfFockMatrices();
    
    for (int32_t i = 0; i < nfock; i++)
    {
        // set up fock matrix type and origin 
        
        auto fcktyp = _aoFock->getFockType(i);
        
        auto idden = _aoFock->getDensityIdentifier(i);
        
        // closed shell restricted Hatree-Fock: 2J + K
        
        if (fcktyp == fockmat::restjk)
        {
            distfock::distRestJK(_fockContainer, i,
                                 _aoDensity->totalDensity(idden),
                                 _aoDensity->getNumberOfColumns(idden),
                                 spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                 nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted Hatree-Fock: J
        
        if (fcktyp == fockmat::restj)
        {
            distfock::distRestJ(_fockContainer, i,
                                _aoDensity->totalDensity(idden),
                                _aoDensity->getNumberOfColumns(idden),
                                spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted Hatree-Fock: K
        
        if (fcktyp == fockmat::restk)
        {
            distfock::distRestK(_fockContainer, i,
                                _aoDensity->totalDensity(idden),
                                _aoDensity->getNumberOfColumns(idden),
                                spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general Coulomb matrix: J
        
        if (fcktyp == fockmat::rgenj)
        {
            distfock::distRestGenJ(_fockContainer, i,
                                   _aoDensity->getDensity(idden),
                                   _aoDensity->getNumberOfColumns(idden),
                                   spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                   nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general exchange matrix: K
        
        if (fcktyp == fockmat::rgenk)
        {
            distfock::distRestGenK(_fockContainer, i,
                                   _aoDensity->getDensity(idden),
                                   _aoDensity->getNumberOfColumns(idden),
                                   spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                   nKetContrPairs, iContrPair);
        }
        
        // closed shell restricted general Fock matrix: 2J - K
        
        if (fcktyp == fockmat::rgenjk)
        {
            distfock::distRestGenJK(_fockContainer, i,
                                    _aoDensity->getDensity(idden),
                                    _aoDensity->getNumberOfColumns(idden),
                                    spherInts, braGtoPairsBlock, ketGtoPairsBlock,
                                    nKetContrPairs, iContrPair);
        }
    }
}

int32_t
CTwoIntsDistribution::_getStartIndexForBatch(const int32_t nShellComponents,
                                             const int32_t iContrPair,
                                             const bool    isBraEqualKet) const
{
    int32_t idx = 0;
    
    if (isBraEqualKet)
    {
        for (int32_t i = 0; i < iContrPair; i++) idx += i + 1; 
    }
    else
    {
        idx = iContrPair * _nColumns;
    }
    
    return idx * nShellComponents;
}

std::ostream&
operator<<(      std::ostream&         output,
           const CTwoIntsDistribution& source)
{
    output << std::endl;
    
    output << "[CTwoIntsDistribution (Object):" << &source << "]" << std::endl;
    
    output << "_distPattern: " << to_string(source._distPattern) << std::endl;
    
    output << "_nRows: " << source._nRows << std::endl;
    
    output << "_nColumns: " << source._nColumns << std::endl;
    
    output << "_idGtoPair: " << source._idGtoPair << std::endl;
    
    output << "_intsData: " << source._intsData << std::endl;
    
    output << "_aoDensity: " << source._aoDensity << std::endl;
    
    output << "_aoFock: " << source._aoFock << std::endl;
    
    output << "_fockContainer: " << source._fockContainer << std::endl;
    
    return output;
}
