//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoContainer.hpp"

#include "AngularMomentum.hpp"

CGtoContainer::CGtoContainer()

    : _maxAngularMomentum(-1)
{

}

CGtoContainer::CGtoContainer(const std::vector<CGtoBlock>& gtoBlocks)

    : _gtoBlocks(gtoBlocks)

    , _maxAngularMomentum(-1)
{
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto angmom = _gtoBlocks[i].getAngularMomentum();

        if (_maxAngularMomentum < angmom) _maxAngularMomentum = angmom;
    }
}

CGtoContainer::CGtoContainer(const CMolecule&       molecule,
                             const CMolecularBasis& basis)

    : _maxAngularMomentum(-1)
{
    _maxAngularMomentum = basis.getMaxAngularMomentum();

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        CGtoBlock gtoblock(molecule, basis, i);

        if (!gtoblock.empty()) _gtoBlocks.push_back(gtoblock);
    }
}

CGtoContainer::CGtoContainer(const CGtoContainer& source)

    : _gtoBlocks(source._gtoBlocks)

    , _maxAngularMomentum(source._maxAngularMomentum)
{

}

CGtoContainer::CGtoContainer(CGtoContainer&& source) noexcept

    : _gtoBlocks(std::move(source._gtoBlocks))

    , _maxAngularMomentum(std::move(source._maxAngularMomentum))
{

}

CGtoContainer::~CGtoContainer()
{

}

CGtoContainer&
CGtoContainer::operator=(const CGtoContainer& source)
{
    if (this == &source) return *this;

    _gtoBlocks = source._gtoBlocks;

    _maxAngularMomentum = source._maxAngularMomentum;

    return *this;
}

CGtoContainer&
CGtoContainer::operator=(CGtoContainer&& source) noexcept
{
    if (this == &source) return *this;

    _gtoBlocks = std::move(source._gtoBlocks);

    _maxAngularMomentum = std::move(source._maxAngularMomentum);

    return *this;
}

bool
CGtoContainer::operator==(const CGtoContainer& other) const
{
     if (_maxAngularMomentum != other._maxAngularMomentum) return false;
    
    if (_gtoBlocks.size() != other._gtoBlocks.size()) return false;

    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        if (_gtoBlocks[i] != other._gtoBlocks[i]) return false;
    }

    return true;
}

bool
CGtoContainer::operator!=(const CGtoContainer& other) const
{
    return !(*this == other);
}

void
CGtoContainer::compress(const CGtoContainer&        source,
                              CMemBlock2D<int32_t>& reducedDimensions,
                        const CVecMemBlock<double>& screeningFactors,
                        const double                screeningThreshold)
{
    // set up pointers to reduced dimensions
    
    auto pidx = reducedDimensions.data(0);
    
    auto cidx = reducedDimensions.data(1);
    
    // loop over GTOs blocks
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto cdim  = _gtoBlocks[i].compress(source._gtoBlocks[i],
                                            screeningFactors[i],
                                            screeningThreshold);
        
        pidx[i] = std::get<0>(cdim);
        
        cidx[i] = std::get<1>(cdim);
    }
}

int32_t
CGtoContainer::getMaxAngularMomentum() const
{
    int32_t maxmom = 0;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto curmom = _gtoBlocks[i].getAngularMomentum();
        
        if (maxmom < curmom) maxmom= curmom;
    }
    
    return maxmom;
}

int32_t
CGtoContainer::getAngularMomentum(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getAngularMomentum();
}

int32_t
CGtoContainer::getNumberOfGtoBlocks() const
{
    return static_cast<int32_t>(_gtoBlocks.size());
}

int32_t
CGtoContainer::getMaxNumberOfPrimGtos() const
{
    int32_t nprimgtos = 0;

    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto cprimgtos = _gtoBlocks[i].getNumberOfPrimGtos();
        
        if (nprimgtos < cprimgtos) nprimgtos = cprimgtos;
    }
    
    return nprimgtos;
}

int32_t
CGtoContainer::getMaxNumberOfContrGtos() const
{
    int32_t ncontrgtos = 0;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto ccontrgtos = _gtoBlocks[i].getNumberOfContrGtos();
        
        if (ncontrgtos < ccontrgtos) ncontrgtos = ccontrgtos;
    }
    
    return ncontrgtos;
}

int32_t
CGtoContainer::getNumberOfPrimGtos(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getNumberOfPrimGtos();
}

int32_t
CGtoContainer::getNumberOfContrGtos(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getNumberOfContrGtos();
}

const int32_t*
CGtoContainer::getStartPositions(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getStartPositions();
}

const int32_t*
CGtoContainer::getEndPositions(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getEndPositions();
}

const int32_t*
CGtoContainer::getIdentifiers(const int32_t iBlock,
                              const int32_t iComponent) const
{
    return _gtoBlocks[iBlock].getIdentifiers(iComponent);
}

const double*
CGtoContainer::getExponents(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getExponents();
}

const double*
CGtoContainer::getNormFactors(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getNormFactors();
}

const double*
CGtoContainer::getCoordinatesX(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getCoordinatesX();
}

const double*
CGtoContainer::getCoordinatesY(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getCoordinatesY();
}

const double*
CGtoContainer::getCoordinatesZ(const int32_t iBlock) const
{
    return _gtoBlocks[iBlock].getCoordinatesZ();
}

CVecMemBlock<double>
CGtoContainer::getPrimBuffer() const
{
    CVecMemBlock<double> mbvec;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        mbvec.push_back(CMemBlock<double>(_gtoBlocks[i].getNumberOfPrimGtos()));
    }
    
    return mbvec;
}

CVecMemBlock2D<double>
CGtoContainer::getPrimAngBuffer(const int32_t nComponents) const
{
    CVecMemBlock2D<double> mbvec;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto cang = _getPrimAngComponents(_gtoBlocks[i].getAngularMomentum());
        
        mbvec.push_back(CMemBlock2D<double>(_gtoBlocks[i].getNumberOfPrimGtos(),
                                            nComponents * cang));
    }
    
    return mbvec;
}

CVecMemBlock2D<double>
CGtoContainer::getCartesianBuffer(const int32_t nComponents) const
{
    CVecMemBlock2D<double> mbvec;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto mang = _gtoBlocks[i].getAngularMomentum();
        
        auto cang = angmom::to_CartesianComponents(mang);
        
        mbvec.push_back(CMemBlock2D<double>(_gtoBlocks[i].getNumberOfContrGtos(),
                                            nComponents * cang));
    }
    
    return mbvec;
}

CVecMemBlock2D<double>
CGtoContainer::getSphericalBuffer(const int32_t nComponents) const
{
    CVecMemBlock2D<double> mbvec;
    
    for (size_t i = 0; i < _gtoBlocks.size(); i++)
    {
        auto mang = _gtoBlocks[i].getAngularMomentum();
        
        auto cang = angmom::to_SphericalComponents(mang);
        
        mbvec.push_back(CMemBlock2D<double>(_gtoBlocks[i].getNumberOfContrGtos(),
                                            nComponents * cang));
    }
    
    return mbvec;
}

int32_t
CGtoContainer::_getPrimAngComponents(const int32_t angularMomentum) const
{
    int32_t ncomp = 0;
    
    for (int32_t i = 0; i <= angularMomentum; i++)
    {
        ncomp += angmom::to_CartesianComponents(i);
    }
    
    return ncomp;
}

std::ostream&
operator<<(      std::ostream&  output,
           const CGtoContainer& source)
{
    output << std::endl;
    
    output << "[CCGtoContainer (Object):" << &source << "]" << std::endl;
    
    output << "_maxAngularMomentum: " << source._maxAngularMomentum << std::endl;
    
    output << "_gtoBlocks: " << std::endl;
    
    for (size_t i = 0; i < source._gtoBlocks.size(); i++)
    {
        output << "_gtoBlocks[" << i << "]: "<< std::endl;

        output << source._gtoBlocks[i] << std::endl;
    }

    return output;
}
