//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoBlock.hpp"

#include <utility>
#include <cmath>

#include "AngularMomentum.hpp"
#include "MatOrder.hpp"

CGtoBlock::CGtoBlock()

    : _angularMomentum(-1)
{

}

CGtoBlock::CGtoBlock(const CMemBlock2D<double>&  gtoPrimitives,
                     const CMemBlock2D<int32_t>& contrPattern,
                     const int32_t               angularMomentum)

    : _gtoPrimitives(gtoPrimitives)

    , _contrPattern(contrPattern)

    , _angularMomentum(angularMomentum)
{

}

CGtoBlock::CGtoBlock(const CMolecule&       molecule,
                     const CMolecularBasis& basis,
                     const int32_t          angularMomentum)

    : _angularMomentum(angularMomentum)
{
    // allocate primitives data

    auto npfuncs = basis.getNumberOfPrimitiveBasisFunctions(molecule,
                                                            _angularMomentum);
    
    if (npfuncs > 0)
    {
        _gtoPrimitives = CMemBlock2D<double>(npfuncs, 5);

        // allocate contraction data

        auto ncfuncs = basis.getNumberOfBasisFunctions(molecule,
                                                       _angularMomentum);
        
        auto angcomp = angmom::to_SphericalComponents(_angularMomentum);

        _contrPattern = CMemBlock2D<int32_t>(ncfuncs, 2 + angcomp);
        
        // determine partial dimensions of AO basis
        
        auto npartdim = basis.getPartialDimensionsOfBasis(molecule,
                                                          _angularMomentum);

        // set up pointers to molecular data
        
        auto molrx = molecule.getCoordinatesX();
        
        auto molry = molecule.getCoordinatesY();
        
        auto molrz = molecule.getCoordinatesZ();
        
        auto idselem = molecule.getIdsElemental();
        
        // contraction pattern data
        
        auto spos = _contrPattern.data(0);
        
        auto epos = _contrPattern.data(1);
        
        // primitives data
        
        auto gtoexps = _gtoPrimitives.data(0);
        
        auto gtonorms = _gtoPrimitives.data(1);
        
        auto coordsx = _gtoPrimitives.data(2);
        
        auto coordsy = _gtoPrimitives.data(3);
        
        auto coordsz = _gtoPrimitives.data(4);

        // loop over atoms in molecule

        int32_t icgto = 0;

        int32_t iprim = 0;

        for (int32_t i = 0; i < molecule.getNumberOfAtoms(); i++)
        {

            // get atom coordinates

            auto rx = molrx[i];

            auto ry = molry[i];

            auto rz = molrz[i];

            // loop over basis functions of i-th atom
            
            auto gtos = basis.getBasisFunctions(idselem[i], _angularMomentum);

            for (size_t j = 0; j < gtos.size(); j++)
            {
                auto nprim = gtos[j].getNumberOfPrimitiveFunctions();

                // set contraction pattern

                spos[icgto] = iprim;

                epos[icgto] = iprim + nprim;

                for (int32_t k = 0; k < angcomp; k++)
                {
                    auto pgtoidx = _contrPattern.data(2 + k);
                    
                    pgtoidx[icgto] = npartdim + k * ncfuncs + icgto;
                }

                // retrieve primitve exponents, norm. factors

                auto pexp  = gtos[j].getExponents();

                auto pnorm = gtos[j].getNormalizationFactors();

                // set up primitives data

                for (int32_t k = 0; k < nprim; k++)
                {
                    // assign exponent, norm. factor

                    gtoexps[iprim + k] = pexp[k];

                    gtonorms[iprim + k] = pnorm[k];

                    // assign atom coordinates

                    coordsx[iprim + k] = rx;

                    coordsy[iprim + k] = ry;

                    coordsz[iprim + k] = rz;
                }

                // update indexes

                iprim += nprim;

                icgto++;
            }
        }
    }
}

CGtoBlock::CGtoBlock(const CGtoBlock& source)

    : _angularMomentum(source._angularMomentum)

    , _contrPattern(source._contrPattern)

    , _gtoPrimitives(source._gtoPrimitives)
{

}

CGtoBlock::CGtoBlock(CGtoBlock&& source) noexcept

    : _angularMomentum(std::move(source._angularMomentum))

    , _contrPattern(std::move(source._contrPattern))

    , _gtoPrimitives(std::move(source._gtoPrimitives))
{

}

CGtoBlock::~CGtoBlock()
{

}

CGtoBlock&
CGtoBlock::operator=(const CGtoBlock& source)
{
    if (this == &source) return *this;

    _angularMomentum = source._angularMomentum;

    _contrPattern = source._contrPattern;

    _gtoPrimitives = source._gtoPrimitives;

    return *this;
}

CGtoBlock&
CGtoBlock::operator=(CGtoBlock&& source) noexcept
{
    if (this == &source) return *this;

    _angularMomentum = std::move(source._angularMomentum);

    _contrPattern = std::move(source._contrPattern);

    _gtoPrimitives = std::move(source._gtoPrimitives);

    return *this;
}

bool
CGtoBlock::operator==(const CGtoBlock& other) const
{
    if (_angularMomentum != other._angularMomentum) return false;

    if (_contrPattern != other._contrPattern) return false;

    if (_gtoPrimitives != other._gtoPrimitives) return false;

    return true;
}

bool
CGtoBlock::operator!=(const CGtoBlock& other) const
{
    return !(*this == other);
}

std::tuple<int32_t, int32_t>
CGtoBlock::compress(const CGtoBlock&         source,
                    const CMemBlock<double>& screeningFactors,
                    const double             screeningThreshold)
{
    if (_angularMomentum != source._angularMomentum)
    {
        return std::make_tuple(0,0);
    }
    
    // zero current data
    
    _gtoPrimitives.zero();
    
    _contrPattern.zero();
    
    // set up pointers to primitives data source
    
    auto srcexps = source.getExponents();
    
    auto srcfacts = source.getNormFactors();
    
    // set up primitive GTOs coordinates data source
    
    auto srcrx = source.getCoordinatesX();
    
    auto srcry = source.getCoordinatesY();
    
    auto srcrz = source.getCoordinatesZ();
    
    // set up pointers to contraction pattern data source 
    
    auto srcspos = source.getStartPositions();
    
    auto srcepos = source.getEndPositions();
    
    // set up pointer to screening factors
    
    auto sfacts = screeningFactors.data();
    
    // determine number of angular components
    
    auto angcomp = angmom::to_SphericalComponents(_angularMomentum);
    
    // set up pointers to primitives
    
    auto pexps = getExponents();
    
    auto pfacts = getNormFactors();
    
    // set up primitive GTOs coordinates
    
    auto prx = getCoordinatesX();
    
    auto pry = getCoordinatesY();
    
    auto prz = getCoordinatesZ();
    
    // set up contraction pattern data
    
    auto cspos = getStartPositions();
    
    auto cepos = getEndPositions();
    
    // primitive and contracted GTOs counters
    
    int32_t npgto = 0;
    
    int32_t ncgto = 0;
    
    // loop over contracted contraction pattern
    
    for (int32_t i = 0; i < _contrPattern.size(0); i++)
    {
        int32_t cprim = 0;
        
        // add primite GTOs data
        
        for (int32_t j = srcspos[i]; j < srcepos[i]; j++)
        {
            if (std::fabs(sfacts[j]) > screeningThreshold)
            {
                auto poff = npgto + cprim;
                
                pexps[poff] = srcexps[j];
                
                pfacts[poff] = srcfacts[j];
                
                prx[poff] = srcrx[j];
                
                pry[poff] = srcry[j];
                
                prz[poff] = srcrz[j];
                
                cprim++;
            }
        }
        
        // add GTOs contraction data
        
        if (cprim > 0)
        {
            cspos[ncgto] = npgto;
            
            cepos[ncgto] = npgto + cprim;
           
            // store GTOs indexes
            
            for (int32_t j = 0; j < angcomp; j++)
            {
                auto srcidx = source.getIdentifiers(j);
            
                auto curidx = getIdentifiers(j);
                
                curidx[ncgto] = srcidx[i];
            }
            
            // update counters
            
            npgto += cprim;
            
            ncgto++;
        }
    }
    
    return std::make_tuple(npgto, ncgto);
}

int32_t
CGtoBlock::getAngularMomentum() const
{
    return _angularMomentum;
}

bool
CGtoBlock::empty() const
{
    return (_contrPattern.size(0) == 0);
}

int32_t
CGtoBlock::getNumberOfPrimGtos() const
{
    return _gtoPrimitives.size(0);
}

int32_t
CGtoBlock::getNumberOfContrGtos() const
{
    return _contrPattern.size(0);
}

const int32_t*
CGtoBlock::getStartPositions() const
{
    return _contrPattern.data(0);
}

int32_t*
CGtoBlock::getStartPositions()
{
    return _contrPattern.data(0);
}

const int32_t*
CGtoBlock::getEndPositions() const
{
    return _contrPattern.data(1);
}

int32_t*
CGtoBlock::getEndPositions()
{
    return _contrPattern.data(1);
}

const int32_t*
CGtoBlock::getIdentifiers(const int32_t iComponent) const
{
    if (iComponent < angmom::to_SphericalComponents(_angularMomentum))
    {
        return _contrPattern.data(2 + iComponent);
    }
    
    return nullptr; 
}

int32_t*
CGtoBlock::getIdentifiers(const int32_t iComponent)
{
    if (iComponent < angmom::to_SphericalComponents(_angularMomentum))
    {
        return _contrPattern.data(2 + iComponent);
    }
    
    return nullptr;
}

const double*
CGtoBlock::getExponents() const
{
    return _gtoPrimitives.data(0);
}

double*
CGtoBlock::getExponents()
{
    return _gtoPrimitives.data(0);
}

const double*
CGtoBlock::getNormFactors() const
{
    return _gtoPrimitives.data(1);
}

double*
CGtoBlock::getNormFactors()
{
    return _gtoPrimitives.data(1);
}

const double*
CGtoBlock::getCoordinatesX() const
{
    return _gtoPrimitives.data(2);
}

double*
CGtoBlock::getCoordinatesX()
{
    return _gtoPrimitives.data(2);
}

const double*
CGtoBlock::getCoordinatesY() const
{
    return _gtoPrimitives.data(3);
}

double*
CGtoBlock::getCoordinatesY()
{
    return _gtoPrimitives.data(3);
}

const double*
CGtoBlock::getCoordinatesZ() const
{
    return _gtoPrimitives.data(4);
}

double*
CGtoBlock::getCoordinatesZ()
{
    return _gtoPrimitives.data(4);
}



std::ostream&
operator<<(      std::ostream& output,
           const CGtoBlock&    source)
{
    output << std::endl;

    output << "[CCGtoBlock (Object):" << &source << "]" << std::endl;

    output << "_angularMomentum: " << source._angularMomentum << std::endl;

    output << "_contrPattern: " << source._contrPattern << std::endl;

    output << "_gtoPrimitives: " << source._gtoPrimitives << std::endl;
    
    return output;
}

