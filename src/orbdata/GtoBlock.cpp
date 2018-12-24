//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
                     const CMemBlock<double>&    gtoNormFactors,
                     const CMemBlock2D<int32_t>& contrPattern,
                     const CMemBlock2D<int32_t>& indexPattern,
                     const int32_t               angularMomentum)

    : _angularMomentum(angularMomentum)

    , _contrPattern(contrPattern)

    , _indexPattern(indexPattern)

    , _gtoPrimitives(gtoPrimitives)

    , _gtoNormFactors(gtoNormFactors)
{

}

CGtoBlock::CGtoBlock(const CMolecule&       molecule,
                     const CMolecularBasis& basis,
                     const int32_t          angularMomentum)

    : CGtoBlock(molecule, basis, 0, molecule.getNumberOfAtoms(), angularMomentum)
{
   
}

CGtoBlock::CGtoBlock(const CMolecule&       molecule,
                     const CMolecularBasis& basis,
                     const int32_t          iAtom,
                     const int32_t          nAtoms,
                     const int32_t          angularMomentum)

    : _angularMomentum(angularMomentum)
{
    // allocate primitives data
    
    auto npfuncs = basis.getNumberOfPrimitiveBasisFunctions(molecule,
                                                            iAtom, nAtoms,
                                                            _angularMomentum);
    
    if (npfuncs > 0)
    {
        _gtoPrimitives = CMemBlock2D<double>(npfuncs, 4);
        
        // allocate normalization factors
        
        auto npfacts = basis.getNumberOfNormalizationFactors(molecule,
                                                             iAtom, nAtoms,
                                                             _angularMomentum);
        
        _gtoNormFactors = CMemBlock<double>(npfacts); 
        
        // allocate contraction data
        
        auto ncfuncs = basis.getNumberOfBasisFunctions(molecule,
                                                       _angularMomentum);
        
        auto angcomp = angmom::to_SphericalComponents(_angularMomentum);
        
        auto ncdim = basis.getNumberOfBasisFunctions(molecule, iAtom, nAtoms,
                                                     _angularMomentum);
        
        auto nrdim = basis.getNumberOfReducedBasisFunctions(molecule,
                                                            iAtom, nAtoms,
                                                            _angularMomentum);
        
        _contrPattern = CMemBlock2D<int32_t>(nrdim, 5);
        
        _indexPattern = CMemBlock2D<int32_t>(ncdim, 2 + angcomp);
        
        // determine partial dimensions of AO basis
        
        auto npartdim = basis.getPartialDimensionsOfBasis(molecule,
                                                          _angularMomentum);
        
        // determine offset in contracted GTOs block
        
        auto ncoff = basis.getNumberOfBasisFunctions(molecule, 0, iAtom,
                                                     _angularMomentum);
        
        // set up pointers to molecular data
        
        auto molrx = molecule.getCoordinatesX();
        
        auto molry = molecule.getCoordinatesY();
        
        auto molrz = molecule.getCoordinatesZ();
        
        auto idselem = molecule.getIdsElemental();
        
        // contraction pattern for primitive GTOs
        
        auto sppos = _contrPattern.data(0);
        
        auto eppos = _contrPattern.data(1);
        
        auto idxatm = _contrPattern.data(2);
        
        auto scpos = _contrPattern.data(3);
        
        auto ecpos = _contrPattern.data(4);
        
        // primitives data
        
        auto gtoexps = _gtoPrimitives.data(0);
        
        auto coordsx = _gtoPrimitives.data(1);
        
        auto coordsy = _gtoPrimitives.data(2);
        
        auto coordsz = _gtoPrimitives.data(3);
        
        auto gtonorms = _gtoNormFactors.data();
        
        // indexing data
        
        auto sfpos = _indexPattern.data(0);
        
        auto efpos = _indexPattern.data(1);
        
        // loop over atoms in molecule
        
        int32_t icgto = 0;
        
        int32_t irgto = 0;
        
        int32_t iprim = 0;
        
        int32_t ifact = 0;
        
        for (int32_t i = iAtom; i < (iAtom + nAtoms); i++)
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
                
                auto ngfuncs = gtos[j].getNumberOfContractedFunctions();
                
                // set contraction pattern
                
                sppos[irgto] = iprim;
                
                eppos[irgto] = iprim + nprim;
                
                idxatm[irgto] = i;
                
                scpos[irgto] = icgto;
                
                ecpos[irgto] = icgto + ngfuncs;
                
                for (int32_t k = 0; k < ngfuncs; k++)
                {
                    for (int32_t l = 0; l < angcomp; l++)
                    {
                        auto pgtoidx = _indexPattern.data(2 + l);
                    
                        pgtoidx[icgto] = npartdim + l * ncfuncs + ncoff + icgto;
                    }
                    
                    sfpos[icgto] = ifact + k * nprim;
                    
                    efpos[icgto] = ifact + (k + 1) * nprim;
                    
                    icgto++;
                }
                
                // retrieve primitve exponents, norm. factors
                
                auto pexp  = gtos[j].getExponents();
                
                auto pnorm = gtos[j].getNormalizationFactors();
                
                // set up primitives data
                
                for (int32_t k = 0; k < nprim; k++)
                {
                    // assign exponent
                    
                    gtoexps[iprim + k] = pexp[k];
                    
                    // assign atom coordinates
                    
                    coordsx[iprim + k] = rx;
                    
                    coordsy[iprim + k] = ry;
                    
                    coordsz[iprim + k] = rz;
                }
                
                // set up normalization factors data
                
                for (int32_t k = 0; k <  ngfuncs * nprim; k++)
                {
                    // assign norm. factor
                    
                    gtonorms[ifact + k] = pnorm[k];
                }
                
                // update indexes
                
                iprim += nprim;
                
                ifact += ngfuncs * nprim;
                
                irgto++;
            }
        }
    }
}

CGtoBlock::CGtoBlock(const CGtoBlock& source)

    : _angularMomentum(source._angularMomentum)

    , _contrPattern(source._contrPattern)

    , _indexPattern(source._indexPattern)

    , _gtoPrimitives(source._gtoPrimitives)

    , _gtoNormFactors(source._gtoNormFactors)
{

}

CGtoBlock::CGtoBlock(CGtoBlock&& source) noexcept

    : _angularMomentum(std::move(source._angularMomentum))

    , _contrPattern(std::move(source._contrPattern))

    , _indexPattern(std::move(source._indexPattern))

    , _gtoPrimitives(std::move(source._gtoPrimitives))

    , _gtoNormFactors(std::move(source._gtoNormFactors))
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
    
    _indexPattern = source._indexPattern;

    _gtoPrimitives = source._gtoPrimitives;
    
    _gtoNormFactors = source._gtoNormFactors;

    return *this;
}

CGtoBlock&
CGtoBlock::operator=(CGtoBlock&& source) noexcept
{
    if (this == &source) return *this;

    _angularMomentum = std::move(source._angularMomentum);

    _contrPattern = std::move(source._contrPattern);
    
    _indexPattern = std::move(source._indexPattern);

    _gtoPrimitives = std::move(source._gtoPrimitives);
    
    _gtoNormFactors = std::move(source._gtoNormFactors);

    return *this;
}

bool
CGtoBlock::operator==(const CGtoBlock& other) const
{
    if (_angularMomentum != other._angularMomentum) return false;

    if (_contrPattern != other._contrPattern) return false;
    
    if (_indexPattern != other._indexPattern) return false;

    if (_gtoPrimitives != other._gtoPrimitives) return false;
    
    if (_gtoNormFactors != other._gtoNormFactors) return false;

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
    
    _gtoNormFactors.zero(); 
    
    _contrPattern.zero();
    
    _indexPattern.zero(); 
    
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
    
    auto sridxatm = source.getAtomicIdentifiers();
    
    auto srcscpos = source.getContrStartPositions();
    
    auto srcecpos = source.getContrEndPositions();
    
    auto srcsfpos = source.getNormFactorsStartPositions();
    
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
    
    auto cidxatm = getAtomicIdentifiers();
    
    auto cscpos = getContrStartPositions();
    
    auto cecpos = getContrEndPositions();
    
    auto csfpos = getNormFactorsStartPositions();
    
    auto cefpos = getNormFactorsEndPositions();
    
    // primitive and contracted GTOs counters
    
    int32_t npgto = 0;
    
    int32_t ngfunc = 0;
    
    int32_t nrgto = 0;
    
    int32_t ncgto = 0;
    
    // loop over contracted contraction pattern
    
    for (int32_t i = 0; i < _contrPattern.size(0); i++)
    {
        auto cgfunc = srcecpos[i] - srcscpos[i];
        
        int32_t cprim = 0;
        
        // add primite GTOs data
        
        for (int32_t j = srcspos[i]; j < srcepos[i]; j++)
        {
            if (sfacts[j] > screeningThreshold)
            {
                auto poff = npgto + cprim;
                
                pexps[poff] = srcexps[j];
                
                prx[poff] = srcrx[j];
                
                pry[poff] = srcry[j];
                
                prz[poff] = srcrz[j];
                
                cprim++;
            }
        }
        
        // add GTOs contraction data
        
        if (cprim > 0)
        {
            
            // contraction pattern
            
            cspos[nrgto] = npgto;
            
            cepos[nrgto] = npgto + cprim;
            
            cidxatm[nrgto] = sridxatm[i];
           
            cscpos[nrgto] = ncgto;
            
            cecpos[nrgto] = ncgto + cgfunc;
            
            // indexing data
            
            for (int32_t j = srcscpos[i]; j < srcecpos[i]; j++)
            {
                for (int32_t k = 0; k < angcomp; k++)
                {
                    auto srcidx = source.getIdentifiers(k);
            
                    auto curidx = getIdentifiers(k);
                
                    curidx[ncgto] = srcidx[j];
                }
                
                csfpos[ncgto] = ngfunc + (j - srcscpos[i]) * cprim;
                
                cefpos[ncgto] = ngfunc + (j - srcscpos[i] + 1) * cprim;
                
                ncgto++;
            }
            
            // normalization factors
            
            int32_t icprim = 0;
            
            for (int32_t j = srcspos[i]; j < srcepos[i]; j++)
            {
                if (sfacts[j] > screeningThreshold)
                {
                    auto joff = j - srcspos[i];
                    
                    for (int32_t k = srcscpos[i]; k < srcecpos[i]; k++)
                    {
                        auto koff = (k - srcscpos[i]) * cprim + icprim;
                        
                        auto soff = srcsfpos[k] + joff;
                        
                        pfacts[ngfunc + koff] = srcfacts[soff];
                    }
                    
                    icprim++;
                }
            }
            
            // update counters
            
            npgto += cprim;
            
            ngfunc += cprim * cgfunc;
            
            nrgto++;
        }
    }
    
    return std::make_tuple(npgto, nrgto);
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
CGtoBlock::getContrStartPositions() const
{
    return _contrPattern.data(3);
}

int32_t*
CGtoBlock::getContrStartPositions()
{
    return _contrPattern.data(3);
}

const int32_t*
CGtoBlock::getContrEndPositions() const
{
    return _contrPattern.data(4);
}

int32_t*
CGtoBlock::getContrEndPositions()
{
    return _contrPattern.data(4);
}

const int32_t*
CGtoBlock::getNormFactorsStartPositions() const
{
    return _indexPattern.data(0);
}

int32_t*
CGtoBlock::getNormFactorsStartPositions()
{
    return _indexPattern.data(0);
}

const int32_t*
CGtoBlock::getNormFactorsEndPositions() const
{
    return _indexPattern.data(1);
}

int32_t*
CGtoBlock::getNormFactorsEndPositions()
{
    return _indexPattern.data(1);
}

const int32_t*
CGtoBlock::getAtomicIdentifiers() const
{
    return _contrPattern.data(2);
}

int32_t*
CGtoBlock::getAtomicIdentifiers()
{
    return _contrPattern.data(2);
}

const int32_t*
CGtoBlock::getIdentifiers(const int32_t iComponent) const
{
    if (iComponent < angmom::to_SphericalComponents(_angularMomentum))
    {
        return _indexPattern.data(2 + iComponent);
    }
    
    return nullptr; 
}

int32_t*
CGtoBlock::getIdentifiers(const int32_t iComponent)
{
    if (iComponent < angmom::to_SphericalComponents(_angularMomentum))
    {
        return _indexPattern.data(2 + iComponent);
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
    return _gtoNormFactors.data();
}

double*
CGtoBlock::getNormFactors()
{
    return _gtoNormFactors.data();
}

const double*
CGtoBlock::getCoordinatesX() const
{
    return _gtoPrimitives.data(1);
}

double*
CGtoBlock::getCoordinatesX()
{
    return _gtoPrimitives.data(1);
}

const double*
CGtoBlock::getCoordinatesY() const
{
    return _gtoPrimitives.data(2);
}

double*
CGtoBlock::getCoordinatesY()
{
    return _gtoPrimitives.data(2);
}

const double*
CGtoBlock::getCoordinatesZ() const
{
    return _gtoPrimitives.data(3);
}

double*
CGtoBlock::getCoordinatesZ()
{
    return _gtoPrimitives.data(3);
}

int32_t
CGtoBlock::getMaxContractionDepth() const
{
    int32_t ndim = 0;
    
    auto spos = getStartPositions();
    
    auto epos = getEndPositions();
    
    for (int32_t i = 0; i < getNumberOfContrGtos(); i++)
    {
        auto cdim = epos[i] - spos[i];
        
        if (cdim > ndim) ndim = cdim;
    }
    
    return ndim; 
}

std::ostream&
operator<<(      std::ostream& output,
           const CGtoBlock&    source)
{
    output << std::endl;

    output << "[CCGtoBlock (Object):" << &source << "]" << std::endl;

    output << "_angularMomentum: " << source._angularMomentum << std::endl;

    output << "_contrPattern: " << source._contrPattern << std::endl;
    
    output << "_indexPattern: " << source._indexPattern << std::endl;

    output << "_gtoPrimitives: " << source._gtoPrimitives << std::endl;
    
    output << "_gtoNormFactors: " << source._gtoNormFactors << std::endl;
    
    return output;
}

