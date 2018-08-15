//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoPairsBlock.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "AngularMomentum.hpp"

CGtoPairsBlock::CGtoPairsBlock()

    : _braAngularMomentum(-1)

    , _ketAngularMomentum(-1)

    , _threshold(1.0e-13)

    , _nOriginalPrimPairs(0)

    , _nScreenedPrimPairs(0)

    , _nOriginalContrPairs(0)

    , _nScreenedContrPairs(0)
{
    
}

CGtoPairsBlock::CGtoPairsBlock(const CMemBlock2D<int32_t>& contrPattern,
                               const CMemBlock2D<double>&  pairFactors,
                               const int32_t               braAngularMomentum,
                               const int32_t               ketAngularMomentum,
                               const double                threshold)

    : _braAngularMomentum(braAngularMomentum)

    , _ketAngularMomentum(ketAngularMomentum)

    , _contrPattern(contrPattern)

    , _pairFactors(pairFactors)

    , _threshold(threshold)

    , _nOriginalPrimPairs(pairFactors.size(0))

    , _nScreenedPrimPairs(pairFactors.size(0))

    , _nOriginalContrPairs(contrPattern.size(0))

    , _nScreenedContrPairs(contrPattern.size(0))
{
    
}

CGtoPairsBlock::CGtoPairsBlock(const CGtoBlock& braGtoBlock,
                               const CGtoBlock& ketGtoBlock,
                               const double     threshold)
    : _braAngularMomentum(braGtoBlock.getAngularMomentum())

    , _ketAngularMomentum(ketGtoBlock.getAngularMomentum())

    , _threshold(threshold)
{
    // set up dimensions of GTOs
    
    auto bpgto = braGtoBlock.getNumberOfPrimGtos();
    
    auto kpgto = ketGtoBlock.getNumberOfPrimGtos();
    
    auto bcgto = braGtoBlock.getNumberOfContrGtos();
    
    auto kcgto = ketGtoBlock.getNumberOfContrGtos();
    
    // set up pointers to GTOs block data on bra side
    
    auto bexp = braGtoBlock.getExponents();
    
    auto bnorm = braGtoBlock.getNormFactors();

    auto brx = braGtoBlock.getCoordinatesX();
    
    auto bry = braGtoBlock.getCoordinatesY();
    
    auto brz = braGtoBlock.getCoordinatesZ();
    
    auto bspos = braGtoBlock.getStartPositions();
    
    auto bepos = braGtoBlock.getEndPositions();
    
    // set up pointers to GTOs block data on ket side
    
    auto kexp = ketGtoBlock.getExponents();
    
    auto knorm = ketGtoBlock.getNormFactors();
    
    auto krx = ketGtoBlock.getCoordinatesX();
    
    auto kry = ketGtoBlock.getCoordinatesY();
    
    auto krz = ketGtoBlock.getCoordinatesZ();
    
    auto kspos = ketGtoBlock.getStartPositions();
    
    auto kepos = ketGtoBlock.getEndPositions();
    
    // determine if bra and ket GTOs block objects are the same
    
    bool symbk = (braGtoBlock == ketGtoBlock);
    
    // set up temporary memory blocks
    
    CMemBlock2D<double> (bpgto * kpgto, 4); 
    
    // fetch pi value
    
    auto fpi = mathconst::getPiValue();
    
    // initialize temporary storage for primitive pairs data
    
    CMemBlock2D<double> ppfacts(bpgto * kpgto, 13);
    
    // set up pointers to Obara-Saika prefactors
    
    auto ppfx = ppfacts.data(0);
    
    auto ppfi = ppfacts.data(1);
    
    auto ppfz = ppfacts.data(2);
    
    auto ppss = ppfacts.data(3);
    
    // set up pointers to P centers
    
    auto pprpx = ppfacts.data(4);
    
    auto pprpy = ppfacts.data(5);

    auto pprpz = ppfacts.data(6);
    
    // set up pointers to R(PA) distances
    
    auto pprpax = ppfacts.data(7);
    
    auto pprpay = ppfacts.data(8);
    
    auto pprpaz = ppfacts.data(9);
    
    // set up pointers to R(PB) distances
    
    auto pprpbx = ppfacts.data(10);
    
    auto pprpby = ppfacts.data(11);
    
    auto pprpbz = ppfacts.data(12);
    
    // initialize temporary storage for contracted pairs data
    
    auto bang = angmom::to_SphericalComponents(_braAngularMomentum);
    
    auto kang = angmom::to_SphericalComponents(_ketAngularMomentum);
    
    CMemBlock2D<int32_t> ppidx(bcgto * kcgto, 2 + bang + kang);
    
    auto ppspos = ppidx.data(0);
    
    auto ppepos = ppidx.data(1);

    // initialize number of pairs
    
    _nOriginalContrPairs = 0;
    
    _nOriginalPrimPairs  = 0;
    
    // loop over contracted GTOs
    
    int32_t idxpgto = 0;
    
    int32_t idxcgto = 0;
    
    for (int32_t i = 0; i < bcgto; i++)
    {
        int32_t joff = (symbk) ? i : 0;
        
        for (int32_t j = joff; j < kcgto; j++)
        {
            // construct pair of contracted GTOs
            
            int32_t nprim = 0;
            
            for (int32_t k = bspos[i]; k < bepos[i]; k++)
            {
                // A center coordinates
                
                auto rax = brx[k];
                
                auto ray = bry[k];
                
                auto raz = brz[k];
                
                // A center primitive GTO
                
                auto fae = bexp[k];
                
                auto fan = bnorm[k];
                
                for (int32_t l = kspos[j];  l < kepos[j]; l++)
                {
                    // R(AB) = A - B
                    
                    auto abx = rax - krx[l];
                    
                    auto aby = ray - kry[l];
                    
                    auto abz = raz - krz[l];
                    
                    // Obara-Saika factors
                    
                    auto fx = fae + kexp[l];
                    
                    auto fi = 1.0 / fx;
                    
                    auto fz = fae * kexp[l] * fi;
                    
                    // pair overlap value
                    
                    auto fovl = fan * knorm[l] * std::pow(fpi * fi, 1.50)
                    
                              * std::exp(-fz * (abx * abx + aby * aby + abz * abz));
                    
                    // update number of orginal primitive pairs
                    
                    _nOriginalPrimPairs++;
                    
                    // add primitive pair to list of primitive pairs
                    
                    if (std::fabs(fovl) > _threshold)
                    {
                        // Obara-Saika factors
                        
                        ppfx[idxpgto + nprim] = fx;
                        
                        ppfi[idxpgto + nprim] = fi;
                        
                        ppfz[idxpgto + nprim] = fz;
                        
                        ppss[idxpgto + nprim] = fovl;
                        
                        // P center
                        
                        pprpx[idxpgto + nprim] = fi * (fae * rax + kexp[l] * krx[l]);
                        
                        pprpy[idxpgto + nprim] = fi * (fae * ray + kexp[l] * kry[l]);
                        
                        pprpz[idxpgto + nprim] = fi * (fae * raz + kexp[l] * krz[l]);
                        
                        // R(PA) distances
                        
                        pprpax[idxpgto + nprim] = pprpx[idxpgto + nprim] - rax;
                        
                        pprpay[idxpgto + nprim] = pprpy[idxpgto + nprim] - ray;
                        
                        pprpaz[idxpgto + nprim] = pprpz[idxpgto + nprim] - raz;
                        
                        // R(PB) distances
                        
                        pprpbx[idxpgto + nprim] = pprpx[idxpgto + nprim] - krx[l];
                        
                        pprpby[idxpgto + nprim] = pprpy[idxpgto + nprim] - kry[l];
                        
                        pprpbz[idxpgto + nprim] = pprpz[idxpgto + nprim] - krz[l];

                        // update local primitives counter 
                        
                        nprim++;
                    }
                }
            }
            
            // update numnber of orginal contracted pairs
            
            _nOriginalContrPairs++;
            
            // add contracted GTO data
            
            if (nprim > 0)
            {
                // pair position in primitive pairs vector
                
                ppspos[idxcgto] = idxpgto;
                
                ppepos[idxcgto] = idxpgto + nprim;
                
                // bra indexes of pair
                
                for (int32_t m = 0; m < bang; m++)
                {
                    auto bidx = braGtoBlock.getIdentifiers(m);
                    
                    auto cidx = ppidx.data(2 +  m);
                    
                    cidx[idxcgto] = bidx[i];
                }
                
                // ket indexes of pair
                
                for (int32_t m = 0; m < kang; m++)
                {
                    auto kidx = ketGtoBlock.getIdentifiers(m);
                    
                    auto cidx = ppidx.data(2 + bang + m);
                    
                    cidx[idxcgto] = kidx[j];
                }
                
                // update indexes
                
                idxpgto += nprim;
                
                idxcgto++;
            }
        }
    }
    
    // set up screened dimensions
    
    _nScreenedPrimPairs = idxpgto;
    
    _nScreenedContrPairs = idxcgto;
    
    // copy pairs data from temporary storage
    
    _pairFactors = ppfacts.slice(0, _nScreenedPrimPairs);
    
    _contrPattern = ppidx.slice(0, _nScreenedContrPairs);
}

CGtoPairsBlock::CGtoPairsBlock(const CGtoBlock& gtoBlock,
                               const double     threshold)

    : CGtoPairsBlock(gtoBlock, gtoBlock, threshold)
{
    
}

CGtoPairsBlock::CGtoPairsBlock(const CGtoPairsBlock& source)

    : _braAngularMomentum(source._braAngularMomentum)

    , _ketAngularMomentum(source._ketAngularMomentum)

    , _contrPattern(source._contrPattern)

    , _pairFactors(source._pairFactors)

    , _threshold(source._threshold)

    , _nOriginalPrimPairs(source._nOriginalPrimPairs)

    , _nScreenedPrimPairs(source._nScreenedPrimPairs)

    , _nOriginalContrPairs(source._nOriginalContrPairs)

    , _nScreenedContrPairs(source._nScreenedContrPairs)
{
    
}

CGtoPairsBlock::CGtoPairsBlock(CGtoPairsBlock&& source) noexcept

    : _braAngularMomentum(std::move(source._braAngularMomentum))

    , _ketAngularMomentum(std::move(source._ketAngularMomentum))

    , _contrPattern(std::move(source._contrPattern))

    , _pairFactors(std::move(source._pairFactors))

    , _threshold(std::move(source._threshold))

    , _nOriginalPrimPairs(std::move(source._nOriginalPrimPairs))

    , _nScreenedPrimPairs(std::move(source._nScreenedPrimPairs))

    , _nOriginalContrPairs(std::move(source._nOriginalContrPairs))

    , _nScreenedContrPairs(std::move(source._nScreenedContrPairs))
{
    
}

CGtoPairsBlock::~CGtoPairsBlock()
{
    
}

CGtoPairsBlock&
CGtoPairsBlock::operator=(const CGtoPairsBlock& source)
{
    if (this == &source) return *this;
    
    _braAngularMomentum = source._braAngularMomentum;
    
    _ketAngularMomentum = source._ketAngularMomentum;
    
    _contrPattern = source._contrPattern;
    
    _pairFactors = source._pairFactors;
    
    _threshold = source._threshold;
    
    _nOriginalPrimPairs = source._nOriginalPrimPairs;
    
    _nScreenedPrimPairs = source._nScreenedPrimPairs;
    
    _nOriginalContrPairs = source._nOriginalContrPairs;
    
    _nScreenedContrPairs = source._nScreenedContrPairs;
    
    return *this;
}

CGtoPairsBlock&
CGtoPairsBlock::operator=(CGtoPairsBlock&& source) noexcept
{
    if (this == &source) return *this;
    
    _braAngularMomentum = std::move(source._braAngularMomentum);
    
    _ketAngularMomentum = std::move(source._ketAngularMomentum);
    
    _contrPattern = std::move(source._contrPattern);
    
    _pairFactors = std::move(source._pairFactors);
    
    _threshold = std::move(source._threshold);
    
    _nOriginalPrimPairs = std::move(source._nOriginalPrimPairs);
    
    _nScreenedPrimPairs = std::move(source._nScreenedPrimPairs);
    
    _nOriginalContrPairs = std::move(source._nOriginalContrPairs);
    
    _nScreenedContrPairs = std::move(source._nScreenedContrPairs);
    
    return *this;
}

bool
CGtoPairsBlock::operator==(const CGtoPairsBlock& other) const
{
    if (_braAngularMomentum != other._braAngularMomentum) return false;
    
    if (_ketAngularMomentum != other._ketAngularMomentum) return false;
    
    if (_contrPattern != other._contrPattern) return false;
    
    if (_pairFactors != other._pairFactors) return false;
    
    if (std::fabs(_threshold - other._threshold) > 1.0e-13) return false;
    
    if (_nOriginalPrimPairs != other._nOriginalPrimPairs) return false;
    
    if (_nScreenedPrimPairs != other._nScreenedPrimPairs) return false;
    
    if (_nOriginalContrPairs != other._nOriginalContrPairs) return false;
    
    if (_nScreenedContrPairs != other._nScreenedContrPairs) return false;
    
    return true;
}

int32_t
CGtoPairsBlock::getBraAngularMomentum() const
{
    return _braAngularMomentum;
}

int32_t
CGtoPairsBlock::getKetAngularMomentum() const
{
    return _ketAngularMomentum;
}

bool
CGtoPairsBlock::empty() const
{
    return (_nScreenedContrPairs == 0);
}

bool
CGtoPairsBlock::operator!=(const CGtoPairsBlock& other) const
{
    return !(*this == other);
}

const double*
CGtoPairsBlock::getFactorsXi() const
{
    return _pairFactors.data(0);
}

const double*
CGtoPairsBlock::getFactorsOneOverXi() const
{
    return _pairFactors.data(1);
}

const double*
CGtoPairsBlock::getFactorsZeta() const
{
    return _pairFactors.data(2);
}

const double*
CGtoPairsBlock::getOverlaps() const
{
    return _pairFactors.data(3); 
}

const double*
CGtoPairsBlock::getCoordinatesPX() const
{
    return _pairFactors.data(4);
}

const double*
CGtoPairsBlock::getCoordinatesPY() const
{
    return _pairFactors.data(5);
}

const double*
CGtoPairsBlock::getCoordinatesPZ() const
{
    return _pairFactors.data(6);
}

const double*
CGtoPairsBlock::getDistancesPAX() const
{
    return _pairFactors.data(7);
}

const double*
CGtoPairsBlock::getDistancesPAY() const
{
    return _pairFactors.data(8);
}

const double*
CGtoPairsBlock::getDistancesPAZ() const
{
    return _pairFactors.data(9);
}

const double*
CGtoPairsBlock::getDistancesPBX() const
{
    return _pairFactors.data(10);
}

const double*
CGtoPairsBlock::getDistancesPBY() const
{
    return _pairFactors.data(11);
}

const double*
CGtoPairsBlock::getDistancesPBZ() const
{
    return _pairFactors.data(12);
}

const int32_t*
CGtoPairsBlock::getStartPositions() const
{
    return _contrPattern.data(0);
}

const int32_t*
CGtoPairsBlock::getEndPositions() const
{
    return _contrPattern.data(1);
}

const int32_t*
CGtoPairsBlock::getBraIdentifiers(const int32_t iComponent) const
{
    return _contrPattern.data(2 + iComponent);
}

const int32_t*
CGtoPairsBlock::getKetIdentifiers(const int32_t iComponent) const
{
    auto bang = angmom::to_SphericalComponents(_braAngularMomentum);
    
    return _contrPattern.data(2 + bang  + iComponent);
}

int32_t
CGtoPairsBlock::getNumberOfOriginalPrimPairs() const
{
    return _nOriginalPrimPairs;
}

int32_t
CGtoPairsBlock::getNumberOfScreenedPrimPairs() const
{
    return _nScreenedPrimPairs;
}

int32_t
CGtoPairsBlock::getNumberOfOriginalContrPairs() const
{
    return _nOriginalContrPairs;
}

int32_t
CGtoPairsBlock::getNumberOfScreenedContrPairs() const
{
    return _nScreenedContrPairs; 
}

std::ostream&
operator<<(      std::ostream&   output,
           const CGtoPairsBlock& source)
{
    output << std::endl;
    
    output << "[CCGtoPairsBlock (Object):" << &source << "]" << std::endl;
    
    output << "_braAngularMomentum: " << source._braAngularMomentum << std::endl;
    
    output << "_ketAngularMomentum: " << source._ketAngularMomentum << std::endl;
    
    output << "_contrPattern: " << source._contrPattern << std::endl;
    
    output << "_pairFactors: " << source._pairFactors << std::endl;
    
    output << "_threshold: " << source._threshold << std::endl;
    
    output << "_nOriginalPrimPairs: " << source._nOriginalPrimPairs << std::endl;
    
    output << "_nScreenedPrimPairs: " << source._nScreenedPrimPairs << std::endl;
    
    output << "_nOriginalContrPairs: " << source._nOriginalContrPairs << std::endl;
    
    output << "_nScreenedContrPairs: " << source._nScreenedContrPairs << std::endl;
    
    return output;
}

