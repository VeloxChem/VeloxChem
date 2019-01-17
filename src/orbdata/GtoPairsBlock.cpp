//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GtoPairsBlock.hpp"

#include <cmath>

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "MathConst.hpp"
#include "AngularMomentum.hpp"
#include "StringFormat.hpp"

CGtoPairsBlock::CGtoPairsBlock()

    : _braAngularMomentum(-1)

    , _ketAngularMomentum(-1)

    , _contrPattern(CMemBlock2D<int32_t>())

    , _indexPattern(CMemBlock2D<int32_t>())

    , _contrFactors(CMemBlock2D<double>())

    , _pairFactors(CMemBlock2D<double>())

    , _normFactors(CMemBlock<double>())

    , _threshold(1.0e-13)

    , _nOriginalPrimPairs(0)

    , _nScreenedPrimPairs(0)

    , _nOriginalRedContrPairs(0)

    , _nScreenedRedContrPairs(0)

    , _nOriginalContrPairs(0)

    , _nScreenedContrPairs(0)
{
    
}

CGtoPairsBlock::CGtoPairsBlock(const CMemBlock2D<int32_t>& contrPattern,
                               const CMemBlock2D<int32_t>& indexPattern,
                               const CMemBlock2D<double>&  contrFactors, 
                               const CMemBlock2D<double>&  pairFactors,
                               const CMemBlock<double>&    normFactors,
                               const int32_t               braAngularMomentum,
                               const int32_t               ketAngularMomentum,
                               const double                threshold)

    : _braAngularMomentum(braAngularMomentum)

    , _ketAngularMomentum(ketAngularMomentum)

    , _contrPattern(contrPattern)

    , _indexPattern(indexPattern)

    , _contrFactors(contrFactors)

    , _pairFactors(pairFactors)

    , _normFactors(normFactors)

    , _threshold(threshold)

    , _nOriginalPrimPairs(pairFactors.size(0))

    , _nScreenedPrimPairs(pairFactors.size(0))

    , _nOriginalRedContrPairs(contrPattern.size(0))

    , _nScreenedRedContrPairs(contrPattern.size(0))

    , _nOriginalContrPairs(indexPattern.size(0))

    , _nScreenedContrPairs(indexPattern.size(0))
{
    
}

CGtoPairsBlock::CGtoPairsBlock(const CGtoBlock& braGtoBlock,
                               const CGtoBlock& ketGtoBlock,
                               const double     threshold)

    : _braAngularMomentum(braGtoBlock.getAngularMomentum())

    , _ketAngularMomentum(ketGtoBlock.getAngularMomentum())

    , _contrPattern(CMemBlock2D<int32_t>())

    , _indexPattern(CMemBlock2D<int32_t>())

    , _contrFactors(CMemBlock2D<double>())

    , _pairFactors(CMemBlock2D<double>())

    , _normFactors(CMemBlock<double>())

    , _threshold(threshold)

    , _nOriginalPrimPairs(0)

    , _nScreenedPrimPairs(0)

    , _nOriginalRedContrPairs(0)

    , _nScreenedRedContrPairs(0)

    , _nOriginalContrPairs(0)

    , _nScreenedContrPairs(0)
{
    // set up dimensions of GTOs
    
    auto bpgto = braGtoBlock.getNumberOfPrimGtos();
    
    auto kpgto = ketGtoBlock.getNumberOfPrimGtos();
    
    auto brgto = braGtoBlock.getNumberOfRedContrGtos();
    
    auto krgto = ketGtoBlock.getNumberOfRedContrGtos();
    
    auto bcgto = braGtoBlock.getNumberOfContrGtos();
    
    auto kcgto = ketGtoBlock.getNumberOfContrGtos();
    
    auto bnfacts = braGtoBlock.getNumberOfNormFactors();
    
    auto knfacts = ketGtoBlock.getNumberOfNormFactors();
    
    // set up pointers to GTOs block data on bra side
    
    auto bexp = braGtoBlock.getExponents();
    
    auto bnorm = braGtoBlock.getNormFactors();

    auto brx = braGtoBlock.getCoordinatesX();
    
    auto bry = braGtoBlock.getCoordinatesY();
    
    auto brz = braGtoBlock.getCoordinatesZ();
    
    auto bspos = braGtoBlock.getStartPositions();
    
    auto bepos = braGtoBlock.getEndPositions();
    
    auto bcspos = braGtoBlock.getContrStartPositions();
    
    auto bcepos = braGtoBlock.getContrEndPositions();
    
    auto bfspos = braGtoBlock.getNormFactorsStartPositions();
    
    // set up pointers to GTOs block data on ket side
    
    auto kexp = ketGtoBlock.getExponents();
    
    auto knorm = ketGtoBlock.getNormFactors();
    
    auto krx = ketGtoBlock.getCoordinatesX();
    
    auto kry = ketGtoBlock.getCoordinatesY();
    
    auto krz = ketGtoBlock.getCoordinatesZ();
    
    auto kspos = ketGtoBlock.getStartPositions();
    
    auto kepos = ketGtoBlock.getEndPositions();
    
    auto kcspos = ketGtoBlock.getContrStartPositions();
    
    auto kcepos = ketGtoBlock.getContrEndPositions();
    
    auto kfspos = ketGtoBlock.getNormFactorsStartPositions();
    
    // determine if bra and ket GTOs block objects are the same
    
    bool symbk = (braGtoBlock == ketGtoBlock);
    
    // fetch pi value
    
    auto fpi = mathconst::getPiValue();
    
    // initialize temporary storage for primitive pairs data
    
    CMemBlock2D<double> ppfacts(bpgto * kpgto, 22);
    
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
    
    // set up pointers to A centers
    
    auto pprax = ppfacts.data(13);
    
    auto ppray = ppfacts.data(14);
    
    auto ppraz = ppfacts.data(15);
    
    // set up pointers to B centers
    
    auto pprbx = ppfacts.data(16);
    
    auto pprby = ppfacts.data(17);
    
    auto pprbz = ppfacts.data(18);
    
    // set up pointers to R(AB) distances
    
    auto pprabx = ppfacts.data(19);
    
    auto ppraby = ppfacts.data(20);
    
    auto pprabz = ppfacts.data(21);
    
    // initialize temporary storage for contracted pairs data
    
    auto bang = angmom::to_SphericalComponents(_braAngularMomentum);
    
    auto kang = angmom::to_SphericalComponents(_ketAngularMomentum);
    
    CMemBlock2D<int32_t> ppidx(brgto * krgto, 4);
    
    CMemBlock2D<int32_t> pcidx(bcgto * kcgto, 2 + bang + kang);
    
    CMemBlock2D<double> pcfacts(brgto * krgto, 3);
    
    CMemBlock<double> pnfacts(bnfacts * knfacts); 
    
    // set up pointers to start and end positions
    
    auto ppspos = ppidx.data(0);
    
    auto ppepos = ppidx.data(1);
    
    auto ppcspos = ppidx.data(2);
    
    auto ppcepos = ppidx.data(3);
    
    auto ppfspos = pcidx.data(0);
    
    auto ppfepos = pcidx.data(1);
    
    // set up pointers to effective P coordinates
    
    auto pcrpx = pcfacts.data(0);
    
    auto pcrpy = pcfacts.data(1);
    
    auto pcrpz = pcfacts.data(2);

    // initialize number of pairs
    
    //_nOriginalContrPairs = 0;
    
    //_nOriginalPrimPairs  = 0;
    
    // loop over contracted GTOs
    
    int32_t idxpgto = 0;
   
    int32_t idxrgto = 0;
    
    int32_t idxcgto = 0;
    
    int32_t idxfact = 0;
    
    for (int32_t i = 0; i < brgto; i++)
    {
        int32_t joff = (symbk) ? i : 0;
        
        auto bcdim = bcepos[i] - bcspos[i];
        
        for (int32_t j = joff; j < krgto; j++)
        {
            auto kcdim = kcepos[j] - kcspos[j];
            
            auto bkcdim = (symbk && (i == j)) ? bcdim * (bcdim + 1) / 2 : bcdim * kcdim;
            
            // effective averaged P coordinates and normalization factor
            
            double repx = 0.0;
            
            double repy = 0.0;
            
            double repz = 0.0;
            
            double epnrm = 0.0;
            
            // construct pair of contracted GTOs
            
            int32_t nprim = 0;
            
            std::vector<double> bkfacts;
            
            for (int32_t k = bspos[i]; k < bepos[i]; k++)
            {
                // A center coordinates
            
                auto rax = brx[k];
            
                auto ray = bry[k];
            
                auto raz = brz[k];
            
                // A center primitive GTO
            
                auto fae = bexp[k];
            
                auto bfan = braGtoBlock.getMaxNormFactor(i, k - bspos[i]);
                
                for (int32_t l = kspos[j]; l < kepos[j]; l++)
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
                    
                    auto kfan = ketGtoBlock.getMaxNormFactor(j, l - kspos[j]);
            
                    auto fovl = std::pow(fpi * fi, 1.50)
                    
                              * std::exp(-fz * (abx * abx + aby * aby + abz * abz));
            
                    // update number of orginal primitive pairs
            
                    _nOriginalPrimPairs++;
            
                    // add primitive pair to list of primitive pairs
            
                    if (std::fabs(bfan * kfan * fovl) > _threshold)
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
            
                        // A Center
            
                        pprax[idxpgto + nprim] = rax;
            
                        ppray[idxpgto + nprim] = ray;
            
                        ppraz[idxpgto + nprim] = raz;
            
                        // B Center
            
                        pprbx[idxpgto + nprim] = krx[l];
            
                        pprby[idxpgto + nprim] = kry[l];
            
                        pprbz[idxpgto + nprim] = krz[l];
            
                        // R(AB) distances
            
                        pprabx[idxpgto + nprim] = abx;
            
                        ppraby[idxpgto + nprim] = aby;
            
                        pprabz[idxpgto + nprim] = abz;
            
                        // effective P coordinates
            
                        auto cab = std::fabs(bfan * kfan);
            
                        repx += cab * pprpx[idxpgto + nprim];
            
                        repy += cab * pprpy[idxpgto + nprim];
            
                        repz += cab * pprpz[idxpgto + nprim];
            
                        epnrm += cab;
                        
                        // combine normalization factors
            
                        for (int32_t m = bcspos[i]; m < bcepos[i]; m++)
                        {
                            int32_t noff = (symbk && (i == j)) ? m : kcspos[j];
                            
                            auto bfact = bnorm[bfspos[m] + (k - bspos[i])];
                            
                            for (int32_t n = noff; n < kcepos[j]; n++)
                            {
                                auto kfact = knorm[kfspos[n] + (l - kspos[j])];
                                
                                bkfacts.push_back(bfact * kfact);
                            }
                        }
                        
                        // update local primitives counter
                        
                        nprim++;
                    }
                }
            }
            
            // update numnber of orginal contracted pairs
            
            _nOriginalRedContrPairs++;
            
            _nOriginalContrPairs += bkcdim;

            // add contracted GTO data

            if (nprim > 0)
            {
                // pair position in primitive pairs vector

                ppspos[idxrgto] = idxpgto;

                ppepos[idxrgto] = idxpgto + nprim;
                
                // pair position in contracted pairs vector
                
                ppcspos[idxrgto] = idxcgto;
                
                ppcepos[idxrgto] = idxcgto + bkcdim;
                
                // effective P coordinates
                
                epnrm = 1.0 / epnrm;
                
                pcrpx[idxrgto] = epnrm * repx;
                
                pcrpy[idxrgto] = epnrm * repy;
                
                pcrpz[idxrgto] = epnrm * repz;
                
                // normalization factors
                
                for (int32_t k = 0; k < nprim; k++)
                {
                    for (int32_t l = 0; l < bkcdim; l++)
                    {
                        pnfacts.at(idxfact + l * nprim + k) = bkfacts[k * bkcdim + l];
                    }
                }

                // indexing information
                
                int32_t iblock = 0;
                
                for (int32_t k = bcspos[i]; k < bcepos[i]; k++)
                {
                    int32_t loff = (symbk && (i == j)) ? k : kcspos[j];
                    
                    for (int32_t l = loff; l < kcepos[j]; l++)
                    {
                        ppfspos[idxcgto + iblock] = idxfact + iblock * nprim;
                        
                        ppfepos[idxcgto + iblock] = idxfact + (iblock + 1) * nprim;
                        
                        // bra indexes of pair
                        
                        for (int32_t m = 0; m < bang; m++)
                        {
                            auto bidx = braGtoBlock.getIdentifiers(m);
                            
                            auto cidx = pcidx.data(2 +  m);
                            
                            cidx[idxcgto + iblock] = bidx[k + iblock];
                        }
                        
                        // ket indexes of pair
                        
                        for (int32_t m = 0; m < kang; m++)
                        {
                            auto kidx = ketGtoBlock.getIdentifiers(m);
                            
                            auto cidx = pcidx.data(2 + bang + m);
                            
                            cidx[idxcgto + iblock] = kidx[l + iblock];
                        }
                        
                        iblock++;
                    }
                }
                
                // update indexes

                idxpgto += nprim;

                idxrgto++;
                
                idxcgto += bkcdim;
                
                idxfact += bkcdim * nprim;
            }
        }
    }
    
    // set up screened dimensions
    
    _nScreenedPrimPairs = idxpgto;
    
    _nScreenedRedContrPairs = idxrgto;
    
    _nScreenedContrPairs = idxcgto;
    
    // copy pairs data from temporary storage
    
    _pairFactors = ppfacts.slice(0, _nScreenedPrimPairs);
    
    _contrPattern = ppidx.slice(0, _nScreenedRedContrPairs);
    
    _indexPattern = pcidx.slice(0, _nScreenedContrPairs);
    
    _contrFactors = pcfacts.slice(0, _nScreenedRedContrPairs);
    
    _normFactors = pnfacts.slice(0, idxfact);
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

    , _indexPattern(source._indexPattern)

    , _contrFactors(source._contrFactors)

    , _pairFactors(source._pairFactors)

    , _normFactors(source._normFactors)

    , _threshold(source._threshold)

    , _nOriginalPrimPairs(source._nOriginalPrimPairs)

    , _nScreenedPrimPairs(source._nScreenedPrimPairs)

    , _nOriginalRedContrPairs(source._nOriginalRedContrPairs)

    , _nScreenedRedContrPairs(source._nScreenedRedContrPairs)

    , _nOriginalContrPairs(source._nOriginalContrPairs)

    , _nScreenedContrPairs(source._nScreenedContrPairs)
{
    
}

CGtoPairsBlock::CGtoPairsBlock(CGtoPairsBlock&& source) noexcept

    : _braAngularMomentum(std::move(source._braAngularMomentum))

    , _ketAngularMomentum(std::move(source._ketAngularMomentum))

    , _contrPattern(std::move(source._contrPattern))

    , _indexPattern(std::move(source._indexPattern))

    , _contrFactors(std::move(source._contrFactors))

    , _pairFactors(std::move(source._pairFactors))

    , _normFactors(std::move(source._normFactors))

    , _threshold(std::move(source._threshold))

    , _nOriginalPrimPairs(std::move(source._nOriginalPrimPairs))

    , _nScreenedPrimPairs(std::move(source._nScreenedPrimPairs))

    , _nOriginalRedContrPairs(std::move(source._nOriginalRedContrPairs))

    , _nScreenedRedContrPairs(std::move(source._nScreenedRedContrPairs))

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
    
    _indexPattern = source._indexPattern;
    
    _contrFactors = source._contrFactors;
    
    _pairFactors = source._pairFactors;
    
    _normFactors = source._normFactors;
    
    _threshold = source._threshold;
    
    _nOriginalPrimPairs = source._nOriginalPrimPairs;
    
    _nScreenedPrimPairs = source._nScreenedPrimPairs;
    
    _nOriginalRedContrPairs = source._nOriginalRedContrPairs;
    
    _nScreenedRedContrPairs = source._nScreenedRedContrPairs;
    
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
    
    _indexPattern = std::move(source._indexPattern);
    
    _contrFactors = std::move(source._contrFactors);
    
    _pairFactors = std::move(source._pairFactors);
    
    _normFactors = std::move(source._normFactors);
    
    _threshold = std::move(source._threshold);
    
    _nOriginalPrimPairs = std::move(source._nOriginalPrimPairs);
    
    _nScreenedPrimPairs = std::move(source._nScreenedPrimPairs);
    
    _nOriginalRedContrPairs = std::move(source._nOriginalRedContrPairs);
    
    _nScreenedRedContrPairs = std::move(source._nScreenedRedContrPairs);
    
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
    
    if (_indexPattern != other._indexPattern) return false;
    
    if (_contrFactors != other._contrFactors) return false;
    
    if (_pairFactors != other._pairFactors) return false;
    
     if (_normFactors != other._normFactors) return false;
    
    if (std::fabs(_threshold - other._threshold) > 1.0e-13) return false;
    
    if (_nOriginalPrimPairs != other._nOriginalPrimPairs) return false;
    
    if (_nScreenedPrimPairs != other._nScreenedPrimPairs) return false;
    
    if (_nOriginalRedContrPairs != other._nOriginalRedContrPairs) return false;
    
    if (_nScreenedRedContrPairs != other._nScreenedRedContrPairs) return false;
    
    if (_nOriginalContrPairs != other._nOriginalContrPairs) return false;
    
    if (_nScreenedContrPairs != other._nScreenedContrPairs) return false;
    
    return true;
}

bool
CGtoPairsBlock::operator!=(const CGtoPairsBlock& other) const
{
    return !(*this == other);
}

std::vector<CGtoPairsBlock>
CGtoPairsBlock::split(const int32_t nodes) const
{
    auto batchSize = _getBlockDimensions();
    
    // determine number of batches
    
    auto nbtch = _nScreenedRedContrPairs / batchSize;
    
    if ((_nScreenedRedContrPairs % batchSize) != 0) nbtch++;
  
    auto ntsk = omp_get_max_threads(); 
     
    if (nbtch > ntsk) nbtch = ntsk; 

    // rescale number of batches over compute nodes

    nbtch = static_cast<int32_t>(nbtch * std::sqrt(nodes));
    
    batchSize = _getBlockDimensions();
    
    if ((_nScreenedRedContrPairs / nbtch) < batchSize)
    {
        nbtch = _nScreenedRedContrPairs / batchSize;
        
        if ((_nScreenedRedContrPairs % batchSize) != 0) nbtch++;
    }

    // set up batches distribution pattern
    
    CMemBlock2D<int32_t> bblk(nbtch, 2);
    
    mpi::batches_pattern(bblk.data(0), _nScreenedRedContrPairs, nbtch);
    
    mathfunc::indexes(bblk.data(1), bblk.data(0), nbtch);
    
    // set up pointers to distribution pattern
    
    auto boff = bblk.data(1);
    
    auto bdim = bblk.data(0);
    
    // primitive space start and end positions
    
    auto spos = getStartPositions();
    
    auto epos = getEndPositions();
    
    // contracted space start and end positions
    
    auto cspos = getContrStartPositions();
    
    auto cepos = getContrEndPositions();
    
    // normalization factors space start and end positions
    
    auto fnspos = getNormFactorsStartPositions();
    
    auto fnepos = getNormFactorsEndPositions();
    
    // split GTOs pairs block
    
    std::vector<CGtoPairsBlock> ppvec;
    
    for (int32_t i = 0; i < nbtch; i++)
    {
        // slice reduced contracted GTOs pairs
        
        auto cidx = boff[i];
        
        auto cdim = bdim[i];
        
        auto cdat = _contrPattern.slice(cidx, cdim);
        
        auto cfac = _contrFactors.slice(cidx, cdim);
        
        // slice primitive GTOs pairs
        
        auto pidx = spos[cidx];
        
        auto pdim = epos[cidx + cdim - 1] - pidx;
        
        auto pdat = _pairFactors.slice(pidx, pdim);
        
        // slice indexing pattern for GTOs pairs
        
        auto cgidx = cspos[cidx];
        
        auto cgdim = cepos[cidx + cdim - 1] - cgidx;
        
        auto cgdat = _indexPattern.slice(cgidx, cgdim);
        
        // adjust contraction pattern for GTOs pairs
        
        auto scurpos = cdat.data(0);
        
        auto ecurpos = cdat.data(1);
        
        auto cscurpos = cdat.data(2);
        
        auto cecurpos = cdat.data(3);
        
        for (int32_t j = 0; j < cdim; j++)
        {
            scurpos[j] = spos[cidx + j] - pidx;
            
            ecurpos[j] = epos[cidx + j] - pidx;
            
            cscurpos[j] = cspos[cidx + j] - cgidx;
            
            cecurpos[j] = cepos[cidx + j] - cgidx;
        }
        
        // slice normalization factors for GTOs pairs
        
        auto cfnidx = fnspos[cspos[cidx]];
        
        auto cfndim = fnepos[cepos[cidx + cdim - 1] - 1] - cfnidx;
        
        auto cfndat = _normFactors.slice(cfnidx, cfndim);
        
        // adjust normalization factors indexing
        
        auto fnscurpos = cgdat.data(0);
        
        auto fnecurpos = cgdat.data(1);
        
        for (int32_t j = 0; j < cgdim; j++)
        {
            fnscurpos[j] = fnspos[cspos[cidx] + j] - cfnidx;
            
            fnecurpos[j] = fnepos[cspos[cidx] + j] - cfnidx;
        }
        
        // slice normalization factors for GTOs pairs
        
        ppvec.push_back(CGtoPairsBlock(cdat, cgdat, cfac, pdat, cfndat,
                                       _braAngularMomentum, _ketAngularMomentum,
                                       _threshold));
    }
    
    return ppvec;
}

CGtoPairsBlock
CGtoPairsBlock::pick(const int32_t iGtoPair) const
{
    if (iGtoPair < _nScreenedContrPairs)
    {
        // primitive space start and end positions
        
        auto spos = getStartPositions();
        
        auto epos = getEndPositions();
        
        // contracted space start and end positions
        
        auto cspos = getContrStartPositions();
        
        auto cepos = getContrEndPositions();
        
        // normalization factors space start and end positions
        
        auto fnspos = getNormFactorsStartPositions();
        
        auto fnepos = getNormFactorsEndPositions();
        
        // slice contracted GTOs pairs
        
        auto cdat = _contrPattern.slice(iGtoPair, 1);
        
        auto cfac = _contrFactors.slice(iGtoPair, 1);
        
        // slice primitive GTOs pairs
        
        auto pidx = spos[iGtoPair];
        
        auto pdim = epos[iGtoPair] - pidx;
        
        auto pdat = _pairFactors.slice(pidx, pdim);
        
        // slice indexing pattern for GTOs pairs
        
        auto cgidx = cspos[iGtoPair];
        
        auto cgdim = cepos[iGtoPair] - cgidx;
        
        auto cgdat = _indexPattern.slice(cgidx, cgdim);
        
        // adjust contraction pattern for GTOs pairs
        
        auto scurpos = cdat.data(0);
        
        auto ecurpos = cdat.data(1);
        
        auto cscurpos = cdat.data(2);
        
        auto cecurpos = cdat.data(3);
        
        scurpos[0] = 0;
            
        ecurpos[0] = pdim;
        
        cscurpos[0] = 0;
        
        cecurpos[0] = cgdim;
        
        // slice normalization factors for GTOs pairs
        
        auto cfnidx = fnspos[cspos[iGtoPair]];
        
        auto cfndim = fnepos[cepos[iGtoPair] - 1] - cfnidx;
        
        auto cfndat = _normFactors.slice(cfnidx, cfndim);
        
        // adjust normalization factors indexing
        
        auto fnscurpos = cgdat.data(0);
        
        auto fnecurpos = cgdat.data(1);
        
        for (int32_t j = 0; j < cgdim; j++)
        {
            fnscurpos[j] = fnspos[cspos[iGtoPair] + j] - cfnidx;
            
            fnecurpos[j] = fnepos[cspos[iGtoPair] + j] - cfnidx;
        }
        
        return CGtoPairsBlock(cdat, cgdat, cfac, pdat, cfndat,
                              _braAngularMomentum, _ketAngularMomentum,
                              _threshold);
    }
    
    return CGtoPairsBlock(); 
}

int32_t
CGtoPairsBlock::compress(const CGtoPairsBlock&     source,
                         const CMemBlock<int32_t>& screeningPattern,
                         const int32_t             nElements)
{
    // clear contrated and primitive pairs data
    
    _contrPattern.zero();
    
    _indexPattern.zero();
    
    _contrFactors.zero();
    
    _pairFactors.zero();
    
    _normFactors.zero();
    
    // set up pointers to source of primitive data
    
    auto sfx   = source.getFactorsXi();
    
    auto sfi   = source.getFactorsOneOverXi();
    
    auto sfz   = source.getFactorsZeta();
    
    auto sss   = source.getOverlaps();
    
    auto srpx  = source.getCoordinatesPX();
    
    auto srpy  = source.getCoordinatesPY();
    
    auto srpz  = source.getCoordinatesPZ();
    
    auto srpax = source.getDistancesPAX();
    
    auto srpay = source.getDistancesPAY();
    
    auto srpaz = source.getDistancesPAZ();
    
    auto srpbx = source.getDistancesPBX();
    
    auto srpby = source.getDistancesPBY();
    
    auto srpbz = source.getDistancesPBZ();
    
    auto srax  = source.getCoordinatesAX();
    
    auto sray  = source.getCoordinatesAY();
    
    auto sraz  = source.getCoordinatesAZ();
    
    auto srbx  = source.getCoordinatesBX();
    
    auto srby  = source.getCoordinatesBY();
    
    auto srbz  = source.getCoordinatesBZ();
    
    auto srabx = source.getDistancesABX();
    
    auto sraby = source.getDistancesABY();
    
    auto srabz = source.getDistancesABZ();
    
    // set up pointers to destination of primitive data
    
    auto dfx   = _pairFactors.data(0);
    
    auto dfi   = _pairFactors.data(1);
    
    auto dfz   = _pairFactors.data(2);
    
    auto dss   = _pairFactors.data(3);
    
    auto drpx  = _pairFactors.data(4);
    
    auto drpy  = _pairFactors.data(5);
    
    auto drpz  = _pairFactors.data(6);
    
    auto drpax = _pairFactors.data(7);
    
    auto drpay = _pairFactors.data(8);
    
    auto drpaz = _pairFactors.data(9);
    
    auto drpbx = _pairFactors.data(10);
    
    auto drpby = _pairFactors.data(11);
    
    auto drpbz = _pairFactors.data(12);
    
    auto drax  = _pairFactors.data(13);
    
    auto dray  = _pairFactors.data(14);
    
    auto draz  = _pairFactors.data(15);
    
    auto drbx  = _pairFactors.data(16);
    
    auto drby  = _pairFactors.data(17);
    
    auto drbz  = _pairFactors.data(18);
    
    auto drabx = _pairFactors.data(19);
    
    auto draby = _pairFactors.data(20);
    
    auto drabz = _pairFactors.data(21);
    
    // set up pointers to source contraction pattern
    
    auto sspos = source.getStartPositions();
    
    auto sepos = source.getEndPositions();
    
    auto scspos = source.getContrStartPositions();
    
    auto scepos = source.getContrEndPositions();
    
    auto sfspos = source.getNormFactorsStartPositions();
    
    auto sfepos = source.getNormFactorsEndPositions();
    
    // set up pointers to destination contraction pattern
    
    auto dspos = _contrPattern.data(0);
    
    auto depos = _contrPattern.data(1);
    
    auto dcspos = _contrPattern.data(2);
    
    auto dcepos = _contrPattern.data(3);
    
    auto dfspos = _indexPattern.data(0);
    
    auto dfepos = _indexPattern.data(1);
    
    // set up pointers to source effective P coordinates
    
    auto scrpx = source.getEffectiveCoordinatesPX();
    
    auto scrpy = source.getEffectiveCoordinatesPY();
    
    auto scrpz = source.getEffectiveCoordinatesPZ();
    
    // set up pointers to destination effective P coordinates
    
    auto dcrpx = _contrFactors.data(0);
    
    auto dcrpy = _contrFactors.data(1);
    
    auto dcrpz = _contrFactors.data(2);
    
    // set up pointers to normalization factors
    
    auto snfacts = source.getNormFactors(); 
    
    auto dnfacts = _normFactors.data();
    
    // set up angular momentum data
    
    auto bang = angmom::to_SphericalComponents(source.getBraAngularMomentum());
    
    auto kang = angmom::to_SphericalComponents(source.getKetAngularMomentum());
    
    // set up pairs counters
    
    int32_t nppairs = 0;
    
    int32_t nrpairs = 0;
    
    int32_t ncpairs = 0;
    
    int32_t nfpairs = 0;
    
    // generate set of screened pairs
    
    for (int32_t i = 0; i < nElements; i++)
    {
        if (screeningPattern.at(i) == 1)
        {
            // set start position
            
            dspos[nrpairs] = nppairs;
            
            dcspos[nrpairs] = ncpairs;
            
            // copy primitive data
            
            for (int32_t j = sspos[i]; j < sepos[i]; j++)
            {
                // Obara-Saika factors
                
                dfx[nppairs] = sfx[j];
                
                dfi[nppairs] = sfi[j];
                
                dfz[nppairs] = sfz[j];
                
                dss[nppairs] = sss[j];
                
                // P center
                
                drpx[nppairs] = srpx[j];
                
                drpy[nppairs] = srpy[j];
                
                drpz[nppairs] = srpz[j];
                
                // R(PA) distances
                
                drpax[nppairs] = srpax[j];
                
                drpay[nppairs] = srpay[j];
                
                drpaz[nppairs] = srpaz[j];
                
                // R(PB) distances
                
                drpbx[nppairs] = srpbx[j];
                
                drpby[nppairs] = srpby[j];
                
                drpbz[nppairs] = srpbz[j];
                
                // A Center
                
                drax[nppairs] = srax[j];
                
                dray[nppairs] = sray[j];
                
                draz[nppairs] = sraz[j];
                
                // B Center
                
                drbx[nppairs] = srbx[j];
                
                drby[nppairs] = srby[j];
                
                drbz[nppairs] = srbz[j];
                
                // R(AB) distances
                
                drabx[nppairs] = srabx[j];
                
                draby[nppairs] = sraby[j];
                
                drabz[nppairs] = srabz[j];
                
                // update primitive pairs counter
                
                nppairs++;
            }
            
            // copy indexing and normalization factors
            
            for (int32_t j = scspos[i]; j < scepos[i]; j++)
            {
                dfspos[ncpairs] = nfpairs;
                
                for (int32_t k = sfspos[j]; k < sfepos[j]; k++)
                {
                    dnfacts[nfpairs] = snfacts[k];
                    
                    nfpairs++;
                }
                
                dfepos[ncpairs] = nfpairs;
                
                // bra indexes of pair
                
                for (int32_t k = 0; k < bang; k++)
                {
                    auto bidx = source.getBraIdentifiers(k);
                    
                    auto cidx = _indexPattern.data(2 +  k);
                    
                    cidx[ncpairs] = bidx[j];
                }
                
                // ket indexes of pair
                
                for (int32_t k = 0; k < kang; k++)
                {
                    auto kidx = source.getKetIdentifiers(k);
                    
                    auto cidx = _indexPattern.data(2 + bang + k);
                    
                    cidx[ncpairs] = kidx[j];
                }
                
                ncpairs++;
            }
            
            // set end position
            
            depos[nrpairs] = nppairs;
            
            dcepos[nrpairs] = ncpairs;
            
            // effective P center
            
            dcrpx[nrpairs] = scrpx[i];
            
            dcrpy[nrpairs] = scrpy[i];
            
            dcrpz[nrpairs] = scrpz[i];
            
            // update contracted pairs counters
            
            nrpairs++;
        }
    }
    
    return nrpairs;
}

int32_t
CGtoPairsBlock::compress(const CGtoPairsBlock&     source,
                         const CMemBlock<int32_t>& screeningPattern)
{
    return compress(source, screeningPattern, screeningPattern.size());
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

const double*
CGtoPairsBlock::getCoordinatesAX() const
{
    return _pairFactors.data(13);
}

const double*
CGtoPairsBlock::getCoordinatesAY() const
{
    return _pairFactors.data(14);
}

const double*
CGtoPairsBlock::getCoordinatesAZ() const
{
    return _pairFactors.data(15);
}

const double*
CGtoPairsBlock::getCoordinatesBX() const
{
    return _pairFactors.data(16);
}

const double*
CGtoPairsBlock::getCoordinatesBY() const
{
    return _pairFactors.data(17);
}

const double*
CGtoPairsBlock::getCoordinatesBZ() const
{
    return _pairFactors.data(18);
}

const double*
CGtoPairsBlock::getDistancesABX() const
{
    return _pairFactors.data(19);
}

const double*
CGtoPairsBlock::getDistancesABY() const
{
    return _pairFactors.data(20);
}

const double*
CGtoPairsBlock::getDistancesABZ() const
{
    return _pairFactors.data(21);
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
CGtoPairsBlock::getContrStartPositions() const
{
    return _contrPattern.data(2);
}

const int32_t*
CGtoPairsBlock::getContrEndPositions() const
{
    return _contrPattern.data(3);
}

const int32_t*
CGtoPairsBlock::getNormFactorsStartPositions() const
{
    return _indexPattern.data(0);
}

const int32_t*
CGtoPairsBlock::getNormFactorsEndPositions() const
{
    return _indexPattern.data(1);
}

const int32_t*
CGtoPairsBlock::getBraIdentifiers(const int32_t iComponent) const
{
    return _indexPattern.data(2 + iComponent);
}

const int32_t*
CGtoPairsBlock::getKetIdentifiers(const int32_t iComponent) const
{
    auto bang = angmom::to_SphericalComponents(_braAngularMomentum);
    
    return _indexPattern.data(2 + bang  + iComponent);
}

const double*
CGtoPairsBlock::getNormFactors() const
{
    return _normFactors.data();
}

CMemBlock2D<double>
CGtoPairsBlock::getDistancesAB() const
{
    CMemBlock2D<double> rab(_nScreenedRedContrPairs, 3);
    
    // set up pointers to R(AB) distances in primitives data
    
    auto abx = getDistancesABX();
    
    auto aby = getDistancesABY();
    
    auto abz = getDistancesABZ();
    
    // set up pointers to R(AB) distances
    
    auto rabx = rab.data(0);
    
    auto raby = rab.data(1);
    
    auto rabz = rab.data(2);
    
    // set up pointer to starting positions
    
    auto spos = getStartPositions();
    
    for (int32_t i = 0; i < _nScreenedRedContrPairs; i++)
    {
        auto ppidx = spos[i];
        
        rabx[i] = abx[ppidx];
        
        raby[i] = aby[ppidx];
        
        rabz[i] = abz[ppidx];
    }
    
    return rab;
}

void
CGtoPairsBlock::getDistancesAB(      CMemBlock2D<double>& abDistances,
                               const int32_t              nContrPairs,
                               const int32_t              nReplicas) const
{
    // set up pointers to R(AB) distances in primitives data
    
    auto abx = getDistancesABX();
    
    auto aby = getDistancesABY();
    
    auto abz = getDistancesABZ();
    
    // set up pointers to R(AB) distances
    
    auto rabx = abDistances.data(0);
    
    auto raby = abDistances.data(1);
    
    auto rabz = abDistances.data(2);
    
    // set up pointer to starting positions in primitives space
    
    auto spos = getStartPositions();
    
    // set up pointer to contraction pattern
    
    auto cspos = getContrStartPositions();
    
    auto cepos = getContrEndPositions();
    
    for (int32_t i = 0; i < nContrPairs; i++)
    {
        auto ppidx = spos[i];
        
        auto dabx = abx[ppidx];
        
        auto daby = aby[ppidx];
        
        auto dabz = abz[ppidx];
        
        for (int32_t j = cspos[i]; j < cepos[i]; j++)
        {
            rabx[j] = dabx;
            
            raby[j] = daby;
            
            rabz[j] = dabz;
        }
    }
    
    // get number of contracted GTOs pairs
    
    auto ncdim = getNumberOfScreenedContrPairs(nContrPairs - 1);
    
    // generate replicas
    
    for (int32_t i = 1; i < nReplicas; i++)
    {
        auto vx = &rabx[i * ncdim];
        
        mathfunc::copy(vx, 0, rabx, 0, ncdim);
        
        auto vy = &raby[i * ncdim];
        
        mathfunc::copy(vy, 0, raby, 0, ncdim);
        
        auto vz = &rabz[i * ncdim];
        
        mathfunc::copy(vz, 0, rabz, 0, ncdim);
    }
}

const double*
CGtoPairsBlock::getEffectiveCoordinatesPX() const
{
    return _contrFactors.data(0);
}

const double*
CGtoPairsBlock::getEffectiveCoordinatesPY() const
{
    return _contrFactors.data(1);
}

const double*
CGtoPairsBlock::getEffectiveCoordinatesPZ() const
{
    return _contrFactors.data(2);
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
CGtoPairsBlock::getMaxContractionDepth() const
{
    // set up pointers to positions data
    
    auto spos = getStartPositions();
    
    auto epos = getEndPositions();
    
    // loop over contracted GTOs pairs
    
    int32_t mpdim = 0;
    
    for (int32_t i = 0; i < _contrPattern.size(0); i++)
    {
        auto cpdim = epos[i] - spos[i];
        
        if (cpdim > mpdim) mpdim = cpdim;
    }
    
    return mpdim;
}

int32_t
CGtoPairsBlock::getNumberOfPrimPairs(const int32_t iContrPair) const
{
    auto epos = getEndPositions();
    
    return epos[iContrPair]; 
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

int32_t
CGtoPairsBlock::getNumberOfScreenedContrPairs(const int32_t iContrPair) const
{
    auto epos = getContrEndPositions();
 
    return epos[iContrPair]; 
}

int32_t
CGtoPairsBlock::getNumberOfOriginalRedContrPairs() const
{
    return _nOriginalRedContrPairs;
}

int32_t
CGtoPairsBlock::getNumberOfScreenedRedContrPairs() const
{
    return _nScreenedRedContrPairs;
}

int32_t
CGtoPairsBlock::getBraMatrixPosition(const int32_t iComponent) const
{
    if (_nScreenedContrPairs > 0)
    {
        auto bidx = getBraIdentifiers(iComponent);
    
        return bidx[0];
    }
    
    return -1;
}

int32_t
CGtoPairsBlock::getKetMatrixPosition(const int32_t iComponent) const
{
    auto kidx = getKetIdentifiers(iComponent);
    
    if (_nScreenedContrPairs > 0)
    {
        auto spos = kidx[0];
    
        for (int32_t i = 1; i < _nScreenedContrPairs; i++)
        {
            auto cpos = kidx[i];
            
            if (spos > cpos) spos = cpos;
        }
        
        return spos;
    }
    
    return -1;
}

int32_t
CGtoPairsBlock::getNumberOfRowsInBraMatrix() const
{
    if (_nScreenedContrPairs > 0)
    {
        auto bidx = getBraIdentifiers(0);
    
        return bidx[_nScreenedContrPairs - 1] - bidx[0] + 1;
    }
    
    return 0;
}

int32_t
CGtoPairsBlock::getNumberOfRowsInKetMatrix() const
{
    auto kidx = getKetIdentifiers(0);
    
    if (_nScreenedContrPairs > 0)
    {
        auto spos = kidx[0];
        
        auto epos = kidx[0];
        
        for (int32_t i = 0; i < _nScreenedContrPairs; i++)
        {
            auto cpos = kidx[i];
            
            if (spos > cpos) spos = cpos;
            
            if (epos < cpos) epos = cpos;
        }
        
        return (epos - spos) + 1;
    }
    
    return 0;
}

int32_t
CGtoPairsBlock::getMaxNumberContrPairs() const
{
    auto spos = getContrStartPositions();
    
    auto epos = getContrEndPositions();
    
    int32_t mdim = 0;
    
    for (int32_t i = 0; i < _nScreenedRedContrPairs; i++)
    {
        auto cdim = epos[i] - spos[i];
        
        if (cdim > mdim) mdim = cdim;
    }
    
    return mdim;
}

std::string
CGtoPairsBlock::getPairType() const
{
    std::string str("(");
    
    str.append(fstr::to_AngularMomentum(_braAngularMomentum));
    
    str.append(",");
    
    str.append(fstr::to_AngularMomentum(_ketAngularMomentum));
    
    str.append(")");
    
    return str; 
}

std::string
CGtoPairsBlock::getRawSizeString() const
{
    std::string str("Contr.: ");
                    
    str.append(fstr::to_string(_nOriginalContrPairs, 8, fmt::left));
    
    str.append(" Prim.: ");
    
    str.append(std::to_string(_nOriginalPrimPairs));
    
    return fstr::format(str, 32, fmt::left);
}

std::string
CGtoPairsBlock::getScreenedSizeString() const
{
    std::string str("Contr.: ");
    
    str.append(fstr::to_string(_nScreenedContrPairs, 8, fmt::left));
    
    str.append(" Prim.: ");
    
    str.append(std::to_string(_nScreenedPrimPairs));
    
    return fstr::format(str, 32, fmt::left);
}

int32_t
CGtoPairsBlock::_getBlockDimensions() const
{
    auto angab = _braAngularMomentum + _ketAngularMomentum ;

    int32_t ndim = 300;  
     
    if (angab >  4) ndim = 50;
    
    if (angab == 4) ndim = 100;
    
    if (angab == 3) ndim = 150;
    
    if (angab == 2) ndim = 200;
    
    if (angab == 1) ndim = 250;

    return ndim;
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
    
    output << "_indexPattern: " << source._indexPattern << std::endl;
    
    output << "_contrFactors: " << source._contrFactors << std::endl;
    
    output << "_pairFactors: " << source._pairFactors << std::endl;
    
    output << "_normFactors: " << source._normFactors << std::endl;
    
    output << "_threshold: " << source._threshold << std::endl;
    
    output << "_nOriginalPrimPairs: " << source._nOriginalPrimPairs << std::endl;
    
    output << "_nScreenedPrimPairs: " << source._nScreenedPrimPairs << std::endl;
    
    output << "_nOriginalRedContrPairs: " << source._nOriginalRedContrPairs << std::endl;
    
    output << "_nScreenedRedContrPairs: " << source._nScreenedRedContrPairs << std::endl;
    
    output << "_nOriginalContrPairs: " << source._nOriginalContrPairs << std::endl;
    
    output << "_nScreenedContrPairs: " << source._nScreenedContrPairs << std::endl;
    
    return output;
}

