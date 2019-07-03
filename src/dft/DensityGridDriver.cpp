//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGridDriver.hpp"

#include <cmath>

#include "GenFunc.hpp"
#include "GtoRecFunc.hpp"
#include "OMPTasks.hpp"
#include "AngularMomentum.hpp"

CDensityGridDriver::CDensityGridDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);

    _thresholdOfDensity = 1.0e-15;

    _thresholdOfPrimGTOs = 1.0e-15;

    _runMode = execmode::cpu;
}

CDensityGridDriver::~CDensityGridDriver()
{
}

CDensityGrid
CDensityGridDriver::generate(const CAODensityMatrix& density,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CMolecularGrid&   molGrid,
                             const xcfun             xcFunctional)
{
    CDensityGrid dgrid(molGrid.getNumberOfGridPoints(), density.getNumberOfDensityMatrices(), xcFunctional, dengrid::ab); 
    
    // execution mode: CPU

    if (_runMode == execmode::cpu)
    {
        _genDensityGridOnCPU(density, molecule, basis, molGrid, xcFunctional);
    }

    // execution mode: CPU/GPU

    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }

    // printf("Density Grid: rank %i points: %i\n",
    //       _locRank, molGrid.getNumberOfGridPoints());
    
    return dgrid;
}

void
CDensityGridDriver::_genDensityGridOnCPU(const CAODensityMatrix& density,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CMolecularGrid&   molGrid,
                                         const xcfun             xcFunctional)
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // set up molecular grid data

    auto mgx = molGrid.getCoordinatesX();

    auto mgy = molGrid.getCoordinatesY();

    auto mgz = molGrid.getCoordinatesZ();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
    
    // set up pointer to density matrix
    
    auto denptr = &density;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfDensityGridPoints(denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize, xcFunctional);
                }
            }
        }
    }

    // delete GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genBatchOfDensityGridPoints(const CAODensityMatrix* aoDensityMatrix,
                                                 const CGtoContainer*    gtoContainer,
                                                 const double*           gridCoordinatesX,
                                                 const double*           gridCoordinatesY,
                                                 const double*           gridCoordinatesZ,
                                                 const int32_t           gridOffset,
                                                 const int32_t           nGridPoints,
                                                 const xcfun             xcFunctional)
{
    // local copy of GTOs containers

    auto gtovec = CGtoContainer(*gtoContainer);
    
    // set up data for density matrices
    
    auto dmattyp = aoDensityMatrix->getDensityType();
    
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        for (int32_t j = i; j < gtovec.getNumberOfGtoBlocks(); j++)
        {
            auto kgtos = gtovec.getGtoBlock(j); 
        }
    }
    
    // set up primitive recursion buffers
//
//    auto nvcomp = _getNumberOfXCComponents(xcFunctional);
//
//    auto pbuffers = gtovec.getPrimAngBuffer(nvcomp);
//
//    // set up contracted Cartesian GTOs buffers
//
//    auto cartbuffers = gtovec.getCartesianBuffer(nvcomp);
//
//    // set up contracted Spherical GTOs buffers
//
//    auto spherbuffer = gtovec.getSphericalBuffer(nvcomp);
//
//    // set up spherical momentum vectors
//
//    auto smomvec = gtovec.getSphericalMomentumVector();
//
//    // loop over batch of grid points
//
//    for (int32_t i = 0; i < nGridPoints; i++)
//    {
//        // grid point coordinates
//
//        auto gx = gridCoordinatesX[gridOffset + i];
//
//        auto gy = gridCoordinatesY[gridOffset + i];
//
//        auto gz = gridCoordinatesZ[gridOffset + i];
//
//        // compute screening factors
//
//        _compScreeningFactors(scrdata, gtovec, gx, gy, gz);
//
//        // update screened GTOs container
//
//        cmpvec.compress(gtovec, redidx, scrdata, _thresholdOfPrimGTOs);
//
//        // loop over GTOs blocks in GTOs container
//
//        for (int32_t j = 0; j < gtovec.getNumberOfGtoBlocks(); j++)
//        {
//            // compute distances
//
//            _compDistances(rdist, cmpvec, redidx, j, gx, gy, gz);
//
//            // compute primitive GTOs values at grid point
//
//            _compPrimGtoValues(pbuffers[j], rdist, cmpvec, redidx, j, xcFunctional);
//
//            // contract Cartesian GTOs values
//
//            _contrPrimGtoValues(cartbuffers[j], pbuffers[j], cmpvec, redidx, j);
//
//            // transform to spherical GTOs values
//
//            _transContrGtoValues(spherbuffer[j], cartbuffers[j], smomvec[j], redidx, j, xcFunctional);
//        }
//
//        // compute density values for grid point
//
//        _compDensityValues(spherbuffer, cmpvec, xcFunctional);
//    }
}

void
CDensityGridDriver::_compScreeningFactors(CVecMemBlock<double>& screenFactors,
                                          const CGtoContainer&  gtoContainer,
                                          const double          gridCoordinateX,
                                          const double          gridCoordinateY,
                                          const double          gridCoordinateZ)
{
    for (int32_t i = 0; i < gtoContainer.getNumberOfGtoBlocks(); i++)
    {
        // set up primitives GTOs data

        auto bang = gtoContainer.getAngularMomentum(i);

        auto pexps = gtoContainer.getExponents(i);

        auto pfacts = gtoContainer.getNormFactors(i);

        // set up primitive GTOs coordinates

        auto coordsx = gtoContainer.getCoordinatesX(i);

        auto coordsy = gtoContainer.getCoordinatesY(i);

        auto coordsz = gtoContainer.getCoordinatesZ(i);

        // set up screening factors

        auto sfacts = screenFactors[i].data();

        // loop over GTOs in GTOs block

        for (int32_t j = 0; j < gtoContainer.getNumberOfPrimGtos(i); j++)
        {
            // compute distances

            auto rx = gridCoordinateX - coordsx[j];

            auto ry = gridCoordinateY - coordsy[j];

            auto rz = gridCoordinateZ - coordsz[j];

            // determine max distance

            auto rmax = (rx > ry) ? rx : ry;

            if (rz > rmax) rmax = rz;

            // compute zero order overlap

            sfacts[j] = pfacts[j] * std::exp(-pexps[j] * (rx * rx + ry * ry + rz * rz));

            // compose screening factor

            sfacts[j] *= _getScaleFactor(bang, rmax);

            // absolute value of screening factor

            sfacts[j] = std::fabs(sfacts[j]);
        }
    }
}

double
CDensityGridDriver::_getScaleFactor(const int32_t angularMomentum, const double maxRadius) const
{
    if (angularMomentum == 0) return 1.0;

    if (angularMomentum == 1) return maxRadius;

    auto r2max = maxRadius * maxRadius;

    if (angularMomentum == 2) return r2max;

    if (angularMomentum == 3) return r2max * maxRadius;

    auto r4max = r2max * r2max;

    if (angularMomentum == 4) return r4max;

    if (angularMomentum == 5) return r4max * maxRadius;

    if (angularMomentum == 6) return r4max * r2max;

    // TO DO: implement higher order if needed

    return 0.0;
}

void
CDensityGridDriver::_compDistances(CMemBlock2D<double>&        distances,
                                   const CGtoContainer&        gtoContainer,
                                   const CMemBlock2D<int32_t>& redDimensions,
                                   const int32_t               iGtoBlock,
                                   const double                gridCoordinateX,
                                   const double                gridCoordinateY,
                                   const double                gridCoordinateZ) const
{
    // set up pointer to reduced dimensions

    auto reddim = redDimensions.data(0);

    // compute distances

    mathfunc::distances(distances.data(0),
                        distances.data(1),
                        distances.data(2),
                        gridCoordinateX,
                        gridCoordinateY,
                        gridCoordinateZ,
                        gtoContainer.getCoordinatesX(iGtoBlock),
                        gtoContainer.getCoordinatesY(iGtoBlock),
                        gtoContainer.getCoordinatesZ(iGtoBlock),
                        reddim[iGtoBlock]);
}

int32_t
CDensityGridDriver::_getNumberOfXCComponents(const xcfun xcFunctional) const
{
    if (xcFunctional == xcfun::lda) return 1;

    if (xcFunctional == xcfun::gga) return 4;

    if (xcFunctional == xcfun::mgga) return 5;

    return 0;
}

void
CDensityGridDriver::_compPrimGtoValues(CMemBlock2D<double>&        gtoValues,
                                       const CMemBlock2D<double>&  distances,
                                       const CGtoContainer&        gtoContainer,
                                       const CMemBlock2D<int32_t>& redDimensions,
                                       const int32_t               iGtoBlock,
                                       const xcfun                 xcFunctional) const
{
    // local density approximation

    if (xcFunctional == xcfun::lda)
    {
        _compPrimGtoValuesForLDA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

        return;
    }

    // generalized gradient approximation

    if (xcFunctional == xcfun::gga)
    {
        _compPrimGtoValuesForGGA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

        return;
    }

    // meta generalized gradient approximation

    if (xcFunctional == xcfun::mgga)
    {
        _compPrimGtoValuesForMGGA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

        return;
    }
}

void
CDensityGridDriver::_compPrimGtoValuesForLDA(CMemBlock2D<double>&        gtoValues,
                                             const CMemBlock2D<double>&  distances,
                                             const CGtoContainer&        gtoContainer,
                                             const CMemBlock2D<int32_t>& redDimensions,
                                             const int32_t               iGtoBlock) const
{
    auto mang = gtoContainer.getAngularMomentum(iGtoBlock);

    // s-type functions

    gtorec::compGtoTypeSForLDA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

    if (mang == 0) return;

    // p-type functions

    gtorec::compGtoTypePForLDA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 1) return;

    // d-type functions

    gtorec::compGtoTypeDForLDA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 2) return;

    // f-type functions

    gtorec::compGtoTypeFForLDA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 3) return;

    // g-type functions

    gtorec::compGtoTypeGForLDA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 4) return;

    // TODO: Implement higher angular momentum (l > 4)
}

void
CDensityGridDriver::_compPrimGtoValuesForGGA(CMemBlock2D<double>&        gtoValues,
                                             const CMemBlock2D<double>&  distances,
                                             const CGtoContainer&        gtoContainer,
                                             const CMemBlock2D<int32_t>& redDimensions,
                                             const int32_t               iGtoBlock) const
{
    auto mang = gtoContainer.getAngularMomentum(iGtoBlock);

    // s-type functions

    gtorec::compGtoTypeSForGGA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

    if (mang == 0) return;

    // p-type functions

    gtorec::compGtoTypePForGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 1) return;

    // d-type functions

    gtorec::compGtoTypeDForGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 2) return;

    // f-type functions

    gtorec::compGtoTypeFForGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 3) return;

    // g-type functions

    gtorec::compGtoTypeGForGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 4) return;

    // TODO: Implement higher angular momentum (l > 4)
}

void
CDensityGridDriver::_compPrimGtoValuesForMGGA(CMemBlock2D<double>&        gtoValues,
                                              const CMemBlock2D<double>&  distances,
                                              const CGtoContainer&        gtoContainer,
                                              const CMemBlock2D<int32_t>& redDimensions,
                                              const int32_t               iGtoBlock) const
{
    auto mang = gtoContainer.getAngularMomentum(iGtoBlock);

    // s-type functions

    gtorec::compGtoTypeSForMGGA(gtoValues, distances, gtoContainer, redDimensions, iGtoBlock);

    if (mang == 0) return;

    // p-type functions

    gtorec::compGtoTypePForMGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 1) return;

    // d-type functions

    gtorec::compGtoTypeDForMGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 2) return;

    // f-type functions

    gtorec::compGtoTypeFForMGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 3) return;

    // g-type functions

    gtorec::compGtoTypeGForMGGA(gtoValues, distances, redDimensions, iGtoBlock);

    if (mang == 4) return;

    // TODO: Implement higher angular momentum (l > 4)
}

void
CDensityGridDriver::_contrPrimGtoValues(CMemBlock2D<double>&        cartGtoValues,
                                        const CMemBlock2D<double>&  primGtoValues,
                                        const CGtoContainer&        gtoContainer,
                                        const CMemBlock2D<int32_t>& redDimensions,
                                        const int32_t               iGtoBlock) const
{
    // contraction pattern (start, end positions)

    auto spos = gtoContainer.getStartPositions(iGtoBlock);

    auto epos = gtoContainer.getEndPositions(iGtoBlock);

    // set up number of contracted GTOs

    auto reddim = redDimensions.data(1);

    auto ngto = reddim[iGtoBlock];

    // set up contracted vectors dimensions, indexes

    auto nvec = cartGtoValues.blocks();

    auto pidx = primGtoValues.blocks() - nvec;

    // contract GTOs

    genfunc::contract(cartGtoValues, primGtoValues, 0, pidx, spos, epos, ngto, nvec);
}

void
CDensityGridDriver::_transContrGtoValues(CMemBlock2D<double>&        spherGtoValues,
                                         const CMemBlock2D<double>&  cartGtoValues,
                                         const CSphericalMomentum&   spherMomentum,
                                         const CMemBlock2D<int32_t>& redDimensions,
                                         const int32_t               iGtoBlock,
                                         const xcfun                 xcFunctional) const
{
    // set up number of contracted GTOs

    auto reddim = redDimensions.data(1);

    auto ngto = reddim[iGtoBlock];

    // set up number of functional contributions

    auto nvcomp = _getNumberOfXCComponents(xcFunctional);

    // transform GTOs

    genfunc::transform(spherGtoValues, cartGtoValues, spherMomentum, 0, 0, ngto, nvcomp);
}

void
CDensityGridDriver::_compDensityValues(const CVecMemBlock2D<double>& spherGtoValues,
                                       const CGtoContainer&          gtoContainer,
                                       const xcfun                   xcFunctional) const
{
    auto ngblk = gtoContainer.getNumberOfGtoBlocks();
    
    for (int32_t i = 0; i < ngblk; i++)
    {
        auto bcomp = angmom::to_SphericalComponents(gtoContainer.getAngularMomentum(i));
        
        auto bdim = gtoContainer.getNumberOfContrGtos(i);
        
        for (int32_t j = 0; j < ngblk; j++)
        {
            auto kcomp = angmom::to_SphericalComponents(gtoContainer.getAngularMomentum(j));
            
            auto kdim = gtoContainer.getNumberOfContrGtos(j);
            
            for (int32_t k = 0; k < bcomp; k++)
            {
                auto bidx = gtoContainer.getIdentifiers(i, k);
                
                for (int32_t l = 0; l < kcomp; l++)
                {
                    auto kidx = gtoContainer.getIdentifiers(j, l);
                    
                    
                }
            }
        }
    }
}
