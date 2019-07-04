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
    // initialize density grid
    
    CDensityGrid dgrid(molGrid.getNumberOfGridPoints(), density.getNumberOfDensityMatrices(), xcFunctional, dengrid::ab);
    
    dgrid.zero(); 
    
    // execution mode: CPU

    if (_runMode == execmode::cpu)
    {
        _genDensityGridOnCPU(dgrid, density, molecule, basis, molGrid, xcFunctional);
    }

    // execution mode: CPU/GPU

    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }
    
    return dgrid;
}

void
CDensityGridDriver::_genDensityGridOnCPU(      CDensityGrid&     denGrid, 
                                         const CAODensityMatrix& density,
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
    
    // set up poinet to density grid
    
    auto dgridptr = &denGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, denptr, dgridptr)
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
                    _genBatchOfDensityGridPoints(dgridptr, denptr, gtovec, mgx, mgy, mgz, tbposition, tbsize, xcFunctional);
                }
            }
        }
    }

    // delete GTOs container

    delete gtovec;
}

void
CDensityGridDriver::_genBatchOfDensityGridPoints(      CDensityGrid*     densityGrid, 
                                                 const CAODensityMatrix* aoDensityMatrix,
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
    
    // loop over GTOs container data
    
    for (int32_t i = 0; i < gtovec.getNumberOfGtoBlocks(); i++)
    {
        auto bgtos = gtovec.getGtoBlock(i);
        
        for (int32_t j = i; j < gtovec.getNumberOfGtoBlocks(); j++)
        {
            _compDensityForGtoBlocks(densityGrid, aoDensityMatrix, bgtos, gtovec.getGtoBlock(j),
                                     gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ,
                                     gridOffset, nGridPoints, xcFunctional);
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
CDensityGridDriver::_compDensityForGtoBlocks(      CDensityGrid*     densityGrid,
                                             const CAODensityMatrix* aoDensityMatrix,
                                             const CGtoBlock&        braGtoBlock,
                                             const CGtoBlock&        ketGtoBlock,
                                             const double*           gridCoordinatesX,
                                             const double*           gridCoordinatesY,
                                             const double*           gridCoordinatesZ,
                                             const int32_t           gridOffset,
                                             const int32_t           nGridPoints,
                                             const xcfun             xcFunctional)
{
    // determine symmetry of bra and ket sides
    
    auto symbk = (braGtoBlock == ketGtoBlock);
    
    // set up max density matrix elements for (i,j) pair of density indexes
    
    auto ndmat = aoDensityMatrix->getNumberOfDensityMatrices();
    
    // angular momentum data for bra and ket
    
    auto bang = braGtoBlock.getAngularMomentum();
    
    auto kang = ketGtoBlock.getAngularMomentum();
    
    // set up Cartesian GTOs buffers
    
    auto nvcomp = _getNumberOfXCComponents(xcFunctional);
    
    auto bncart = angmom::to_CartesianComponents(bang);
    
    auto kncart = angmom::to_CartesianComponents(kang);
    
    CMemBlock2D<double> bcartbuff(nGridPoints, nvcomp * bncart);
    
    CMemBlock2D<double> kcartbuff(nGridPoints, nvcomp * kncart);
    
    // set up spherical GTOs buffers
    
    auto bnspher = angmom::to_SphericalComponents(bang);
    
    auto knspher = angmom::to_SphericalComponents(kang);
    
    CMemBlock2D<double> bspherbuff(nGridPoints, nvcomp * bnspher);
    
    CMemBlock2D<double> kspherbuff(nGridPoints, nvcomp * knspher);
    
    CMemBlock<double> maxdq(ndmat);
    
    for (int32_t i = 0; i < braGtoBlock.getNumberOfContrGtos(); i++)
    {
        _compGtoValuesOnGrid(bspherbuff, bcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                             braGtoBlock, i, xcFunctional);
        
        auto jstart = symbk ? i : 0;
        
        for (int32_t j = jstart; j < ketGtoBlock.getNumberOfContrGtos(); j++)
        {
            _compGtoValuesOnGrid(kspherbuff, kcartbuff, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset,
                                 ketGtoBlock, j, xcFunctional);
        }
    }
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
CDensityGridDriver::_compGtoValuesOnGrid(      CMemBlock2D<double>& cartGtoGridBuffer,
                                               CMemBlock2D<double>& spherGtoGridBuffer,
                                         const double*              gridCoordinatesX,
                                         const double*              gridCoordinatesY,
                                         const double*              gridCoordinatesZ,
                                         const int32_t              gridOffset,
                                         const CGtoBlock&           gtoBlock,
                                         const int32_t              iContrGto,
                                         const xcfun                xcFunctional) const
{
    
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
