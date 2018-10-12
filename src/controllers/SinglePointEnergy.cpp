//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "SinglePointEnergy.hpp"

#include "MpiFunc.hpp"
#include "MolXYZReader.hpp"
#include "BasisReader.hpp"
//#include "AtomicDensityReader.hpp"

#include "GridDriver.hpp"
#include "DensityGridDriver.hpp"
#include "XCFuncType.hpp"

#include "OverlapIntegralsDriver.hpp"
#include "KineticEnergyIntegralsDriver.hpp"
#include "NuclearPotentialIntegralsDriver.hpp"
#include "ElectronicPotentialIntegralsDriver.hpp"
#include "ThreeCenterElectronRepulsionIntegralsDriver.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"
#include "SADGuessDriver.hpp"
#include "AOFockMatrix.hpp"
#include "MemBlock2D.hpp"

CSinglePointEnergy::CSinglePointEnergy(const int32_t  globRank,
                                       const int32_t  globNodes,
                                       const execmode runMode)

    : CBaseJob(globRank, globNodes, runMode)
{

}

void
CSinglePointEnergy::set(const std::string&   pathToBasisSets,
                        const CInputData&    inputData,
                              COutputStream& oStream)
{
    if (_globRank == mpi::master()) _startHeader(oStream);
    
    // read molecular geometry
    
    if (_globRank == mpi::master())
    {
        CMolXYZReader rdrmolxyz;

        rdrmolxyz.parse(_molecule, inputData, oStream);

        _state = rdrmolxyz.getState();
    }

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
    
    if (!_state) return;

    // broadcast molecular geometry

    _molecule.broadcast(_globRank, MPI_COMM_WORLD);

    // print molecular geometry

    if (_globRank == mpi::master())
    {
        _molecule.printGeometry(oStream);

        _state = _molecule.checkProximity(0.1, oStream);
    }

    oStream.flush();
    
    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

    if (!_state) return;
    
    // read AO basis from basis set library

    if (_globRank == mpi::master())
    {
        CBasisReader rdraobasis;

        rdraobasis.parse(inputData, oStream);

         _state = rdraobasis.getState();

        // read AO basis
        
        if (_state)
        {
            _aoBasis = rdraobasis.getAOBasis(pathToBasisSets, _molecule,
                                             oStream);
        }

        _state = rdraobasis.getState();
        
        // read RI-J basis
        
        //if (_state)
        //{
        //    _riBasis = rdraobasis.getRIJBasis(pathToBasisSets, _molecule,
        //                                      oStream);
        // }

        //_state = rdraobasis.getState();
        
        // read minimal AO basis
        
        if (_state)
        {
            _minBasis = rdraobasis.getMinBasis(pathToBasisSets, _molecule,
                                               oStream);
        }
        
        _state = rdraobasis.getState();
    }

    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);

    if (!_state) return;

    // broadcast AO basis

    _aoBasis.broadcast(_globRank, MPI_COMM_WORLD);
    
    // broadcast RI basis
    
    //_riBasis.broadcast(_globRank, MPI_COMM_WORLD);
    
    // print atomic orbitals i.e. AO basis

    if (_globRank == mpi::master()) _aoBasis.printBasis("Atomic Orbitals",
                                                        _molecule, oStream);
    
    // print RI basis i.e. Coulomb fitting
    
    //if (_globRank == mpi::master()) _riBasis.printBasis("Coulomb Fitting Orbitals",
    //                                                    _molecule, oStream);

//    if (_globRank == mpi::master())
//    {
//        // TODO: move to proper initial guess object
//
//        //CAtomicDensityReader adrdr;
//
//        //_density = adrdr.getAtomicDensities(pathToBasisSets, _molecule, _aoBasis, oStream);
//
//        //_state = adrdr.getState();
//    }
//
//    mpi::bcast_bool(_state, _globRank, MPI_COMM_WORLD);
//
//    if (!_state) return;

    // TODO: add other keywords...
}

void
CSinglePointEnergy::run(COutputStream& oStream,
                        MPI_Comm       comm)
{
    // generate molecular grid

    //CGridDriver drvgrid(_globRank, _globNodes, _runMode, comm);

    //auto molgrid = drvgrid.generate(_molecule, oStream, comm);

    //molgrid.distribute(_globRank, _globNodes, comm);
    
    // generate density grid
    
    //CDensityGridDriver drvDenGrid(_globRank, _globNodes, _runMode, comm);
    
    //drvDenGrid.generate(_molecule, _aoBasis, molgrid, xcfun::mgga, oStream, comm);
    
    // compute overlap integrals
    
    COverlapIntegralsDriver ovldrv(_globRank, _globNodes, comm);
    
    auto ovlmat = ovldrv.compute(_molecule, _aoBasis, oStream, comm);

    // compute kinetic energy integrals
    
    CKineticEnergyIntegralsDriver kindrv(_globRank, _globNodes, comm);
    
    auto kinmat = kindrv.compute(_molecule, _aoBasis, oStream, comm);
    
    // compute nuclear potential integrals
    
    CNuclearPotentialIntegralsDriver npotdrv(_globRank, _globNodes, comm);
    
    auto npotmat = npotdrv.compute(_molecule, _aoBasis, oStream, comm);

    // compute electronic potential integrals
    
    //CElectronicPotentialIntegralsDriver epotdrv(_globRank, _globNodes, comm);
    
    //auto epotmat = epotdrv.compute(_molecule, _aoBasis, oStream, comm);
    
    // compute three center electron repulsion integrals
    
    // CThreeCenterElectronRepulsionIntegralsDriver ridrv(_globRank, _globNodes, comm);
    
    // ridrv.compute(_molecule, _aoBasis, _riBasis, 1.0e-13, oStream, comm);
    
    // compute SAD initial guess
    
    CSADGuessDriver saddrv(_globRank, _globNodes, comm);

    auto ovlmat_min_ao = ovldrv.compute(_molecule, _minBasis, _aoBasis, oStream, comm);
    
    auto dsad = saddrv.compute(_molecule, _minBasis, _aoBasis,
                               ovlmat_min_ao, ovlmat, oStream, comm);

    // compute electron repulsion integrals
    
    CElectronRepulsionIntegralsDriver eridrv(_globRank, _globNodes, comm);
    
    //CGtoBlock sgtos(_molecule, _aoBasis, 0);
    
    //CGtoBlock pgtos(_molecule, _aoBasis, 1);
    
    //CGtoBlock bgtos(_molecule, _aoBasis, 4);
    
    //CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    //CGtoPairsBlock kpairs(sgtos, pgtos, 1.0e-13);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-12, _molecule, _aoBasis,
                                 oStream, comm);

    CAOFockMatrix fock(dsad);
    
    eridrv.compute(fock, dsad, _molecule, _aoBasis, qqdata, oStream, comm);
    
    //fock.addCoreHamiltonian(kinmat, npotmat, 0);
    
    //std::cout << kinmat.getString();
    
    //std::cout << npotmat.getString();
    
    std::cout << dsad.getString();
    
    std::cout << fock.getString();
}

void
CSinglePointEnergy::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header;

    oStream << "==============================================";

    oStream << fmt::end;

    oStream << "=   Single Point Energy For Single Molecule  =";

    oStream << fmt::end;

    oStream << "==============================================" ;

    oStream << fmt::end  << fmt::blank;
}
