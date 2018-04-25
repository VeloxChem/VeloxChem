//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GridDriver.hpp"

#include <array>

#include "MpiFunc.hpp"
#include "StringFormat.hpp"
#include "ChemicalElement.hpp"
#include "Log3Quadrature.hpp"
#include "LebedevLaikovQuadrature.hpp"
#include "MathConst.hpp"
#include "PartitionFunc.hpp"

CGridDriver::CGridDriver(const int32_t  globRank,
                         const int32_t  globNodes,
                         const execmode runMode,
                               MPI_Comm comm)

    : _gridLevel(5)

    , _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)

    , _thresholdOfWeight(1.0e-15)

    , _runMode(runMode)
{
    _locRank  = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CGridDriver::~CGridDriver()
{

}

void
CGridDriver::setLevel(const int32_t  gridLevel,
                            MPI_Comm comm)
{
    if (_globRank == mpi::master())
    {
        if ((gridLevel > 0) && (gridLevel < 7)) _gridLevel = gridLevel;
    }

    if (_isLocalMode)
    {
        // FIX ME: global master to local master transfer.
    }

    mpi::bcast(_gridLevel, comm);
}

CMolecularGrid
CGridDriver::generate(const CMolecule&     molecule,
                            COutputStream& oStream,
                            MPI_Comm       comm) const
{
    CSystemClock time;
    
    if (_globRank == mpi::master()) _startHeader(molecule, oStream);
    
    // initialize molecular grid
    
    CMolecularGrid molgrid;
    
    // execution mode: CPU
    
    if (_runMode == execmode::cpu)
    {
        molgrid = _genGridPointsOnCPU(molecule, comm);
    }
    
    // execution mode: CPU/GPU
    
    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }
    
    if (_isLocalMode)
    {
        // FIX ME: add tranfer from local master to global master
    }
    
    if (_globRank == mpi::master()) _finishHeader(time, molgrid, oStream);
    
    return molgrid;
}

int32_t
CGridDriver::_getNumberOfRadialPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 6> nPoints{{20, 25, 30, 35, 45, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Ne atoms

    if ((idElemental > 2)  && (idElemental < 11))
    {
        const std::array<int32_t, 6> nPoints{{25, 30, 35, 40, 50, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Na-Ar atoms

    if ((idElemental > 10)  && (idElemental < 19))
    {
        const std::array<int32_t, 6> nPoints{{30, 35, 40, 45, 55, 200}};

        return nPoints[_gridLevel - 1];
    }

    // K-Kr atoms

    if ((idElemental > 18)  && (idElemental < 37))
    {
        const std::array<int32_t, 6> nPoints{{35, 40, 45, 50, 60, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36)  && (idElemental < 87))
    {
        const std::array<int32_t, 6> nPoints{{40, 45, 50, 55, 65, 200}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

int32_t
CGridDriver::_getNumberOfAngularPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 6> nPoints{{50, 110, 194, 302, 434, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Kr atoms

    if ((idElemental > 2)  && (idElemental < 37))
    {
        const std::array<int32_t, 6> nPoints{{110, 194, 302, 434, 590, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36)  && (idElemental < 87))
    {
        const std::array<int32_t, 6> nPoints{{194, 302, 434, 590, 770, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

void
CGridDriver::_startHeader(const CMolecule&     molecule,
                                COutputStream& oStream) const
{
    oStream << fmt::header << "Numerical Grid For DFT Integration" << fmt::end;
    
    oStream << std::string(36, '=') << fmt::end << fmt::blank;
    
    std::string str("Grid Level          : ");
    
    str.append(std::to_string(_gridLevel));
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end;
    
    str.assign("Radial Quadrature   : Log3");
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end;
    
    str.assign("Angular Quadrature  : Lebedev-Laikov");
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end;
    
    str.assign("Partitioning Scheme : SSF");
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end << fmt::blank;
    
    oStream << "  Atom ";
    
    oStream << fstr::format(std::string("Angular Points"), 15, fmt::center);
    
    oStream << fstr::format(std::string("Radial Points"), 15, fmt::center);
    
    oStream << fmt::end << fmt::blank;
    
    auto molcomp = molecule.getElementalComposition();
    
    int32_t npoints = 0;
    
    for (auto i = molcomp.cbegin(); i != molcomp.cend(); ++i)
    {
        std::string label("  ");
        
        CChemicalElement elem;
        
        elem.setAtomType(*i);
        
        label.append(elem.getName());
        
        oStream << fstr::format(label, 6, fmt::left);
        
        auto apoints = _getNumberOfAngularPoints(*i);
        
        auto rpoints = _getNumberOfRadialPoints(*i);
        
        oStream << fstr::format(std::to_string(apoints), 15, fmt::center);
        
        oStream << fstr::format(std::to_string(rpoints), 15, fmt::center);
        
        oStream << fmt::end;
        
        npoints += apoints * rpoints * molecule.getNumberOfAtoms(*i);
    }
    
    oStream << fmt::blank;
    
    str.assign("Full Grid       : ");
    
    str.append(std::to_string(npoints));
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end;
}

void
CGridDriver::_finishHeader(const CSystemClock&   time,
                           const CMolecularGrid& molecularGrid,
                                 COutputStream&  oStream) const
{
    std::string str("Pruned Grid     : ");
    
    str.append(std::to_string(molecularGrid.getNumberOfGridPoints()));
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end;
    
    str.assign("Generation Time : ");
    
    str.append(time.getElapsedTime());
    
    oStream << fstr::format(str, 54, fmt::left) << fmt::end << fmt::blank;
}

CMolecularGrid
CGridDriver::_genGridPointsOnCPU(const CMolecule& molecule,
                                       MPI_Comm   comm) const
{
    // molecular data
    
    auto natoms = molecule.getNumberOfAtoms();
    
    auto idselem = molecule.getIdsElemental();
    
    auto molrx = molecule.getCoordinatesX();
    
    auto molry = molecule.getCoordinatesY();
    
    auto molrz = molecule.getCoordinatesZ();
    
    auto mdist = molecule.getMinDistances();
    
    auto molrm = mdist.data();
    
    // determine dimensions of atoms batch for each MPI node
    
    auto nodatm = mpi::batch_size(natoms, _locRank, _locNodes);
    
    auto nodoff = mpi::batch_offset(natoms, _locRank, _locNodes);
    
    // allocate raw grid
    
    int32_t bpoints = _getBatchSize(idselem, nodoff, nodatm);
    
    CMemBlock2D<double>* rawgrid = new CMemBlock2D<double>(bpoints, 4);
    
    // generate atomic grid for each atom in molecule
    
    #pragma omp parallel shared(rawgrid, idselem, molrx, molry, molrz, molrm,\
                                natoms, nodatm, nodoff)
    {
        #pragma omp single nowait
        {
            int32_t gridoff = 0;
            
            for (int32_t i = 0; i < nodatm; i++)
            {
                auto iatom = nodoff + i;
                
                auto ielem = idselem[iatom];
                
                auto minrad = molrm[iatom];
                
                #pragma omp task firstprivate(ielem, iatom, minrad, gridoff)
                {
                    _genAtomGridPoints(rawgrid, minrad, gridoff,
                                       molrx, molry, molrz,
                                       natoms, ielem, iatom);
                }
                
                gridoff += _getNumberOfRadialPoints(ielem)
                
                         * _getNumberOfAngularPoints(ielem);
            }
        }
    }
    
    // screen raw grid points & create prunned grid
    
    bpoints = _screenRawGridPoints(rawgrid);
    
    auto prngrid = rawgrid->slice(0, bpoints);
    
    delete rawgrid;
    
    // create molecular grid on master node 
    
    return CMolecularGrid(prngrid.gather(_locRank, _locNodes, comm));
}

int32_t
CGridDriver::_getBatchSize(const int32_t* idsElemental,
                           const int32_t  offset,
                           const int32_t  nAtoms) const
{
    int32_t npoints = 0;
    
    for (int32_t i = 0; i < nAtoms; i++)
    {
        auto idx = idsElemental[offset + i];
        
        npoints += _getNumberOfRadialPoints(idx)
        
                * _getNumberOfAngularPoints(idx);
    }
    
    return npoints;
}

void
CGridDriver::_genAtomGridPoints(      CMemBlock2D<double>* rawGridPoints,
                                const double               minDistance, 
                                const int32_t              gridOffset,
                                const double*              atomCoordinatesX,
                                const double*              atomCoordinatesY,
                                const double*              atomCoordinatesZ,
                                const int32_t              nAtoms,
                                const int32_t              idElemental,
                                const int32_t              idAtomic) const
{
    // determine number of grid points
    
    auto nrpoints = _getNumberOfRadialPoints(idElemental);
    
    auto napoints = _getNumberOfAngularPoints(idElemental);
    
    // generate radial grid points
    
    CLog3Quadrature rquad(nrpoints, idElemental);
    
    auto rpoints = rquad.generate();
    
    auto rrx = rpoints.data(0);
    
    auto rrw = rpoints.data(1);
    
    auto rrf = 4.0 * mathconst::getPiValue();
    
    // generate angular grid points
    
    CLebedevLaikovQuadrature aquad(napoints);
    
    auto apoints = aquad.generate();
    
    auto rax = apoints.data(0);
    
    auto ray = apoints.data(1);
    
    auto raz = apoints.data(2);
    
    auto raw = apoints.data(3);
    
    // atom coordinates
    
    auto atmx = atomCoordinatesX[idAtomic];
    
    auto atmy = atomCoordinatesY[idAtomic];
    
    auto atmz = atomCoordinatesZ[idAtomic];
    
    // assemble atom grid points from quadratures
    
    for (int32_t i = 0; i < nrpoints; i++)
    {
        // radial quadrature point
        
        auto crx = rrx[i];
        
        auto crw = rrf * rrw[i];
        
        // set up pointers to raw grid
        
        auto gx = rawGridPoints->data(0, gridOffset + i * napoints);
        
        auto gy = rawGridPoints->data(1, gridOffset + i * napoints);
        
        auto gz = rawGridPoints->data(2, gridOffset + i * napoints);
        
        auto gw = rawGridPoints->data(3, gridOffset + i * napoints);
        
        // contract with angular quadrature
        
        #pragma omp simd
        for (int32_t j = 0; j < napoints; j++)
        {
            gx[j] = crx * rax[j] + atmx;
            
            gy[j] = crx * ray[j] + atmy;
            
            gz[j] = crx * raz[j] + atmz;
            
            gw[j] = crw * raw[j];
        }
    }
    
    // apply partitioning function
    
    partfunc::ssf(rawGridPoints, minDistance, gridOffset, nrpoints * napoints,
                  atomCoordinatesX, atomCoordinatesY, atomCoordinatesZ,
                  nAtoms, idAtomic);
}

int32_t
CGridDriver::_screenRawGridPoints(CMemBlock2D<double>* rawGridPoints) const
{
    // set up pointers to grid points data
    
    auto gx = rawGridPoints->data(0);
    
    auto gy = rawGridPoints->data(1);
    
    auto gz = rawGridPoints->data(2);
    
    auto gw = rawGridPoints->data(3);
    
    // loop over grid points
    
    int32_t cpoints = 0;

    for (int32_t i = 0; i < rawGridPoints->size(0); i++)
    {
        if (gw[i] > _thresholdOfWeight)
        {
            gx[cpoints] = gx[i];
            
            gy[cpoints] = gy[i];
            
            gz[cpoints] = gz[i];
            
            gw[cpoints] = gw[i];
            
            cpoints++;
        }
    }
    
    return cpoints;
}

