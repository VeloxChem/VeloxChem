//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "PropertyPlasmon.hpp"

#include "MpiFunc.hpp"
#include "MolXYZReader.hpp"
#include "CMMParamReader.hpp"
#include "StringFormat.hpp"
#include "ChemicalElement.hpp"

CPropertyPlasmon::CPropertyPlasmon(const int32_t  globRank,
                                   const int32_t  globNodes,
                                   const execmode runMode)

    : CBaseJob(globRank, globNodes, runMode)
{

}

void
CPropertyPlasmon::set(const std::string&   pathToBasisSets, 
                      const std::string&   pathToForceFields,
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
    
    // read force field data from force fields library
    
    if (_globRank == mpi::master())
    {
        CCMMParamReader rdrparams;
        
        rdrparams.parse(inputData, oStream);
        
        _state = rdrparams.getState();
        
        if (_state)
        {
            _cmmParameters = rdrparams.getCMMParameters(pathToForceFields,
                                                        _molecule, oStream);
        }
        
        _state = rdrparams.getState();
    }
    
    mpi::bcast(_state, _globRank, MPI_COMM_WORLD);
    
    if (!_state) return;
}

void
CPropertyPlasmon::run(COutputStream& oStream,
                      MPI_Comm       comm)
{
    if (!_state) return;
    
    if (_globRank == mpi::master()) _printSetup(oStream);
    
    
}

void
CPropertyPlasmon::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header;

    oStream << "=====================================";

    oStream << fmt::end;

    oStream << "= Classical CMM Plasmon Calculation =";

    oStream << fmt::end;

    oStream << "=====================================";

    oStream << fmt::end << fmt::blank;
}

void
CPropertyPlasmon::_printSetup(COutputStream& oStream) const
{
    std::string str("CMM Force Field Setup");
    
    oStream << fmt::header;
    
    oStream << str << fmt::end;
    
    oStream << std::string(str.size() + 2, '=') << fmt::end << fmt::blank;
    
    str.assign("Model          : ");
    
    str.append(to_string(_cmmParameters[0].getForceFieldModel()));
    
    oStream << fstr::format(str, 46, fmt::left) << fmt::end;
    
    str.assign("Composition    : ");
    
    auto elmlist = _molecule.getElementalComposition();
    
    for (auto i = elmlist.cbegin(); i != elmlist.cend(); ++i)
    {
        CChemicalElement ce;
        
        ce.setAtomType(*i);
        
        str.append(ce.getName());
        
        str.append(" ");
    }
    
    oStream << fstr::format(str, 46, fmt::left) << fmt::end;
    
    str.assign("Atoms          : ");
    
    str.append(std::to_string(_molecule.getNumberOfAtoms()));
    
    oStream << fstr::format(str, 46, fmt::left) << fmt::end;
    
    str.assign("Polarizability : ");
    
    str.append("Static");
    
    oStream << fstr::format(str, 46, fmt::left) << fmt::end << fmt::blank;
    
    if (_molecule.getNumberOfAtoms() < 2000)
    {
        _molecule.printGeometry(oStream);
    }
}

