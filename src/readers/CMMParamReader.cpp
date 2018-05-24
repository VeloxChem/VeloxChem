//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "CMMParamReader.hpp"

#include <fstream>

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"

CCMMParamReader::CCMMParamReader()

    : _state(true)
{

}

CCMMParamReader::~CCMMParamReader()
{

}

bool
CCMMParamReader::getState() const
{
    return _state;
}

void
CCMMParamReader::parse(const CInputData&    inputData,
                             COutputStream& oStream)
{
    auto ngroups = inputData.getNumberOfControlGroups("cmmparam");

    if (ngroups > 0)
    {
        oStream << fmt::info << "Parsing @cmmparam group..." << fmt::end;

        if (ngroups == 1)
        {
            auto cgroup = inputData.getControlGroup(0, "cmmparam");

            auto ncommands = cgroup.getNumberOfCommands();

            if (ncommands == 1)
            {
                auto iline = cgroup.getCommand(0);

                if (iline.getNumberOfKeywords() == 1)
                {
                    _label = iline.getUpcasedKeyword(0);
                }
                else
                {
                    _syntaxForceField(oStream);
                }

                if (!_state) return;

                oStream << fmt::info << "...done." << fmt::end << fmt::blank;
            }
            else
            {
                _syntaxForceField(oStream);
            }
        }
        else
        {
            _errorUniqueGroup(oStream);
        }
    }
    else
    {
        oStream << fmt::info << "Setting default CMM force field parameters.";

        oStream << fmt::end;

        _label = std::string("CMM-FF-JENSEN");
    }
}

std::vector<CCMMParameters>
CCMMParamReader::getCMMParameters(const std::string&   pathToForceFields,
                                  const CMolecule&     molecule,
                                        COutputStream& oStream)
{
    std::vector<CCMMParameters> ffparams;

    auto fname = pathToForceFields;

    fname.append(_label);
    
    auto elmlist = molecule.getElementalComposition();

    for (auto i = elmlist.cbegin(); i != elmlist.cend(); ++i)
    {
        auto atmparam = _readAtomForceField(*i, fname, oStream);

        if (_state)
        {
            ffparams.push_back(atmparam);
        }
        else
        {
            break;
        }
    }

    return ffparams;
}

CCMMParameters
CCMMParamReader::_readAtomForceField(const int32_t        idElemental,
                                     const std::string&   fileName,
                                          COutputStream& oStream)
{
    CCMMParameters atmparam;

    std::ifstream istream;

    istream.open(fileName.c_str(), std::ifstream::in);

    if (istream.good())
    {
        bool fb = false;

        CChemicalElement ce;

        ce.setAtomType(idElemental);

        auto lbl = ce.getName();

        while (!istream.eof())
        {
            std::string str;

            std::getline(istream, str);

            CInputLine iline(str);

            if (iline.isEmpty()) continue;

            // additional check if it's correct force field

            if (iline.isControlKeyword("force_field"))
            {
                if (iline.isKeyword(1, _label))
                {
                    continue;
                }
                else
                {
                    _errorCorruptedForceField(oStream);
                }
            }

            // find and read atom basis set

            if (iline.isControlKeyword("atomff"))
            {
                if (!iline.isKeyword(1, lbl)) continue;

                fb = true;

                continue;
            }

            // scroll over force fields of other atoms

            if (!fb) continue;

            // stop reading atom force field

            if (iline.isControlKeyword("end")) break;

            // read header of single shell

            _readForceFieldHeader(atmparam, iline, oStream);
            
            if (_state)
            {
                // read number of Lorentzians
            
                std::getline(istream, str);
        
                iline = CInputLine(str);
            
                auto nfreq = _readNumberOfLorentzians(iline, oStream);
                
                if (!_state) break; 
            
                // add each Lorentzian
            
                for (int32_t i = 0; i < nfreq; i++)
                {
                    std::getline(istream, str);
                
                    iline = CInputLine(str);
                
                    _readLorentzian(atmparam, iline, oStream);
                    
                    if (!_state) break;
                }
                
                if (!_state) break; 
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        _errorForceField(oStream);
    }

    istream.close();

    return atmparam;
}

void
CCMMParamReader::_syntaxForceField(COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::error << "A @cmmparam group must contain a single lines";
    
    oStream << fmt::end <<  "[Force Field Label]" << fmt::end;
    
    oStream << "...." << fmt::end << fmt::blank;
}

void
CCMMParamReader::_errorUniqueGroup(COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::cerror;
    
    oStream << "Multiple definitions of @cmmparam group are encountered!";
    
    oStream << fmt::end;
    
    oStream << "Please use a single @cmmparam group!" << fmt::end;
    
    oStream << fmt::blank;
}

void
CCMMParamReader::_errorForceField(COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::cerror << "Unsupported force field ";
    
    oStream << _label << " is requested!" << fmt::end;
    
    oStream << "Please consult manual for list of supported force fields and ";
    
    oStream << "adjust @cmmparam input group!" << fmt::end << fmt::blank;
}

void
CCMMParamReader::_errorCorruptedForceField(COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::cerror << "Data file of requested force field ";
    
    oStream << _label << " is corrupted in force fields library!" << fmt::end;
    
    oStream << "Please consult manual and update force fields library!";
    
    oStream << fmt::end << fmt::blank;
}

void
CCMMParamReader::_readForceFieldHeader(      CCMMParameters& ffParameters,
                                       const CInputLine&     inputLine,
                                             COutputStream&  oStream)
{
    if ((inputLine.getNumberOfKeywords() != 3) &&
        (inputLine.getNumberOfKeywords() != 5))
    {
        _errorCorruptedForceField(oStream);
        
        return;
    }
    
    if (inputLine.isKeyword(0, "orig"))
    {
        ffParameters.setForceFieldModel(cmmtyp::original);
        
        ffParameters.setPolarizability(inputLine.getRealNumber(1));
        
        ffParameters.setExponent(inputLine.getRealNumber(2));
        
        return;
    }
    
    if (inputLine.isKeyword(0, "ext"))
    {
        ffParameters.setForceFieldModel(cmmtyp::enhanced);
        
        ffParameters.setPolarizability(inputLine.getRealNumber(1));
        
        ffParameters.setExponent(inputLine.getRealNumber(2));
        
        ffParameters.setCoordinationNumber(inputLine.getRealNumber(3));
        
        ffParameters.setCoordinationFactor(inputLine.getRealNumber(4));
        
        return;
    }
}

int32_t
CCMMParamReader::_readNumberOfLorentzians(const CInputLine&    inputLine,
                                                COutputStream& oStream)
{
    if (inputLine.getNumberOfKeywords() == 1)
    {
        return inputLine.getIntegerNumber(0);
    }
    
    _errorCorruptedForceField(oStream);
    
    return 0;
}

void
CCMMParamReader::_readLorentzian(      CCMMParameters& ffParameters,
                                 const CInputLine&     inputLine,
                                       COutputStream&  oStream)
{
    if (inputLine.getNumberOfKeywords() != 3)
    {
        _errorCorruptedForceField(oStream);
        
        return;
    }
    
    ffParameters.addFrequency(inputLine.getRealNumber(0),
                              inputLine.getRealNumber(1),
                              inputLine.getRealNumber(2));
}
