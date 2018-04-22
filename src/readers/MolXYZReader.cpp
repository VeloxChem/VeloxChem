//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MolXYZReader.hpp"

#include "ChemicalElement.hpp"
#include "Codata.hpp"

CMolXYZReader::CMolXYZReader()

    : _state(true)

    , _needUnitConversion(true)

    , _molCharge(0.0)

    , _molMultiplicity(1)
{

}

CMolXYZReader::~CMolXYZReader()
{

}

bool CMolXYZReader::getState() const
{
    return _state;
}

void CMolXYZReader::parse(CMolecule& molecule, const CInputData& inputData,
                          COutputStream& oStream)
{
    auto nGroups = inputData.getNumberOfControlGroups("molxyz");

    if (nGroups > 0)
    {
        if (nGroups == 1)
        {
            parse(molecule, inputData, 0, oStream);
        }
        else
        {
            _errorUniqueGroup(oStream);
        }
    }
}

void CMolXYZReader::parse(CMolecule& molecule, const CInputData& inputData,
                          const int32_t iGroup, COutputStream &oStream)
{
    auto nGroups = inputData.getNumberOfControlGroups("molxyz");
    
    if (nGroups > 0)
    {
       if (iGroup < nGroups)
       {
           oStream << fmt::info << "Parsing @molxyz group..." << fmt::end;
           
           auto controlGroup = inputData.getControlGroup(iGroup, "molxyz");
           
           auto nCommands = controlGroup.getNumberOfCommands();
           
           if (nCommands > 1)
           {
               _setDimensions(nCommands - 1);
               
               for (int32_t i = 0; i < nCommands; i++)
               {
                   if (i == 0)
                   {
                       _parseHeader(controlGroup.getCommand(i), oStream);
                   }
                   else
                   {
                       _parseAtom(controlGroup.getCommand(i), (i - 1), oStream);
                   }
                   
                   if (!_state) return;
               }
               
               // create molecule
               
               molecule = CMolecule(_coordinates, _charges, _masses, _labels,
                                    _idsElemental);
               
               // set up charge and multiplicity
               
               molecule.setCharge(_molCharge);
               
               molecule.setMultiplicity(_molMultiplicity);
               
               _checkMultiplicity(molecule, controlGroup.getCommand(0), oStream);
               
               if (!_state) return;
               
               oStream << fmt::info << "...done." << fmt::end << fmt::blank;
           }
           else
           {
               _syntaxIncompleteGroup(oStream);
           }
       }
       else
       {
           molecule = CMolecule();
       }
    }
}

void CMolXYZReader::_parseHeader(const CInputLine& inputLine,
                                 COutputStream& oStream)
{
    auto nkeys = inputLine.getNumberOfKeywords();

    if ((nkeys == 2) || (nkeys == 3))
    {
        if (inputLine.isIntegerNumber(0) && inputLine.isIntegerNumber(1))
        {
            auto mvalue = inputLine.getIntegerNumber(0);

            _molCharge = static_cast<double>(mvalue);

            mvalue = inputLine.getIntegerNumber(1);

            if (mvalue > 0)
            {
                _molMultiplicity = mvalue;
            }
            else
            {
                _syntaxHeader(inputLine, oStream);
            }

            if (nkeys == 3) _parseUnits(inputLine, oStream);
        }
        else
        {
            _syntaxHeader(inputLine, oStream);
        }
    }
    else
    {
        _syntaxHeader(inputLine, oStream);
    }
}

void CMolXYZReader::_parseAtom(const CInputLine& inputLine, const int32_t iAtom,
                               COutputStream& oStream)
{
    auto nkeys = inputLine.getNumberOfKeywords();

    if ((nkeys == 4) || (nkeys == 5))
    {
        if (inputLine.isRealNumber(1) &&
            inputLine.isRealNumber(2) &&
            inputLine.isRealNumber(3))
        {
            CChemicalElement chemelm;

            if (chemelm.setAtomType(inputLine.getUpcasedKeyword(0)))
            {
                if (nkeys == 5)
                {
                    if (inputLine.isIntegerNumber(4))
                    {
                        if (!chemelm.setIsotope(inputLine.getIntegerNumber(4)))
                        {
                            _errorIsotope(inputLine, oStream);
                        }

                        if (!_state) return;
                    }
                    else
                    {
                        _syntaxAtom(inputLine, oStream);
                    }
                }
                
                // determine total number of atoms
                
                int32_t natoms = static_cast<int32_t>(_charges.size());
                
                // set coordinates of i-th atom
                
                auto rx = inputLine.getRealNumber(1);
                
                auto ry = inputLine.getRealNumber(2);
                
                auto rz = inputLine.getRealNumber(3);
                
                if (_needUnitConversion)
                {
                    auto factor = 1.0 / units::getBohrValueInAngstroms();

                    rx *= factor;
                    
                    ry *= factor;
                    
                    rz *= factor;
                }
                
                _coordinates[iAtom] = rx;
                
                _coordinates[natoms + iAtom] = ry;
                
                _coordinates[2 * natoms + iAtom] = rz;
                
                // set up properties of i-th atom
                
                _charges[iAtom] = chemelm.getAtomicCharge();
                
                _masses[iAtom] = chemelm.getAtomicMass();
                
                _labels[iAtom] = chemelm.getName();
                
                _idsElemental[iAtom] = chemelm.getIdentifier();
            }
            else
            {
                _errorChemicalElement(inputLine, oStream);
            }
        }
        else
        {
            _syntaxAtom(inputLine, oStream);
        }
    }
    else
    {
        _syntaxAtom(inputLine, oStream);
    }
}

void CMolXYZReader::_parseUnits(const CInputLine& inputLine,
                                COutputStream& oStream)
{
    _needUnitConversion = true;

    if (inputLine.isKeyword(2, "Angstroms")) return;

    if (inputLine.isKeyword(2, "Angs")) return;

    if (inputLine.isKeyword(2, "Bohrs"))
    {
        _needUnitConversion = false;

        return;
    }

    if (inputLine.isKeyword(2, "Au"))
    {
        _needUnitConversion = false;

        return;
    }

    _errorUnits(inputLine, oStream);
}

void CMolXYZReader::_errorUniqueGroup(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror;

    oStream << "Multiple definitions of @molxyz group are not allowed ";
    
    oStream << "for this type of calculation!" << fmt::end;

    oStream << "Please combine @molxyz groups into single @molxyz group!";

    oStream << fmt::end << fmt::blank;
}

void CMolXYZReader::_errorMultiplicity(const CInputLine& inputLine,
                                       COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror;

    oStream << "Multiplicity defined in @molxyz group is incompatble with ";

    oStream << "number of electrons in molecule!" << fmt::end;

    oStream << "Please adjust charge or multiplicity in input line: ";

    oStream << inputLine.getOriginalString() << fmt::blank;
}

void CMolXYZReader::_errorUnits(const CInputLine& inputLine,
                                COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported units type ";

    oStream << inputLine.getKeyword(2) << " is requested!" << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

void CMolXYZReader::_errorChemicalElement(const CInputLine& inputLine,
                                          COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported chemical element ";

    oStream << inputLine.getKeyword(0) << " is requested!" << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

void CMolXYZReader::_errorIsotope(const CInputLine& inputLine,
                                  COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported isotope ";

    oStream << inputLine.getKeyword(4) << " of chemical element ";

    oStream << inputLine.getKeyword(0) << " is requested!" << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

void CMolXYZReader::_syntaxIncompleteGroup(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "A @molxyz group must contain at least two lines:";

    oStream << fmt::end;

    oStream << "[Charge] [Multiplicity] (Units)" << fmt::end;

    oStream << "[Label] [Coord X] [Coord Y] [Coord Z] (Isotope)" << fmt::end;

    oStream << "...." << fmt::end << fmt::blank;
}

void CMolXYZReader::_syntaxHeader(const CInputLine& inputLine,
                                  COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "The header of @molxyz group is defined as:";

    oStream << fmt::end;

    oStream << "[Charge] [Multiplicity] (Units)" << fmt::end;

    oStream << "a) Parameter [Charge] must be an integer number" << fmt::end;

    oStream << "b) Parameter [Multiplicity] must be an integer number greater";

    oStream << " than zero" << fmt::end;

    oStream << "c) Optional parameter (Units) can be Angstroms, Angs, Bohrs,";

    oStream << " or Au" << fmt::end;

    oStream << fmt::end << "Please correct input line: ";

    oStream << inputLine.getOriginalString() << fmt::end << fmt::blank;
}

void CMolXYZReader::_syntaxAtom(const CInputLine& inputLine,
                                COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "An atom in @molxyz group is defined as:";

    oStream << fmt::end;

    oStream << "[Label] [Coord X] [Coord Y] [Coord Z] (Isotope)";

    oStream << fmt::end;

    oStream << "a) Parameter [Label] must be label of chemical element";

    oStream << fmt::end;

    oStream << "b) Parameter [Coord X] must be a real number";

    oStream << fmt::end;

    oStream << "c) Parameter [Coord Y] must be a real number";

    oStream << fmt::end;

    oStream << "d) Parameter [Coord Z] must be a real number";

    oStream << fmt::end;

    oStream << "e) Optional parameter (Isotope) can be an integer number, ";

    oStream << "which is isotope label" << fmt::end;

    oStream << fmt::end << "Please correct input line: ";

    oStream << inputLine.getOriginalString() << fmt::end << fmt::blank;
}

void CMolXYZReader::_checkMultiplicity(const CMolecule& molecule,
                                       const CInputLine& inputLine,
                                       COutputStream& oStream)
{
    auto multip = molecule.getMultiplicity() % 2;

    auto nelec = molecule.getNumberOfElectrons() % 2;

    if ((multip == 0) && (nelec != 1)) _state = false;

    if ((multip == 1) && (nelec != 0)) _state = false;

    if (!_state) _errorMultiplicity(inputLine, oStream);
}

void CMolXYZReader::_setDimensions(const int32_t nAtoms)
{
    _coordinates = std::vector<double>(3 * nAtoms, 0.0);
    
    _charges = std::vector<double>(nAtoms, 0.0);
    
    _masses = std::vector<double>(nAtoms, 0.0);
    
    _labels = std::vector<std::string>(nAtoms, std::string(" "));
    
    _idsElemental = std::vector<int32_t>(nAtoms, 0);
}
