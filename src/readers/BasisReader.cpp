//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "BasisReader.hpp"

#include <fstream>

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"

CBasisReader::CBasisReader()

    : _state(true)
{

}

CBasisReader::~CBasisReader()
{

}

bool
CBasisReader::getState() const
{
    return _state;
}

void
CBasisReader::parse(const CInputData&    inputData,
                          COutputStream& oStream)
{
    auto ngroups = inputData.getNumberOfControlGroups("basis");

    if (ngroups > 0)
    {
        oStream << fmt::info << "Parsing @basis group..." << fmt::end;

        if (ngroups == 1)
        {
            auto cgroup = inputData.getControlGroup(0, "basis");

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
                    _syntaxBasisSet(oStream);
                }

                if (!_state) return;

                oStream << fmt::info << "...done." << fmt::end << fmt::blank;
            }
            else
            {
                _syntaxBasisSet(oStream);
            }
        }
        else
        {
            _errorUniqueGroup(oStream);
        }
    }
    else
    {
        oStream << fmt::info << "Setting default AO basis set to def2-SVP.";

        oStream << fmt::end;

        _label = std::string("DEF2-SVP");
    }
}

CMolecularBasis
CBasisReader::getAOBasis(const std::string&   pathToBasisSets,
                         const CMolecule&     molecule,
                               COutputStream& oStream)
{
    CMolecularBasis molbasis;

    auto fname = pathToBasisSets;

    fname.append(_label);

    auto elmlist = molecule.getElementalComposition();

    for (auto i = elmlist.cbegin(); i != elmlist.cend(); ++i)
    {
        CAtomBasis atmbasis = _readAtomBasis(*i, fname, oStream);

        if (_state)
        {
            molbasis.addAtomBasis(atmbasis);
        }
        else
        {
            break;
        }
    }

    molbasis.setLabel(_label);

    return molbasis;
}

CAtomBasis
CBasisReader::_readAtomBasis(const int32_t        idElemental,
                             const std::string&   fileName,
                                   COutputStream& oStream)
{
    CAtomBasis atmbasis;

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

            // additional check if it's correct basis set

            if (iline.isControlKeyword("basis_set"))
            {
                if (iline.isKeyword(1, _label))
                {
                    continue;
                }
                else
                {
                    _errorCorruptedBasisSet(oStream);
                }
            }

            // find and read atom basis set

            if (iline.isControlKeyword("atombasis"))
            {
                if (!iline.isKeyword(1, lbl)) continue;

                fb = true;

                continue;
            }

            // scroll over atomic basis sets of other atoms

            if (!fb) continue;

            // stop reading atom basis set

            if (iline.isControlKeyword("end")) break;

            // read header of single shell

            auto bpar = _readShellHeader(iline, oStream);

            // read primitive basis functions

            if (_state)
            {
                CBasisFunction bf;

                for (int32_t i = 0; i < std::get<1>(bpar); i++)
                {
                    std::getline(istream, str);

                    iline = CInputLine(str);

                    auto pfun = _readPrimitveBasisFuction(iline, oStream);

                    if (_state)
                    {
                        bf.add(std::get<0>(pfun), std::get<1>(pfun));
                    }
                    else
                    {
                        break;
                    }
                }

                if (_state)
                {
                    bf.setAngularMomentum(std::get<0>(bpar));

                    bf.normalize();

                    atmbasis.addBasisFunction(bf);
                }
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        _errorBasisSet(oStream);
    }

    istream.close();

    atmbasis.setIdElemental(idElemental);

    return atmbasis;
}

std::tuple<int32_t, int32_t>
CBasisReader::_readShellHeader(const CInputLine&    inputLine,
                                     COutputStream& oStream)
{
    if (inputLine.getNumberOfKeywords() != 3)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1);
    }

    if (inputLine.getIntegerNumber(2) != 1)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1);
    }

    auto nbfuncs = inputLine.getIntegerNumber(1);

    if (nbfuncs < 1)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1);
    }

    auto mang = fstr::to_AngularMomentum(inputLine.getKeyword(0));

    if (mang < 0)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1);
    }

    return std::make_tuple(mang, nbfuncs);
}

std::tuple<double, double>
CBasisReader::_readPrimitveBasisFuction(const CInputLine&    inputLine,
                                              COutputStream& oStream)
{
    if (inputLine.getNumberOfKeywords() != 2)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(0.0, 0.0);
    }

    auto rexp = inputLine.getRealNumber(0);

    if (rexp <= 1.0e-6)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(0.0, 0.0);
    }

    auto rcoef = inputLine.getRealNumber(1);

    return  std::make_tuple(rexp, rcoef);
}

void
CBasisReader::_errorUniqueGroup(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror;

    oStream << "Multiple definitions of @basis group are encountered!";

    oStream << fmt::end;

    oStream << "Please use a single @basis group!" << fmt::end;

    oStream << fmt::blank;
}

void
CBasisReader::_errorBasisSet(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported basis set ";

    oStream << _label << " is requested!" << fmt::end;

    oStream << "Please consult manual for list of supported basis sets and ";

    oStream << "adjust @basis input group!" << fmt::end << fmt::blank;
}

void
CBasisReader::_errorCorruptedBasisSet(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Data file of requested basis set ";

    oStream << _label << " is corrupted in basis sets library!" << fmt::end;

    oStream << "Please consult manual and update basis sets library!";

    oStream << fmt::end << fmt::blank;
}

void
CBasisReader::_syntaxBasisSet(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "A @basis group must contain a single lines";

    oStream << fmt::end <<  "[Basis Set Label]" << fmt::end;

    oStream << "...." << fmt::end << fmt::blank;
}
