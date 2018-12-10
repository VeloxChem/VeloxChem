//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "BasisReader.hpp"

#include <fstream>

#include "ChemicalElement.hpp"
#include "StringFormat.hpp"
#include "ErrorHandler.hpp"

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
CBasisReader::setLabel(const std::string& label)
{
    _label = fstr::upcase(label);
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

        oStream << fmt::end << fmt::blank;

        _label = std::string("DEF2-SVP");
    }
}

CMolecularBasis
CBasisReader::getAOBasis(const std::string&   pathToBasisSets,
                         const CMolecule&     molecule,
                               COutputStream& oStream)
{
    return _getBasis(_label, pathToBasisSets, molecule, oStream);
}

CMolecularBasis
CBasisReader::getRIJBasis(const std::string&   pathToBasisSets,
                          const CMolecule&     molecule,
                                COutputStream& oStream)
{
    auto rilbl = _getLabelForRIJBasis(oStream);
    
    if (!_state) return CMolecularBasis();
    
    return _getBasis(rilbl, pathToBasisSets, molecule, oStream);
}

CMolecularBasis
CBasisReader::getMinBasis(const std::string&   pathToBasisSets,
                          const CMolecule&     molecule,
                                COutputStream& oStream)
{
    auto minlbl = std::string("MIN-CC-PVDZ");
    
    return _getBasis(minlbl, pathToBasisSets, molecule, oStream);
}

CAtomBasis
CBasisReader::_readAtomBasis(const int32_t        idElemental,
                             const std::string&   basisLabel,
                             const std::string&   fileName,
                                   COutputStream& oStream)
{
    CAtomBasis atmbasis;

    std::ifstream istream;

    istream.open(fileName.c_str(), std::ifstream::in);

    std::string errbas("readAtomBasis: Cannot find basis set ");

    errors::assertMsgCritical(istream.good(), errbas + fileName);

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
                if (iline.isKeyword(1, basisLabel))
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
                auto npgto = std::get<1>(bpar);
                
                auto ncgto = std::get<2>(bpar);
                
                std::vector<double> pexps(npgto, 0.0);
                
                std::vector<double> pcoefs(npgto * ncgto, 0.0);

                for (int32_t i = 0; i < npgto; i++)
                {
                    std::getline(istream, str);

                    iline = CInputLine(str);

                    _readPrimitveBasisFuction(iline, pexps, pcoefs, i, npgto,
                                              ncgto, oStream);

                    if (!_state) break;
                }

                if (_state)
                {
                    CBasisFunction bf(pexps, pcoefs, ncgto, std::get<0>(bpar));

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

std::tuple<int32_t, int32_t, int32_t>
CBasisReader::_readShellHeader(const CInputLine&    inputLine,
                                     COutputStream& oStream)
{
    if (inputLine.getNumberOfKeywords() != 3)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1, -1);
    }

    auto mang = fstr::to_AngularMomentum(inputLine.getKeyword(0));
    
    if (mang < 0)
    {
        _errorCorruptedBasisSet(oStream);
        
        return std::make_tuple(-1, -1, -1);
    }
    
    auto nbfuncs = inputLine.getIntegerNumber(1);

    if (nbfuncs < 1)
    {
        _errorCorruptedBasisSet(oStream);

        return std::make_tuple(-1, -1, -1);
    }

    auto ncfuncs = inputLine.getIntegerNumber(2);
    
    if (ncfuncs < 1)
    {
        _errorCorruptedBasisSet(oStream);
        
        return std::make_tuple(-1, -1, -1);
    }
    
    return std::make_tuple(mang, nbfuncs, ncfuncs);
}

void
CBasisReader::_readPrimitveBasisFuction(const CInputLine&          inputLine,
                                              std::vector<double>& exponents,
                                              std::vector<double>& normFactors,
                                        const int32_t              iExponent,
                                        const int32_t              nExponents,
                                        const int32_t              nContrVectors,
                                              COutputStream&       oStream)
{
    if (inputLine.getNumberOfKeywords() != (nContrVectors + 1))
    {
        _errorCorruptedBasisSet(oStream);
        
        exponents.clear();
        
        normFactors.clear();

        return;
    }
    
    // read exponent

    auto rexp = inputLine.getRealNumber(0);

    if (rexp <= 1.0e-6)
    {
        _errorCorruptedBasisSet(oStream);
        
        exponents.clear();
        
        normFactors.clear();

        return;
    }

    exponents[iExponent] = rexp;
    
    // read normalization factors
    
    for (int32_t i = 0; i < nContrVectors; i++)
    {
        auto rcoef = inputLine.getRealNumber(i + 1);
                
        normFactors[i * nExponents + iExponent] = rcoef;
    }
    
    return;
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

CMolecularBasis
CBasisReader::_getBasis(const std::string&   basisLabel,
                        const std::string&   pathToBasisSets,
                        const CMolecule&     molecule,
                              COutputStream& oStream)
{
    CMolecularBasis molbasis;
    
    auto fname = pathToBasisSets;

    fname.append("/");
    
    fname.append(basisLabel);
    
    auto elmlist = molecule.getElementalComposition();
    
    for (auto i = elmlist.cbegin(); i != elmlist.cend(); ++i)
    {
        CAtomBasis atmbasis = _readAtomBasis(*i, basisLabel, fname, oStream);
        
        if (_state)
        {
            molbasis.addAtomBasis(atmbasis);
        }
        else
        {
            break;
        }
    }
    
    molbasis.setLabel(basisLabel);
    
    return molbasis;
}

std::string
CBasisReader::_getLabelForRIJBasis(COutputStream& oStream)
{
    // DEF2 basis sets family
    
    std::vector<std::string> d2lbls({"DEF2-SV(P)", "DEF2-SVP", "DEF2-TZVP",
                                     "DEF2-TZVPP", "DEF2-QZVP", "DEF2-QZVPP"});
    
    for (size_t i = 0; i < d2lbls.size(); i++)
    {
        if (d2lbls[i] == _label) return std::string("DEF2-RI-J");
    }
    
    // ADD handling of other basis sets here...
    
    _errorRIJBasisSet(oStream);
    
    return std::string();
}

void
CBasisReader::_errorRIJBasisSet(COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::cerror << "Coulomb fitting basis is not supported for  ";
    
    oStream << _label << " basis set!" << fmt::end;
    
    oStream << "Please consult manual for list of supported basis sets and ";
    
    oStream << "adjust @basis input group!" << fmt::end << fmt::blank;
}

