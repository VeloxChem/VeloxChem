//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "EnvironmentReader.hpp"

CEnvironmentReader::CEnvironmentReader()

    : _state(true)
{

}

CEnvironmentReader::~CEnvironmentReader()
{

}

bool
CEnvironmentReader::getState() const
{
    return _state;
}

void CEnvironmentReader::parse(const CInputData&    inputData,
                                COutputStream& oStream)
{
    oStream << fmt::info << "Parsing mandatory @progenv group..." << fmt::end;

    auto ngroups = inputData.getNumberOfControlGroups("progenv");

    if (ngroups == 1)
    {
        auto cgroup = inputData.getControlGroup(0, "progenv");

        auto ncommands = cgroup.getNumberOfCommands();

        for (int32_t i = 0; i < ncommands; i++)
        {
            auto iline = cgroup.getCommand(i);

            if (_addPathToBasisSets(iline, oStream)) continue;

            // TODO: add other types of jobs

            if (!_state) return;

            _errorUnknownVariable(iline, oStream);

            return;
        }

        oStream << fmt::info << "...done." << fmt::end << fmt::blank;
    }
    else
    {
        _errorUniqueGroup(ngroups, oStream);
    }
}

std::string CEnvironmentReader::getPathToBasisSets() const
{
    return _pathToBasisSets;
}

void CEnvironmentReader::_errorUniqueGroup(const size_t nGroups,
                                           COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror;

    if (nGroups > 1)
    {
        oStream << "Multiple definitions of @progenv group are encountered!";

        oStream << fmt::end;

        oStream << "Please combine @progenv groups into single @progenv group!";

        oStream << fmt::end;
    }
    else
    {
        oStream << "Unable to find a mandatory @progenv group!" << fmt::end;

        oStream << "Please define environmental variables(s) in ";

        oStream << "@progenv group!" << fmt::end;
    }

    oStream << fmt::blank;
}

void CEnvironmentReader::_errorUnknownVariable(const CInputLine& inputLine,
                                               COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "An unknown environmental variable is defined!";

    oStream << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

bool CEnvironmentReader::_addPathToBasisSets(const CInputLine& inputLine,
                                             COutputStream& oStream)
{
    if (inputLine.isKeyword(0, "BasisLibrary:"))
    {
        if (inputLine.getNumberOfKeywords() == 2)
        {
            _pathToBasisSets = inputLine.getKeyword(1);

            if (_pathToBasisSets[_pathToBasisSets.size()] != '/')
            {
                _pathToBasisSets.append("/");
            }

            return true;
        }
        else
        {
            _syntaxBasisLibrary(inputLine, oStream);
        }
    }

    return false;
}

void CEnvironmentReader::_syntaxBasisLibrary(const CInputLine& inputLine,
                                             COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "An environmental variable BasisLibrary must be ";

    oStream << "defined as: " << fmt::end;

    oStream << fmt::error << "BasisLibrary: [Full Path To Basis Sets Library]";

    oStream << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}
