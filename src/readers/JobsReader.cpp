//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "JobsReader.hpp"

#include "JobType.hpp"

CJobsReader::CJobsReader()

    : _state(true)

    , _runMode(execmode::cpu)
{

}

CJobsReader::~CJobsReader()
{

}

bool
CJobsReader::getState() const
{
    return _state;
}

execmode
CJobsReader::getRunMode() const
{
    return _runMode;
}

void
CJobsReader::parse(      std::vector<int32_t>& listOfJobIds,
                   const CInputData&           inputData,
                         COutputStream&        oStream)
{
    oStream << fmt::info << "Parsing mandatory @jobs group..." << fmt::end;

    auto ngroups = inputData.getNumberOfControlGroups("jobs");

    if (ngroups == 1)
    {
        auto cgroup = inputData.getControlGroup(0, "jobs");

        auto ncommands = cgroup.getNumberOfCommands();

        for (int32_t i = 0; i < ncommands; i++)
        {
            auto iline = cgroup.getCommand(i);

            if (_addExecutionMode(iline, oStream))
            {
                continue;
            }
            
            if (_addSinglePoint(listOfJobIds, iline, oStream))
            {
                continue;
            }

            if (_addOptimization(listOfJobIds, iline, oStream))
            {
                continue;
            }
            
            // TODO: add other types of jobs

            if (!_state) return;

            _errorUnknownJobType(iline, oStream);

            return;
        }

        oStream << fmt::info << "...done." << fmt::end << fmt::blank;
    }
    else
    {
        _errorUniqueGroup(ngroups, oStream);
    }
}

void
CJobsReader::_errorUniqueGroup(const size_t         nGroups,
                                     COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror;

    if (nGroups > 1)
    {
        oStream << "Multiple definitions of @jobs group are encountered!";

        oStream << fmt::end;

        oStream << "Please combine @jobs groups into single @jobs group!";

        oStream << fmt::end;
    }
    else
    {
        oStream << "Unable to find a mandatory @jobs group!" << fmt::end;

        oStream << "Please define calculation type(s) in @jobs group!";

        oStream << fmt::end;
    }

    oStream << fmt::blank;
}

void
CJobsReader::_errorUnknownJobType(const CInputLine&    inputLine,
                                        COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported job type requested!" << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}
                
bool
CJobsReader::_addExecutionMode(const CInputLine&    inputLine,
                                     COutputStream& oStream)
{
    if (inputLine.isKeyword(0, "RunMode:"))
    {
        if (inputLine.getNumberOfKeywords() == 2)
        {
            if (inputLine.isKeyword(1, "CPU/GPU"))
            {
                _runMode = execmode::cpu_gpu;

                return true;
            }
            
            if (inputLine.isKeyword(1, "CPU"))
            {
                _runMode = execmode::cpu; 
                
                return true;
            }

            // TODO: Other types of single point calculations

            _errorUnknownCalculationType("RunMode", inputLine, oStream);
        }
        else
        {
            _syntaxRunMode(inputLine, oStream);
        }
    }

    return false;
}

bool
CJobsReader::_addSinglePoint(      std::vector<int32_t>& listOfJobIds,
                             const CInputLine&           inputLine,
                                   COutputStream&        oStream)
{
    if (inputLine.isKeyword(0, "SinglePoint:"))
    {
        if (inputLine.getNumberOfKeywords() == 2)
        {
            if (inputLine.isKeyword(1, "Energy"))
            {
                listOfJobIds.push_back(to_int(job::sp_energy));

                return true;
            }

            // TODO: Other types of single point calculations

            _errorUnknownCalculationType("Single Point", inputLine, oStream);
        }
        else
        {
            _syntaxSinglePoint(inputLine, oStream);
        }
    }

    return false;
}

bool
CJobsReader::_addOptimization(      std::vector<int32_t>& listOfJobIds,
                              const CInputLine&           inputLine,
                                    COutputStream&        oStream)
{
    if (inputLine.isKeyword(0, "Optimization:"))
    {
        if (inputLine.getNumberOfKeywords() == 2)
        {
            if (inputLine.isKeyword(1, "Geometry"))
            {
                listOfJobIds.push_back(to_int(job::opt_geometry));

                return true;
            }

            // TODO: Other types of optimization calculations

            _errorUnknownCalculationType("Optimization", inputLine, oStream);
        }
        else
        {
            _syntaxOptimization(inputLine, oStream);
        }
    }

    return false;
}

void
CJobsReader::_syntaxRunMode(const CInputLine&    inputLine,
                                  COutputStream& oStream)
{
    _state = false;
    
    oStream << fmt::error << "A Run mode is selected as: ";
    
    oStream << fmt::end;
    
    oStream << fmt::error << "RunMode:: [execution mode] " << fmt::end;
    
    oStream << "Please correct input line: " << inputLine.getOriginalString();
    
    oStream << fmt::end << fmt::blank;
}

void
CJobsReader::_syntaxSinglePoint(const CInputLine&    inputLine,
                                      COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "A Single Point calculation is defined as: ";

    oStream << fmt::end;

    oStream << fmt::error << "SinglePoint: [calculation type] " << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

void
CJobsReader::_errorUnknownCalculationType(const char*          calcType,
                                          const CInputLine&    inputLine,
                                                COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Unsupported type of " << calcType;

    oStream <<  " is requested!" << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

void
CJobsReader::_syntaxOptimization(const CInputLine&    inputLine,
                                       COutputStream& oStream)
{
    _state = false;

    oStream << fmt::error << "An Optimization calculation is defined as: ";

    oStream << fmt::end;

    oStream << fmt::error << "Optimization: [calculation type] " << fmt::end;

    oStream << "Please correct input line: " << inputLine.getOriginalString();

    oStream << fmt::end << fmt::blank;
}

