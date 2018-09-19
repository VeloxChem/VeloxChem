//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "InputStream.hpp"

CInputStream::CInputStream(const std::string&   iFilename,
                                 COutputStream& oStream)

    : _state(true)

    , _iFilename(iFilename)
{
    if (!_iFilename.empty())
    {
        std::ifstream fstream;

        fstream.open(_iFilename.c_str(), std::ifstream::in);

        if (!fstream.good()) _errorFileOpen(oStream);

        fstream.close();
    }
}

CInputStream::~CInputStream()
{

}

void
CInputStream::read(CInputData&    inputData,
                   COutputStream& oStream)
{
    if (!_iFilename.empty())
    {
        std::ifstream fstream;

        fstream.open(_iFilename.c_str(), std::ifstream::in);

        if (fstream.good())
        {
            _startMessage(oStream);

            bool ingroup = false;

            CControlGroup cgroup;

            size_t egroups = 0;

            while (!fstream.eof())
            {
                std::string str;

                std::getline(fstream, str);

                CInputLine iline(str);

                if (!iline.isEmpty())
                {
                    if (iline.isControlLine())
                    {
                        if (iline.isControlKeyword("end"))
                        {
                            if (cgroup.isEmpty())
                            {
                                egroups++;
                            }
                            else
                            {
                                inputData.addControlGroup(cgroup);
                            }

                            cgroup.clear();

                            ingroup = false;

                            continue;
                        }

                        if (ingroup)
                        {
                            _errorControlGroup(iline, oStream);

                            fstream.close();

                            return;
                        }
                        else
                        {
                            cgroup.setHeader(iline);

                            ingroup = true;

                            continue;
                        }
                    }

                    if (ingroup) cgroup.addCommand(iline);
                }
            }

            if (_state) _finishMessage(inputData, egroups, oStream);
        }
        else
        {
            _errorFileOpen(oStream);
        }

        fstream.close();
    }
}

bool
CInputStream::getState() const
{
    return _state;
}

void
CInputStream::_errorFileOpen(COutputStream& oStream)
{
    _state = false;

    oStream << fmt::cerror << "Failed to open input file ";

    oStream << _iFilename << "!" << fmt::end << fmt::blank;
}

void
CInputStream::_errorControlGroup(const CInputLine&    inputLine,
                                       COutputStream& oStream)
{
    _state = false;

    oStream << fmt::blank << fmt::cerror;

    oStream << "Nesting of control groups is not allowed!" << fmt::end;

    oStream << "Please correct input line: ";

    oStream << inputLine.getOriginalString() << fmt::end;
}

void
CInputStream::_startMessage(COutputStream& oStream) const
{
    oStream << fmt::info << "Reading input file ";

    oStream << _iFilename << "..." << fmt::end;
}

void
CInputStream::_finishMessage(const CInputData&    inpuData,
                             const size_t         nEmptyGroups,
                                   COutputStream& oStream) const
{
    auto ngroups = inpuData.getNumberOfControlGroups() + nEmptyGroups;

    if (ngroups > 0)
    {
        oStream << fmt::info << "Found ";

        if (ngroups == 1)
        {
            oStream << "control group." << fmt::end;
        }
        else
        {
            oStream << std::to_string(ngroups) << " control groups.";

            oStream << fmt::end;
        }
    }

    if (nEmptyGroups > 0)
    {
        oStream << "Discarded ";

        if (nEmptyGroups == 1)
        {
            oStream << "empty control group." << fmt::end;
        }
        else
        {
            oStream << std::to_string(nEmptyGroups);

            oStream << " empty control groups." << fmt::end;
        }
    }

    oStream << "...done." << fmt::end << fmt::blank;
}
