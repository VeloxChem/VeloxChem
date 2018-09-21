//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CommandLineReader.hpp"

#include <iostream>

CCommandLineReader::CCommandLineReader(int    argc,
                                       char** argv)

    : _state(false)
{
    if (argc == 3)
    {
        _state = true;

        _iFilename.assign(argv[1]);

        _oFilename.assign(argv[2]);
    }
}

CCommandLineReader::~CCommandLineReader()
{

}

std::string
CCommandLineReader::getInputFilename() const
{
    return _iFilename;
}

std::string
CCommandLineReader::getOutputFilename() const
{
    return _oFilename;
}

bool
CCommandLineReader::checkParameters() const
{
    if (!_state) _printSyntaxHelp();

    return _state;
}

void
CCommandLineReader::_printSyntaxHelp() const
{
    std::cout << "@ Usage: VeloxChemMP.x [input file] [output file]";

    std::cout << std::endl;

    std::cout << "@ [input  file] - name of the input  file." << std::endl;

    std::cout << "@ [output file] - name of the output file." << std::endl;
}
