//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef CommandLineReader_hpp
#define CommandLineReader_hpp

#include <string>

/**
 Class CCommandLineReader handles command line parameters reading and parsing.
 
 @author Z. Rinkevicius
 */
class CCommandLineReader
{
    /**
     The state of command line reader: true - no errors, false - otherwise.
     */
    bool _state;

    /**
     The name of input file.
     */
    std::string _iFilename;

    /**
     The name of output file.
     */
    std::string _oFilename;
    
    /**
     Prints a description of command line input syntax to standard output
     stream.
     */
    void _printSyntaxHelp() const;

public:

    /**
     Creates a command line reader object.
     
     @param argc the number of command line arguments.
     @param argv the array of command line arguments.
     */
    CCommandLineReader(int    argc,
                       char** argv);

    /**
     Destroys a command line reader object.
     */
    ~CCommandLineReader();

    /**
     Gets a name of input file parsed from command line input.

     @return the name of input file.
     */
    std::string getInputFilename() const;

    /**
     Gets a name of output file parsed from command line input.

     @return the name of output file.
     */
    std::string getOutputFilename() const;

    /**
     Determines if command line parameters obey input syntax requirements.

     @return true if command line parameters satisfy syntax requirements, false
            otherwise.
     */
    bool checkParameters() const;
};

#endif /* CommandLineReader_hpp */
