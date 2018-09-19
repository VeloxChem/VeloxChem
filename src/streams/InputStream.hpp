//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef InputStream_hpp
#define InputStream_hpp

#include <string>

#include "OutputStream.hpp"
#include "InputData.hpp"

/**
 Class CInputStream implements input stream, which provides text input handling
 functionality.
 
 @author Z. Rinkevicius
 */
class CInputStream
{
    /**
     The state of input stream: true  - no errors, false - otherwise.
     */
    bool _state;
    
    /**
     The name of input file.
     */
    std::string _iFilename;

    /**
     Prints file opening error message to output stream and sets input stream
     object state to abnormal.

     @param oStream the output stream.
     */
    void _errorFileOpen(COutputStream& oStream);

    /**
     Prints control group parsing error mesage to output stream and sets input
     stream object state to abnormal.

     @param inputLine the input line object, which can not be parsed.
     @param oStream the output stream.
     */
    void _errorControlGroup(const CInputLine&    inputLine,
                                  COutputStream& oStream);

    /**
     Prints input parsing start message to output stream.

     @param oStream the output stream.
     */
    void _startMessage(COutputStream& oStream) const;

    /**
     Prints input parsing finish message alongside with input parsing
     statistics to output stream.

     @param inpuData the parsed input stored in input data.
     @param nEmptyGroups the number of empty control group objects.
     @param oStream the output stream.
     */
    void _finishMessage(const CInputData&    inpuData,
                        const size_t         nEmptyGroups,
                              COutputStream& oStream) const;

public:

    /**
     Creates an input stream object from name of output file. Errors
     are printed to output stream.
     
     @param iFilename the name of input file.
     @param oStream the output stream.
     */
    CInputStream(const std::string&   iFilename,
                       COutputStream& oStream);

    /**
     Destroys an input stream object.
     */
    ~CInputStream();

    /**
     Reads input file and stores parsed input into input data object. Errors
     are printed to output stream.

     @param inputData the input data object.
     @param oStream the output stream.
     */
    void read(CInputData&    inputData,
              COutputStream& oStream);

    /**
     Get a state of input stream object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;
};

#endif /* InputStream_hpp */
