//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OutputStream_hpp
#define OutputStream_hpp

#include <string>
#include <vector>

#include "FmtType.hpp"
#include "OutputLine.hpp"

/**
 Class COutputStream implements output stream, which provides text output
 handling functionality.
 
 @author Z. Rinkevicius
 */
class COutputStream
{
    /**
     The state of output stream: true  - no errors, false - otherwise.
     */
    bool _state;

    /**
     The name of output file.
     */
    std::string _filename;

    /**
     The temporary output line.
     */
    std::string _line;

    /**
     The maixum width of formatted output line.
     */
    size_t _maxWidth;

    /**
     The start symbol of formatted output line.
     */
    char _leftSymbol;

    /**
     The termination symbol of formatted output line.
     */
    char _rightSymbol;

    /**
     The background symbol of formatted output line.
     */
    char _fillSymbol;

    /**
     The formatting key of text.
     */
    fmt _alignment;

    /**
     The prefix added to start of formatted output line.
     */
    std::string _prefix;

    /**
     The buffer of formatted output lines stored in packed form.
     */
    std::vector<COutputLine> _buffer;

    /**
     Prints a file open error message to standard error stream.
     */
    void _errorFileOpen() const;
    
    /**
     Adds a output line to stream's buffer in packed format, and empties
     temporary output line.
     */
    void _addLineToBuffer();

    /**
     Gets a position of temporary output line in formatted output line.

     @return the position in formatted output line.
     */
    size_t _getLinePosition() const;

    /**
     Appends a string to temporary output line.

     @param line the string.
     */
    void _appendToLine(const std::string& line);

public:

    /**
     Creates an output stream object from name of output file.
     
     @param filename the  name of output file.
     */
    COutputStream(const std::string& filename);

    /**
     Destroys an output stream object.
     */
    ~COutputStream();

    /**
     Gets a state of output stream object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Writes buffer of formatted output lines to output file, and empties buffer
     afterwords.
     */
    void flush();
    
    /**
     Inserts formatting key into temporary output line. If formatting key is
     fmt::end, the temporary output line is formatted, packed, and added to
     buffer of packed output lines; output line is emptied afterwards.

     @param output the output stream.
     @param source the formatting key.
     */
    friend COutputStream& operator<<(COutputStream& output, const fmt& source);

    /**
     Inserts string into temporary output line.

     @param output the output stream.
     @param source the string.
     */
    friend COutputStream& operator<<(COutputStream& output,
                                     const std::string& source);

    /**
     Inserts C string into temporary output line.

     @param output the output stream.
     @param source the C string.
     */
    friend COutputStream& operator<<(COutputStream& output, const char* source);
};

#endif /* OutputStream_hpp */
