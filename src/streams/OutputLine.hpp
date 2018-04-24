//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OutputLine_hpp
#define OutputLine_hpp

#include <string>
#include <fstream>
#include <ostream>

/**
 Class COutputLine stores information about output line and provides set of
 methods for formatted writing of output linet to various streams.
 
 @author Z. Rinkevicius
 */
class COutputLine
{
    /**
      The non-redundant string in output line.
     */
    std::string _line; //

    /**
     The offset of non-redundant string from start of formatted output line.
     */
    size_t _offset;

    /**
     The first symbol in output line.
     */
    char _leftSymbol;

    /**
     The last symbol in output line.
     */
    char _rightSymbol;

    /**
     The background fill symbol of formatted output line.
     */
    char _fillSymbol;

    /**
     The total width of formatted output line.
     */
    size_t _width;

public:

    /**
     Creates an empty output line object.
     */
    COutputLine();

    /**
     Creates an output line object from string and formatting data.
     
     @param line the string.
     @param offset the offset of string in formatted output line.
     @param leftSymol the starting symbol of formatted output line.
     @param rightSymbol the end symbol of formatted output line.
     @param fillSymbol the background symbol in formatted output line.
     @param width the width of formatted output line.
     */
    COutputLine(const std::string& line,
                const size_t       offset,
                const char         leftSymol,
                const char         rightSymbol,
                const char         fillSymbol,
                const size_t       width);

    /**
     Creates an output line object by copying other output line object.
     
     @param source the output line object.
     */
    COutputLine(const COutputLine& source);

    /**
     Creates an output line object by by moving other output line object.
     
     @param source the output line object.
     */
    COutputLine(COutputLine&& source) noexcept;

    /**
     Destroys an output line object.
     */
    ~COutputLine();

    /**
     Assigns an output line object by copying other output line object.
     
     @param source the output line object.
     */
    COutputLine& operator=(const COutputLine& source);

    /**
     Assigns an output line object by moving other output line object.
     
     @param source the output line object.
     */
    COutputLine& operator=(COutputLine&& source) noexcept;

    /**
     Compares output line object with other input line object.
     
     @param other the output line object.
     @return true if output line objects are equal, false otherwise.
     */
    bool operator==(const COutputLine& other) const;

    /**
     Compares output line object with other output line object.
     
     @param other the output line object.
     @return true if output line objects are not equal, false otherwise.
     */
    bool operator!=(const COutputLine& other) const;

    /**
     Converts output line object to formatted text line and insert it into
     output file stream.
     
     @param output the output file stream.
     @param source the output line object.
     */
    friend std::ofstream& operator<<(      std::ofstream& output,
                                     const COutputLine&   source);

    /**
     Converts output line object to formatted text line and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the output line object.
     */
    friend std::ostream& operator<<(      std::ostream& output,
                                    const COutputLine&  source);
};

#endif /* OutputLine_hpp */
