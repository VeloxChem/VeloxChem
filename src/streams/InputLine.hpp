//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef InputLine_hpp
#define InputLine_hpp

#include <string>
#include <ostream>

/**
 Class CInputLine stores information about input line and provides set of
 methods for extracting data from input line.
 
 @author Z. Rinkevicius
 */
class CInputLine
{
    /**
     The original input line.
     */
    std::string _originalString; 

    /**
     The parsed input line.
     */
    std::string _parsedString;

    /**
     Converts original input line to parsed input line by removing spaces at the
     begining of original input line and part of original line behind any
     comment symbol.
     */
    void _trimParsedString();

public:

    /**
     Creates an empty input line object.
     */
    CInputLine();

    /**
     Creates an input line object from string.
     
     @param str the string.
     */
    CInputLine(const std::string& str);

    /**
     Creates an input line object by copying other input line object.
     
     @param source the input line object.
     */
    CInputLine(const CInputLine& source);

    /**
     Creates an input line object by by moving other input line object.
     
     @param source the input line object.
     */
    CInputLine(CInputLine&& source) noexcept;

    /**
     Destroys an input line object.
     */
    ~CInputLine();

    /**
     Assigns an input line object by copying other input line object.

     @param source the input line object.
     */
    CInputLine& operator=(const CInputLine& source);

    /**
     Assigns an input line object by moving other input line object.
     
     @param source the input line object.
     */
    CInputLine& operator=(CInputLine&& source) noexcept;

    /**
     Compares input line object with other input line object.

     @param other the input line object.
     @return true if input line objects are equal, false otherwise.
     */
    bool operator==(const CInputLine& other) const;

    /**
     Compares input line object with other input line object.
     
     @param other the input line object.
     @return true if input line objects are not equal, false otherwise.
     */
    bool operator!=(const CInputLine& other) const;

    /**
     Gets requested keyword from parsed input line.

     @param iKeyword the index of keyword in parsed input line.
     @return the keyword.
     */
    std::string getKeyword(const size_t iKeyword) const;

    /**
     Gets requested uppercased keyword from parsed input line.

     @param iKeyword the index of keyword in parsed input line.
     @return the uppercased keyword.
     */
    std::string getUpcasedKeyword(const size_t iKeyword) const;

    /**
     Gets number of keywords in parsed input line.

     @return the number of keywords.
     */
    size_t getNumberOfKeywords() const;

    /**
     Determines if requested keyword in parsed input line is real number.

     @param iKeyword the index of keyword in parsed input line.
     @return true if requested keyword is real number, false otherwise.
     */
    bool isRealNumber(const size_t iKeyword) const;

    /**
     Converts request keyword to real number.

     @param iKeyword the index of keyword in parsed input line.
     @return the real number.
     */
    double getRealNumber(const size_t iKeyword) const;

    /**
     Determines if requested keyword in parsed input line is integer number.

     @param iKeyword the index of keyword in parsed input line.
     @return true if requested keyword is integer number, false otherwise.
     */
    bool isIntegerNumber(const size_t iKeyword) const;

    /**
     Converts request keyword to integer number.

     @param iKeyword the index of keyword in parsed input line.
     @return the integer number.
     */
    int32_t getIntegerNumber(const size_t iKeyword) const;

    /**
     Determines if requested keyword is equal to label.

     @param iKeyword the index of keyword in parsed input line.
     @param keyLabel the label.
     @return true if requested keyword is equal label, false otherwise.
     */
    bool isKeyword(const size_t iKeyword, const std::string& keyLabel) const;

    /**
     Determines if requested keyword is equal to label.

     @param iKeyword the index of keyword in parsed input line.
     @param keyLabel the label.
     @return true if requested keyword is equal label, false otherwise.
     */
    bool isKeyword(const size_t iKeyword, const char* keyLabel) const;

    /**
     Determines if parsed input line starts from conrol keyword equal to label.

     @param keyLabel the label.
     @return true if parsed input line starts from control keyword equal to
             label, false otherwise.
     */
    bool isControlKeyword(const std::string& keyLabel) const;

    /**
     Determines if parsed input line starts from conrol keyword equal to label.

     @param keyLabel the label.
     @return true if parsed input line starts from control keyword equal to
             label, false otherwise.
     */
    bool isControlKeyword(const char* keyLabel) const;

    /**
     Determines if parsed input line starts from control keyword.

     @return true if parsed input line starts from control keyword, false
             otherwise.
     */
    bool isControlLine() const;

    /**
     Gets parsed input line from input line object.

     @return the parsed input line.
     */
    std::string getParsedString() const;

    /**
     Gets original input line from input line object.

     @return the original input line.
     */
    std::string getOriginalString() const;

    /**
     Determines if parsed input line is empty.

     @return true if parsed input line is empty, false otherwise.
     */
    bool isEmpty() const;

    /**
     Resets original and parsed input lines to empty strings.
     */
    void clear();

    /**
     Converts input line object to text line and insert it into output text
     stream.

     @param output the output text stream.
     @param source the input line object.
     */
    friend std::ostream& operator<<(std::ostream& output,
                                    const CInputLine& source);
};

#endif /* InputLine_hpp */
