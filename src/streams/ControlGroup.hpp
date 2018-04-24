//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ControlGroup_hpp
#define ControlGroup_hpp

#include <ostream>
#include <vector>

#include "InputLine.hpp"

/**
 Class CControlGroup stores information about single control group and provides
 set of methods for manipulating control group data.
 
 @author Z. Rinkevicius
 */
class CControlGroup
{
    /**
     The header line of control group
     */
    CInputLine _header;

    /**
     The vector of command lines in control group.
     */
    std::vector<CInputLine> _commands;

public:

    /**
     Creates an empty control group object.
     */
    CControlGroup();

    /**
     Creates a control group object by copying other control group object.
     
     @param source the control group object.
     */
    CControlGroup(const CControlGroup& source);
    
    /**
     Creates a control group object by by moving other control group object.
     
     @param source the control group object.
     */
    CControlGroup(CControlGroup&& source) noexcept;

    /**
     Destroys a control group object.
     */
    ~CControlGroup();

    /**
     Assigns a control group object by copying other control group object.
     
     @param source the control group object.
     */
    CControlGroup& operator=(const CControlGroup& source);

    /**
     Assigns a control group object by moving other control group object.
     
     @param source the control group object.
     */
    CControlGroup& operator=(CControlGroup&& source) noexcept;

    /**
     Compares control group object with other control group object.
     
     @param other the control group object.
     @return true if control group objects are equal, false otherwise.
     */
    bool operator==(const CControlGroup& other) const;

    /**
     Compares control group object with other control group object.
     
     @param other the control group object.
     @return true if control group objects are not equal, false otherwise.
     */
    bool operator!=(const CControlGroup& other) const;

    /**
     Sets the header of control group.

     @param header the input line object.
     */
    void setHeader(const CInputLine& header);

    /**
     Adds input line object to vector of commands.

     @param command the input line object.
     */
    void addCommand(const CInputLine& command);

    /**
     Sets header to empty input line object, and empties vector of command
     lines.
     */
    void clear();

    /**
     Determines if vector of command lines is empty.

     @return true if vector of command lines is empty, false otherwise.
     */
    bool isEmpty() const;

    /**
     Determines if header starts from control keyword.

     @param name the name of control keyword.
     @return true if header starts from control keyword with requested name,
             false otherwise.
     */
    bool isNameOfControlGroup(const std::string& name) const;

    /**
     Determines if header starts from control keyword.
     
     @param name the name of control keyword.
     @return true if header starts from control keyword with requested name,
     false otherwise.
     */
    bool isNameOfControlGroup(const char* name) const;

    /**
     Determines number of command lines in vector of command lines.

     @return the number of command lines in vector of command lines.
     */
    int32_t getNumberOfCommands() const;

    /**
     Determines number of input lines starting with requested keyword
     in vector of command lines.

     @param keyword the keyword.
     @return the number of command lines starting with requested keyword in
             vector of command lines.
     */
    int32_t getNumberOfCommands(const std::string& keyword) const;

    /**
     Determines number of input lines starting with requested keyword
     in vector of command lines.
     
     @param keyword the keyword.
     @return the number of command lines starting with requested keyword in
     vector of command lines.
     */
    int32_t getNumberOfCommands(const char* keyword) const;

    /**
     Gets requested command line from vector of command lines.

     @param index the index of requested command line in vector of command
            lines.
     @return the command line.
     */
    CInputLine getCommand(const int32_t index) const;

    /**
     Gets command line starting with request keyword from vector of command
     lines.

     @param keyword the keyword.
     @return the command line.
     */
    CInputLine getCommand(const std::string& keyword) const;

    /**
     Gets command line starting with request keyword from vector of command
     lines.
     
     @param keyword the keyword.
     @return the command line.
     */
    CInputLine getCommand(const char* keyword) const;

    /**
     Converts output line object to formatted text and insert it into
     output text stream.
     
     @param output the output text stream.
     @param source the output line object.
     */
    friend std::ostream& operator<<(std::ostream& output,
                                    const CControlGroup& source);
};

#endif /* ControlGroup_hpp */
