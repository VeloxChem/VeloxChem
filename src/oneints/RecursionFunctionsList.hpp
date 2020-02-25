//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef RecursionFunctionsList_hpp
#define RecursionFunctionsList_hpp

#include "RecursionFunction.hpp"

/**
 Class CRecursionFunctionsList stores information about list of recursion function
 objects and provides methods for manipulating with list data.

 @author Z. Rinkevicius
 */
class CRecursionFunctionsList
{
    /**
     The vector of recursion function objects.
     */
    std::vector<CRecursionFunction> _recursionFunctions;

   public:
    /**
     Creates an empty recursion functions list object.
     */
    CRecursionFunctionsList();

    /**
     Creates a recursion functions list object.

     @param recursionFunctions the vector of recursion function objects.
     */
    CRecursionFunctionsList(const std::vector<CRecursionFunction>& recursionFunctions);

    /**
     Creates a recursion functions list object by copying other recursion functions list object.

     @param source the recursion functions list object.
     */
    CRecursionFunctionsList(const CRecursionFunctionsList& source);

    /**
     Creates a recursion functions list object by moving other recursion functions list object.

     @param source the recursion functions list object.
     */
    CRecursionFunctionsList(CRecursionFunctionsList&& source) noexcept;

    /**
     Destroys a recursion functions list object.
     */
    ~CRecursionFunctionsList();

    /**
     Assigns a recursion functions list object by copying other recursion functions list object.

     @param source the recursion functions list object.
     */
    CRecursionFunctionsList& operator=(const CRecursionFunctionsList& source);

    /**
     Assigns a recursion functions list object by moving other recursion functions list object.

     @param source the recursion functions list object.
     */
    CRecursionFunctionsList& operator=(CRecursionFunctionsList&& source) noexcept;

    /**
     Compares recursion functions list object with other recursion functions list object.

     @param other the recursion functions list object.
     @return true if recursion functions list objects are equal, false otherwise.
     */
    bool operator==(const CRecursionFunctionsList& other) const;

    /**
     Compares recursion functions list object with other recursion functions list object.

     @param other the recursion functions list object.
     @return true if recursion functions list objects are not equal, false otherwise.
     */
    bool operator!=(const CRecursionFunctionsList& other) const;

    /**
     Adds unique recursion function to recursion functions list.

     @param recursionFunction the recursion function.
     */
    void add(const CRecursionFunction& recursionFunction);

    /**
     Applies selected recursion function to recursion term object.

     @param recursionTerm the recursion term object.
     @return the vector of recursion term objects.
     */
    std::vector<CRecursionTerm> compute(const CRecursionTerm& recursionTerm) const;

    /**
     Finds recursion function index by it's label and recursion targed in
     recursion functions list.

     @param label the factor label.
     @return the index of recursion function in recursion functions list object.
     */
    int32_t find(const std::string& label) const;

    /**
     Converts recursion functions list object to text format and insert it into
     output text stream.

     @param output the output text stream.
     @param source the factors list object.
     */
    friend std::ostream& operator<<(std::ostream& output, const CRecursionFunctionsList& source);
};

#endif /* RecursionFunctionsList_hpp */
