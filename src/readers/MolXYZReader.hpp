//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolXYZReader_hpp
#define MolXYZReader_hpp

#include <cstdint>
#include <vector>

#include "InputLine.hpp"
#include "InputData.hpp"
#include "OutputStream.hpp"
#include "Molecule.hpp"

/**
 Class CMolXYZReader handles parsing of @molxyz control group, which
 contains molecular geometry.
 
 @author Z. Rinkevicius
 */
class CMolXYZReader
{
    /**
     The state of molxyz reader object : true - no errors, false - otherwise.
     */
    bool _state;

    /**
     The units conversion flag: if true units conversion of molecular geometry
     is required.
     */
    bool _needUnitConversion;
    
    /**
     The Cartesian coordinates of atoms.
     */
    std::vector<double> _coordinates;
    
    /**
     The charges of atoms.
     */
    std::vector<double> _charges;
    
    /**
     The masses of atoms.
     */
    std::vector<double> _masses;
    
    /**
     The names of atoms.
     */
    std::vector<std::string> _labels;

    /**
     The chemical element identifiers of atoms.
     */
    std::vector<int32_t> _idsElemental;
    
    /**
     The charge of molecule.
     */
    double _molCharge;
    
    /**
     The multiplicity of molecule.
     */
    int32_t _molMultiplicity;
    
    /**
     Reads and parses header line of @molxyz control group, determones charge
     and multiplicity of molecule. Parsing errors are printed to output stream.

     @param inputLine the input line object.
     @param oStream the output stream.
     */
    void _parseHeader(const CInputLine&    inputLine,
                            COutputStream& oStream);
    
    /**
     Reads and parses atomic data line from @molxyz control group, adds it to
     correcpnding atomic property vector. Parsing errors are printed to output
     stream.

     @param inputLine the input line object.
     @param iAtom the index of parsed atom.
     @param oStream the output stream.
     */
    void _parseAtom(const CInputLine&    inputLine,
                    const int32_t        iAtom,
                          COutputStream& oStream);
   
    /**
     Reads and parses @molxyz line in @molxyz control group, sets up flag for
     units conversion. Parsing errors are printed to output stream.

     @param inputLine the input line object.
     @param oStream the output stream.
     */
    void _parseUnits(const CInputLine&    inputLine,
                           COutputStream& oStream);

    /**
     Prints error message for multiple definitions of @molxyz control group to
     output stream and sets molxyz reader object state to abnormal.

     @param oStream the output stream.
     */
    void _errorUniqueGroup(COutputStream& oStream);

    /**
     Prints error message for wrong multiplicity of molecule to output stream
     // and sets molxyz reader object state to abnormal.

     @param inputLine the header line of @molxyz control group.
     @param oStream the output stream.
     */
    void _errorMultiplicity(const CInputLine&    inputLine,
                                  COutputStream& oStream);

    /**
     Prints error message for wrong definition of units in @molxyz line of
     @molxyz control group to output stream and sets molxyz reader object state
     to abnormal.

     @param inputLine the @molxyz line of @molxyz control group.
     @param oStream the output stream.
     */
    void _errorUnits(const CInputLine&    inputLine,
                           COutputStream& oStream);

    /**
     Prints error message for encounter of unsupported chemical element to
     output stream and sets molxyz reader object state to abnormal.

     @param inputLine the input line object with unsuppoprted chemical element.
     @param oStream the output stream.
     */
    void _errorChemicalElement(const CInputLine&    inputLine,
                                     COutputStream& oStream);

    /**
     Prints error message for encounter of unsupported isotope of chemical
     element to output stream and sets molxyz reader object state to abnormal.

     @param inputLine the input line object with unsuppoprted isotope of
            chemical element.
     @param oStream the output stream.
     */
    void _errorIsotope(const CInputLine&    inputLine,
                             COutputStream& oStream);

    /**
     Prints error message for syntax error of imcomplete @molxyz control group
     to output stream and sets molxyz reader object state to abnormal.

     @param oStream the output stream.
     */
    void _syntaxIncompleteGroup(COutputStream& oStream);

    /**
     Prints error message for syntax error in header of @molxyz control group
     to output stream and sets molxyz reader object state to abnormal.

     @param inputLine the header line of @molxyz control group.
     @param oStream the output stream.
     */
    void _syntaxHeader(const CInputLine&    inputLine,
                             COutputStream& oStream);

    /**
     Prints error message for syntax error in atom input to output stream and
     sets molxyz reader object state to abnormal.

     @param inputLine  the input line object with atom input.
     @param oStream the output stream.
     */
    void _syntaxAtom(const CInputLine&    inputLine,
                           COutputStream& oStream);

    /**
     Determines if multiplicity of molecule is correctly defined in header of
     @molxyz. If multiplicity is set incorrectly, prints error message to output
     stream and sets molxyz reader object state to abnormal.

     @param molecule the molecule.
     @param inputLine the header line of @molxyz control group.
     @param oStream the output stream.
     */
    void _checkMultiplicity(const CMolecule&     molecule,
                            const CInputLine&    inputLine,
                                  COutputStream& oStream);
    
    /**
     Sets up dimensions of atomic data vectors.

     @param nAtoms the number of atoms.
     */
    void _setDimensions(const int32_t nAtoms);

public:

    /**
     Creates a molxyz reader object.
     */
    CMolXYZReader();

    /**
     Destroys a molxyz reader object.
     */
    ~CMolXYZReader();

    /**
     Get a state of molxyz reader object.
     
     @return true if no errors, false otherwise.
     */
    bool getState() const;

    /**
     Creates molecule by reading and parsing first @molxyz control group from
     input data object. Parsing errors are printed to output stream.

     @param molecule the molecule object.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void parse(      CMolecule&     molecule,
               const CInputData&    inputData,
                     COutputStream& oStream);
    
    /**
     Creates molecule by reading and parsing specific @molxyz control group from
     input data object. Parsing errors are printed to output stream.

     @param molecule the molecule object.
     @param inputData the input data object.
     @param iGroup the index of @molxyz control group.
     @param oStream the output stream.
     */
    void parse(      CMolecule&     molecule,
               const CInputData&    inputData,
               const int32_t        iGroup,
                     COutputStream& oStream);
};

#endif /* MolXYZReader_hpp */
