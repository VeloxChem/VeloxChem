//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ChemicalElement_hpp
#define ChemicalElement_hpp

#include <cstdint>
#include <string>
#include <ostream>

/**
 Class CChemicalElement stores information about chemical elements H-Rn and
 provides methods for retrieving various properties of chemical elements H-Rn.
 
 @author Z. Rinkevicius
 */
class CChemicalElement
{
    /**
     The name of chemical element.
     */
    std::string _atomicLabel;

    /**
     The charge of chemical element.
     */
    double _atomicCharge;

    /**
     The mass of chemical element.
     */
    double _atomicMass;

    /**
     The number of chemical element i.e. index in periodic table.
     */
    int32_t _atomicNumber;

    /**
     Sets chemical element properties to properties of fictonal "dummy"
     element.
     */
    void _selectDummyAtom();

    /**
     Sets chemical element properties to properties of hydrogen most abudant
     isotope.
     */
    void _selectHydrogenAtom();

    /**
     Sets chemical element properties to properties of helium most abudant
     isotope.
     */
    void _selectHeliumAtom();

    /**
     Sets chemical element properties to properties of lithium most abudant
     isotope.
     */
    void _selectLithiumAtom();

    /**
     Sets chemical element properties to properties of beryllium most abudant
     isotope.
     */
    void _selectBerylliumAtom();

    /**
     Sets chemical element properties to properties of boron most abudant
     isotope.
     */
    void _selectBoronAtom();

    /**
     Sets chemical element properties to properties of carbon most abudant
     isotope.
     */
    void _selectCarbonAtom();

    /**
     Sets chemical element properties to properties of nitrogen most abudant
     isotope.
     */
    void _selectNitrogenAtom();

    /**
     Sets chemical element properties to properties of oxygen most abudant
     isotope.
     */
    void _selectOxygenAtom();

    /**
     Sets chemical element properties to properties of flourine most abudant
     isotope.
     */
    void _selectFlourineAtom();

    /**
     Sets chemical element properties to properties of neon most abudant
     isotope.
     */
    void _selectNeonAtom();

    /**
     Sets chemical element properties to properties of sodium most abudant
     isotope.
     */
    void _selectSodiumAtom();

    /**
     Sets chemical element properties to properties of magnesium most abudant
     isotope.
     */
    void _selectMagnesiumAtom();

    /**
     Sets chemical element properties to properties of aluminium most abudant
     isotope.
     */
    void _selectAluminiumAtom();

    /**
     Sets chemical element properties to properties of silicon most abudant
     isotope.
     */
    void _selectSiliconAtom();

    /**
     Sets chemical element properties to properties of phosphorus most abudant
     isotope.
     */
    void _selectPhosphorusAtom();

    /**
     Sets chemical element properties to properties of sulfur most abudant
     isotope.
     */
    void _selectSulfurAtom();

    /**
     Sets chemical element properties to properties of chlorine most abudant
     isotope.
     */
    void _selectChlorineAtom();

    /**
     Sets chemical element properties to properties of argon most abudant
     isotope.
     */
    void _selectArgonAtom();

    /**
     Sets chemical element properties to properties of potasium most abudant
     isotope.
     */
    void _selectPotasiumAtom();

    /**
     Sets chemical element properties to properties of calcium most abudant
     isotope.
     */
    void _selectCalciumAtom();

    /**
     Sets chemical element properties to properties of scandium most abudant
     isotope.
     */
    void _selectScandiumAtom();

    /**
     Sets chemical element properties to properties of titanium most abudant
     isotope.
     */
    void _selectTitaniumAtom();

    /**
     Sets chemical element properties to properties of vanadium most abudant
     isotope.
     */
    void _selectVanadiumAtom();

    /**
     Sets chemical element properties to properties of chromium most abudant
     isotope.
     */
    void _selectChromiumAtom();

    /**
     Sets chemical element properties to properties of manganese most abudant
     isotope.
     */
    void _selectManganeseAtom();

    /**
     Sets chemical element properties to properties of iron most abudant
     isotope.
     */
    void _selectIronAtom();

    /**
     Sets chemical element properties to properties of cobalt most abudant
     isotope.
     */
    void _selectCobaltAtom();

    /**
     Sets chemical element properties to properties of nickel most abudant
     isotope.
     */
    void _selectNickelAtom();

    /**
     Sets chemical element properties to properties of copper most abudant
     isotope.
     */
    void _selectCopperAtom();

    /**
     Sets chemical element properties to properties of zinc most abudant
     isotope.
     */
    void _selectZincAtom();

    /**
     Sets chemical element properties to properties of galium most abudant
     isotope.
     */
    void _selectGaliumAtom();

    /**
     Sets chemical element properties to properties of germanium most abudant
     isotope.
     */
    void _selectGermaniumAtom();

    /**
     Sets chemical element properties to properties of arsenic most abudant
     isotope.
     */
    void _selectArsenicAtom();

    /**
     Sets chemical element properties to properties of selenium most abudant
     isotope.
     */
    void _selectSeleniumAtom();

    /**
     Sets chemical element properties to properties of bromine most abudant
     isotope.
     */
    void _selectBromineAtom();

    /**
     Sets chemical element properties to properties of krypton most abudant
     isotope.
     */
    void _selectKryptonAtom();

    /**
     Sets chemical element properties to properties of rubidium most abudant
     isotope.
     */
    void _selectRubidiumAtom();

    /**
     Sets chemical element properties to properties of strontium most abudant
     isotope.
     */
    void _selectStrontiumAtom();

    /**
     Sets chemical element properties to properties of yttrium most abudant
     isotope.
     */
    void _selectYttriumAtom();

    /**
     Sets chemical element properties to properties of zirconium most abudant
     isotope.
     */
    void _selectZirconiumAtom();

    /**
     Sets chemical element properties to properties of niobium most abudant
     isotope.
     */
    void _selectNiobiumAtom();

    /**
     Sets chemical element properties to properties of molybdenum most abudant
     isotope.
     */
    void _selectMolybdenumAtom();

    /**
     Sets chemical element properties to properties of technetium most abudant
     isotope.
     */
    void _selectTechnetiumAtom();

    /**
     Sets chemical element properties to properties of ruthenium most abudant
     isotope.
     */
    void _selectRutheniumAtom();

    /**
     Sets chemical element properties to properties of rhodium most abudant
     isotope.
     */
    void _selectRhodiumAtom();

    /**
     Sets chemical element properties to properties of paladium most abudant
     isotope.
     */
    void _selectPaladiumAtom();

    /**
     Sets chemical element properties to properties of silver most abudant
     isotope.
     */
    void _selectSilverAtom();
    
    /**
     Sets chemical element properties to properties of cadmium most abudant
     isotope.
     */
    void _selectCadmiumAtom();

    /**
     Sets chemical element properties to properties of indium most abudant
     isotope.
     */
    void _selectIndiumAtom();

    // void _selectTinAtom():
    //
    // Sets chemical element properties to properties of tin most abudant
    // isotope.

    /**
     Sets chemical element properties to properties of tin most abudant
     isotope.
     */
    void _selectTinAtom();

    /**
     Sets chemical element properties to properties of antimony most abudant
     isotope.
     */
    void _selectAntimonyAtom();

    /**
     Sets chemical element properties to properties of tellurium most abudant
     isotope.
     */
    void _selectTelluriumAtom();

    /**
     Sets chemical element properties to properties of iodine most abudant
     isotope.
     */
    void _selectIodineAtom();

    /**
     Sets chemical element properties to properties of xeon most abudant
     isotope.
     */
    void _selectXenonAtom();

    /**
     Sets chemical element properties to properties of cesium most abudant
     isotope.
     */
    void _selectCesiumAtom();

    /**
     Sets chemical element properties to properties of barium most abudant
     isotope.
     */
    void _selectBariumAtom();

    /**
     Sets chemical element properties to properties of lanthanum most abudant
     isotope.
     */
    void _selectLanthanumAtom();

    // void _selectCeriumAtom():
    //
    // Sets chemical element properties to properties of cerium most abudant
    // isotope.

    /**
     Sets chemical element properties to properties of cerium most abudant
     isotope.
     */
    void _selectCeriumAtom();

    /**
     Sets chemical element properties to properties of praseodymium most
     abudant isotope.
     */
    void _selectPraseodymiumAtom();

    /**
     Sets chemical element properties to properties of neodymium most abudant
     isotope.
     */
    void _selectNeodymiumAtom();

    /**
     Sets chemical element properties to properties of promethium most abudant
     isotope.
     */
    void _selectPromethiumAtom();

    /**
     Sets chemical element properties to properties of samarium most abudant
     isotope.
     */
    void _selectSamariumAtom();

    /**
     Sets chemical element properties to properties of europium most abudant
     isotope.
     */
    void _selectEuropiumAtom();

    /**
     Sets chemical element properties to properties of gadolinium most abudant
     isotope.
     */
    void _selectGadoliniumAtom();

    /**
     Sets chemical element properties to properties of terbium most abudant
     isotope.
     */
    void _selectTerbiumAtom();

    /**
     Sets chemical element properties to properties of dysprosium most abudant
     isotope.
     */
    void _selectDysprosiumAtom();

    /**
     Sets chemical element properties to properties of holmium most abudant
     isotope.
     */
    void _selectHolmiumAtom();

    /**
     Sets chemical element properties to properties of erbium most abudant
     isotope.
     */
    void _selectErbiumAtom();

    /**
     Sets chemical element properties to properties of thulium most abudant
     isotope.
     */
    void _selectThuliumAtom();

    /**
     Sets chemical element properties to properties of ytterbium most abudant
     isotope.
     */
    void _selectYtterbiumAtom();

    /**
     Sets chemical element properties to properties of luthenium most abudant
     isotope.
     */
    void _selectLutheniumAtom();

    /**
     Sets chemical element properties to properties of hafnium most abudant
     isotope.
     */
    void _selectHafniumAtom();

    /**
     Sets chemical element properties to properties of tantalum most abudant
     isotope.
     */
    void _selectTantalumAtom();

    /**
     Sets chemical element properties to properties of tungsten most abudant
     isotope.
     */
    void _selectTungstenAtom();

    /**
     Sets chemical element properties to properties of rhenium most abudant
     isotope.
     */
    void _selectRheniumAtom();

    /**
     Sets chemical element properties to properties of osmium most abudant
     isotope.
     */
    void _selectOsmiumAtom();

    /**
     Sets chemical element properties to properties of iridium most abudant
     isotope.
     */
    void _selectIridiumAtom();

    /**
     Sets chemical element properties to properties of platinum most abudant
     isotope.
     */
    void _selectPlatinumAtom();

    /**
     Sets chemical element properties to properties of gold most abudant
     isotope.
     */
    void _selectGoldAtom();

    /**
     Sets chemical element properties to properties of mercury most abudant
     isotope.
     */
    void _selectMercuryAtom();

    /**
     Sets chemical element properties to properties of thallium most abudant
     isotope.
     */
    void _selectThalliumAtom();

    /**
     Sets chemical element properties to properties of lead most abudant
     isotope.
     */
    void _selectLeadAtom();

    /**
     Sets chemical element properties to properties of bismuth most abudant
     isotope.
     */
    void _selectBismuthAtom();

    /**
     Sets chemical element properties to properties of polonium most abudant
     isotope.
     */
    void _selectPoloniumAtom();

    /**
     Sets chemical element properties to properties of astatine most abudant
     isotope.
     */
    void _selectAstatineAtom();

    /**
     Sets chemical element properties to properties of radon most abudant
     isotope.
     */
    void _selectRadonAtom();

    /**
     Sets chemical element mass according to given isotope number for hydrogen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectHydrogenIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for helium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectHeliumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for lithium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectLithiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for beryllium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectBerylliumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for baron.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectBoronIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for carbon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCarbonIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for nitrogen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectNitrogenIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for oxygen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectOxygenIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for flourine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectFlourineIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for neon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectNeonIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for sodium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSodiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for magnesium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectMagnesiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for aluminium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectAluminiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for silicon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSiliconIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for phosphorus.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPhosphorusIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for sulfur.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSulfurIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for chlorine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectChlorineIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for argon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectArgonIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for potasium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPotasiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for calcium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCalciumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for scandium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectScandiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for titanium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTitaniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for vanadium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectVanadiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for chromium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectChromiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for manganese.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectManganeseIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for iron.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectIronIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for cobalt.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCobaltIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for nickel.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectNickelIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for copper.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCopperIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for zinc.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectZincIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for gallium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectGalliumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for germanium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectGermaniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for arsenic.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectArsenicIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for selenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSeleniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for bromine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectBromineIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for krypton.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectKryptonIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for rubidium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectRubidiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for strontium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectStrontiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for yttrium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectYttriumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for zirconium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectZirconiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for niobium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    
    bool _selectNiobiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for molybdenum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectMolybdenumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for technetium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTechnetiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for ruthenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectRutheniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for rhodium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectRhodiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for paladium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPaladiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for silver.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSilverIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for cadmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCadmiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for indium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectIndiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for tin.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTinIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for antimony.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectAntimonyIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for tellurium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTelluriumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for iodine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectIodineIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for xenon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectXenonIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for cesium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCesiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for barium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectBariumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for lanthanum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectLanthanumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for cerium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectCeriumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for praseodymium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPraseodymiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for neodymium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectNeodymiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for promethium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPromethiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for samarium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectSamariumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for europium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectEuropiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for gadolinium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectGadoliniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for terbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTerbiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for dysprosium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectDysprosiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for holmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectHolmiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for erbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectErbiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for thulium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectThuliumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for ytterbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectYtterbiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for luthenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectLutheniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for hafnium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectHafniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for tantalum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTantalumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for tungsten.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectTungstenIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for rhenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectRheniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for osmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectOsmiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for iridium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectIridiumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for platinum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPlatinumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for gold.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectGoldIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for mercury.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectMercuryIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for thallium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectThalliumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for lead.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectLeadIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for bismuth.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectBismuthIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for polonium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectPoloniumIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for astatine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectAstatineIsotopeMass(const int32_t isotopeLabel);

    /**
     Sets chemical element mass according to given isotope number for radon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.
     
     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool _selectRadonIsotopeMass(const int32_t isotopeLabel);

public:

    /**
     Creates an empty chemical element object.
     */
    CChemicalElement();

    /**
     Constructs a chemical element with specific atomic properties.
     
     @param atomicLabel the name of chemical element as string.
     @param atomicCharge the charge.
     @param atomicMass the mass.
     @paeam atomicNumber the chemical element number.
     */
    CChemicalElement(const std::string& atomicLabel, const double atomicCharge,
                     const double atomicMass, const int32_t atomicNumber);

    /**
     Destroys a chemical element object.
     */
    ~CChemicalElement();

    /**
     Compares chemical element object with other chemical element object.
     
     @param other the chemical element object.
     @return true if chemical element objects are equal, false otherwise.
     */
    bool operator==(const CChemicalElement& other) const;

    /**
     Compares chemical element object with other chemical element object.
     
     @param other the chemical element object.
     @return true if chemical element objects are not equal, false otherwise.
     */
    bool operator!=(const CChemicalElement& other) const;

    /**
     Sets chemical element properties using name of chemical element.
     Special cases:
     a) if name is "BQ", sets chemical element properties to "dummy" atom.
     b) if unsupported name of chemical element is given, chemical element
     properties are not updated and function returns false.

     @param atomLabel the name of chemical element.
     @return true if success, false otherwise.
     */
    bool setAtomType(const std::string& atomLabel);

    /**
     Sets chemical element object properties using chemical element number.
     Special cases:
     a) if chemical element number is zero, sets chemical element properties
     to "dummy" atom.
     b) if unsupported chemical element number is given, chemical element
     properties are not updated and function returns false.

     @param idElemental the chemical element number.
     @return true if success, false otherwise.
     */
    bool setAtomType(const int32_t idElemental);

    /**
     Sets chemical element object properties dependent on isotope using isotope
     number.
     Special cases:
     a) if isotope number is zero, sets chemical element properties to most
     abudant isotope.
     b) if unsupported isotope number is given, chemical element properties
     are not updated and function returns false.

     @param isotopeLabel the isotope number.
     @return true if success, false otherwise.
     */
    bool setIsotope(const int32_t isotopeLabel);

    /**
     Gets name of chemical element.

     @return the name of chemical element.
     */
    std::string getName() const;

    /**
     Gets chemical element number.

     @return the chemical element number.
     */
    int32_t getIdentifier() const;

    /**
     Gets mass of chemical element (isotope specific).

     @return the mass of chemical element.
     */
    double getAtomicMass() const;

    /**
     Gets charge of chemical element.

     @return the charge of chemical element.
     */
    double getAtomicCharge() const;

    /**
     Converts chemical element object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the chemical element object.
     */
    friend std::ostream& operator<<(std::ostream& output,
                                    const CChemicalElement& source);
};

#endif /* ChemicalElement_hpp */
