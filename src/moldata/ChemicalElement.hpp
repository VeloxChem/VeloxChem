#ifndef ChemicalElement_hpp
#define ChemicalElement_hpp

#include <cstdint>
#include <string>

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
    std::string _label;

    /**
     The charge of chemical element.
     */
    double _charge;

    /**
     The mass of chemical element.
     */
    double _mass;

    /**
     The number of chemical element i.e. index in periodic table.
     */
    int64_t _number;

    /**
     Sets chemical element properties to properties of fictonal "dummy"
     element.
     */
    auto _selectDummyAtom() -> void;

    /**
     Sets chemical element properties to properties of hydrogen most abudant
     isotope.
     */
    auto _selectHydrogenAtom() -> void;

    /**
     Sets chemical element properties to properties of helium most abudant
     isotope.
     */
    auto _selectHeliumAtom() -> void;

    /**
     Sets chemical element properties to properties of lithium most abudant
     isotope.
     */
    auto _selectLithiumAtom() -> void;

    /**
     Sets chemical element properties to properties of beryllium most abudant
     isotope.
     */
    auto _selectBerylliumAtom() -> void;

    /**
     Sets chemical element properties to properties of boron most abudant
     isotope.
     */
    auto _selectBoronAtom() -> void;

    /**
     Sets chemical element properties to properties of carbon most abudant
     isotope.
     */
    auto _selectCarbonAtom() -> void;

    /**
     Sets chemical element properties to properties of nitrogen most abudant
     isotope.
     */
    auto _selectNitrogenAtom() -> void;

    /**
     Sets chemical element properties to properties of oxygen most abudant
     isotope.
     */
    auto _selectOxygenAtom() -> void;

    /**
     Sets chemical element properties to properties of flourine most abudant
     isotope.
     */
    auto _selectFlourineAtom() -> void;

    /**
     Sets chemical element properties to properties of neon most abudant
     isotope.
     */
    auto _selectNeonAtom() -> void;

    /**
     Sets chemical element properties to properties of sodium most abudant
     isotope.
     */
    auto _selectSodiumAtom() -> void;

    /**
     Sets chemical element properties to properties of magnesium most abudant
     isotope.
     */
    auto _selectMagnesiumAtom() -> void;

    /**
     Sets chemical element properties to properties of aluminium most abudant
     isotope.
     */
    auto _selectAluminiumAtom() -> void;

    /**
     Sets chemical element properties to properties of silicon most abudant
     isotope.
     */
    auto _selectSiliconAtom() -> void;

    /**
     Sets chemical element properties to properties of phosphorus most abudant
     isotope.
     */
    auto _selectPhosphorusAtom() -> void;

    /**
     Sets chemical element properties to properties of sulfur most abudant
     isotope.
     */
    auto _selectSulfurAtom() -> void;

    /**
     Sets chemical element properties to properties of chlorine most abudant
     isotope.
     */
    auto _selectChlorineAtom() -> void;

    /**
     Sets chemical element properties to properties of argon most abudant
     isotope.
     */
    auto _selectArgonAtom() -> void;

    /**
     Sets chemical element properties to properties of potasium most abudant
     isotope.
     */
    auto _selectPotasiumAtom() -> void;

    /**
     Sets chemical element properties to properties of calcium most abudant
     isotope.
     */
    auto _selectCalciumAtom() -> void;

    /**
     Sets chemical element properties to properties of scandium most abudant
     isotope.
     */
    auto _selectScandiumAtom() -> void;

    /**
     Sets chemical element properties to properties of titanium most abudant
     isotope.
     */
    auto _selectTitaniumAtom() -> void;

    /**
     Sets chemical element properties to properties of vanadium most abudant
     isotope.
     */
    auto _selectVanadiumAtom() -> void;

    /**
     Sets chemical element properties to properties of chromium most abudant
     isotope.
     */
    auto _selectChromiumAtom() -> void;

    /**
     Sets chemical element properties to properties of manganese most abudant
     isotope.
     */
    auto _selectManganeseAtom() -> void;

    /**
     Sets chemical element properties to properties of iron most abudant
     isotope.
     */
    auto _selectIronAtom() -> void;

    /**
     Sets chemical element properties to properties of cobalt most abudant
     isotope.
     */
    auto _selectCobaltAtom() -> void;

    /**
     Sets chemical element properties to properties of nickel most abudant
     isotope.
     */
    auto _selectNickelAtom() -> void;

    /**
     Sets chemical element properties to properties of copper most abudant
     isotope.
     */
    auto _selectCopperAtom() -> void;

    /**
     Sets chemical element properties to properties of zinc most abudant
     isotope.
     */
    auto _selectZincAtom() -> void;

    /**
     Sets chemical element properties to properties of galium most abudant
     isotope.
     */
    auto _selectGaliumAtom() -> void;

    /**
     Sets chemical element properties to properties of germanium most abudant
     isotope.
     */
    auto _selectGermaniumAtom() -> void;

    /**
     Sets chemical element properties to properties of arsenic most abudant
     isotope.
     */
    auto _selectArsenicAtom() -> void;

    /**
     Sets chemical element properties to properties of selenium most abudant
     isotope.
     */
    auto _selectSeleniumAtom() -> void;

    /**
     Sets chemical element properties to properties of bromine most abudant
     isotope.
     */
    auto _selectBromineAtom() -> void;

    /**
     Sets chemical element properties to properties of krypton most abudant
     isotope.
     */
    auto _selectKryptonAtom() -> void;

    /**
     Sets chemical element properties to properties of rubidium most abudant
     isotope.
     */
    auto _selectRubidiumAtom() -> void;

    /**
     Sets chemical element properties to properties of strontium most abudant
     isotope.
     */
    auto _selectStrontiumAtom() -> void;

    /**
     Sets chemical element properties to properties of yttrium most abudant
     isotope.
     */
    auto _selectYttriumAtom() -> void;

    /**
     Sets chemical element properties to properties of zirconium most abudant
     isotope.
     */
    auto _selectZirconiumAtom() -> void;

    /**
     Sets chemical element properties to properties of niobium most abudant
     isotope.
     */
    auto _selectNiobiumAtom() -> void;

    /**
     Sets chemical element properties to properties of molybdenum most abudant
     isotope.
     */
    auto _selectMolybdenumAtom() -> void;

    /**
     Sets chemical element properties to properties of technetium most abudant
     isotope.
     */
    auto _selectTechnetiumAtom() -> void;

    /**
     Sets chemical element properties to properties of ruthenium most abudant
     isotope.
     */
    auto _selectRutheniumAtom() -> void;

    /**
     Sets chemical element properties to properties of rhodium most abudant
     isotope.
     */
    auto _selectRhodiumAtom() -> void;

    /**
     Sets chemical element properties to properties of paladium most abudant
     isotope.
     */
    auto _selectPaladiumAtom() -> void;

    /**
     Sets chemical element properties to properties of silver most abudant
     isotope.
     */
    auto _selectSilverAtom() -> void;

    /**
     Sets chemical element properties to properties of cadmium most abudant
     isotope.
     */
    auto _selectCadmiumAtom() -> void;

    /**
     Sets chemical element properties to properties of indium most abudant
     isotope.
     */
    auto _selectIndiumAtom() -> void;

    /**
     Sets chemical element properties to properties of tin most abudant
     isotope.
     */
    auto _selectTinAtom() -> void;

    /**
     Sets chemical element properties to properties of antimony most abudant
     isotope.
     */
    auto _selectAntimonyAtom() -> void;

    /**
     Sets chemical element properties to properties of tellurium most abudant
     isotope.
     */
    auto _selectTelluriumAtom() -> void;

    /**
     Sets chemical element properties to properties of iodine most abudant
     isotope.
     */
    auto _selectIodineAtom() -> void;

    /**
     Sets chemical element properties to properties of xeon most abudant
     isotope.
     */
    auto _selectXenonAtom() -> void;

    /**
     Sets chemical element properties to properties of cesium most abudant
     isotope.
     */
    auto _selectCesiumAtom() -> void;

    /**
     Sets chemical element properties to properties of barium most abudant
     isotope.
     */
    auto _selectBariumAtom() -> void;

    /**
     Sets chemical element properties to properties of lanthanum most abudant
     isotope.
     */
    auto _selectLanthanumAtom() -> void;

    /**
     Sets chemical element properties to properties of cerium most abudant
     isotope.
     */
    auto _selectCeriumAtom() -> void;

    /**
     Sets chemical element properties to properties of praseodymium most
     abudant isotope.
     */
    auto _selectPraseodymiumAtom() -> void;

    /**
     Sets chemical element properties to properties of neodymium most abudant
     isotope.
     */
    auto _selectNeodymiumAtom() -> void;

    /**
     Sets chemical element properties to properties of promethium most abudant
     isotope.
     */
    auto _selectPromethiumAtom() -> void;

    /**
     Sets chemical element properties to properties of samarium most abudant
     isotope.
     */
    auto _selectSamariumAtom() -> void;

    /**
     Sets chemical element properties to properties of europium most abudant
     isotope.
     */
    auto _selectEuropiumAtom() -> void;

    /**
     Sets chemical element properties to properties of gadolinium most abudant
     isotope.
     */
    auto _selectGadoliniumAtom() -> void;

    /**
     Sets chemical element properties to properties of terbium most abudant
     isotope.
     */
    auto _selectTerbiumAtom() -> void;

    /**
     Sets chemical element properties to properties of dysprosium most abudant
     isotope.
     */
    auto _selectDysprosiumAtom() -> void;

    /**
     Sets chemical element properties to properties of holmium most abudant
     isotope.
     */
    auto _selectHolmiumAtom() -> void;

    /**
     Sets chemical element properties to properties of erbium most abudant
     isotope.
     */
    auto _selectErbiumAtom() -> void;

    /**
     Sets chemical element properties to properties of thulium most abudant
     isotope.
     */
    auto _selectThuliumAtom() -> void;

    /**
     Sets chemical element properties to properties of ytterbium most abudant
     isotope.
     */
    auto _selectYtterbiumAtom() -> void;

    /**
     Sets chemical element properties to properties of luthenium most abudant
     isotope.
     */
    auto _selectLutheniumAtom() -> void;

    /**
     Sets chemical element properties to properties of hafnium most abudant
     isotope.
     */
    auto _selectHafniumAtom() -> void;

    /**
     Sets chemical element properties to properties of tantalum most abudant
     isotope.
     */
    auto _selectTantalumAtom() -> void;

    /**
     Sets chemical element properties to properties of tungsten most abudant
     isotope.
     */
    auto _selectTungstenAtom() -> void;

    /**
     Sets chemical element properties to properties of rhenium most abudant
     isotope.
     */
    auto _selectRheniumAtom() -> void;

    /**
     Sets chemical element properties to properties of osmium most abudant
     isotope.
     */
    auto _selectOsmiumAtom() -> void;

    /**
     Sets chemical element properties to properties of iridium most abudant
     isotope.
     */
    auto _selectIridiumAtom() -> void;

    /**
     Sets chemical element properties to properties of platinum most abudant
     isotope.
     */
    auto _selectPlatinumAtom() -> void;

    /**
     Sets chemical element properties to properties of gold most abudant
     isotope.
     */
    auto _selectGoldAtom() -> void;

    /**
     Sets chemical element properties to properties of mercury most abudant
     isotope.
     */
    auto _selectMercuryAtom() -> void;

    /**
     Sets chemical element properties to properties of thallium most abudant
     isotope.
     */
    auto _selectThalliumAtom() -> void;

    /**
     Sets chemical element properties to properties of lead most abudant
     isotope.
     */
    auto _selectLeadAtom() -> void;

    /**
     Sets chemical element properties to properties of bismuth most abudant
     isotope.
     */
    auto _selectBismuthAtom() -> void;

    /**
     Sets chemical element properties to properties of polonium most abudant
     isotope.
     */
    auto _selectPoloniumAtom() -> void;

    /**
     Sets chemical element properties to properties of astatine most abudant
     isotope.
     */
    auto _selectAstatineAtom() -> void;

    /**
     Sets chemical element properties to properties of radon most abudant
     isotope.
     */
    auto _selectRadonAtom() -> void;

    /**
     Sets chemical element mass according to given isotope number for hydrogen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectHydrogenIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for helium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectHeliumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for lithium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectLithiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for beryllium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectBerylliumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for baron.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectBoronIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for carbon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCarbonIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for nitrogen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectNitrogenIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for oxygen.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectOxygenIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for flourine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectFlourineIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for neon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectNeonIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for sodium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSodiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for magnesium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectMagnesiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for aluminium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectAluminiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for silicon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSiliconIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for phosphorus.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPhosphorusIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for sulfur.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSulfurIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for chlorine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectChlorineIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for argon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectArgonIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for potasium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPotasiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for calcium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCalciumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for scandium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectScandiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for titanium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTitaniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for vanadium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectVanadiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for chromium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectChromiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for manganese.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectManganeseIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for iron.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectIronIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for cobalt.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCobaltIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for nickel.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectNickelIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for copper.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCopperIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for zinc.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectZincIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for gallium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectGalliumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for germanium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectGermaniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for arsenic.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectArsenicIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for selenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSeleniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for bromine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectBromineIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for krypton.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectKryptonIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for rubidium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectRubidiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for strontium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectStrontiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for yttrium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectYttriumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for zirconium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectZirconiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for niobium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */

    auto _selectNiobiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for molybdenum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectMolybdenumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for technetium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTechnetiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for ruthenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectRutheniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for rhodium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectRhodiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for paladium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPaladiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for silver.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSilverIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for cadmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCadmiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for indium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectIndiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for tin.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTinIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for antimony.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectAntimonyIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for tellurium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTelluriumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for iodine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectIodineIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for xenon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectXenonIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for cesium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCesiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for barium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectBariumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for lanthanum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectLanthanumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for cerium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectCeriumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for praseodymium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPraseodymiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for neodymium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectNeodymiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for promethium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPromethiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for samarium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectSamariumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for europium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectEuropiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for gadolinium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectGadoliniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for terbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTerbiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for dysprosium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectDysprosiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for holmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectHolmiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for erbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectErbiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for thulium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectThuliumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for ytterbium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectYtterbiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for luthenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectLutheniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for hafnium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectHafniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for tantalum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTantalumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for tungsten.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectTungstenIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for rhenium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectRheniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for osmium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectOsmiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for iridium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectIridiumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for platinum.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPlatinumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for gold.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectGoldIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for mercury.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectMercuryIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for thallium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectThalliumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for lead.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectLeadIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for bismuth.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectBismuthIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for polonium.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectPoloniumIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for astatine.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectAstatineIsotopeMass(const int64_t label) -> bool;

    /**
     Sets chemical element mass according to given isotope number for radon.
     Special cases:
     a) if isotope number is zero, sets mass to mass of most abudant isotope
     of chemical element.
     b) if unsupported isotope number is given, mass of chemical element is not
     updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto _selectRadonIsotopeMass(const int64_t label) -> bool;

   public:
    /**
     Creates an empty chemical element object.
     */
    CChemicalElement() = default;

    /**
     Constructs a chemical element with specific atomic properties.

     @param label the name of chemical element as string.
     @param charge the charge.
     @param mass the mass.
     @param number the chemical element number.
     */
    CChemicalElement(const std::string& label, const double charge, const double mass, const int64_t number);

    /**
     Sets chemical element properties using name of chemical element.
     Special cases:
     a) if name is "BQ", sets chemical element properties to "dummy" atom.
     b) if unsupported name of chemical element is given, chemical element
     properties are not updated and function returns false.

     @param label the name of chemical element.
     @return true if success, false otherwise.
     */
    auto setAtomType(const std::string& label) -> bool;

    /**
     Sets chemical element object properties using chemical element number.
     Special cases:
     a) if chemical element number is zero, sets chemical element properties
     to "dummy" atom.
     b) if unsupported chemical element number is given, chemical element
     properties are not updated and function returns false.

     @param identifier the chemical element number.
     @return true if success, false otherwise.
     */
    auto setAtomType(const int64_t identifier) -> bool;

    /**
     Sets chemical element object properties dependent on isotope using isotope
     number.
     Special cases:
     a) if isotope number is zero, sets chemical element properties to most
     abudant isotope.
     b) if unsupported isotope number is given, chemical element properties
     are not updated and function returns false.

     @param label the isotope number.
     @return true if success, false otherwise.
     */
    auto setIsotope(const int64_t label) -> bool;

    /**
     Gets name of chemical element.

     @return the name of chemical element.
     */
    auto getName() const -> std::string;

    /**
     Gets chemical element number.

     @return the chemical element number.
     */
    auto getIdentifier() const -> int64_t;

    /**
     Gets mass of chemical element (isotope specific).

     @return the mass of chemical element.
     */
    auto getAtomicMass() const -> double;

    /**
     Gets charge of chemical element.

     @return the charge of chemical element.
     */
    auto getAtomicCharge() const -> double;

    /**
     Gets maximum angular momentum of occupied atomic shells in chemical
     element.

     @return the maximum angular momentum.
     */
    auto getMaxAngularMomentum() const -> int64_t;

    /**
     Gets maximum number of supported chemical elements.

     @return the charge of chemical element.
     */
    auto getMaxIdentifier() const -> int64_t;
};

#endif /* ChemicalElement_hpp */
