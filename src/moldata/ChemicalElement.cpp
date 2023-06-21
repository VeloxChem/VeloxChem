#include "ChemicalElement.hpp"

CChemicalElement::CChemicalElement(const std::string& label,
                                   const double       charge,
                                   const double       mass,
                                   const int64_t      number)

    : _label(label)

    , _charge(charge)

    , _mass(mass)

    , _number(number)
{
}

auto
CChemicalElement::setAtomType(const std::string& label) -> bool
{
    if (label.compare("BQ") == 0)
    {
        _selectDummyAtom();

        return true;
    }

    if (label.compare("H") == 0)
    {
        _selectHydrogenAtom();

        return true;
    }

    if (label.compare("HE") == 0)
    {
        _selectHeliumAtom();

        return true;
    }

    if (label.compare("LI") == 0)
    {
        _selectLithiumAtom();

        return true;
    }

    if (label.compare("BE") == 0)
    {
        _selectBerylliumAtom();

        return true;
    }

    if (label.compare("B") == 0)
    {
        _selectBoronAtom();

        return true;
    }

    if (label.compare("C") == 0)
    {
        _selectCarbonAtom();

        return true;
    }

    if (label.compare("N") == 0)
    {
        _selectNitrogenAtom();

        return true;
    }

    if (label.compare("O") == 0)
    {
        _selectOxygenAtom();

        return true;
    }

    if (label.compare("F") == 0)
    {
        _selectFlourineAtom();

        return true;
    }

    if (label.compare("NE") == 0)
    {
        _selectNeonAtom();

        return true;
    }

    if (label.compare("NA") == 0)
    {
        _selectSodiumAtom();

        return true;
    }

    if (label.compare("MG") == 0)
    {
        _selectMagnesiumAtom();

        return true;
    }

    if (label.compare("AL") == 0)
    {
        _selectAluminiumAtom();

        return true;
    }

    if (label.compare("SI") == 0)
    {
        _selectSiliconAtom();

        return true;
    }

    if (label.compare("P") == 0)
    {
        _selectPhosphorusAtom();

        return true;
    }

    if (label.compare("S") == 0)
    {
        _selectSulfurAtom();

        return true;
    }

    if (label.compare("CL") == 0)
    {
        _selectChlorineAtom();

        return true;
    }

    if (label.compare("AR") == 0)
    {
        _selectArgonAtom();

        return true;
    }

    if (label.compare("K") == 0)
    {
        _selectPotasiumAtom();

        return true;
    }

    if (label.compare("CA") == 0)
    {
        _selectCalciumAtom();

        return true;
    }

    if (label.compare("SC") == 0)
    {
        _selectScandiumAtom();

        return true;
    }

    if (label.compare("TI") == 0)
    {
        _selectTitaniumAtom();

        return true;
    }

    if (label.compare("V") == 0)
    {
        _selectVanadiumAtom();

        return true;
    }

    if (label.compare("CR") == 0)
    {
        _selectChromiumAtom();

        return true;
    }

    if (label.compare("MN") == 0)
    {
        _selectManganeseAtom();

        return true;
    }

    if (label.compare("FE") == 0)
    {
        _selectIronAtom();

        return true;
    }

    if (label.compare("CO") == 0)
    {
        _selectCobaltAtom();

        return true;
    }

    if (label.compare("NI") == 0)
    {
        _selectNickelAtom();

        return true;
    }

    if (label.compare("CU") == 0)
    {
        _selectCopperAtom();

        return true;
    }

    if (label.compare("ZN") == 0)
    {
        _selectZincAtom();

        return true;
    }

    if (label.compare("GA") == 0)
    {
        _selectGaliumAtom();

        return true;
    }

    if (label.compare("GE") == 0)
    {
        _selectGermaniumAtom();

        return true;
    }

    if (label.compare("AS") == 0)
    {
        _selectArsenicAtom();

        return true;
    }

    if (label.compare("SE") == 0)
    {
        _selectSeleniumAtom();

        return true;
    }

    if (label.compare("BR") == 0)
    {
        _selectBromineAtom();

        return true;
    }

    if (label.compare("KR") == 0)
    {
        _selectKryptonAtom();

        return true;
    }

    if (label.compare("RB") == 0)
    {
        _selectRubidiumAtom();

        return true;
    }

    if (label.compare("SR") == 0)
    {
        _selectStrontiumAtom();

        return true;
    }

    if (label.compare("Y") == 0)
    {
        _selectYttriumAtom();

        return true;
    }

    if (label.compare("ZR") == 0)
    {
        _selectZirconiumAtom();

        return true;
    }

    if (label.compare("NB") == 0)
    {
        _selectNiobiumAtom();

        return true;
    }

    if (label.compare("MO") == 0)
    {
        _selectMolybdenumAtom();
        return true;
    }

    if (label.compare("TC") == 0)
    {
        _selectTechnetiumAtom();

        return true;
    }

    if (label.compare("RU") == 0)
    {
        _selectRutheniumAtom();

        return true;
    }

    if (label.compare("RH") == 0)
    {
        _selectRhodiumAtom();

        return true;
    }

    if (label.compare("PD") == 0)
    {
        _selectPaladiumAtom();

        return true;
    }

    if (label.compare("AG") == 0)
    {
        _selectSilverAtom();

        return true;
    }

    if (label.compare("CD") == 0)
    {
        _selectCadmiumAtom();

        return true;
    }

    if (label.compare("IN") == 0)
    {
        _selectIndiumAtom();

        return true;
    }

    if (label.compare("SN") == 0)
    {
        _selectTinAtom();

        return true;
    }

    if (label.compare("SB") == 0)
    {
        _selectAntimonyAtom();

        return true;
    }

    if (label.compare("TE") == 0)
    {
        _selectTelluriumAtom();

        return true;
    }

    if (label.compare("I") == 0)
    {
        _selectIodineAtom();

        return true;
    }

    if (label.compare("XE") == 0)
    {
        _selectXenonAtom();

        return true;
    }

    if (label.compare("CS") == 0)
    {
        _selectCesiumAtom();

        return true;
    }

    if (label.compare("BA") == 0)
    {
        _selectBariumAtom();

        return true;
    }

    if (label.compare("LA") == 0)
    {
        _selectLanthanumAtom();

        return true;
    }

    if (label.compare("CE") == 0)
    {
        _selectCeriumAtom();

        return true;
    }

    if (label.compare("PR") == 0)
    {
        _selectPraseodymiumAtom();

        return true;
    }

    if (label.compare("ND") == 0)
    {
        _selectNeodymiumAtom();

        return true;
    }

    if (label.compare("PM") == 0)
    {
        _selectPromethiumAtom();

        return true;
    }

    if (label.compare("SM") == 0)
    {
        _selectSamariumAtom();

        return true;
    }

    if (label.compare("EU") == 0)
    {
        _selectEuropiumAtom();

        return true;
    }

    if (label.compare("GD") == 0)
    {
        _selectGadoliniumAtom();

        return true;
    }

    if (label.compare("TB") == 0)
    {
        _selectTerbiumAtom();

        return true;
    }

    if (label.compare("DY") == 0)
    {
        _selectDysprosiumAtom();

        return true;
    }

    if (label.compare("HO") == 0)
    {
        _selectHolmiumAtom();

        return true;
    }

    if (label.compare("ER") == 0)
    {
        _selectErbiumAtom();

        return true;
    }

    if (label.compare("TM") == 0)
    {
        _selectThuliumAtom();

        return true;
    }

    if (label.compare("YB") == 0)
    {
        _selectYtterbiumAtom();

        return true;
    }

    if (label.compare("LU") == 0)
    {
        _selectLutheniumAtom();

        return true;
    }

    if (label.compare("HF") == 0)
    {
        _selectHafniumAtom();

        return true;
    }

    if (label.compare("TA") == 0)
    {
        _selectTantalumAtom();

        return true;
    }

    if (label.compare("W") == 0)
    {
        _selectTungstenAtom();

        return true;
    }

    if (label.compare("RE") == 0)
    {
        _selectRheniumAtom();

        return true;
    }

    if (label.compare("OS") == 0)
    {
        _selectOsmiumAtom();

        return true;
    }

    if (label.compare("IR") == 0)
    {
        _selectIridiumAtom();

        return true;
    }

    if (label.compare("PT") == 0)
    {
        _selectPlatinumAtom();

        return true;
    }

    if (label.compare("AU") == 0)
    {
        _selectGoldAtom();

        return true;
    }

    if (label.compare("HG") == 0)
    {
        _selectMercuryAtom();

        return true;
    }

    if (label.compare("TL") == 0)
    {
        _selectThalliumAtom();

        return true;
    }

    if (label.compare("PB") == 0)
    {
        _selectLeadAtom();

        return true;
    }

    if (label.compare("BI") == 0)
    {
        _selectBismuthAtom();

        return true;
    }

    if (label.compare("PO") == 0)
    {
        _selectPoloniumAtom();

        return true;
    }

    if (label.compare("AT") == 0)
    {
        _selectAstatineAtom();

        return true;
    }

    if (label.compare("RN") == 0)
    {
        _selectRadonAtom();

        return true;
    }

    return false;
}

auto
CChemicalElement::setAtomType(const int64_t identifier) -> bool
{
    bool fg = true;

    switch (identifier)
    {
        case 0:

            _selectDummyAtom();

            break;
            
        case 1:

            _selectHydrogenAtom();

            break;

        case 2:

            _selectHeliumAtom();

            break;

        case 3:

            _selectLithiumAtom();

            break;

        case 4:

            _selectBerylliumAtom();

            break;

        case 5:

            _selectBoronAtom();

            break;

        case 6:

            _selectCarbonAtom();

            break;

        case 7:

            _selectNitrogenAtom();

            break;

        case 8:

            _selectOxygenAtom();

            break;

        case 9:

            _selectFlourineAtom();

            break;

        case 10:

            _selectNeonAtom();

            break;

        case 11:

            _selectSodiumAtom();

            break;

        case 12:

            _selectMagnesiumAtom();

            break;

        case 13:

            _selectAluminiumAtom();

            break;

        case 14:

            _selectSiliconAtom();

            break;

        case 15:

            _selectPhosphorusAtom();

            break;

        case 16:

            _selectSulfurAtom();

            break;

        case 17:

            _selectChlorineAtom();

            break;

        case 18:

            _selectArgonAtom();

            break;

        case 19:

            _selectPotasiumAtom();

            break;

        case 20:

            _selectCalciumAtom();

            break;

        case 21:

            _selectScandiumAtom();

            break;

        case 22:

            _selectTitaniumAtom();

            break;

        case 23:

            _selectVanadiumAtom();

            break;

        case 24:

            _selectChromiumAtom();

            break;

        case 25:

            _selectManganeseAtom();

            break;

        case 26:

            _selectIronAtom();

            break;

        case 27:

            _selectCobaltAtom();

            break;

        case 28:

            _selectNickelAtom();

            break;

        case 29:

            _selectCopperAtom();

            break;

        case 30:

            _selectZincAtom();

            break;

        case 31:

            _selectGaliumAtom();

            break;

        case 32:

            _selectGermaniumAtom();

            break;

        case 33:

            _selectArsenicAtom();

            break;

        case 34:

            _selectSeleniumAtom();

            break;

        case 35:

            _selectBromineAtom();

            break;

        case 36:

            _selectKryptonAtom();

            break;

        case 37:

            _selectRubidiumAtom();

            break;

        case 38:

            _selectStrontiumAtom();

            break;

        case 39:

            _selectYttriumAtom();

            break;

        case 40:

            _selectZirconiumAtom();

            break;

        case 41:

            _selectNiobiumAtom();

            break;

        case 42:

            _selectMolybdenumAtom();

            break;

        case 43:

            _selectTechnetiumAtom();

            break;

        case 44:

            _selectRutheniumAtom();

            break;

        case 45:

            _selectRhodiumAtom();

            break;

        case 46:

            _selectPaladiumAtom();

            break;

        case 47:

            _selectSilverAtom();

            break;

        case 48:

            _selectCadmiumAtom();

            break;

        case 49:

            _selectIndiumAtom();

            break;

        case 50:

            _selectTinAtom();

            break;

        case 51:

            _selectAntimonyAtom();

            break;

        case 52:

            _selectTelluriumAtom();

            break;

        case 53:

            _selectIodineAtom();

            break;

        case 54:

            _selectXenonAtom();

            break;

        case 55:

            _selectCesiumAtom();

            break;

        case 56:

            _selectBariumAtom();

            break;

        case 57:

            _selectLanthanumAtom();

            break;

        case 58:

            _selectCeriumAtom();

            break;

        case 59:

            _selectPraseodymiumAtom();

            break;

        case 60:

            _selectNeodymiumAtom();

            break;

        case 61:

            _selectPromethiumAtom();

            break;

        case 62:

            _selectSamariumAtom();

            break;

        case 63:

            _selectEuropiumAtom();

            break;

        case 64:

            _selectGadoliniumAtom();

            break;

        case 65:

            _selectTerbiumAtom();

            break;

        case 66:

            _selectDysprosiumAtom();

            break;

        case 67:

            _selectHolmiumAtom();

            break;

        case 68:

            _selectErbiumAtom();

            break;

        case 69:

            _selectThuliumAtom();

            break;

        case 70:

            _selectYtterbiumAtom();

            break;

        case 71:

            _selectLutheniumAtom();

            break;

        case 72:

            _selectHafniumAtom();

            break;

        case 73:

            _selectTantalumAtom();

            break;

        case 74:

            _selectTungstenAtom();

            break;

        case 75:

            _selectRheniumAtom();

            break;

        case 76:

            _selectOsmiumAtom();

            break;

        case 77:

            _selectIridiumAtom();

            break;

        case 78:

            _selectPlatinumAtom();

            break;

        case 79:

            _selectGoldAtom();

            break;

        case 80:

            _selectMercuryAtom();

            break;

        case 81:

            _selectThalliumAtom();

            break;

        case 82:

            _selectLeadAtom();

            break;

        case 83:

            _selectBismuthAtom();

            break;

        case 84:

            _selectPoloniumAtom();

            break;

        case 85:

            _selectAstatineAtom();

            break;

        case 86:

            _selectRadonAtom();

            break;

        default:

            fg = false;

            break;
    }

    return fg;
}

auto
CChemicalElement::setIsotope(const int64_t label) -> bool
{
    if (_number < 0) return false;

    if (_number == 0) return true;

    if (_number == 1)
    {
        bool flg = _selectHydrogenIsotopeMass(label);

        return flg;
    }

    if (_number == 2)
    {
        bool flg = _selectHeliumIsotopeMass(label);

        return flg;
    }

    if (_number == 3)
    {
        bool flg = _selectLithiumIsotopeMass(label);

        return flg;
    }

    if (_number == 4)
    {
        bool flg = _selectBerylliumIsotopeMass(label);

        return flg;
    }

    if (_number == 5)
    {
        bool flg = _selectBoronIsotopeMass(label);

        return flg;
    }

    if (_number == 6)
    {
        bool flg = _selectCarbonIsotopeMass(label);

        return flg;
    }

    if (_number == 7)
    {
        bool flg = _selectNitrogenIsotopeMass(label);

        return flg;
    }

    if (_number == 8)
    {
        bool flg = _selectOxygenIsotopeMass(label);

        return flg;
    }

    if (_number == 9)
    {
        bool flg = _selectFlourineIsotopeMass(label);

        return flg;
    }

    if (_number == 10)
    {
        bool flg = _selectNeonIsotopeMass(label);

        return flg;
    }

    if (_number == 11)
    {
        bool flg = _selectSodiumIsotopeMass(label);

        return flg;
    }

    if (_number == 12)
    {
        bool flg = _selectMagnesiumIsotopeMass(label);

        return flg;
    }

    if (_number == 13)
    {
        bool flg = _selectAluminiumIsotopeMass(label);

        return flg;
    }

    if (_number == 14)
    {
        bool flg = _selectSiliconIsotopeMass(label);

        return flg;
    }

    if (_number == 15)
    {
        bool flg = _selectPhosphorusIsotopeMass(label);

        return flg;
    }

    if (_number == 16)
    {
        bool flg = _selectSulfurIsotopeMass(label);

        return flg;
    }

    if (_number == 17)
    {
        bool flg = _selectChlorineIsotopeMass(label);

        return flg;
    }

    if (_number == 18)
    {
        bool flg = _selectArgonIsotopeMass(label);

        return flg;
    }

    if (_number == 19)
    {
        bool flg = _selectPotasiumIsotopeMass(label);

        return flg;
    }

    if (_number == 20)
    {
        bool flg = _selectCalciumIsotopeMass(label);

        return flg;
    }

    if (_number == 21)
    {
        bool flg = _selectScandiumIsotopeMass(label);

        return flg;
    }

    if (_number == 22)
    {
        bool flg = _selectTitaniumIsotopeMass(label);

        return flg;
    }

    if (_number == 23)
    {
        bool flg = _selectVanadiumIsotopeMass(label);

        return flg;
    }

    if (_number == 24)
    {
        bool flg = _selectChromiumIsotopeMass(label);

        return flg;
    }

    if (_number == 25)
    {
        bool flg = _selectManganeseIsotopeMass(label);

        return flg;
    }

    if (_number == 26)
    {
        bool flg = _selectIronIsotopeMass(label);

        return flg;
    }

    if (_number == 27)
    {
        bool flg = _selectCobaltIsotopeMass(label);

        return flg;
    }

    if (_number == 28)
    {
        bool flg = _selectNickelIsotopeMass(label);

        return flg;
    }

    if (_number == 29)
    {
        bool flg = _selectCopperIsotopeMass(label);

        return flg;
    }

    if (_number == 30)
    {
        bool flg = _selectZincIsotopeMass(label);

        return flg;
    }

    if (_number == 31)
    {
        bool flg = _selectGalliumIsotopeMass(label);

        return flg;
    }

    if (_number == 32)
    {
        bool flg = _selectGermaniumIsotopeMass(label);

        return flg;
    }

    if (_number == 33)
    {
        bool flg = _selectArsenicIsotopeMass(label);

        return flg;
    }

    if (_number == 34)
    {
        bool flg = _selectSeleniumIsotopeMass(label);

        return flg;
    }

    if (_number == 35)
    {
        bool flg = _selectBromineIsotopeMass(label);

        return flg;
    }

    if (_number == 36)
    {
        bool flg = _selectKryptonIsotopeMass(label);

        return flg;
    }

    if (_number == 37)
    {
        bool flg = _selectRubidiumIsotopeMass(label);

        return flg;
    }

    if (_number == 38)
    {
        bool flg = _selectStrontiumIsotopeMass(label);

        return flg;
    }

    if (_number == 39)
    {
        bool flg = _selectYttriumIsotopeMass(label);

        return flg;
    }

    if (_number == 40)
    {
        bool flg = _selectZirconiumIsotopeMass(label);

        return flg;
    }

    if (_number == 41)
    {
        bool flg = _selectNiobiumIsotopeMass(label);

        return flg;
    }

    if (_number == 42)
    {
        bool flg = _selectMolybdenumIsotopeMass(label);

        return flg;
    }

    if (_number == 43)
    {
        bool flg = _selectTechnetiumIsotopeMass(label);

        return flg;
    }

    if (_number == 44)
    {
        bool flg = _selectRutheniumIsotopeMass(label);

        return flg;
    }

    if (_number == 45)
    {
        bool flg = _selectRhodiumIsotopeMass(label);

        return flg;
    }

    if (_number == 46)
    {
        bool flg = _selectPaladiumIsotopeMass(label);

        return flg;
    }

    if (_number == 47)
    {
        bool flg = _selectSilverIsotopeMass(label);

        return flg;
    }

    if (_number == 48)
    {
        bool flg = _selectCadmiumIsotopeMass(label);

        return flg;
    }

    if (_number == 49)
    {
        bool flg = _selectIndiumIsotopeMass(label);

        return flg;
    }

    if (_number == 50)
    {
        bool flg = _selectTinIsotopeMass(label);

        return flg;
    }

    if (_number == 51)
    {
        bool flg = _selectAntimonyIsotopeMass(label);

        return flg;
    }

    if (_number == 52)
    {
        bool flg = _selectTelluriumIsotopeMass(label);

        return flg;
    }

    if (_number == 53)
    {
        bool flg = _selectIodineIsotopeMass(label);

        return flg;
    }

    if (_number == 54)
    {
        bool flg = _selectXenonIsotopeMass(label);

        return flg;
    }

    if (_number == 55)
    {
        bool flg = _selectCesiumIsotopeMass(label);

        return flg;
    }

    if (_number == 56)
    {
        bool flg = _selectBariumIsotopeMass(label);

        return flg;
    }

    if (_number == 57)
    {
        bool flg = _selectLanthanumIsotopeMass(label);

        return flg;
    }

    if (_number == 58)
    {
        bool flg = _selectCeriumIsotopeMass(label);

        return flg;
    }

    if (_number == 59)
    {
        bool flg = _selectPraseodymiumIsotopeMass(label);

        return flg;
    }

    if (_number == 60)
    {
        bool flg = _selectNeodymiumIsotopeMass(label);

        return flg;
    }

    if (_number == 61)
    {
        bool flg = _selectPromethiumIsotopeMass(label);

        return flg;
    }

    if (_number == 62)
    {
        bool flg = _selectSamariumIsotopeMass(label);

        return flg;
    }

    if (_number == 63)
    {
        bool flg = _selectEuropiumIsotopeMass(label);

        return flg;
    }

    if (_number == 64)
    {
        bool flg = _selectGadoliniumIsotopeMass(label);

        return flg;
    }

    if (_number == 65)
    {
        bool flg = _selectTerbiumIsotopeMass(label);

        return flg;
    }

    if (_number == 66)
    {
        bool flg = _selectDysprosiumIsotopeMass(label);

        return flg;
    }

    if (_number == 67)
    {
        bool flg = _selectHolmiumIsotopeMass(label);

        return flg;
    }

    if (_number == 68)
    {
        bool flg = _selectErbiumIsotopeMass(label);

        return flg;
    }

    if (_number == 69)
    {
        bool flg = _selectThuliumIsotopeMass(label);

        return flg;
    }

    if (_number == 70)
    {
        bool flg = _selectYtterbiumIsotopeMass(label);

        return flg;
    }

    if (_number == 71)
    {
        bool flg = _selectLutheniumIsotopeMass(label);

        return flg;
    }

    if (_number == 72)
    {
        bool flg = _selectHafniumIsotopeMass(label);

        return flg;
    }

    if (_number == 73)
    {
        bool flg = _selectTantalumIsotopeMass(label);

        return flg;
    }

    if (_number == 74)
    {
        bool flg = _selectTungstenIsotopeMass(label);

        return flg;
    }

    if (_number == 75)
    {
        bool flg = _selectRheniumIsotopeMass(label);

        return flg;
    }

    if (_number == 76)
    {
        bool flg = _selectOsmiumIsotopeMass(label);

        return flg;
    }

    if (_number == 77)
    {
        bool flg = _selectIridiumIsotopeMass(label);

        return flg;
    }

    if (_number == 78)
    {
        bool flg = _selectPlatinumIsotopeMass(label);

        return flg;
    }

    if (_number == 79)
    {
        bool flg = _selectGoldIsotopeMass(label);

        return flg;
    }

    if (_number == 80)
    {
        bool flg = _selectMercuryIsotopeMass(label);

        return flg;
    }

    if (_number == 81)
    {
        bool flg = _selectThalliumIsotopeMass(label);

        return flg;
    }

    if (_number == 82)
    {
        bool flg = _selectLeadIsotopeMass(label);

        return flg;
    }

    if (_number == 83)
    {
        bool flg = _selectBismuthIsotopeMass(label);

        return flg;
    }

    if (_number == 84)
    {
        bool flg = _selectPoloniumIsotopeMass(label);

        return flg;
    }

    if (_number == 85)
    {
        bool flg = _selectAstatineIsotopeMass(label);

        return flg;
    }

    if (_number == 86)
    {
        bool flg = _selectRadonIsotopeMass(label);

        return flg;
    }
    return false;
}

auto
CChemicalElement::getName() const -> std::string
{
    return _label;
}

auto
CChemicalElement::getIdentifier() const -> int64_t
{
    return _number;
}

auto
CChemicalElement::getAtomicCharge() const -> double
{
    return _charge;
}

auto
CChemicalElement::getAtomicMass() const -> double
{
    return _mass;
}

auto
CChemicalElement::getMaxAngularMomentum() const -> int64_t
{
    if ((_number > 0) && (_number < 5)) return 0;

    if ((_number > 4) && (_number < 21)) return 1;

    if ((_number > 20) && (_number < 57)) return 2;

    if ((_number > 56) && (_number < 87)) return 3;

    return -1;
}

auto
CChemicalElement::getMaxIdentifier() const -> int64_t
{
    return 86; 
}

auto
CChemicalElement::_selectDummyAtom() -> void
{
    _label.assign("Bq");

    _number = 0;

    _charge = 0.0;

    _mass = 0.0;
}

auto
CChemicalElement::_selectHydrogenAtom() -> void
{
    _label.assign("H");

    _number = 1;

    _charge = 1.0;

    _mass = 1.007825;
}

auto
CChemicalElement::_selectHeliumAtom() -> void
{
    _label.assign("He");

    _number = 2;

    _charge = 2.0;

    _mass = 4.002603;
}

auto
CChemicalElement::_selectLithiumAtom() -> void
{
    _label.assign("Li");

    _number = 3;

    _charge = 3.0;

    _mass = 7.016005;
}

auto
CChemicalElement::_selectBerylliumAtom() -> void
{
    _label.assign("Be");

    _number = 4;

    _charge = 4.0;

    _mass = 9.012182;
}

auto
CChemicalElement::_selectBoronAtom() -> void
{
    _label.assign("B");

    _number = 5;

    _charge = 5.0;

    _mass = 11.009305;
}

auto
CChemicalElement::_selectCarbonAtom() -> void
{
    _label.assign("C");

    _number = 6;

    _charge = 6.0;

    _mass = 12.000000;
}

auto
CChemicalElement::_selectNitrogenAtom() -> void
{
    _label.assign("N");

    _number = 7;

    _charge = 7.0;

    _mass = 14.003074;
}

auto
CChemicalElement::_selectOxygenAtom() -> void
{
    _label.assign("O");

    _number = 8;

    _charge = 8.0;

    _mass = 15.994915;
}

auto
CChemicalElement::_selectFlourineAtom() -> void
{
    _label.assign("F");

    _number = 9;

    _charge = 9.0;

    _mass = 18.998403;
}

auto
CChemicalElement::_selectNeonAtom() -> void
{
    _label.assign("Ne");

    _number = 10;

    _charge = 10.0;

    _mass = 19.992440;
}

auto
CChemicalElement::_selectSodiumAtom() -> void
{
    _label.assign("Na");

    _number = 11;

    _charge = 11.0;

    _mass = 22.989769;
}

auto
CChemicalElement::_selectMagnesiumAtom() -> void
{
    _label.assign("Mg");

    _number = 12;

    _charge = 12.0;

    _mass = 23.985042;
}

auto
CChemicalElement::_selectAluminiumAtom() -> void
{
    _label.assign("Al");

    _number = 13;

    _charge = 13.0;

    _mass = 26.981539;
}

auto
CChemicalElement::_selectSiliconAtom() -> void
{
    _label.assign("Si");

    _number = 14;

    _charge = 14.0;

    _mass = 27.976927;
}

auto
CChemicalElement::_selectPhosphorusAtom() -> void
{
    _label.assign("P");

    _number = 15;

    _charge = 15.0;

    _mass = 30.973762;
}

auto
CChemicalElement::_selectSulfurAtom() -> void
{
    _label.assign("S");

    _number = 16;

    _charge = 16.0;

    _mass = 31.972071;
}

auto
CChemicalElement::_selectChlorineAtom() -> void
{
    _label.assign("Cl");

    _number = 17;

    _charge = 17.0;

    _mass = 34.968853;
}

auto
CChemicalElement::_selectArgonAtom() -> void
{
    _label.assign("Ar");

    _number = 18;

    _charge = 18.0;

    _mass = 39.962383;
}

auto
CChemicalElement::_selectPotasiumAtom() -> void
{
    _label.assign("K");

    _number = 19;

    _charge = 19.0;

    _mass = 38.963707;
}

auto
CChemicalElement::_selectCalciumAtom() -> void
{
    _label.assign("Ca");

    _number = 20;

    _charge = 20.0;

    _mass = 39.962591;
}

auto
CChemicalElement::_selectScandiumAtom() -> void
{
    _label.assign("Sc");

    _number = 21;

    _charge = 21.0;

    _mass = 44.955912;
}

auto
CChemicalElement::_selectTitaniumAtom() -> void
{
    _label.assign("Ti");

    _number = 22;

    _charge = 22.0;

    _mass = 47.947946;
}

auto
CChemicalElement::_selectVanadiumAtom() -> void
{
    _label.assign("V");

    _number = 23;

    _charge = 23.0;

    _mass = 50.943960;
}

auto
CChemicalElement::_selectChromiumAtom() -> void
{
    _label.assign("Cr");

    _number = 24;

    _charge = 24.0;

    _mass = 51.940508;
}

auto
CChemicalElement::_selectManganeseAtom() -> void
{
    _label.assign("Mn");

    _number = 25;

    _charge = 25.0;

    _mass = 54.938045;
}

auto
CChemicalElement::_selectIronAtom() -> void
{
    _label.assign("Fe");

    _number = 26;

    _charge = 26.0;

    _mass = 55.934938;
}

auto
CChemicalElement::_selectCobaltAtom() -> void
{
    _label.assign("Co");

    _number = 27;

    _charge = 27.0;

    _mass = 58.933195;
}

auto
CChemicalElement::_selectNickelAtom() -> void
{
    _label.assign("Ni");

    _number = 28;

    _charge = 28.0;

    _mass = 57.935343;
}

auto
CChemicalElement::_selectCopperAtom() -> void
{
    _label.assign("Cu");

    _number = 29;

    _charge = 29.0;

    _mass = 62.929598;
}

auto
CChemicalElement::_selectZincAtom() -> void
{
    _label.assign("Zn");

    _number = 30;

    _charge = 30.0;

    _mass = 63.929142;
}

auto
CChemicalElement::_selectGaliumAtom() -> void
{
    _label.assign("Ga");

    _number = 31;

    _charge = 31.0;

    _mass = 68.925574;
}

auto
CChemicalElement::_selectGermaniumAtom() -> void
{
    _label.assign("Ge");

    _number = 32;

    _charge = 32.0;

    _mass = 73.921178;
}

auto
CChemicalElement::_selectArsenicAtom() -> void
{
    _label.assign("As");

    _number = 33;

    _charge = 33.0;

    _mass = 74.921597;
}

auto
CChemicalElement::_selectSeleniumAtom() -> void
{
    _label.assign("Se");

    _number = 34;

    _charge = 34.0;

    _mass = 79.916521;
}

auto
CChemicalElement::_selectBromineAtom() -> void
{
    _label.assign("Br");

    _number = 35;

    _charge = 35.0;

    _mass = 78.918337;
}

auto
CChemicalElement::_selectKryptonAtom() -> void
{
    _label.assign("Kr");

    _number = 36;

    _charge = 36.0;

    _mass = 83.911507;
}

auto
CChemicalElement::_selectRubidiumAtom() -> void
{
    _label.assign("Rb");

    _number = 37;

    _charge = 37.0;

    _mass = 84.911790;
}

auto
CChemicalElement::_selectStrontiumAtom() -> void
{
    _label.assign("Sr");

    _number = 38;

    _charge = 38.0;

    _mass = 87.905612;
}

auto
CChemicalElement::_selectYttriumAtom() -> void
{
    _label.assign("Y");

    _number = 39;

    _charge = 39.0;

    _mass = 88.905848;
}

auto
CChemicalElement::_selectZirconiumAtom() -> void
{
    _label.assign("Zr");

    _number = 40;

    _charge = 40.0;

    _mass = 89.904704;
}

auto
CChemicalElement::_selectNiobiumAtom() -> void
{
    _label.assign("Nb");

    _number = 41;

    _charge = 41.0;

    _mass = 92.906378;
}

auto
CChemicalElement::_selectMolybdenumAtom() -> void
{
    _label.assign("Mo");

    _number = 42;

    _charge = 42.0;

    _mass = 97.905408;
}

auto
CChemicalElement::_selectTechnetiumAtom() -> void
{
    _label.assign("Tc");

    _number = 43;

    _charge = 43.0;

    _mass = 97.907216;
}

auto
CChemicalElement::_selectRutheniumAtom() -> void
{
    _label.assign("Ru");

    _number = 44;

    _charge = 44.0;

    _mass = 101.904349;
}

auto
CChemicalElement::_selectRhodiumAtom() -> void
{
    _label.assign("Rh");

    _number = 45;

    _charge = 45.0;

    _mass = 102.905504;
}

auto
CChemicalElement::_selectPaladiumAtom() -> void
{
    _label.assign("Pd");

    _number = 46;

    _charge = 46.0;

    _mass = 105.903486;
}

auto
CChemicalElement::_selectSilverAtom() -> void
{
    _label.assign("Ag");

    _number = 47;

    _charge = 47.0;

    _mass = 106.905097;
}

auto
CChemicalElement::_selectCadmiumAtom() -> void
{
    _label.assign("Cd");

    _number = 48;

    _charge = 48.0;

    _mass = 113.903359;
}

auto
CChemicalElement::_selectIndiumAtom() -> void
{
    _label.assign("In");

    _number = 49;

    _charge = 49.0;

    _mass = 114.903878;
}

auto
CChemicalElement::_selectTinAtom() -> void
{
    _label.assign("Sn");

    _number = 50;

    _charge = 50.0;

    _mass = 119.902195;
}

auto
CChemicalElement::_selectAntimonyAtom() -> void
{
    _label.assign("Sb");

    _number = 51;

    _charge = 51.0;

    _mass = 120.903816;
}

auto
CChemicalElement::_selectTelluriumAtom() -> void
{
    _label.assign("Te");

    _number = 52;

    _charge = 52.0;

    _mass = 129.906224;
}

auto
CChemicalElement::_selectIodineAtom() -> void
{
    _label.assign("I");

    _number = 53;

    _charge = 53.0;

    _mass = 126.904473;
}

auto
CChemicalElement::_selectXenonAtom() -> void
{
    _label.assign("Xe");

    _number = 54;

    _charge = 54.0;

    _mass = 131.904153;
}

auto
CChemicalElement::_selectCesiumAtom() -> void
{
    _label.assign("Cs");

    _number = 55;

    _charge = 55.0;

    _mass = 132.905452;
}

auto
CChemicalElement::_selectBariumAtom() -> void
{
    _label.assign("Ba");

    _number = 56;

    _charge = 56.0;

    _mass = 137.905247;
}

auto
CChemicalElement::_selectLanthanumAtom() -> void
{
    _label.assign("La");

    _number = 57;

    _charge = 57.0;

    _mass = 138.906353;
}

auto
CChemicalElement::_selectCeriumAtom() -> void
{
    _label.assign("Ce");

    _number = 58;

    _charge = 58.0;

    _mass = 139.905439;
}

auto
CChemicalElement::_selectPraseodymiumAtom() -> void
{
    _label.assign("Pr");

    _number = 59;

    _charge = 59.0;

    _mass = 140.907653;
}

auto
CChemicalElement::_selectNeodymiumAtom() -> void
{
    _label.assign("Nd");

    _number = 60;

    _charge = 60.0;

    _mass = 141.907723;
}

auto
CChemicalElement::_selectPromethiumAtom() -> void
{
    _label.assign("Pm");

    _number = 61;

    _charge = 61.0;

    _mass = 146.915139;
}

auto
CChemicalElement::_selectSamariumAtom() -> void
{
    _label.assign("Sm");

    _number = 62;

    _charge = 62.0;

    _mass = 151.919732;
}

auto
CChemicalElement::_selectEuropiumAtom() -> void
{
    _label.assign("Eu");

    _number = 63;

    _charge = 63.0;

    _mass = 152.921230;
}

auto
CChemicalElement::_selectGadoliniumAtom() -> void
{
    _label.assign("Gd");

    _number = 64;

    _charge = 64.0;

    _mass = 157.924104;
}

auto
CChemicalElement::_selectTerbiumAtom() -> void
{
    _label.assign("Tb");

    _number = 65;

    _charge = 65.0;

    _mass = 158.925347;
}

auto
CChemicalElement::_selectDysprosiumAtom() -> void
{
    _label.assign("Dy");

    _number = 66;

    _charge = 66.0;

    _mass = 163.929175;
}

auto
CChemicalElement::_selectHolmiumAtom() -> void
{
    _label.assign("Ho");

    _number = 67;

    _charge = 67.0;

    _mass = 164.930322;
}

auto
CChemicalElement::_selectErbiumAtom() -> void
{
    _label.assign("Er");

    _number = 68;

    _charge = 68.0;

    _mass = 165.930293;
}

auto
CChemicalElement::_selectThuliumAtom() -> void
{
    _label.assign("Tm");

    _number = 69;

    _charge = 69.0;

    _mass = 168.934213;
}

auto
CChemicalElement::_selectYtterbiumAtom() -> void
{
    _label.assign("Yb");

    _number = 70;

    _charge = 70.0;

    _mass = 173.938862;
}

auto
CChemicalElement::_selectLutheniumAtom() -> void
{
    _label.assign("Lu");

    _number = 71;

    _charge = 71.0;

    _mass = 174.940772;
}

auto
CChemicalElement::_selectHafniumAtom() -> void
{
    _label.assign("Hf");

    _number = 72;

    _charge = 72.0;

    _mass = 179.946550;
}

auto
CChemicalElement::_selectTantalumAtom() -> void
{
    _label.assign("Ta");

    _number = 73;

    _charge = 73.0;

    _mass = 180.947996;
}

auto
CChemicalElement::_selectTungstenAtom() -> void
{
    _label.assign("W");

    _number = 74;

    _charge = 74.0;

    _mass = 183.950931;
}

auto
CChemicalElement::_selectRheniumAtom() -> void
{
    _label.assign("Re");

    _number = 75;

    _charge = 75.0;

    _mass = 186.955753;
}

auto
CChemicalElement::_selectOsmiumAtom() -> void
{
    _label.assign("Os");

    _number = 76;

    _charge = 76.0;

    _mass = 191.961481;
}

auto
CChemicalElement::_selectIridiumAtom() -> void
{
    _label.assign("Ir");

    _number = 77;

    _charge = 77.0;

    _mass = 192.962926;
}

auto
CChemicalElement::_selectPlatinumAtom() -> void
{
    _label.assign("Pt");

    _number = 78;

    _charge = 78.0;

    _mass = 194.964791;
}

auto
CChemicalElement::_selectGoldAtom() -> void
{
    _label.assign("Au");

    _number = 79;

    _charge = 79.0;

    _mass = 196.966569;
}

auto
CChemicalElement::_selectMercuryAtom() -> void
{
    _label.assign("Hg");

    _number = 80;

    _charge = 80.0;

    _mass = 201.970643;
}

auto
CChemicalElement::_selectThalliumAtom() -> void
{
    _label.assign("Tl");

    _number = 81;

    _charge = 81.0;

    _mass = 204.974428;
}

auto
CChemicalElement::_selectLeadAtom() -> void
{
    _label.assign("Pb");

    _number = 82;

    _charge = 82.0;

    _mass = 207.976652;
}

auto
CChemicalElement::_selectBismuthAtom() -> void
{
    _label.assign("Bi");

    _number = 83;

    _charge = 83.0;

    _mass = 208.980399;
}

auto
CChemicalElement::_selectPoloniumAtom() -> void
{
    _label.assign("Po");

    _number = 84;

    _charge = 84.0;

    _mass = 208.982430;
}

auto
CChemicalElement::_selectAstatineAtom() -> void
{
    _label.assign("At");

    _number = 85;

    _charge = 85.0;

    _mass = 209.987148;
}

auto
CChemicalElement::_selectRadonAtom() -> void
{
    _label.assign("Rn");

    _number = 86;

    _charge = 86.0;

    _mass = 222.017578;
}

auto
CChemicalElement::_selectHydrogenIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 1.007825;

            return true;

            break;

        case 1:

            _mass = 1.007825;

            return true;

            break;

        case 2:

            _mass = 2.014101;

            return true;

            break;

        case 3:

            _mass = 3.016049;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectHeliumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 4.002603;

            return true;

            break;

        case 3:

            _mass = 3.016029;

            return true;

            break;

        case 4:

            _mass = 4.002603;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectLithiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 7.016005;

            return true;

            break;

        case 6:

            _mass = 6.015123;

            return true;

            break;

        case 7:

            _mass = 7.016005;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectBerylliumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 9.012182;

            return true;

            break;

        case 9:

            _mass = 9.012182;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectBoronIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 11.009305;

            return true;

            break;

        case 10:

            _mass = 10.012937;

            return true;

            break;

        case 11:

            _mass = 11.009305;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCarbonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 12.000000;

            return true;

            break;

        case 12:

            _mass = 12.000000;

            return true;

            break;

        case 13:

            _mass = 13.003355;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectNitrogenIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 14.003074;

            return true;

            break;

        case 14:

            _mass = 14.003074;

            return true;

            break;

        case 15:

            _mass = 15.000109;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectOxygenIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 15.994915;

            return true;

            break;

        case 16:

            _mass = 15.994915;

            return true;

            break;

        case 17:

            _mass = 16.999132;

            return true;

            break;

        case 18:

            _mass = 17.999161;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectFlourineIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 18.998403;

            return true;

            break;

        case 19:

            _mass = 18.998403;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectNeonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 19.992440;

            return true;

            break;

        case 20:

            _mass = 19.992440;

            return true;

            break;

        case 21:

            _mass = 20.993847;

            return true;

            break;

        case 22:

            _mass = 21.991385;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSodiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 22.989769;

            return true;

            break;

        case 23:

            _mass = 22.989769;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectMagnesiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 23.985042;

            return true;

            break;

        case 24:

            _mass = 23.985042;

            return true;

            break;

        case 25:

            _mass = 24.985837;

            return true;

            break;

        case 26:

            _mass = 25.982593;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectAluminiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 26.981539;

            return true;

            break;

        case 27:

            _mass = 26.981539;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSiliconIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 27.976927;

            return true;

            break;

        case 28:

            _mass = 27.976927;

            return true;

            break;

        case 29:

            _mass = 28.976495;

            return true;

            break;

        case 30:

            _mass = 29.973770;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPhosphorusIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 30.973762;

            return true;

            break;

        case 31:

            _mass = 30.973762;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSulfurIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 31.972071;

            return true;

            break;

        case 32:

            _mass = 31.972071;

            return true;

            break;

        case 33:

            _mass = 32.971459;

            return true;

            break;

        case 34:

            _mass = 33.967867;

            return true;

            break;

        case 36:

            _mass = 35.967081;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectChlorineIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 34.968853;

            return true;

            break;

        case 35:

            _mass = 34.968853;

            return true;

            break;

        case 37:

            _mass = 36.965903;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectArgonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 39.962383;

            return true;

            break;

        case 36:

            _mass = 35.967545;

            return true;

            break;

        case 38:

            _mass = 37.962732;

            return true;

            break;

        case 40:

            _mass = 39.962383;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPotasiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 38.963707;

            return true;

            break;

        case 39:

            _mass = 38.963707;

            return true;

            break;

        case 40:

            _mass = 39.963998;

            return true;

            break;

        case 41:

            _mass = 40.961826;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCalciumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 39.962591;

            return true;

            break;

        case 40:

            _mass = 39.962591;

            return true;

            break;

        case 42:

            _mass = 41.958618;

            return true;

            break;

        case 43:

            _mass = 42.958767;

            return true;

            break;

        case 44:

            _mass = 43.955482;

            return true;

            break;

        case 46:

            _mass = 45.953693;

            return true;

            break;

        case 48:

            _mass = 47.952534;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectScandiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 44.955912;

            return true;

            break;

        case 45:

            _mass = 44.955912;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTitaniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 47.947946;

            return true;

            break;

        case 46:

            _mass = 45.952632;

            return true;

            break;

        case 47:

            _mass = 46.951763;

            return true;

            break;

        case 48:

            _mass = 47.947946;

            return true;

            break;

        case 49:

            _mass = 48.947870;

            return true;

            break;

        case 50:

            _mass = 49.944791;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectVanadiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 50.943960;

            return true;

            break;

        case 50:

            _mass = 49.947159;

            return true;

            break;

        case 51:

            _mass = 50.943960;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectChromiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 51.940508;

            return true;

            break;

        case 50:

            _mass = 49.946044;

            return true;

            break;

        case 52:

            _mass = 51.940508;

            return true;

            break;

        case 53:

            _mass = 52.940649;

            return true;

            break;

        case 54:

            _mass = 53.938880;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectManganeseIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 54.938045;

            return true;

            break;

        case 55:

            _mass = 54.938045;

            return true;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectIronIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 55.934938;

            return true;

            break;

        case 54:

            _mass = 53.939611;

            return true;

            break;

        case 56:

            _mass = 55.934938;

            return true;

            break;

        case 57:

            _mass = 56.935394;

            return true;

            break;

        case 58:

            _mass = 57.933276;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCobaltIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 58.933195;

            return true;

            break;

        case 59:

            _mass = 58.933195;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectNickelIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 57.935343;

            return true;

            break;

        case 58:

            _mass = 57.935343;

            return true;

            break;

        case 60:

            _mass = 59.930786;

            return true;

            break;

        case 61:

            _mass = 60.931056;

            return true;

            break;

        case 62:

            _mass = 61.928345;

            return true;

            break;

        case 64:

            _mass = 63.927966;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCopperIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 62.929598;

            return true;

            break;

        case 63:

            _mass = 62.929598;

            return true;

            break;

        case 65:

            _mass = 64.927790;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectZincIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 63.929142;

            return true;

            break;

        case 64:

            _mass = 63.929142;

            return true;

            break;

        case 66:

            _mass = 65.926033;

            return true;

            break;

        case 67:

            _mass = 66.927127;

            return true;

            break;

        case 68:

            _mass = 67.924844;

            return true;

            break;

        case 70:

            _mass = 69.925319;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectGalliumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 68.925574;

            return true;

            break;

        case 69:

            _mass = 68.925574;

            return true;

            break;

        case 71:

            _mass = 70.924701;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectGermaniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 73.921178;

            return true;

            break;

        case 70:

            _mass = 69.924247;

            return true;

            break;

        case 72:

            _mass = 71.922076;

            return true;

            break;

        case 73:

            _mass = 72.923459;

            return true;

            break;

        case 74:

            _mass = 73.921178;

            return true;

            break;

        case 76:

            _mass = 75.921403;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectArsenicIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 74.921597;

            return true;

            break;

        case 75:

            _mass = 74.921597;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSeleniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 79.916521;

            return true;

            break;

        case 74:

            _mass = 73.922476;

            return true;

            break;

        case 76:

            _mass = 75.919214;

            return true;

            break;

        case 77:

            _mass = 76.919914;

            return true;

            break;

        case 78:

            _mass = 77.917309;

            return true;

            break;

        case 80:

            _mass = 79.916521;

            return true;

            break;

        case 82:

            _mass = 81.916699;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectBromineIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 78.918337;

            return true;

            break;

        case 79:

            _mass = 78.918337;

            return true;

            break;

        case 81:

            _mass = 80.916291;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectKryptonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 83.911507;

            return true;

            break;

        case 78:

            _mass = 77.920365;

            return true;

            break;

        case 80:

            _mass = 79.916379;

            return true;

            break;

        case 82:

            _mass = 81.913484;

            return true;

            break;

        case 83:

            _mass = 82.914136;

            return true;

            break;

        case 84:

            _mass = 83.911507;

            return true;

            break;

        case 86:

            _mass = 85.910611;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectRubidiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 84.911790;

            return true;

            break;

        case 85:

            _mass = 84.911790;

            return true;

            break;

        case 87:

            _mass = 86.909181;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectStrontiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 87.905612;

            return true;

            break;

        case 84:

            _mass = 83.913425;

            return true;

            break;

        case 86:

            _mass = 85.909260;

            return true;

            break;

        case 87:

            _mass = 86.908877;

            return true;

            break;

        case 88:

            _mass = 87.905612;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectYttriumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 88.905848;

            return true;

            break;

        case 89:

            _mass = 88.905848;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectZirconiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 89.904704;

            return true;

            break;

        case 90:

            _mass = 89.904704;

            return true;

            break;

        case 91:

            _mass = 90.905646;

            return true;

            break;

        case 92:

            _mass = 91.905041;

            return true;

            break;

        case 94:

            _mass = 93.906315;

            return true;

            break;

        case 96:

            _mass = 95.908273;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectNiobiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 92.906378;

            return true;

            break;

        case 93:

            _mass = 92.906378;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectMolybdenumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 97.905408;

            return true;

            break;

        case 92:

            _mass = 91.906811;

            return true;

            break;

        case 94:

            _mass = 93.905088;

            return true;

            break;

        case 95:

            _mass = 94.905842;

            return true;

            break;

        case 96:

            _mass = 95.904680;

            return true;

            break;

        case 97:

            _mass = 96.906022;

            return true;

            break;

        case 98:

            _mass = 97.905408;

            return true;

            break;

        case 100:

            _mass = 99.907477;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTechnetiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 97.907216;

            return true;

            break;

        case 98:

            _mass = 97.907216;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectRutheniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 101.904349;

            return true;

            break;

        case 96:

            _mass = 95.907598;

            return true;

            break;

        case 98:

            _mass = 97.905287;

            return true;

            break;

        case 99:

            _mass = 98.905939;

            return true;

            break;

        case 100:

            _mass = 99.904220;

            return true;

            break;

        case 101:

            _mass = 100.905582;

            return true;

            break;

        case 102:

            _mass = 101.904349;

            return true;

            break;

        case 104:

            _mass = 103.905433;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectRhodiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 102.905504;

            return true;

            break;

        case 103:

            _mass = 102.905504;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPaladiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 105.903486;

            return true;

            break;

        case 102:

            _mass = 101.905609;

            return true;

            break;

        case 104:

            _mass = 103.904036;

            return true;

            break;

        case 105:

            _mass = 104.905085;

            return true;

            break;

        case 106:

            _mass = 105.903486;

            return true;

            break;

        case 108:

            _mass = 107.903892;

            return true;

            break;

        case 110:

            _mass = 109.905153;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSilverIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 106.905097;

            return true;

            break;

        case 107:

            _mass = 106.905097;

            return true;

            break;

        case 109:

            _mass = 108.904752;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCadmiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 113.903359;

            return true;

            break;

        case 106:

            _mass = 105.906459;

            return true;

            break;

        case 108:

            _mass = 107.904184;

            return true;

            break;

        case 110:

            _mass = 109.903002;

            return true;

            break;

        case 111:

            _mass = 110.904178;

            return true;

            break;

        case 112:

            _mass = 111.902758;

            return true;

            break;

        case 113:

            _mass = 112.904402;

            return true;

            break;

        case 114:

            _mass = 113.903359;

            return true;

            break;

        case 116:

            _mass = 115.904756;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectIndiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 114.903878;

            return true;

            break;

        case 113:

            _mass = 112.904058;

            return true;

            break;

        case 115:

            _mass = 114.903878;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTinIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 119.902195;

            return true;

            break;

        case 112:

            _mass = 111.904818;

            return true;

            break;

        case 114:

            _mass = 113.902779;

            return true;

            break;

        case 115:

            _mass = 114.903342;

            return true;

            break;

        case 116:

            _mass = 115.901741;

            return true;

            break;

        case 117:

            _mass = 116.902952;

            return true;

            break;

        case 118:

            _mass = 117.901603;

            return true;

            break;

        case 119:

            _mass = 118.903308;

            return true;

            break;

        case 120:

            _mass = 119.902195;

            return true;

            break;

        case 122:

            _mass = 121.903439;

            return true;

            break;

        case 124:

            _mass = 123.905274;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectAntimonyIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 120.903816;

            return true;

            break;

        case 121:

            _mass = 120.903816;

            return true;

            break;

        case 123:

            _mass = 122.904214;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTelluriumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 129.906224;

            return true;

            break;

        case 120:

            _mass = 119.904020;

            return true;

            break;

        case 122:

            _mass = 121.903044;

            return true;

            break;

        case 123:

            _mass = 122.904270;

            return true;

            break;

        case 124:

            _mass = 123.902818;

            return true;

            break;

        case 125:

            _mass = 124.904431;

            return true;

            break;

        case 126:

            _mass = 125.903312;

            return true;

            break;

        case 128:

            _mass = 127.904463;

            return true;

            break;

        case 130:

            _mass = 129.906224;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectIodineIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 126.904473;

            return true;

            break;

        case 127:

            _mass = 126.904473;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectXenonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 131.904153;

            return true;

            break;

        case 124:

            _mass = 123.905893;

            return true;

            break;

        case 126:

            _mass = 125.904274;

            return true;

            break;

        case 128:

            _mass = 127.903531;

            return true;

            break;

        case 129:

            _mass = 128.904779;

            return true;

            break;

        case 130:

            _mass = 129.903508;

            return true;

            break;

        case 131:

            _mass = 130.905082;

            return true;

            break;

        case 132:

            _mass = 131.904153;

            return true;

            break;

        case 134:

            _mass = 133.905395;

            return true;

            break;

        case 136:

            _mass = 135.907219;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCesiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 132.905452;

            return true;

            break;

        case 133:

            _mass = 132.905452;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectBariumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 137.905247;

            return true;

            break;

        case 130:

            _mass = 129.906321;

            return true;

            break;

        case 132:

            _mass = 131.905061;

            return true;

            break;

        case 134:

            _mass = 133.904508;

            return true;

            break;

        case 135:

            _mass = 134.905689;

            return true;

            break;

        case 136:

            _mass = 135.904576;

            return true;

            break;

        case 137:

            _mass = 136.905827;

            return true;

            break;

        case 138:

            _mass = 137.905247;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectLanthanumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 138.906353;

            return true;

            break;

        case 138:

            _mass = 137.907112;

            return true;

            break;

        case 139:

            _mass = 138.906353;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectCeriumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 139.905439;

            return true;

            break;

        case 136:

            _mass = 135.907172;

            return true;

            break;

        case 138:

            _mass = 137.905991;

            return true;

            break;

        case 140:

            _mass = 139.905439;

            return true;

            break;

        case 142:

            _mass = 141.909244;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPraseodymiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 140.907653;

            return true;

            break;

        case 141:

            _mass = 140.907653;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectNeodymiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 141.907723;

            return true;

            break;

        case 142:

            _mass = 141.907723;

            return true;

            break;

        case 143:

            _mass = 142.909814;

            return true;

            break;

        case 144:

            _mass = 143.910087;

            return true;

            break;

        case 145:

            _mass = 144.912574;

            return true;

            break;

        case 146:

            _mass = 145.913117;

            return true;

            break;

        case 148:

            _mass = 147.916893;

            return true;

            break;

        case 150:

            _mass = 149.920891;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPromethiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 146.915139;

            return true;

            break;

        case 147:

            _mass = 146.915139;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectSamariumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 151.919732;

            return true;

            break;

        case 144:

            _mass = 143.911999;

            return true;

            break;

        case 147:

            _mass = 146.914898;

            return true;

            break;

        case 148:

            _mass = 147.914823;

            return true;

            break;

        case 149:

            _mass = 148.917185;

            return true;

            break;

        case 150:

            _mass = 149.917276;

            return true;

            break;

        case 152:

            _mass = 151.919732;

            return true;

            break;

        case 154:

            _mass = 153.922209;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectEuropiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 152.921230;

            return true;

            break;

        case 151:

            _mass = 150.919850;

            return true;

            break;

        case 153:

            _mass = 152.921230;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectGadoliniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 157.924104;

            return true;

            break;

        case 152:

            _mass = 151.919791;

            return true;

            break;

        case 154:

            _mass = 153.920866;

            return true;

            break;

        case 155:

            _mass = 154.922622;

            return true;

            break;

        case 156:

            _mass = 155.922123;

            return true;

            break;

        case 157:

            _mass = 156.923960;

            return true;

            break;

        case 158:

            _mass = 157.924104;

            return true;

            break;

        case 160:

            _mass = 159.927054;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTerbiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 158.925347;

            return true;

            break;

        case 159:

            _mass = 158.925347;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectDysprosiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 163.929175;

            return true;

            break;

        case 156:

            _mass = 155.924283;

            return true;

            break;

        case 158:

            _mass = 157.924409;

            return true;

            break;

        case 160:

            _mass = 159.925198;

            return true;

            break;

        case 161:

            _mass = 160.926933;

            return true;

            break;

        case 162:

            _mass = 161.926798;

            return true;

            break;

        case 163:

            _mass = 162.928731;

            return true;

            break;

        case 164:

            _mass = 163.929175;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectHolmiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 164.930322;

            return true;

            break;

        case 165:

            _mass = 164.930322;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectErbiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 165.930293;

            return true;

            break;

        case 162:

            _mass = 161.928778;

            return true;

            break;

        case 164:

            _mass = 163.929200;

            return true;

            break;

        case 166:

            _mass = 165.930293;

            return true;

            break;

        case 167:

            _mass = 166.932048;

            return true;

            break;

        case 168:

            _mass = 167.932370;

            return true;

            break;

        case 170:

            _mass = 169.935464;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectThuliumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 168.934213;

            return true;

            break;

        case 169:

            _mass = 168.934213;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectYtterbiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 173.938862;

            return true;

            break;

        case 168:

            _mass = 167.933897;

            return true;

            break;

        case 170:

            _mass = 169.934762;

            return true;

            break;

        case 171:

            _mass = 170.936326;

            return true;

            break;

        case 172:

            _mass = 171.936382;

            return true;

            break;

        case 173:

            _mass = 172.938211;

            return true;

            break;

        case 174:

            _mass = 173.938862;

            return true;

            break;

        case 176:

            _mass = 175.942572;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectLutheniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 174.940772;

            return true;

            break;

        case 175:

            _mass = 174.940772;

            return true;

            break;

        case 176:

            _mass = 175.942686;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectHafniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 179.946550;

            return true;

            break;

        case 174:

            _mass = 173.940046;

            return true;

            break;

        case 176:

            _mass = 175.941409;

            return true;

            break;

        case 177:

            _mass = 176.943221;

            return true;

            break;

        case 178:

            _mass = 177.943699;

            return true;

            break;

        case 179:

            _mass = 178.945816;

            return true;

            break;

        case 180:

            _mass = 179.946550;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTantalumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 180.947996;

            return true;

            break;

        case 180:

            _mass = 179.947465;

            return true;

            break;

        case 181:

            _mass = 180.947996;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectTungstenIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 183.950931;

            return true;

            break;

        case 180:

            _mass = 179.946704;

            return true;

            break;

        case 182:

            _mass = 181.948204;

            return true;

            break;

        case 183:

            _mass = 182.950223;

            return true;

            break;

        case 184:

            _mass = 183.950931;

            return true;

            break;

        case 186:

            _mass = 185.954364;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectRheniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 186.955753;

            return true;

            break;

        case 185:

            _mass = 184.952955;

            return true;

            break;

        case 187:

            _mass = 186.955753;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectOsmiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 191.961481;

            return true;

            break;

        case 184:

            _mass = 183.952489;

            return true;

            break;

        case 186:

            _mass = 185.953838;

            return true;

            break;

        case 187:

            _mass = 186.955751;

            return true;

            break;

        case 188:

            _mass = 187.955838;

            return true;

            break;

        case 189:

            _mass = 188.958148;

            return true;

            break;

        case 190:

            _mass = 189.958447;

            return true;

            break;

        case 192:

            _mass = 191.961481;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectIridiumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 192.962926;

            return true;

            break;

        case 191:

            _mass = 190.960594;

            return true;

            break;

        case 193:

            _mass = 192.962926;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPlatinumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 194.964791;

            return true;

            break;

        case 190:

            _mass = 189.959932;

            return true;

            break;

        case 192:

            _mass = 191.961038;

            return true;

            break;

        case 194:

            _mass = 193.962680;

            return true;

            break;

        case 195:

            _mass = 194.964791;

            return true;

            break;

        case 196:

            _mass = 195.964952;

            return true;

            break;

        case 198:

            _mass = 197.967893;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectGoldIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 196.966569;

            return true;

            break;

        case 197:

            _mass = 196.966569;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectMercuryIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 201.970643;

            return true;

            break;

        case 196:

            _mass = 195.965833;

            return true;

            break;

        case 198:

            _mass = 197.966769;

            return true;

            break;

        case 199:

            _mass = 198.968280;

            return true;

            break;

        case 200:

            _mass = 199.968326;

            return true;

            break;

        case 201:

            _mass = 200.970302;

            return true;

            break;

        case 202:

            _mass = 201.970643;

            return true;

            break;

        case 204:

            _mass = 203.973494;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectThalliumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 204.974428;

            return true;

            break;

        case 203:

            _mass = 202.972344;

            return true;

            break;

        case 205:

            _mass = 204.974428;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectLeadIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 207.976652;

            return true;

            break;

        case 204:

            _mass = 203.973044;

            return true;

            break;

        case 206:

            _mass = 205.974465;

            return true;

            break;

        case 207:

            _mass = 206.975897;

            return true;

            break;

        case 208:

            _mass = 207.976652;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectBismuthIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 208.980399;

            return true;

            break;

        case 209:

            _mass = 208.980399;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectPoloniumIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 208.982430;

            return true;

            break;

        case 209:

            _mass = 208.982430;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectAstatineIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 209.987148;

            return true;

            break;

        case 210:

            _mass = 209.987148;

            return true;

            break;

        default:

            return false;

            break;
    }
}

auto
CChemicalElement::_selectRadonIsotopeMass(const int64_t label) -> bool
{
    switch (label)
    {
        case 0:

            _mass = 222.017578;

            return true;

            break;

        case 222:

            _mass = 222.017578;

            return true;

            break;

        default:

            return false;

            break;
    }
}

