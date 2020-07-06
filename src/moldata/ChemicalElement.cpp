//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ChemicalElement.hpp"

#include <cmath>

CChemicalElement::CChemicalElement()

    : _atomicCharge(0.0)

    , _atomicMass(0.0)

    , _atomicNumber(-1)
{
}

CChemicalElement::CChemicalElement(const std::string& atomicLabel, const double atomicCharge, const double atomicMass, const int32_t atomicNumber)

    : _atomicLabel(atomicLabel)

    , _atomicCharge(atomicCharge)

    , _atomicMass(atomicMass)

    , _atomicNumber(atomicNumber)
{
}

CChemicalElement::~CChemicalElement()
{
}

bool
CChemicalElement::operator==(const CChemicalElement& other) const
{
    if (_atomicLabel != other._atomicLabel) return false;

    if (std::fabs(_atomicCharge - other._atomicCharge) > 1.0e-13) return false;

    if (std::fabs(_atomicMass - other._atomicMass) > 1.0e-13) return false;

    if (_atomicNumber != other._atomicNumber) return false;

    return true;
}

bool
CChemicalElement::operator!=(const CChemicalElement& other) const
{
    return !(*this == other);
}

bool
CChemicalElement::setAtomType(const std::string& atomLabel)
{
    if (atomLabel.compare("BQ") == 0)
    {
        _selectDummyAtom();

        return true;
    }

    if (atomLabel.compare("H") == 0)
    {
        _selectHydrogenAtom();

        return true;
    }

    if (atomLabel.compare("HE") == 0)
    {
        _selectHeliumAtom();

        return true;
    }

    if (atomLabel.compare("LI") == 0)
    {
        _selectLithiumAtom();

        return true;
    }

    if (atomLabel.compare("BE") == 0)
    {
        _selectBerylliumAtom();

        return true;
    }

    if (atomLabel.compare("B") == 0)
    {
        _selectBoronAtom();

        return true;
    }

    if (atomLabel.compare("C") == 0)
    {
        _selectCarbonAtom();

        return true;
    }

    if (atomLabel.compare("N") == 0)
    {
        _selectNitrogenAtom();

        return true;
    }

    if (atomLabel.compare("O") == 0)
    {
        _selectOxygenAtom();

        return true;
    }

    if (atomLabel.compare("F") == 0)
    {
        _selectFlourineAtom();

        return true;
    }

    if (atomLabel.compare("NE") == 0)
    {
        _selectNeonAtom();

        return true;
    }

    if (atomLabel.compare("NA") == 0)
    {
        _selectSodiumAtom();

        return true;
    }

    if (atomLabel.compare("MG") == 0)
    {
        _selectMagnesiumAtom();

        return true;
    }

    if (atomLabel.compare("AL") == 0)
    {
        _selectAluminiumAtom();

        return true;
    }

    if (atomLabel.compare("SI") == 0)
    {
        _selectSiliconAtom();

        return true;
    }

    if (atomLabel.compare("P") == 0)
    {
        _selectPhosphorusAtom();

        return true;
    }

    if (atomLabel.compare("S") == 0)
    {
        _selectSulfurAtom();

        return true;
    }

    if (atomLabel.compare("CL") == 0)
    {
        _selectChlorineAtom();

        return true;
    }

    if (atomLabel.compare("AR") == 0)
    {
        _selectArgonAtom();

        return true;
    }

    if (atomLabel.compare("K") == 0)
    {
        _selectPotasiumAtom();

        return true;
    }

    if (atomLabel.compare("CA") == 0)
    {
        _selectCalciumAtom();

        return true;
    }

    if (atomLabel.compare("SC") == 0)
    {
        _selectScandiumAtom();

        return true;
    }

    if (atomLabel.compare("TI") == 0)
    {
        _selectTitaniumAtom();

        return true;
    }

    if (atomLabel.compare("V") == 0)
    {
        _selectVanadiumAtom();

        return true;
    }

    if (atomLabel.compare("CR") == 0)
    {
        _selectChromiumAtom();

        return true;
    }

    if (atomLabel.compare("MN") == 0)
    {
        _selectManganeseAtom();

        return true;
    }

    if (atomLabel.compare("FE") == 0)
    {
        _selectIronAtom();

        return true;
    }

    if (atomLabel.compare("CO") == 0)
    {
        _selectCobaltAtom();

        return true;
    }

    if (atomLabel.compare("NI") == 0)
    {
        _selectNickelAtom();

        return true;
    }

    if (atomLabel.compare("CU") == 0)
    {
        _selectCopperAtom();

        return true;
    }

    if (atomLabel.compare("ZN") == 0)
    {
        _selectZincAtom();

        return true;
    }

    if (atomLabel.compare("GA") == 0)
    {
        _selectGaliumAtom();

        return true;
    }

    if (atomLabel.compare("GE") == 0)
    {
        _selectGermaniumAtom();

        return true;
    }

    if (atomLabel.compare("AS") == 0)
    {
        _selectArsenicAtom();

        return true;
    }

    if (atomLabel.compare("SE") == 0)
    {
        _selectSeleniumAtom();

        return true;
    }

    if (atomLabel.compare("BR") == 0)
    {
        _selectBromineAtom();

        return true;
    }

    if (atomLabel.compare("KR") == 0)
    {
        _selectKryptonAtom();

        return true;
    }

    if (atomLabel.compare("RB") == 0)
    {
        _selectRubidiumAtom();

        return true;
    }

    if (atomLabel.compare("SR") == 0)
    {
        _selectStrontiumAtom();

        return true;
    }

    if (atomLabel.compare("Y") == 0)
    {
        _selectYttriumAtom();

        return true;
    }

    if (atomLabel.compare("ZR") == 0)
    {
        _selectZirconiumAtom();

        return true;
    }

    if (atomLabel.compare("NB") == 0)
    {
        _selectNiobiumAtom();

        return true;
    }

    if (atomLabel.compare("MO") == 0)
    {
        _selectMolybdenumAtom();
        return true;
    }

    if (atomLabel.compare("TC") == 0)
    {
        _selectTechnetiumAtom();

        return true;
    }

    if (atomLabel.compare("RU") == 0)
    {
        _selectRutheniumAtom();

        return true;
    }

    if (atomLabel.compare("RH") == 0)
    {
        _selectRhodiumAtom();

        return true;
    }

    if (atomLabel.compare("PD") == 0)
    {
        _selectPaladiumAtom();

        return true;
    }

    if (atomLabel.compare("AG") == 0)
    {
        _selectSilverAtom();

        return true;
    }

    if (atomLabel.compare("CD") == 0)
    {
        _selectCadmiumAtom();

        return true;
    }

    if (atomLabel.compare("IN") == 0)
    {
        _selectIndiumAtom();

        return true;
    }

    if (atomLabel.compare("SN") == 0)
    {
        _selectTinAtom();

        return true;
    }

    if (atomLabel.compare("SB") == 0)
    {
        _selectAntimonyAtom();

        return true;
    }

    if (atomLabel.compare("TE") == 0)
    {
        _selectTelluriumAtom();

        return true;
    }

    if (atomLabel.compare("I") == 0)
    {
        _selectIodineAtom();

        return true;
    }

    if (atomLabel.compare("XE") == 0)
    {
        _selectXenonAtom();

        return true;
    }

    if (atomLabel.compare("CS") == 0)
    {
        _selectCesiumAtom();

        return true;
    }

    if (atomLabel.compare("BA") == 0)
    {
        _selectBariumAtom();

        return true;
    }

    if (atomLabel.compare("LA") == 0)
    {
        _selectLanthanumAtom();

        return true;
    }

    if (atomLabel.compare("CE") == 0)
    {
        _selectCeriumAtom();

        return true;
    }

    if (atomLabel.compare("PR") == 0)
    {
        _selectPraseodymiumAtom();

        return true;
    }

    if (atomLabel.compare("ND") == 0)
    {
        _selectNeodymiumAtom();

        return true;
    }

    if (atomLabel.compare("PM") == 0)
    {
        _selectPromethiumAtom();

        return true;
    }

    if (atomLabel.compare("SM") == 0)
    {
        _selectSamariumAtom();

        return true;
    }

    if (atomLabel.compare("EU") == 0)
    {
        _selectEuropiumAtom();

        return true;
    }

    if (atomLabel.compare("GD") == 0)
    {
        _selectGadoliniumAtom();

        return true;
    }

    if (atomLabel.compare("TB") == 0)
    {
        _selectTerbiumAtom();

        return true;
    }

    if (atomLabel.compare("DY") == 0)
    {
        _selectDysprosiumAtom();

        return true;
    }

    if (atomLabel.compare("HO") == 0)
    {
        _selectHolmiumAtom();

        return true;
    }

    if (atomLabel.compare("ER") == 0)
    {
        _selectErbiumAtom();

        return true;
    }

    if (atomLabel.compare("TM") == 0)
    {
        _selectThuliumAtom();

        return true;
    }

    if (atomLabel.compare("YB") == 0)
    {
        _selectYtterbiumAtom();

        return true;
    }

    if (atomLabel.compare("LU") == 0)
    {
        _selectLutheniumAtom();

        return true;
    }

    if (atomLabel.compare("HF") == 0)
    {
        _selectHafniumAtom();

        return true;
    }

    if (atomLabel.compare("TA") == 0)
    {
        _selectTantalumAtom();

        return true;
    }

    if (atomLabel.compare("W") == 0)
    {
        _selectTungstenAtom();

        return true;
    }

    if (atomLabel.compare("RE") == 0)
    {
        _selectRheniumAtom();

        return true;
    }

    if (atomLabel.compare("OS") == 0)
    {
        _selectOsmiumAtom();

        return true;
    }

    if (atomLabel.compare("IR") == 0)
    {
        _selectIridiumAtom();

        return true;
    }

    if (atomLabel.compare("PT") == 0)
    {
        _selectPlatinumAtom();

        return true;
    }

    if (atomLabel.compare("AU") == 0)
    {
        _selectGoldAtom();

        return true;
    }

    if (atomLabel.compare("HG") == 0)
    {
        _selectMercuryAtom();

        return true;
    }

    if (atomLabel.compare("TL") == 0)
    {
        _selectThalliumAtom();

        return true;
    }

    if (atomLabel.compare("PB") == 0)
    {
        _selectLeadAtom();

        return true;
    }

    if (atomLabel.compare("BI") == 0)
    {
        _selectBismuthAtom();

        return true;
    }

    if (atomLabel.compare("PO") == 0)
    {
        _selectPoloniumAtom();

        return true;
    }

    if (atomLabel.compare("AT") == 0)
    {
        _selectAstatineAtom();

        return true;
    }

    if (atomLabel.compare("RN") == 0)
    {
        _selectRadonAtom();

        return true;
    }

    return false;
}

bool
CChemicalElement::setAtomType(const int32_t idElemental)
{
    bool fg = true;

    switch (idElemental)
    {
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

bool
CChemicalElement::setIsotope(const int32_t isotopeLabel)
{
    if (_atomicNumber < 0) return false;

    if (_atomicNumber == 0) return true;

    if (_atomicNumber == 1)
    {
        bool flg = _selectHydrogenIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 2)
    {
        bool flg = _selectHeliumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 3)
    {
        bool flg = _selectLithiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 4)
    {
        bool flg = _selectBerylliumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 5)
    {
        bool flg = _selectBoronIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 6)
    {
        bool flg = _selectCarbonIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 7)
    {
        bool flg = _selectNitrogenIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 8)
    {
        bool flg = _selectOxygenIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 9)
    {
        bool flg = _selectFlourineIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 10)
    {
        bool flg = _selectNeonIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 11)
    {
        bool flg = _selectSodiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 12)
    {
        bool flg = _selectMagnesiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 13)
    {
        bool flg = _selectAluminiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 14)
    {
        bool flg = _selectSiliconIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 15)
    {
        bool flg = _selectPhosphorusIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 16)
    {
        bool flg = _selectSulfurIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 17)
    {
        bool flg = _selectChlorineIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 18)
    {
        bool flg = _selectArgonIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 19)
    {
        bool flg = _selectPotasiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 20)
    {
        bool flg = _selectCalciumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 21)
    {
        bool flg = _selectScandiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 22)
    {
        bool flg = _selectTitaniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 23)
    {
        bool flg = _selectVanadiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 24)
    {
        bool flg = _selectChromiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 25)
    {
        bool flg = _selectManganeseIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 26)
    {
        bool flg = _selectIronIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 27)
    {
        bool flg = _selectCobaltIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 28)
    {
        bool flg = _selectNickelIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 29)
    {
        bool flg = _selectCopperIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 30)
    {
        bool flg = _selectZincIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 31)
    {
        bool flg = _selectGalliumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 32)
    {
        bool flg = _selectGermaniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 33)
    {
        bool flg = _selectArsenicIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 34)
    {
        bool flg = _selectSeleniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 35)
    {
        bool flg = _selectBromineIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 36)
    {
        bool flg = _selectKryptonIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 37)
    {
        bool flg = _selectRubidiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 38)
    {
        bool flg = _selectStrontiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 39)
    {
        bool flg = _selectYttriumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 40)
    {
        bool flg = _selectZirconiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 41)
    {
        bool flg = _selectNiobiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 42)
    {
        bool flg = _selectMolybdenumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 43)
    {
        bool flg = _selectTechnetiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 44)
    {
        bool flg = _selectRutheniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 45)
    {
        bool flg = _selectRhodiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 46)
    {
        bool flg = _selectPaladiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 47)
    {
        bool flg = _selectSilverIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 48)
    {
        bool flg = _selectCadmiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 49)
    {
        bool flg = _selectIndiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 50)
    {
        bool flg = _selectTinIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 51)
    {
        bool flg = _selectAntimonyIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 52)
    {
        bool flg = _selectTelluriumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 53)
    {
        bool flg = _selectIodineIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 54)
    {
        bool flg = _selectXenonIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 55)
    {
        bool flg = _selectCesiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 56)
    {
        bool flg = _selectBariumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 57)
    {
        bool flg = _selectLanthanumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 58)
    {
        bool flg = _selectCeriumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 59)
    {
        bool flg = _selectPraseodymiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 60)
    {
        bool flg = _selectNeodymiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 61)
    {
        bool flg = _selectPromethiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 62)
    {
        bool flg = _selectSamariumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 63)
    {
        bool flg = _selectEuropiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 64)
    {
        bool flg = _selectGadoliniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 65)
    {
        bool flg = _selectTerbiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 66)
    {
        bool flg = _selectDysprosiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 67)
    {
        bool flg = _selectHolmiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 68)
    {
        bool flg = _selectErbiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 69)
    {
        bool flg = _selectThuliumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 70)
    {
        bool flg = _selectYtterbiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 71)
    {
        bool flg = _selectLutheniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 72)
    {
        bool flg = _selectHafniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 73)
    {
        bool flg = _selectTantalumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 74)
    {
        bool flg = _selectTungstenIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 75)
    {
        bool flg = _selectRheniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 76)
    {
        bool flg = _selectOsmiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 77)
    {
        bool flg = _selectIridiumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 78)
    {
        bool flg = _selectPlatinumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 79)
    {
        bool flg = _selectGoldIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 80)
    {
        bool flg = _selectMercuryIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 81)
    {
        bool flg = _selectThalliumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 82)
    {
        bool flg = _selectLeadIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 83)
    {
        bool flg = _selectBismuthIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 84)
    {
        bool flg = _selectPoloniumIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 85)
    {
        bool flg = _selectAstatineIsotopeMass(isotopeLabel);

        return flg;
    }

    if (_atomicNumber == 86)
    {
        bool flg = _selectRadonIsotopeMass(isotopeLabel);

        return flg;
    }
    return false;
}

std::string
CChemicalElement::getName() const
{
    return _atomicLabel;
}

int32_t
CChemicalElement::getIdentifier() const
{
    return _atomicNumber;
}

double
CChemicalElement::getAtomicCharge() const
{
    return _atomicCharge;
}

double
CChemicalElement::getAtomicMass() const
{
    return _atomicMass;
}

int32_t
CChemicalElement::getMaxAngularMomentum() const
{
    if ((_atomicNumber > 0) && (_atomicNumber < 5)) return 0;

    if ((_atomicNumber > 4) && (_atomicNumber < 21)) return 1;

    if ((_atomicNumber > 20) && (_atomicNumber < 57)) return 2;

    if ((_atomicNumber > 56) && (_atomicNumber < 87)) return 3;

    return -1;
}

void
CChemicalElement::_selectDummyAtom()
{
    _atomicLabel.assign("Bq");

    _atomicNumber = 0;

    _atomicCharge = 0.0;

    _atomicMass = 0.0;
}

void
CChemicalElement::_selectHydrogenAtom()
{
    _atomicLabel.assign("H");

    _atomicNumber = 1;

    _atomicCharge = 1.0;

    _atomicMass = 1.007825;
}

void
CChemicalElement::_selectHeliumAtom()
{
    _atomicLabel.assign("He");

    _atomicNumber = 2;

    _atomicCharge = 2.0;

    _atomicMass = 4.002603;
}

void
CChemicalElement::_selectLithiumAtom()
{
    _atomicLabel.assign("Li");

    _atomicNumber = 3;

    _atomicCharge = 3.0;

    _atomicMass = 7.016005;
}

void
CChemicalElement::_selectBerylliumAtom()
{
    _atomicLabel.assign("Be");

    _atomicNumber = 4;

    _atomicCharge = 4.0;

    _atomicMass = 9.012182;
}

void
CChemicalElement::_selectBoronAtom()
{
    _atomicLabel.assign("B");

    _atomicNumber = 5;

    _atomicCharge = 5.0;

    _atomicMass = 11.009305;
}

void
CChemicalElement::_selectCarbonAtom()
{
    _atomicLabel.assign("C");

    _atomicNumber = 6;

    _atomicCharge = 6.0;

    _atomicMass = 12.000000;
}

void
CChemicalElement::_selectNitrogenAtom()
{
    _atomicLabel.assign("N");

    _atomicNumber = 7;

    _atomicCharge = 7.0;

    _atomicMass = 14.003074;
}

void
CChemicalElement::_selectOxygenAtom()
{
    _atomicLabel.assign("O");

    _atomicNumber = 8;

    _atomicCharge = 8.0;

    _atomicMass = 15.994915;
}

void
CChemicalElement::_selectFlourineAtom()
{
    _atomicLabel.assign("F");

    _atomicNumber = 9;

    _atomicCharge = 9.0;

    _atomicMass = 18.998403;
}

void
CChemicalElement::_selectNeonAtom()
{
    _atomicLabel.assign("Ne");

    _atomicNumber = 10;

    _atomicCharge = 10.0;

    _atomicMass = 19.992440;
}

void
CChemicalElement::_selectSodiumAtom()
{
    _atomicLabel.assign("Na");

    _atomicNumber = 11;

    _atomicCharge = 11.0;

    _atomicMass = 22.989769;
}

void
CChemicalElement::_selectMagnesiumAtom()
{
    _atomicLabel.assign("Mg");

    _atomicNumber = 12;

    _atomicCharge = 12.0;

    _atomicMass = 23.985042;
}

void
CChemicalElement::_selectAluminiumAtom()
{
    _atomicLabel.assign("Al");

    _atomicNumber = 13;

    _atomicCharge = 13.0;

    _atomicMass = 26.981539;
}

void
CChemicalElement::_selectSiliconAtom()
{
    _atomicLabel.assign("Si");

    _atomicNumber = 14;

    _atomicCharge = 14.0;

    _atomicMass = 27.976927;
}

void
CChemicalElement::_selectPhosphorusAtom()
{
    _atomicLabel.assign("P");

    _atomicNumber = 15;

    _atomicCharge = 15.0;

    _atomicMass = 30.973762;
}

void
CChemicalElement::_selectSulfurAtom()
{
    _atomicLabel.assign("S");

    _atomicNumber = 16;

    _atomicCharge = 16.0;

    _atomicMass = 31.972071;
}

void
CChemicalElement::_selectChlorineAtom()
{
    _atomicLabel.assign("Cl");

    _atomicNumber = 17;

    _atomicCharge = 17.0;

    _atomicMass = 34.968853;
}

void
CChemicalElement::_selectArgonAtom()
{
    _atomicLabel.assign("Ar");

    _atomicNumber = 18;

    _atomicCharge = 18.0;

    _atomicMass = 39.962383;
}

void
CChemicalElement::_selectPotasiumAtom()
{
    _atomicLabel.assign("K");

    _atomicNumber = 19;

    _atomicCharge = 19.0;

    _atomicMass = 38.963707;
}

void
CChemicalElement::_selectCalciumAtom()
{
    _atomicLabel.assign("Ca");

    _atomicNumber = 20;

    _atomicCharge = 20.0;

    _atomicMass = 39.962591;
}

void
CChemicalElement::_selectScandiumAtom()
{
    _atomicLabel.assign("Sc");

    _atomicNumber = 21;

    _atomicCharge = 21.0;

    _atomicMass = 44.955912;
}

void
CChemicalElement::_selectTitaniumAtom()
{
    _atomicLabel.assign("Ti");

    _atomicNumber = 22;

    _atomicCharge = 22.0;

    _atomicMass = 47.947946;
}

void
CChemicalElement::_selectVanadiumAtom()
{
    _atomicLabel.assign("V");

    _atomicNumber = 23;

    _atomicCharge = 23.0;

    _atomicMass = 50.943960;
}

void
CChemicalElement::_selectChromiumAtom()
{
    _atomicLabel.assign("Cr");

    _atomicNumber = 24;

    _atomicCharge = 24.0;

    _atomicMass = 51.940508;
}

void
CChemicalElement::_selectManganeseAtom()
{
    _atomicLabel.assign("Mn");

    _atomicNumber = 25;

    _atomicCharge = 25.0;

    _atomicMass = 54.938045;
}

void
CChemicalElement::_selectIronAtom()
{
    _atomicLabel.assign("Fe");

    _atomicNumber = 26;

    _atomicCharge = 26.0;

    _atomicMass = 55.934938;
}

void
CChemicalElement::_selectCobaltAtom()
{
    _atomicLabel.assign("Co");

    _atomicNumber = 27;

    _atomicCharge = 27.0;

    _atomicMass = 58.933195;
}

void
CChemicalElement::_selectNickelAtom()
{
    _atomicLabel.assign("Ni");

    _atomicNumber = 28;

    _atomicCharge = 28.0;

    _atomicMass = 57.935343;
}

void
CChemicalElement::_selectCopperAtom()
{
    _atomicLabel.assign("Cu");

    _atomicNumber = 29;

    _atomicCharge = 29.0;

    _atomicMass = 62.929598;
}

void
CChemicalElement::_selectZincAtom()
{
    _atomicLabel.assign("Zn");

    _atomicNumber = 30;

    _atomicCharge = 30.0;

    _atomicMass = 63.929142;
}

void
CChemicalElement::_selectGaliumAtom()
{
    _atomicLabel.assign("Ga");

    _atomicNumber = 31;

    _atomicCharge = 31.0;

    _atomicMass = 68.925574;
}

void
CChemicalElement::_selectGermaniumAtom()
{
    _atomicLabel.assign("Ge");

    _atomicNumber = 32;

    _atomicCharge = 32.0;

    _atomicMass = 73.921178;
}

void
CChemicalElement::_selectArsenicAtom()
{
    _atomicLabel.assign("As");

    _atomicNumber = 33;

    _atomicCharge = 33.0;

    _atomicMass = 74.921597;
}

void
CChemicalElement::_selectSeleniumAtom()
{
    _atomicLabel.assign("Se");

    _atomicNumber = 34;

    _atomicCharge = 34.0;

    _atomicMass = 79.916521;
}

void
CChemicalElement::_selectBromineAtom()
{
    _atomicLabel.assign("Br");

    _atomicNumber = 35;

    _atomicCharge = 35.0;

    _atomicMass = 78.918337;
}

void
CChemicalElement::_selectKryptonAtom()
{
    _atomicLabel.assign("Kr");

    _atomicNumber = 36;

    _atomicCharge = 36.0;

    _atomicMass = 83.911507;
}

void
CChemicalElement::_selectRubidiumAtom()
{
    _atomicLabel.assign("Rb");

    _atomicNumber = 37;

    _atomicCharge = 37.0;

    _atomicMass = 84.911790;
}

void
CChemicalElement::_selectStrontiumAtom()
{
    _atomicLabel.assign("Sr");

    _atomicNumber = 38;

    _atomicCharge = 38.0;

    _atomicMass = 87.905612;
}

void
CChemicalElement::_selectYttriumAtom()
{
    _atomicLabel.assign("Y");

    _atomicNumber = 39;

    _atomicCharge = 39.0;

    _atomicMass = 88.905848;
}

void
CChemicalElement::_selectZirconiumAtom()
{
    _atomicLabel.assign("Zr");

    _atomicNumber = 40;

    _atomicCharge = 40.0;

    _atomicMass = 89.904704;
}

void
CChemicalElement::_selectNiobiumAtom()
{
    _atomicLabel.assign("Nb");

    _atomicNumber = 41;

    _atomicCharge = 41.0;

    _atomicMass = 92.906378;
}

void
CChemicalElement::_selectMolybdenumAtom()
{
    _atomicLabel.assign("Mo");

    _atomicNumber = 42;

    _atomicCharge = 42.0;

    _atomicMass = 97.905408;
}

void
CChemicalElement::_selectTechnetiumAtom()
{
    _atomicLabel.assign("Tc");

    _atomicNumber = 43;

    _atomicCharge = 43.0;

    _atomicMass = 97.907216;
}

void
CChemicalElement::_selectRutheniumAtom()
{
    _atomicLabel.assign("Ru");

    _atomicNumber = 44;

    _atomicCharge = 44.0;

    _atomicMass = 101.904349;
}

void
CChemicalElement::_selectRhodiumAtom()
{
    _atomicLabel.assign("Rh");

    _atomicNumber = 45;

    _atomicCharge = 45.0;

    _atomicMass = 102.905504;
}

void
CChemicalElement::_selectPaladiumAtom()
{
    _atomicLabel.assign("Pd");

    _atomicNumber = 46;

    _atomicCharge = 46.0;

    _atomicMass = 105.903486;
}

void
CChemicalElement::_selectSilverAtom()
{
    _atomicLabel.assign("Ag");

    _atomicNumber = 47;

    _atomicCharge = 47.0;

    _atomicMass = 106.905097;
}

void
CChemicalElement::_selectCadmiumAtom()
{
    _atomicLabel.assign("Cd");

    _atomicNumber = 48;

    _atomicCharge = 48.0;

    _atomicMass = 113.903359;
}

void
CChemicalElement::_selectIndiumAtom()
{
    _atomicLabel.assign("In");

    _atomicNumber = 49;

    _atomicCharge = 49.0;

    _atomicMass = 114.903878;
}

void
CChemicalElement::_selectTinAtom()
{
    _atomicLabel.assign("Sn");

    _atomicNumber = 50;

    _atomicCharge = 50.0;

    _atomicMass = 119.902195;
}

void
CChemicalElement::_selectAntimonyAtom()
{
    _atomicLabel.assign("Sb");

    _atomicNumber = 51;

    _atomicCharge = 51.0;

    _atomicMass = 120.903816;
}

void
CChemicalElement::_selectTelluriumAtom()
{
    _atomicLabel.assign("Te");

    _atomicNumber = 52;

    _atomicCharge = 52.0;

    _atomicMass = 129.906224;
}

void
CChemicalElement::_selectIodineAtom()
{
    _atomicLabel.assign("I");

    _atomicNumber = 53;

    _atomicCharge = 53.0;

    _atomicMass = 126.904473;
}

void
CChemicalElement::_selectXenonAtom()
{
    _atomicLabel.assign("Xe");

    _atomicNumber = 54;

    _atomicCharge = 54.0;

    _atomicMass = 131.904153;
}

void
CChemicalElement::_selectCesiumAtom()
{
    _atomicLabel.assign("Cs");

    _atomicNumber = 55;

    _atomicCharge = 55.0;

    _atomicMass = 132.905452;
}

void
CChemicalElement::_selectBariumAtom()
{
    _atomicLabel.assign("Ba");

    _atomicNumber = 56;

    _atomicCharge = 56.0;

    _atomicMass = 137.905247;
}

void
CChemicalElement::_selectLanthanumAtom()
{
    _atomicLabel.assign("La");

    _atomicNumber = 57;

    _atomicCharge = 57.0;

    _atomicMass = 138.906353;
}

void
CChemicalElement::_selectCeriumAtom()
{
    _atomicLabel.assign("Ce");

    _atomicNumber = 58;

    _atomicCharge = 58.0;

    _atomicMass = 139.905439;
}

void
CChemicalElement::_selectPraseodymiumAtom()
{
    _atomicLabel.assign("Pr");

    _atomicNumber = 59;

    _atomicCharge = 59.0;

    _atomicMass = 140.907653;
}

void
CChemicalElement::_selectNeodymiumAtom()
{
    _atomicLabel.assign("Nd");

    _atomicNumber = 60;

    _atomicCharge = 60.0;

    _atomicMass = 141.907723;
}

void
CChemicalElement::_selectPromethiumAtom()
{
    _atomicLabel.assign("Pm");

    _atomicNumber = 61;

    _atomicCharge = 61.0;

    _atomicMass = 146.915139;
}

void
CChemicalElement::_selectSamariumAtom()
{
    _atomicLabel.assign("Sm");

    _atomicNumber = 62;

    _atomicCharge = 62.0;

    _atomicMass = 151.919732;
}

void
CChemicalElement::_selectEuropiumAtom()
{
    _atomicLabel.assign("Eu");

    _atomicNumber = 63;

    _atomicCharge = 63.0;

    _atomicMass = 152.921230;
}

void
CChemicalElement::_selectGadoliniumAtom()
{
    _atomicLabel.assign("Gd");

    _atomicNumber = 64;

    _atomicCharge = 64.0;

    _atomicMass = 157.924104;
}

void
CChemicalElement::_selectTerbiumAtom()
{
    _atomicLabel.assign("Tb");

    _atomicNumber = 65;

    _atomicCharge = 65.0;

    _atomicMass = 158.925347;
}

void
CChemicalElement::_selectDysprosiumAtom()
{
    _atomicLabel.assign("Dy");

    _atomicNumber = 66;

    _atomicCharge = 66.0;

    _atomicMass = 163.929175;
}

void
CChemicalElement::_selectHolmiumAtom()
{
    _atomicLabel.assign("Ho");

    _atomicNumber = 67;

    _atomicCharge = 67.0;

    _atomicMass = 164.930322;
}

void
CChemicalElement::_selectErbiumAtom()
{
    _atomicLabel.assign("Er");

    _atomicNumber = 68;

    _atomicCharge = 68.0;

    _atomicMass = 165.930293;
}

void
CChemicalElement::_selectThuliumAtom()
{
    _atomicLabel.assign("Tm");

    _atomicNumber = 69;

    _atomicCharge = 69.0;

    _atomicMass = 168.934213;
}

void
CChemicalElement::_selectYtterbiumAtom()
{
    _atomicLabel.assign("Tm");

    _atomicNumber = 70;

    _atomicCharge = 70.0;

    _atomicMass = 173.938862;
}

void
CChemicalElement::_selectLutheniumAtom()
{
    _atomicLabel.assign("Lu");

    _atomicNumber = 71;

    _atomicCharge = 71.0;

    _atomicMass = 174.940772;
}

void
CChemicalElement::_selectHafniumAtom()
{
    _atomicLabel.assign("Hf");

    _atomicNumber = 72;

    _atomicCharge = 72.0;

    _atomicMass = 179.946550;
}

void
CChemicalElement::_selectTantalumAtom()
{
    _atomicLabel.assign("Ta");

    _atomicNumber = 73;

    _atomicCharge = 73.0;

    _atomicMass = 180.947996;
}

void
CChemicalElement::_selectTungstenAtom()
{
    _atomicLabel.assign("W");

    _atomicNumber = 74;

    _atomicCharge = 74.0;

    _atomicMass = 183.950931;
}

void
CChemicalElement::_selectRheniumAtom()
{
    _atomicLabel.assign("Re");

    _atomicNumber = 75;

    _atomicCharge = 75.0;

    _atomicMass = 186.955753;
}

void
CChemicalElement::_selectOsmiumAtom()
{
    _atomicLabel.assign("Os");

    _atomicNumber = 76;

    _atomicCharge = 76.0;

    _atomicMass = 191.961481;
}

void
CChemicalElement::_selectIridiumAtom()
{
    _atomicLabel.assign("Ir");

    _atomicNumber = 77;

    _atomicCharge = 77.0;

    _atomicMass = 192.962926;
}

void
CChemicalElement::_selectPlatinumAtom()
{
    _atomicLabel.assign("Pt");

    _atomicNumber = 78;

    _atomicCharge = 78.0;

    _atomicMass = 194.964791;
}

void
CChemicalElement::_selectGoldAtom()
{
    _atomicLabel.assign("Au");

    _atomicNumber = 79;

    _atomicCharge = 79.0;

    _atomicMass = 196.966569;
}

void
CChemicalElement::_selectMercuryAtom()
{
    _atomicLabel.assign("Hg");

    _atomicNumber = 80;

    _atomicCharge = 80.0;

    _atomicMass = 201.970643;
}

void
CChemicalElement::_selectThalliumAtom()
{
    _atomicLabel.assign("Tl");

    _atomicNumber = 81;

    _atomicCharge = 81.0;

    _atomicMass = 204.974428;
}

void
CChemicalElement::_selectLeadAtom()
{
    _atomicLabel.assign("Pb");

    _atomicNumber = 82;

    _atomicCharge = 82.0;

    _atomicMass = 207.976652;
}

void
CChemicalElement::_selectBismuthAtom()
{
    _atomicLabel.assign("Bi");

    _atomicNumber = 83;

    _atomicCharge = 83.0;

    _atomicMass = 208.980399;
}

void
CChemicalElement::_selectPoloniumAtom()
{
    _atomicLabel.assign("Po");

    _atomicNumber = 84;

    _atomicCharge = 84.0;

    _atomicMass = 208.982430;
}

void
CChemicalElement::_selectAstatineAtom()
{
    _atomicLabel.assign("At");

    _atomicNumber = 85;

    _atomicCharge = 85.0;

    _atomicMass = 209.987148;
}

void
CChemicalElement::_selectRadonAtom()
{
    _atomicLabel.assign("Rn");

    _atomicNumber = 86;

    _atomicCharge = 86.0;

    _atomicMass = 222.017578;
}

bool
CChemicalElement::_selectHydrogenIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 1.007825;

            return true;

            break;

        case 1:

            _atomicMass = 1.007825;

            return true;

            break;

        case 2:

            _atomicMass = 2.014101;

            return true;

            break;

        case 3:

            _atomicMass = 3.016049;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectHeliumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 4.002603;

            return true;

            break;

        case 3:

            _atomicMass = 3.016029;

            return true;

            break;

        case 4:

            _atomicMass = 4.002603;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectLithiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 7.016005;

            return true;

            break;

        case 6:

            _atomicMass = 6.015123;

            return true;

            break;

        case 7:

            _atomicMass = 7.016005;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectBerylliumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 9.012182;

            return true;

            break;

        case 9:

            _atomicMass = 9.012182;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectBoronIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 11.009305;

            return true;

            break;

        case 10:

            _atomicMass = 10.012937;

            return true;

            break;

        case 11:

            _atomicMass = 11.009305;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCarbonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 12.000000;

            return true;

            break;

        case 12:

            _atomicMass = 12.000000;

            return true;

            break;

        case 13:

            _atomicMass = 13.003355;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectNitrogenIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 14.003074;

            return true;

            break;

        case 14:

            _atomicMass = 14.003074;

            return true;

            break;

        case 15:

            _atomicMass = 15.000109;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectOxygenIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 15.994915;

            return true;

            break;

        case 16:

            _atomicMass = 15.994915;

            return true;

            break;

        case 17:

            _atomicMass = 16.999132;

            return true;

            break;

        case 18:

            _atomicMass = 17.999161;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectFlourineIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 18.998403;

            return true;

            break;

        case 19:

            _atomicMass = 18.998403;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectNeonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 19.992440;

            return true;

            break;

        case 20:

            _atomicMass = 19.992440;

            return true;

            break;

        case 21:

            _atomicMass = 20.993847;

            return true;

            break;

        case 22:

            _atomicMass = 21.991385;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSodiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 22.989769;

            return true;

            break;

        case 23:

            _atomicMass = 22.989769;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectMagnesiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 23.985042;

            return true;

            break;

        case 24:

            _atomicMass = 23.985042;

            return true;

            break;

        case 25:

            _atomicMass = 24.985837;

            return true;

            break;

        case 26:

            _atomicMass = 25.982593;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectAluminiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 26.981539;

            return true;

            break;

        case 27:

            _atomicMass = 26.981539;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSiliconIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 27.976927;

            return true;

            break;

        case 28:

            _atomicMass = 27.976927;

            return true;

            break;

        case 29:

            _atomicMass = 28.976495;

            return true;

            break;

        case 30:

            _atomicMass = 29.973770;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPhosphorusIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 30.973762;

            return true;

            break;

        case 31:

            _atomicMass = 30.973762;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSulfurIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 31.972071;

            return true;

            break;

        case 32:

            _atomicMass = 31.972071;

            return true;

            break;

        case 33:

            _atomicMass = 32.971459;

            return true;

            break;

        case 34:

            _atomicMass = 33.967867;

            return true;

            break;

        case 36:

            _atomicMass = 35.967081;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectChlorineIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 34.968853;

            return true;

            break;

        case 35:

            _atomicMass = 34.968853;

            return true;

            break;

        case 37:

            _atomicMass = 36.965903;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectArgonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 39.962383;

            return true;

            break;

        case 36:

            _atomicMass = 35.967545;

            return true;

            break;

        case 38:

            _atomicMass = 37.962732;

            return true;

            break;

        case 40:

            _atomicMass = 39.962383;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPotasiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 38.963707;

            return true;

            break;

        case 39:

            _atomicMass = 38.963707;

            return true;

            break;

        case 40:

            _atomicMass = 39.963998;

            return true;

            break;

        case 41:

            _atomicMass = 40.961826;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCalciumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 39.962591;

            return true;

            break;

        case 40:

            _atomicMass = 39.962591;

            return true;

            break;

        case 42:

            _atomicMass = 41.958618;

            return true;

            break;

        case 43:

            _atomicMass = 42.958767;

            return true;

            break;

        case 44:

            _atomicMass = 43.955482;

            return true;

            break;

        case 46:

            _atomicMass = 45.953693;

            return true;

            break;

        case 48:

            _atomicMass = 47.952534;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectScandiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 44.955912;

            return true;

            break;

        case 45:

            _atomicMass = 44.955912;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTitaniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 47.947946;

            return true;

            break;

        case 46:

            _atomicMass = 45.952632;

            return true;

            break;

        case 47:

            _atomicMass = 46.951763;

            return true;

            break;

        case 48:

            _atomicMass = 47.947946;

            return true;

            break;

        case 49:

            _atomicMass = 48.947870;

            return true;

            break;

        case 50:

            _atomicMass = 49.944791;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectVanadiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 50.943960;

            return true;

            break;

        case 50:

            _atomicMass = 49.947159;

            return true;

            break;

        case 51:

            _atomicMass = 50.943960;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectChromiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 51.940508;

            return true;

            break;

        case 50:

            _atomicMass = 49.946044;

            return true;

            break;

        case 52:

            _atomicMass = 51.940508;

            return true;

            break;

        case 53:

            _atomicMass = 52.940649;

            return true;

            break;

        case 54:

            _atomicMass = 53.938880;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectManganeseIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 54.938045;

            return true;

            break;

        case 55:

            _atomicMass = 54.938045;

            return true;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectIronIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 55.934938;

            return true;

            break;

        case 54:

            _atomicMass = 53.939611;

            return true;

            break;

        case 56:

            _atomicMass = 55.934938;

            return true;

            break;

        case 57:

            _atomicMass = 56.935394;

            return true;

            break;

        case 58:

            _atomicMass = 57.933276;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCobaltIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 58.933195;

            return true;

            break;

        case 59:

            _atomicMass = 58.933195;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectNickelIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 57.935343;

            return true;

            break;

        case 58:

            _atomicMass = 57.935343;

            return true;

            break;

        case 60:

            _atomicMass = 59.930786;

            return true;

            break;

        case 61:

            _atomicMass = 60.931056;

            return true;

            break;

        case 62:

            _atomicMass = 61.928345;

            return true;

            break;

        case 64:

            _atomicMass = 63.927966;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCopperIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 62.929598;

            return true;

            break;

        case 63:

            _atomicMass = 62.929598;

            return true;

            break;

        case 65:

            _atomicMass = 64.927790;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectZincIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 63.929142;

            return true;

            break;

        case 64:

            _atomicMass = 63.929142;

            return true;

            break;

        case 66:

            _atomicMass = 65.926033;

            return true;

            break;

        case 67:

            _atomicMass = 66.927127;

            return true;

            break;

        case 68:

            _atomicMass = 67.924844;

            return true;

            break;

        case 70:

            _atomicMass = 69.925319;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectGalliumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 68.925574;

            return true;

            break;

        case 69:

            _atomicMass = 68.925574;

            return true;

            break;

        case 71:

            _atomicMass = 70.924701;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectGermaniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 73.921178;

            return true;

            break;

        case 70:

            _atomicMass = 69.924247;

            return true;

            break;

        case 72:

            _atomicMass = 71.922076;

            return true;

            break;

        case 73:

            _atomicMass = 72.923459;

            return true;

            break;

        case 74:

            _atomicMass = 73.921178;

            return true;

            break;

        case 76:

            _atomicMass = 75.921403;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectArsenicIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 74.921597;

            return true;

            break;

        case 75:

            _atomicMass = 74.921597;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSeleniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 79.916521;

            return true;

            break;

        case 74:

            _atomicMass = 73.922476;

            return true;

            break;

        case 76:

            _atomicMass = 75.919214;

            return true;

            break;

        case 77:

            _atomicMass = 76.919914;

            return true;

            break;

        case 78:

            _atomicMass = 77.917309;

            return true;

            break;

        case 80:

            _atomicMass = 79.916521;

            return true;

            break;

        case 82:

            _atomicMass = 81.916699;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectBromineIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 78.918337;

            return true;

            break;

        case 79:

            _atomicMass = 78.918337;

            return true;

            break;

        case 81:

            _atomicMass = 80.916291;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectKryptonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 83.911507;

            return true;

            break;

        case 78:

            _atomicMass = 77.920365;

            return true;

            break;

        case 80:

            _atomicMass = 79.916379;

            return true;

            break;

        case 82:

            _atomicMass = 81.913484;

            return true;

            break;

        case 83:

            _atomicMass = 82.914136;

            return true;

            break;

        case 84:

            _atomicMass = 83.911507;

            return true;

            break;

        case 86:

            _atomicMass = 85.910611;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectRubidiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 84.911790;

            return true;

            break;

        case 85:

            _atomicMass = 84.911790;

            return true;

            break;

        case 87:

            _atomicMass = 86.909181;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectStrontiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 87.905612;

            return true;

            break;

        case 84:

            _atomicMass = 83.913425;

            return true;

            break;

        case 86:

            _atomicMass = 85.909260;

            return true;

            break;

        case 87:

            _atomicMass = 86.908877;

            return true;

            break;

        case 88:

            _atomicMass = 87.905612;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectYttriumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 88.905848;

            return true;

            break;

        case 89:

            _atomicMass = 88.905848;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectZirconiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 89.904704;

            return true;

            break;

        case 90:

            _atomicMass = 89.904704;

            return true;

            break;

        case 91:

            _atomicMass = 90.905646;

            return true;

            break;

        case 92:

            _atomicMass = 91.905041;

            return true;

            break;

        case 94:

            _atomicMass = 93.906315;

            return true;

            break;

        case 96:

            _atomicMass = 95.908273;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectNiobiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 92.906378;

            return true;

            break;

        case 93:

            _atomicMass = 92.906378;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectMolybdenumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 97.905408;

            return true;

            break;

        case 92:

            _atomicMass = 91.906811;

            return true;

            break;

        case 94:

            _atomicMass = 93.905088;

            return true;

            break;

        case 95:

            _atomicMass = 94.905842;

            return true;

            break;

        case 96:

            _atomicMass = 95.904680;

            return true;

            break;

        case 97:

            _atomicMass = 96.906022;

            return true;

            break;

        case 98:

            _atomicMass = 97.905408;

            return true;

            break;

        case 100:

            _atomicMass = 99.907477;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTechnetiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 97.907216;

            return true;

            break;

        case 98:

            _atomicMass = 97.907216;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectRutheniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 101.904349;

            return true;

            break;

        case 96:

            _atomicMass = 95.907598;

            return true;

            break;

        case 98:

            _atomicMass = 97.905287;

            return true;

            break;

        case 99:

            _atomicMass = 98.905939;

            return true;

            break;

        case 100:

            _atomicMass = 99.904220;

            return true;

            break;

        case 101:

            _atomicMass = 100.905582;

            return true;

            break;

        case 102:

            _atomicMass = 101.904349;

            return true;

            break;

        case 104:

            _atomicMass = 103.905433;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectRhodiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 102.905504;

            return true;

            break;

        case 103:

            _atomicMass = 102.905504;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPaladiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 105.903486;

            return true;

            break;

        case 102:

            _atomicMass = 101.905609;

            return true;

            break;

        case 104:

            _atomicMass = 103.904036;

            return true;

            break;

        case 105:

            _atomicMass = 104.905085;

            return true;

            break;

        case 106:

            _atomicMass = 105.903486;

            return true;

            break;

        case 108:

            _atomicMass = 107.903892;

            return true;

            break;

        case 110:

            _atomicMass = 109.905153;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSilverIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 106.905097;

            return true;

            break;

        case 107:

            _atomicMass = 106.905097;

            return true;

            break;

        case 109:

            _atomicMass = 108.904752;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCadmiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 113.903359;

            return true;

            break;

        case 106:

            _atomicMass = 105.906459;

            return true;

            break;

        case 108:

            _atomicMass = 107.904184;

            return true;

            break;

        case 110:

            _atomicMass = 109.903002;

            return true;

            break;

        case 111:

            _atomicMass = 110.904178;

            return true;

            break;

        case 112:

            _atomicMass = 111.902758;

            return true;

            break;

        case 113:

            _atomicMass = 112.904402;

            return true;

            break;

        case 114:

            _atomicMass = 113.903359;

            return true;

            break;

        case 116:

            _atomicMass = 115.904756;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectIndiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 114.903878;

            return true;

            break;

        case 113:

            _atomicMass = 112.904058;

            return true;

            break;

        case 115:

            _atomicMass = 114.903878;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTinIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 119.902195;

            return true;

            break;

        case 112:

            _atomicMass = 111.904818;

            return true;

            break;

        case 114:

            _atomicMass = 113.902779;

            return true;

            break;

        case 115:

            _atomicMass = 114.903342;

            return true;

            break;

        case 116:

            _atomicMass = 115.901741;

            return true;

            break;

        case 117:

            _atomicMass = 116.902952;

            return true;

            break;

        case 118:

            _atomicMass = 117.901603;

            return true;

            break;

        case 119:

            _atomicMass = 118.903308;

            return true;

            break;

        case 120:

            _atomicMass = 119.902195;

            return true;

            break;

        case 122:

            _atomicMass = 121.903439;

            return true;

            break;

        case 124:

            _atomicMass = 123.905274;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectAntimonyIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 120.903816;

            return true;

            break;

        case 121:

            _atomicMass = 120.903816;

            return true;

            break;

        case 123:

            _atomicMass = 122.904214;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTelluriumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 129.906224;

            return true;

            break;

        case 120:

            _atomicMass = 119.904020;

            return true;

            break;

        case 122:

            _atomicMass = 121.903044;

            return true;

            break;

        case 123:

            _atomicMass = 122.904270;

            return true;

            break;

        case 124:

            _atomicMass = 123.902818;

            return true;

            break;

        case 125:

            _atomicMass = 124.904431;

            return true;

            break;

        case 126:

            _atomicMass = 125.903312;

            return true;

            break;

        case 128:

            _atomicMass = 127.904463;

            return true;

            break;

        case 130:

            _atomicMass = 129.906224;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectIodineIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 126.904473;

            return true;

            break;

        case 127:

            _atomicMass = 126.904473;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectXenonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 131.904153;

            return true;

            break;

        case 124:

            _atomicMass = 123.905893;

            return true;

            break;

        case 126:

            _atomicMass = 125.904274;

            return true;

            break;

        case 128:

            _atomicMass = 127.903531;

            return true;

            break;

        case 129:

            _atomicMass = 128.904779;

            return true;

            break;

        case 130:

            _atomicMass = 129.903508;

            return true;

            break;

        case 131:

            _atomicMass = 130.905082;

            return true;

            break;

        case 132:

            _atomicMass = 131.904153;

            return true;

            break;

        case 134:

            _atomicMass = 133.905395;

            return true;

            break;

        case 136:

            _atomicMass = 135.907219;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCesiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 132.905452;

            return true;

            break;

        case 133:

            _atomicMass = 132.905452;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectBariumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 137.905247;

            return true;

            break;

        case 130:

            _atomicMass = 129.906321;

            return true;

            break;

        case 132:

            _atomicMass = 131.905061;

            return true;

            break;

        case 134:

            _atomicMass = 133.904508;

            return true;

            break;

        case 135:

            _atomicMass = 134.905689;

            return true;

            break;

        case 136:

            _atomicMass = 135.904576;

            return true;

            break;

        case 137:

            _atomicMass = 136.905827;

            return true;

            break;

        case 138:

            _atomicMass = 137.905247;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectLanthanumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 138.906353;

            return true;

            break;

        case 138:

            _atomicMass = 137.907112;

            return true;

            break;

        case 139:

            _atomicMass = 138.906353;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectCeriumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 139.905439;

            return true;

            break;

        case 136:

            _atomicMass = 135.907172;

            return true;

            break;

        case 138:

            _atomicMass = 137.905991;

            return true;

            break;

        case 140:

            _atomicMass = 139.905439;

            return true;

            break;

        case 142:

            _atomicMass = 141.909244;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPraseodymiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 140.907653;

            return true;

            break;

        case 141:

            _atomicMass = 140.907653;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectNeodymiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 141.907723;

            return true;

            break;

        case 142:

            _atomicMass = 141.907723;

            return true;

            break;

        case 143:

            _atomicMass = 142.909814;

            return true;

            break;

        case 144:

            _atomicMass = 143.910087;

            return true;

            break;

        case 145:

            _atomicMass = 144.912574;

            return true;

            break;

        case 146:

            _atomicMass = 145.913117;

            return true;

            break;

        case 148:

            _atomicMass = 147.916893;

            return true;

            break;

        case 150:

            _atomicMass = 149.920891;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPromethiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 146.915139;

            return true;

            break;

        case 147:

            _atomicMass = 146.915139;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectSamariumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 151.919732;

            return true;

            break;

        case 144:

            _atomicMass = 143.911999;

            return true;

            break;

        case 147:

            _atomicMass = 146.914898;

            return true;

            break;

        case 148:

            _atomicMass = 147.914823;

            return true;

            break;

        case 149:

            _atomicMass = 148.917185;

            return true;

            break;

        case 150:

            _atomicMass = 149.917276;

            return true;

            break;

        case 152:

            _atomicMass = 151.919732;

            return true;

            break;

        case 154:

            _atomicMass = 153.922209;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectEuropiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 152.921230;

            return true;

            break;

        case 151:

            _atomicMass = 150.919850;

            return true;

            break;

        case 153:

            _atomicMass = 152.921230;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectGadoliniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 157.924104;

            return true;

            break;

        case 152:

            _atomicMass = 151.919791;

            return true;

            break;

        case 154:

            _atomicMass = 153.920866;

            return true;

            break;

        case 155:

            _atomicMass = 154.922622;

            return true;

            break;

        case 156:

            _atomicMass = 155.922123;

            return true;

            break;

        case 157:

            _atomicMass = 156.923960;

            return true;

            break;

        case 158:

            _atomicMass = 157.924104;

            return true;

            break;

        case 160:

            _atomicMass = 159.927054;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTerbiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 158.925347;

            return true;

            break;

        case 159:

            _atomicMass = 158.925347;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectDysprosiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 163.929175;

            return true;

            break;

        case 156:

            _atomicMass = 155.924283;

            return true;

            break;

        case 158:

            _atomicMass = 157.924409;

            return true;

            break;

        case 160:

            _atomicMass = 159.925198;

            return true;

            break;

        case 161:

            _atomicMass = 160.926933;

            return true;

            break;

        case 162:

            _atomicMass = 161.926798;

            return true;

            break;

        case 163:

            _atomicMass = 162.928731;

            return true;

            break;

        case 164:

            _atomicMass = 163.929175;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectHolmiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 164.930322;

            return true;

            break;

        case 165:

            _atomicMass = 164.930322;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectErbiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 165.930293;

            return true;

            break;

        case 162:

            _atomicMass = 161.928778;

            return true;

            break;

        case 164:

            _atomicMass = 163.929200;

            return true;

            break;

        case 166:

            _atomicMass = 165.930293;

            return true;

            break;

        case 167:

            _atomicMass = 166.932048;

            return true;

            break;

        case 168:

            _atomicMass = 167.932370;

            return true;

            break;

        case 170:

            _atomicMass = 169.935464;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectThuliumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 168.934213;

            return true;

            break;

        case 169:

            _atomicMass = 168.934213;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectYtterbiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 173.938862;

            return true;

            break;

        case 168:

            _atomicMass = 167.933897;

            return true;

            break;

        case 170:

            _atomicMass = 169.934762;

            return true;

            break;

        case 171:

            _atomicMass = 170.936326;

            return true;

            break;

        case 172:

            _atomicMass = 171.936382;

            return true;

            break;

        case 173:

            _atomicMass = 172.938211;

            return true;

            break;

        case 174:

            _atomicMass = 173.938862;

            return true;

            break;

        case 176:

            _atomicMass = 175.942572;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectLutheniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 174.940772;

            return true;

            break;

        case 175:

            _atomicMass = 174.940772;

            return true;

            break;

        case 176:

            _atomicMass = 175.942686;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectHafniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 179.946550;

            return true;

            break;

        case 174:

            _atomicMass = 173.940046;

            return true;

            break;

        case 176:

            _atomicMass = 175.941409;

            return true;

            break;

        case 177:

            _atomicMass = 176.943221;

            return true;

            break;

        case 178:

            _atomicMass = 177.943699;

            return true;

            break;

        case 179:

            _atomicMass = 178.945816;

            return true;

            break;

        case 180:

            _atomicMass = 179.946550;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTantalumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 180.947996;

            return true;

            break;

        case 180:

            _atomicMass = 179.947465;

            return true;

            break;

        case 181:

            _atomicMass = 180.947996;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectTungstenIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 183.950931;

            return true;

            break;

        case 180:

            _atomicMass = 179.946704;

            return true;

            break;

        case 182:

            _atomicMass = 181.948204;

            return true;

            break;

        case 183:

            _atomicMass = 182.950223;

            return true;

            break;

        case 184:

            _atomicMass = 183.950931;

            return true;

            break;

        case 186:

            _atomicMass = 185.954364;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectRheniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 186.955753;

            return true;

            break;

        case 185:

            _atomicMass = 184.952955;

            return true;

            break;

        case 187:

            _atomicMass = 186.955753;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectOsmiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 191.961481;

            return true;

            break;

        case 184:

            _atomicMass = 183.952489;

            return true;

            break;

        case 186:

            _atomicMass = 185.953838;

            return true;

            break;

        case 187:

            _atomicMass = 186.955751;

            return true;

            break;

        case 188:

            _atomicMass = 187.955838;

            return true;

            break;

        case 189:

            _atomicMass = 188.958148;

            return true;

            break;

        case 190:

            _atomicMass = 189.958447;

            return true;

            break;

        case 192:

            _atomicMass = 191.961481;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectIridiumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 192.962926;

            return true;

            break;

        case 191:

            _atomicMass = 190.960594;

            return true;

            break;

        case 193:

            _atomicMass = 192.962926;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPlatinumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 194.964791;

            return true;

            break;

        case 190:

            _atomicMass = 189.959932;

            return true;

            break;

        case 192:

            _atomicMass = 191.961038;

            return true;

            break;

        case 194:

            _atomicMass = 193.962680;

            return true;

            break;

        case 195:

            _atomicMass = 194.964791;

            return true;

            break;

        case 196:

            _atomicMass = 195.964952;

            return true;

            break;

        case 198:

            _atomicMass = 197.967893;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectGoldIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 196.966569;

            return true;

            break;

        case 197:

            _atomicMass = 196.966569;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectMercuryIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 201.970643;

            return true;

            break;

        case 196:

            _atomicMass = 195.965833;

            return true;

            break;

        case 198:

            _atomicMass = 197.966769;

            return true;

            break;

        case 199:

            _atomicMass = 198.968280;

            return true;

            break;

        case 200:

            _atomicMass = 199.968326;

            return true;

            break;

        case 201:

            _atomicMass = 200.970302;

            return true;

            break;

        case 202:

            _atomicMass = 201.970643;

            return true;

            break;

        case 204:

            _atomicMass = 203.973494;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectThalliumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 204.974428;

            return true;

            break;

        case 203:

            _atomicMass = 202.972344;

            return true;

            break;

        case 205:

            _atomicMass = 204.974428;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectLeadIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 207.976652;

            return true;

            break;

        case 204:

            _atomicMass = 203.973044;

            return true;

            break;

        case 206:

            _atomicMass = 205.974465;

            return true;

            break;

        case 207:

            _atomicMass = 206.975897;

            return true;

            break;

        case 208:

            _atomicMass = 207.976652;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectBismuthIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 208.980399;

            return true;

            break;

        case 209:

            _atomicMass = 208.980399;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectPoloniumIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 208.982430;

            return true;

            break;

        case 209:

            _atomicMass = 208.982430;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectAstatineIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 209.987148;

            return true;

            break;

        case 210:

            _atomicMass = 209.987148;

            return true;

            break;

        default:

            return false;

            break;
    }
}

bool
CChemicalElement::_selectRadonIsotopeMass(const int32_t isotopeLabel)
{
    switch (isotopeLabel)
    {
        case 0:

            _atomicMass = 222.017578;

            return true;

            break;

        case 222:

            _atomicMass = 222.017578;

            return true;

            break;

        default:

            return false;

            break;
    }
}

std::ostream&
operator<<(std::ostream& output, const CChemicalElement& source)
{
    output << std::endl;

    output << "[CChemicalElement (Object):" << &source << "]" << std::endl;

    output << "_atomicLabel: " << source._atomicLabel << std::endl;

    output << "_atomicCharge: " << source._atomicCharge << std::endl;

    output << "_atomicMass: " << source._atomicMass << std::endl;

    output << "_atomicNumber: " << source._atomicNumber << std::endl;

    return output;
}
