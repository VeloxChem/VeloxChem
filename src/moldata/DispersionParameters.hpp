//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef DispersionParameters_hpp
#define DispersionParameters_hpp

#include <string>

/**
 Class CDispersionParameters implements the density functional parameters for
 the D4 dispersion correction (reference: dftd4 v2.4.0).

 @author X. Li
 */
class CDispersionParameters
{
    /**
     The s6 parameter (scaling of dipole-dipole dispersion).
     */
    double _s6;

    /**
     The s8 parameter (scaling of dipole-quadrupole dispersion).
     */
    double _s8;

    /**
     The s10 parameter (scaling of higher order dispersion).
     */
    double _s10;

    /**
     The a1 parameter (scaling of vdW-Radius in finite damping).
     */
    double _a1;

    /**
     The a2 parameter (constant offset off vdW-Radius in finite damping).
     */
    double _a2;

    /**
     The s9 parameter (scaling of non-addititive dispersion).
     */
    double _s9;

    /**
     The alp parameter (exponent of zero damping).
     */
    int32_t _alp;

    /**
     The beta parameter (range separation parameter for Fermi-damping).
     */
    double _beta;

    /**
     Sets default parameters.
     */
    void _setDefaultParameters();

    /**
     Sets parameters for a given density functional.

     @param xcLabel the label of the density functional.
     */
    void _setFunctionalParameters(const std::string& xcLabel);

    /**
     Sets the s6, s8, a1, a2 parameters.

     @param s6 the s6 parameter.
     @param s8 the s8 parameter.
     @param a1 the a1 parameter.
     @param a2 the a2 parameter.
     */
    void _setFourParameters(const double s6, const double s8, const double a1, const double a2);

   public:
    /**
     Creates a dispersion parameters object.
     */
    CDispersionParameters();

    /**
     Creates a dispersion parameters object for a given density functional.

     @param xcLabel the label of the density functional.
     */
    CDispersionParameters(const std::string& xcLabel);

    /**
     Destroys a dispersion parameters object.
     */
    ~CDispersionParameters();

    /**
     Gets the s6 parameter.

     @return the s6 parameter.
     */
    double getS6() const;

    /**
     Gets the s8 parameter.

     @return the s8 parameter.
     */
    double getS8() const;

    /**
     Gets the s10 parameter.

     @return the s10 parameter.
     */
    double getS10() const;

    /**
     Gets the a1 parameter.

     @return the a1 parameter.
     */
    double getA1() const;

    /**
     Gets the a2 parameter.

     @return the a2 parameter.
     */
    double getA2() const;

    /**
     Gets the alp parameter.

     @return the alp parameter.
     */
    int32_t getAlp() const;
};

#endif /* DispersionParameters_hpp */
