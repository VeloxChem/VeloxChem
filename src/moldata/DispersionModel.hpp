//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DispersionModel_hpp
#define DispersionModel_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "MemBlock.hpp"
#include "Molecule.hpp"

/**
 Class CDispersionModel implements the D4 dispersion correction
 (reference: dftd4 v2.4.0).

 @author X. Li
 */
class CDispersionModel
{
    /**
     The dispersion energy.
     */
    double _energy;

    /**
     The dispersion gradient.
     */
    CDenseMatrix _gradient;

    /**
     The wf parameter (weighting factor).
     */
    double _wf;

    /**
     The g_a parameter (charge scale height).
     */
    double _g_a;

    /**
     The g_c parameter (charge scale steepness).
     */
    double _g_c;

    /**
     The atoms data of the dispersion model.
     */
    std::vector<int32_t> _atoms;

    /**
     The nref data of the dispersion model.
     */
    std::vector<int32_t> _nref;

    /**
     The ncount data of the dispersion model.
     */
    std::vector<std::vector<int32_t>> _ncount;

    /**
     The cn data of the dispersion model.
     */
    std::vector<std::vector<double>> _cn;

    /**
     The q data of the dispersion model.
     */
    std::vector<std::vector<double>> _q;

    /**
     The alpha data of the dispersion model.
     */
    std::vector<std::vector<std::vector<double>>> _alpha;

    /**
     The c6 data of the dispersion model.
     */
    std::vector<double> _c6;

    /**
     The ndim parameter (calculation dimension) of the dispersion model.
     */
    int32_t _ndim;

    /**
     The gw data (Gaussian weights) of the dispersion model.
     */
    std::vector<double> _gw;

    /**
     The dgw data of the dispersion model.
     */
    std::vector<double> _dgw;

    /**
     The refc6 data (reference C6 coefficients) of the dispersion model.
     */
    std::vector<double> _refc6;

    /**
     The charge scaling function.

     @param a parameter a.
     @param c parameter c.
     @param qref parameter qref.
     @param qmod parameter qmod.
     @return zeta.
     */
    double _zeta(const double a, const double c, const double qref, const double qmod);

    /**
     Derivative of charge scaling function w.r.t. charge.

     @param a parameter a.
     @param c parameter c.
     @param qref parameter qref.
     @param qmod parameter qmod.
     @return dzeta.
     */
    double _dzeta(const double a, const double c, const double qref, const double qmod);

    /**
     The coordination number Gaussian weight.

     @param wf parameter wf.
     @param cn parameter cn.
     @param cnref parameter cnref.
     @return the coordination number Gaussian weight.
     */
    double _cngw(const double wf, const double cn, const double cnref);

    /**
     Gets the weights for numerical Casimir-Polder integration.

     @return the weights.
     */
    CMemBlock<double> _getWeights();

    /**
     Gets the maximum value in a 2d vector.

     @param data the 2d vector.
     @return the maximum value.
     */
    double _getMaximum(const std::vector<std::vector<double>>& data);

    /**
     Gets the nearest integer value.

     @param value the input value.
     @return the nearest integer value.
     */
    int32_t _getRound(double value);

    /**
     Initializes the dispersion model.

     @param molecule the molecule.
     */
    void _initialize(const CMolecule& molecule);

    /**
     Computes Gaussian weights and reference c6 coefficients.

     @param molecule the molecule.
     */
    void _compWeightsAndCoefficients(const CMolecule& molecule);

    /**
     Computes two-body contribution to dispersion energy and gradient.

     @param molecule the molecule.
     @param xcLabel the label of the density functional.
     @param gradient the derivative matrix of dimension (3,N).
     @return the dispersion energy.
     */
    double _compTwoBodyContribution(const CMolecule& molecule, const std::string& xcLabel, CDenseMatrix& gradient);

    /**
     Computes three-body contribution to dispersion energy and gradient.

     @param molecule the molecule.
     @param xcLabel the label of the density functional.
     @param gradient the derivative matrix of dimension (3,N).
     @return the three-body contribution to the dispersion energy.
     */
    double _compThreeBodyContribution(const CMolecule& molecule, const std::string& xcLabel, CDenseMatrix& gradient);

   public:
    /**
     Creates a dispersion model object.
     */
    CDispersionModel();

    /**
     Destroys a dispersion model object.
     */
    ~CDispersionModel();

    /**
     Computes dispersion energy and gradient for a given molecule and a given
     density functional.

     @param molecule the molecule.
     @param xcLabel the label of the density functional.
     @return the dispersion energy.
     */
    void compute(const CMolecule& molecule, const std::string& xcLabel);

    /**
     Gets dispersion energy.

     @return the dispersion energy.
     */
    double getEnergy() const;

    /**
     Gets dispersion gradient.

     @return the dispersion gradient.
     */
    CDenseMatrix getGradient() const;
};

#endif /* DispersionModel_hpp */
