//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AtomBasisSetter.hpp"

namespace vlxbas { // vlxbas namespace
    
CAtomBasis
getAtomBasisEmpty()
{
    CAtomBasis bas;

    bas.setMaxAngularMomentum(-1);

    bas.setIdElemental(-1);

    return bas;
}

CAtomBasis
getAtomBasisForH()
{
    CAtomBasis bas;

    bas.setIdElemental(1);

    bas.addBasisFunction(CBasisFunction({1.301070100000e+01, 1.962257200000e+00,
                                         4.445379600000e-01},
                                        {1.968215800000e-02, 1.379652400000e-01,
                                         4.783193500000e-01}, 0));

    bas.addBasisFunction(CBasisFunction({1.219496200000e-01},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({8.000000000000e-01},
                                        {1.000000000000e+00}, 1));

    return bas;
}

CAtomBasis
getAtomBasisForLi()
{
    CAtomBasis bas;

    bas.setIdElemental(3);

    bas.addBasisFunction(CBasisFunction({2.662778551600e+02, 4.006978344700e+01,
                                         9.055994438900e+00, 2.450300905100e+00,
                                         7.220957185500e-01},
                                        {6.492015032500e-03, 4.774786321500e-02,
                                         2.026879611100e-01, 4.860657481700e-01,
                                         4.362697795500e-01}, 0));

    bas.addBasisFunction(CBasisFunction({1.450000000000e+00,
                                         3.000000000000e-01},
                                        {2.586000000000e-01,
                                         1.000000000000e+00}, 1));

    bas.addBasisFunction(CBasisFunction({5.281088472100e-02},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({2.096094879800e-02},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({8.200000000000e-02},
                                        {1.000000000000e+00}, 1));

    return bas;
}
    
CAtomBasis
getAtomBasisForLiX()
{
    CAtomBasis bas = vlxbas::getAtomBasisForLi();
    
    bas.addBasisFunction(CBasisFunction({3.000000000000e-02},
                                        {1.000000000000e+00}, 0));
    
    return bas;
}

CAtomBasis
getNormalizedAtomBasisForH()
{
    CAtomBasis bas;

    bas.setIdElemental(1);

    std::vector<CBasisFunction> bf;

    bf.push_back(CBasisFunction ({1.301070100000e+01,
                                  1.962257200000e+00,
                                  4.445379600000e-01},
                                 {1.968215800000e-02,
                                  1.379652400000e-01,
                                  4.783193500000e-01}, 0));

    bf.push_back(CBasisFunction ({1.219496200000e-01},
                                 {1.000000000000e+00}, 0));

    bf.push_back(CBasisFunction ({8.000000000000e-01},
                                 {1.000000000000e+00}, 1));

    for (size_t i = 0; i < bf.size(); i++)
    {
        bf[i].normalize();

        bas.addBasisFunction(bf[i]);
    }

    return bas;
}

CAtomBasis
getNormalizedAtomBasisForHe()
{
    CAtomBasis bas;

    bas.setIdElemental(2);

    std::vector<CBasisFunction> bf;

    bf.push_back(CBasisFunction({3.835493673700e+01,
                                 5.768908147900e+00,
                                 1.239940703500e+00},
                                {2.381428890500e-02,
                                 1.549090677700e-01,
                                 4.699809663300e-01}, 0));

    bf.push_back(CBasisFunction({2.975781595300e-01},
                                {1.000000000000e+00}, 0));

    bf.push_back(CBasisFunction({1.000000000000e+00},
                                {1.000000000000e+00}, 1));

    for (size_t i = 0; i < bf.size(); i++)
    {
        bf[i].normalize();

        bas.addBasisFunction(bf[i]);
    }

    return bas;
}
    
CAtomBasis
getAtomBasisSPDForHe()
{
    CAtomBasis bas;

    bas.setIdElemental(2);

    std::vector<CBasisFunction> bf;

    bf.push_back(CBasisFunction({2.2000000000000000},
                                {1.0000000000000000},
                                0));
    
    bf.push_back(CBasisFunction({1.3000000000000000},
                                {1.0000000000000000},
                                1));
    
    bf.push_back(CBasisFunction({0.8000000000000000},
                                {1.0000000000000000},
                                2));

    for (size_t i = 0; i < bf.size(); i++)
    {
        bf[i].normalize();

        bas.addBasisFunction(bf[i]);
    }

    return bas;
}
  
CAtomBasis
getNormalizedAtomBasisForO()
{
    CAtomBasis bas;

    bas.setIdElemental(8);

    std::vector<CBasisFunction> bf;

    bf.push_back(CBasisFunction ({ 2.266176778500e+03,
                                   3.408701019100e+02,
                                   7.736313516700e+01,
                                   2.147964494000e+01,
                                   6.658943312400e+00},
                                 {-5.343180992600e-03,
                                  -3.989003923000e-02,
                                  -1.785391198500e-01,
                                  -4.642768495900e-01,
                                  -4.430974517200e-01}, 0));
                       
    bf.push_back(CBasisFunction ({ 8.097597566800e-01},
                                 { 1.000000000000e+00}, 0));
                       
    bf.push_back(CBasisFunction ({ 2.553077223400e-01},
                                 { 1.000000000000e+00}, 0));
                       
    bf.push_back(CBasisFunction ({ 1.772150431700e+01,
                                   3.863550544000e+00,
                                   1.048092088300e+00},
                                 { 4.339457319300e-02,
                                   2.309412076500e-01,
                                   5.137531106400e-01}, 1));
                       
    bf.push_back(CBasisFunction ({ 2.764154441100e-01},
                                 { 1.000000000000e+00}, 1));
                       
    bf.push_back(CBasisFunction ({ 1.200000000000e+00},
                                 { 1.000000000000e+00}, 2));

    for (size_t i = 0; i < bf.size(); i++)
    {
        bf[i].normalize();

        bas.addBasisFunction(bf[i]);
    }

    return bas;
}

CAtomBasis
getNormalizedAtomBasisForSe()
{
    CAtomBasis bas;

    bas.setIdElemental(34);

    std::vector<CBasisFunction> bf;

    bf.push_back(CBasisFunction ({ 1.066122002700e+05, 1.600760470100e+04,
                                   3.642169970700e+03, 1.028591299300e+03,
                                   3.353029888800e+02, 1.158012915400e+02},
                                 { 1.427488911300e-03, 1.094352511400e-02,
                                   5.437417159600e-02, 1.901809294700e-01,
                                   3.891302169600e-01, 3.062020708800e-01}, 0));

    bf.push_back(CBasisFunction ({ 2.229332502000e+02, 2.672625793400e+01,
                                   1.112450192300e+01},
                                 {-1.119880800900e-01, 6.475212420700e-01,
                                   4.424197665100e-01}, 0));

    bf.push_back(CBasisFunction ({ 1.988852006100e+01, 3.366847380300e+00,
                                   1.524927783900e+00},
                                 {-2.285722776200e-01, 7.359135995100e-01,
                                   4.433019957700e-01}, 0));

    bf.push_back(CBasisFunction ({ 2.963003791200e-01},
                                 { 1.000000000000e+00}, 0));

    bf.push_back(CBasisFunction ({ 1.100928897400e-01},
                                 { 1.000000000000e+00}, 0));

    bf.push_back(CBasisFunction ({ 1.455906812000e+03, 3.437510183100e+02,
                                   1.092955496400e+02, 3.970771102200e+01,
                                   1.495018523200e+01},
                                 { 8.820359704300e-03, 6.687585196700e-02,
                                   2.664057808000e-01, 5.383492842200e-01,
                                   3.630328199300e-01}, 1));

    bf.push_back(CBasisFunction ({ 8.020896209400e+00, 3.293464975600e+00,
                                   1.405860243800e+00},
                                 { 3.415380702500e-01, 5.725790658300e-01,
                                   2.554981322200e-01}, 1));

    bf.push_back(CBasisFunction ({ 4.507612322600e-01},
                                 { 1.000000000000e+00}, 1));

    bf.push_back(CBasisFunction ({ 1.335341332500e-01},
                                 { 1.000000000000e+00}, 1));

    bf.push_back(CBasisFunction ({ 9.449402404400e+01, 2.718818526000e+01,
                                   9.509156735200e+00, 3.417051685300e+00},
                                 { 2.349010109800e-02, 1.374773576700e-01,
                                   3.664992912400e-01, 4.875098988400e-01}, 2));

    bf.push_back(CBasisFunction ({ 1.147959008300e+00},
                                 { 1.000000000000e+00}, 2));

    bf.push_back(CBasisFunction ({ 3.380000000000e-01},
                                 { 1.000000000000e+00}, 2));

    for (size_t i = 0; i < bf.size(); i++)
    {
        bf[i].normalize();

        bas.addBasisFunction(bf[i]);
    }

    return bas;
}

CAtomBasis
getMinimalBasisForH()
{
    CAtomBasis bas;

    bas.setIdElemental(1);

    CBasisFunction bf_1 ({1.301000000000e+01, 1.962000000000e+00,
                          4.446000000000e-01, 1.220000000000e-01},
                         {1.968500000000e-02, 1.379770000000e-01,
                          4.781480000000e-01, 5.012400000000e-01}, 0);
    bf_1.normalize();

    bas.addBasisFunction(bf_1);

    return bas;
}

CAtomBasis
getMinimalBasisForHe()
{
    CAtomBasis bas;

    bas.setIdElemental(2);

    CBasisFunction bf_1 ({3.836000000000e+01, 5.770000000000e+00,
                          1.240000000000e+00, 2.976000000000e-01},
                         {2.380900000000e-02, 1.548910000000e-01,
                          4.699870000000e-01, 5.130270000000e-01}, 0);

    bf_1.normalize();

    bas.addBasisFunction(bf_1);

    return bas;
}

CAtomBasis
getMinimalBasisForC()
{
    CAtomBasis bas;

    bas.setIdElemental(6);

    CBasisFunction bf_1 ({ 6.665000000000e+03,  1.000000000000e+03,
                           2.280000000000e+02,  6.471000000000e+01,
                           2.106000000000e+01,  7.495000000000e+00,
                           2.797000000000e+00,  5.215000000000e-01,
                           1.596000000000e-01},
                         { 6.920000000000e-04,  5.329000000000e-03,
                           2.707700000000e-02,  1.017180000000e-01,
                           2.747400000000e-01,  4.485640000000e-01,
                           2.850740000000e-01,  1.520400000000e-02,
                          -3.191000000000e-03}, 0);

    CBasisFunction bf_2 ({ 6.665000000000e+03,  1.000000000000e+03,
                           2.280000000000e+02,  6.471000000000e+01,
                           2.106000000000e+01,  7.495000000000e+00,
                           2.797000000000e+00,  5.215000000000e-01,
                           1.596000000000e-01},
                         {-1.460000000000e-04, -1.154000000000e-03,
                          -5.725000000000e-03, -2.331200000000e-02,
                          -6.395500000000e-02, -1.499810000000e-01,
                          -1.272620000000e-01,  5.445290000000e-01,
                           5.804960000000e-01}, 0);

    CBasisFunction bf_3 ({ 9.439000000000e+00,  2.002000000000e+00,
                           5.456000000000e-01,  1.517000000000e-01},
                         { 3.810900000000e-02,  2.094800000000e-01,
                           5.085570000000e-01,  4.688420000000e-01}, 1);

    bf_1.normalize();

    bf_2.normalize();

    bf_3.normalize();

    bas.addBasisFunction(bf_1);

    bas.addBasisFunction(bf_2);

    bas.addBasisFunction(bf_3);

    return bas;
}

CAtomBasis
getMinimalBasisForN()
{
    CAtomBasis bas;

    bas.setIdElemental(7);


    CBasisFunction bf_1 ({ 9.046000000000e+03,  1.357000000000e+03,
                           3.093000000000e+02,  8.773000000000e+01,
                           2.856000000000e+01,  1.021000000000e+01,
                           3.838000000000e+00,  7.466000000000e-01,
                           2.248000000000e-01},
                         { 7.000000000000e-04,  5.389000000000e-03,
                           2.740600000000e-02,  1.032070000000e-01,
                           2.787230000000e-01,  4.485400000000e-01,
                           2.782380000000e-01,  1.544000000000e-02,
                          -2.864000000000e-03}, 0);

    CBasisFunction bf_2 ({ 9.046000000000e+03,  1.357000000000e+03,
                           3.093000000000e+02,  8.773000000000e+01,
                           2.856000000000e+01,  1.021000000000e+01,
                           3.838000000000e+00,  7.466000000000e-01,
                           2.248000000000e-01},
                         {-1.530000000000e-04, -1.208000000000e-03,
                          -5.992000000000e-03, -2.454400000000e-02,
                          -6.745900000000e-02, -1.580780000000e-01,
                          -1.218310000000e-01,  5.490030000000e-01,
                           5.788150000000e-01}, 0);

    CBasisFunction bf_3 ({ 1.355000000000e+01,  2.917000000000e+00,
                           7.973000000000e-01,  2.185000000000e-01},
                         { 3.991900000000e-02,  2.171690000000e-01,
                           5.103190000000e-01,  4.622140000000e-01}, 1);

    bf_1.normalize();

    bf_2.normalize();

    bf_3.normalize();

    bas.addBasisFunction(bf_1);

    bas.addBasisFunction(bf_2);

    bas.addBasisFunction(bf_3);

    return bas;
}

CAtomBasis
getMinimalBasisForO()
{
    CAtomBasis bas;

    bas.setIdElemental(8);

    CBasisFunction bf_1 ({ 1.172000000000e+04,  1.759000000000e+03,
                           4.008000000000e+02,  1.137000000000e+02,
                           3.703000000000e+01,  1.327000000000e+01,
                           5.025000000000e+00,  1.013000000000e+00,
                           3.023000000000e-01},
                         { 7.100000000000e-04,  5.470000000000e-03,
                           2.783700000000e-02,  1.048000000000e-01,
                           2.830620000000e-01,  4.487190000000e-01,
                           2.709520000000e-01,  1.545800000000e-02,
                          -2.585000000000e-03}, 0);

    CBasisFunction bf_2 ({ 1.172000000000e+04,  1.759000000000e+03,
                           4.008000000000e+02,  1.137000000000e+02,
                           3.703000000000e+01,  1.327000000000e+01,
                           5.025000000000e+00,  1.013000000000e+00,
                           3.023000000000e-01},
                         {-1.600000000000e-04, -1.263000000000e-03,
                          -6.267000000000e-03, -2.571600000000e-02,
                          -7.092400000000e-02, -1.654110000000e-01,
                          -1.169550000000e-01,  5.573680000000e-01,
                           5.727590000000e-01}, 0);

    CBasisFunction bf_3 ({ 1.770000000000e+01,  3.854000000000e+00,
                           1.046000000000e+00,  2.753000000000e-01},
                         { 4.301800000000e-02,  2.289130000000e-01,
                           5.087280000000e-01,  4.605310000000e-01}, 1);

    bf_1.normalize();

    bf_2.normalize();

    bf_3.normalize();

    bas.addBasisFunction(bf_1);

    bas.addBasisFunction(bf_2);

    bas.addBasisFunction(bf_3);

    return bas;
}

CAtomBasis
getGenContrBasisForLi()
{
    CAtomBasis bas;
    
    bas.setIdElemental(3);
    
    CBasisFunction bf1({1.400, 1.200, 0.300},
                       {0.800, 0.200, 0.300, 0.400, -0.200, 0.900},
                       2, 0);
    
    bas.addBasisFunction(bf1);
    
    CBasisFunction bf2({1.500}, {1.000}, 0);
    
    bas.addBasisFunction(bf2);
    
    CBasisFunction bf3({0.640}, {1.000}, 1);
    
    bas.addBasisFunction(bf3);
    
    return bas;
}
    
CAtomBasis
getTestBasisForH()
{
    CAtomBasis bas;
    
    bas.setIdElemental(1);
    
    for (int32_t i = 0; i < 5; i++)
    {
        CBasisFunction bfone({3.0, 2.0}, {0.5, 0.5}, i);
    
        bfone.normalize();
    
        CBasisFunction bftwo({0.8}, {1.0}, i);
    
        bftwo.normalize();
    
        bas.addBasisFunction(bfone);
    
        bas.addBasisFunction(bftwo);
    }
    
    return bas; 
}
    
CAtomBasis
getTestBasisForLi()
{
    CAtomBasis bas;
    
    bas.setIdElemental(3);
    
    for (int32_t i = 0; i < 5; i++)
    {
        CBasisFunction bfone({2.8, 1.5}, {0.5, 0.5}, i);
        
        bfone.normalize();
        
        CBasisFunction bftwo({1.2}, {1.0}, i);
        
        bftwo.normalize();
        
        bas.addBasisFunction(bfone);
        
        bas.addBasisFunction(bftwo);
    }
        
    return bas;
}
    
} // vlxbas namespace
