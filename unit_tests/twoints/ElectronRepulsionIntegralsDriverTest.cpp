//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionIntegralsDriverTest.hpp"

#include "mpi.h"

#include "ElectronRepulsionIntegralsDriver.hpp"
#include "AODensityMatrixSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"
#include "EriScreenerType.hpp"

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(1);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));
    
    std::vector<double> r0000vals{ 1.629232440467330};
    
    vlxtest::compare(r0000vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));
    
    std::vector<double> r0101vals{ 1.296490921648330};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));
    
    std::vector<double> r0123vals{ 0.569573705357336};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));
    
    std::vector<double> r0312vals{ 0.091957062417462};
    
    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));
    
    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0000vals,fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));
    
    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));
    
    std::vector<double> r0123vals{-0.048499014664650, -0.088914860218526, -0.032332676443100};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));
    
    std::vector<double> r0312vals{-0.065342004894742, -0.119793675640361, -0.043561336596495};
    
    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));
    
    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0000vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));
    
    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(11), kpairs.pick(1));
    
    std::vector<double> r2301vals{-0.048499014664650, -0.088914860218526, -0.032332676443100};
    
    vlxtest::compare(r2301vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(6), kpairs.pick(3));
    
    std::vector<double> r1203vals{-0.065342004894742, -0.119793675640361, -0.043561336596495};
    
    vlxtest::compare(r1203vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 0.271300324109232,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.271300324109232,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.271300324109232};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.158211411411533,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.158211411411533,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.158211411411533};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(11));

    std::vector<double> r0123vals{ 0.014876488624304, -0.023320174063975, -0.008480063295991,
                                  -0.023320174063975, -0.015157068882331, -0.015546782709317,
                                  -0.008480063295991, -0.015546782709317,  0.021943208037630};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(6));

    std::vector<double> r0312vals{ 0.056736244044378,  0.081977329675453,  0.029809938063801,
                                   0.081977329675453,  0.162313108020340,  0.054651553116969,
                                   0.029809938063801,  0.054651553116969,  0.031894628991210};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.357416654829900,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.357416654829900,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.357416654829900};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 1.043481505026390,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.043481505026390,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.043481505026390};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.472688274973429,  0.033593694100454,  0.012215888763801,
                                   0.033593694100454,  0.515952881011892,  0.022395796066969,
                                   0.012215888763801,  0.022395796066969,  0.462508367670261};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.023339095978464, -0.095784363951802, -0.034830677800655,
                                  -0.095784363951802, -0.100019554565523, -0.063856242634535,
                                  -0.034830677800655, -0.063856242634535,  0.052364660812343};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.357416654829900,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.357416654829900,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.357416654829900};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 1.043481505026390,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.043481505026390,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.043481505026390};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.527902875740139,  0.037278021172033,  0.013555644062557,
                                   0.037278021172033,  0.575912448461696,  0.024852014114689,
                                   0.013555644062557,  0.024852014114689,  0.516606505688008};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.032590929327399, -0.067864620590098, -0.024678043850945,
                                  -0.067864620590098, -0.054810475978031, -0.045243080393399,
                                  -0.024678043850945, -0.045243080393399,  0.053155965869853};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(27);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.024740456719223, -0.079866951735308, -0.029042527903748,
                                   0.006120732901895, -0.011000649215887, -0.007422660648759,
                                   0.002225721055234, -0.007422660648759,  0.006712518241378,
                                   0.006120732901895, -0.011000649215887, -0.007422660648759,
                                  -0.057903022654361, -0.071646093782930, -0.038602015102907,
                                  -0.007422660648759, -0.007333766143925,  0.012306283442527,
                                   0.002225721055234, -0.007422660648759,  0.006712518241378,
                                  -0.007422660648759, -0.007333766143925,  0.012306283442527,
                                  -0.040189855197096, -0.073681401194676, -0.014244346707131};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.025187805950051, -0.039734523702261, -0.014448917709913,
                                   0.027147577114811,  0.097741908029784,  0.032527520363146,
                                   0.009871846223568,  0.032527520363146,  0.020119416254094,
                                   0.063991863435128,  0.077645024582338,  0.032527520363146,
                                   0.041163878682117,  0.069023990377715,  0.027442585788078,
                                   0.032527520363146,  0.051763349721559,  0.036885596465839,
                                   0.023269768521865,  0.032527520363146,  0.000022532806649,
                                   0.032527520363146,  0.065161272019856,  0.000041310145522,
                                  -0.036458613093572, -0.066840790671549, -0.026648694985836};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(27);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.027005659425976, -0.006710881773728, -0.002440320644992,
                                  -0.006710881773728,  0.063368692728713,  0.008139179236825,
                                  -0.002440320644992,  0.008139179236825,  0.043945651368109,
                                   0.087349676872220,  0.012062933467287,  0.008139179236825,
                                   0.012062933467287,  0.078336635411377,  0.008041955644858,
                                   0.008139179236825,  0.008041955644858,  0.080567027508199,
                                   0.031763518862625,  0.008139179236825, -0.007360107893318,
                                   0.008139179236825,  0.042245795152476, -0.013493531137749,
                                  -0.007360107893318, -0.013493531137749,  0.015537355000097};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.014522101231219,  0.053494871659473,  0.019452680603445,
                                   0.082466044299379,  0.067897378051834,  0.043882287079406,
                                   0.029987652472501,  0.043882287079406, -0.036821716114931,
                                  -0.030937906977868,  0.129754078748514,  0.043882287079406,
                                   0.113951620944929,  0.128792581148996,  0.075967747296619,
                                   0.043882287079406,  0.086502719165676, -0.067506479544040,
                                  -0.011250147991952,  0.043882287079406,  0.025034984581749,
                                   0.043882287079406,  0.045264918701223,  0.045897471733207,
                                   0.009232526778164,  0.016926299093301, -0.022979063572148};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CMemBlock<double> fints(81);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 1.324120036848930,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.165219477049000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.165219477049000,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.165219477049000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.324120036848930,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.165219477049000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   1.165219477049000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.165219477049000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.324120036848930};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.985817209017402,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.866952226906758,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.866952226906758,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.866952226906758,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.985817209017402,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.866952226906758,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.866952226906758,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.866952226906758,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.985817209017402};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));

    std::vector<double> r0123vals{ 0.434400987263440,  0.012306427508151,  0.004475064548418,
                                   0.012306427508151,  0.477839565430902,  0.020216204162471,
                                   0.004475064548418,  0.020216204162471,  0.429596350952278,
                                   0.012099176064835, -0.001117051089860, -0.001395907024805,
                                  -0.001117051089860,  0.028426445525904,  0.005175785063223,
                                  -0.001395907024805,  0.005175785063223,  0.026275344338962,
                                   0.004399700387213, -0.001395907024805,  0.002214090673878,
                                  -0.001395907024805,  0.016888845185986, -0.005935570998635,
                                   0.002214090673878, -0.005935570998635,  0.003002714764875,
                                   0.012099176064835, -0.001117051089860, -0.001395907024805,
                                  -0.001117051089860,  0.028426445525904,  0.005175785063223,
                                  -0.001395907024805,  0.005175785063223,  0.026275344338962,
                                   0.477572650693298,  0.028633696969220,  0.016964209347191,
                                   0.028633696969220,  0.486859990413360,  0.019089131312813,
                                   0.016964209347191,  0.019089131312813,  0.463435809570639,
                                   0.020078036533594,  0.005175785063223, -0.005935570998635,
                                   0.005175785063223,  0.018950963683936, -0.005430205309213,
                                  -0.005935570998635, -0.005430205309213,  0.005504977068938,
                                   0.004399700387213, -0.001395907024805,  0.002214090673878,
                                  -0.001395907024805,  0.016888845185986, -0.005935570998635,
                                   0.002214090673878, -0.005935570998635,  0.003002714764875,
                                   0.020078036533594,  0.005175785063223, -0.005935570998635,
                                   0.005175785063223,  0.018950963683936, -0.005430205309213,
                                  -0.005935570998635, -0.005430205309213,  0.005504977068938,
                                   0.429659154419950,  0.026482595782278,  0.003078078926081,
                                   0.026482595782278,  0.463765527775914,  0.005643144697815,
                                   0.003078078926081,  0.005643144697815,  0.428169504502362};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));

    std::vector<double> r0312vals{ 0.016566258811693, -0.029593027022679, -0.010761100735520,
                                  -0.038062728927181, -0.032163055222041, -0.022143397895054,
                                  -0.013840992337157, -0.022143397895054,  0.020679144300248,
                                  -0.020486981860483,  0.102859587179471,  0.033939425136731,
                                   0.039799609976381,  0.077734864017329,  0.044178971788063,
                                   0.014971833782310,  0.051445039107019, -0.037952465807687,
                                  -0.007449811585630,  0.033939425136731,  0.021867777194091,
                                   0.014971833782310,  0.027245589113936,  0.028596021673293,
                                   0.004071370268595,  0.008614336546164, -0.012779262492248,
                                  -0.013642431687775,  0.049661060177614,  0.020076584474713,
                                   0.101813087958240,  0.070890313844621,  0.051740521240073,
                                   0.033397707892799,  0.041394562319480, -0.037952465807687,
                                  -0.010574122614948,  0.053314566777922,  0.021149890495462,
                                   0.061784268682424,  0.068705659285387,  0.041189512454949,
                                   0.021149890495462,  0.035543044518615, -0.028199031361167,
                                  -0.010967950427853,  0.041394562319480,  0.028596021673293,
                                   0.051740521240073,  0.047260209229747,  0.058695986924846,
                                   0.008614336546164,  0.015165591578048, -0.023428647902454,
                                  -0.004960884250100,  0.020076584474713,  0.001751029044775,
                                   0.033397707892799,  0.027245589113936,  0.000957210507559,
                                   0.022114012304969,  0.029408597539190, -0.015268189827778,
                                  -0.010967950427853,  0.051445039107019,  0.000957210507559,
                                   0.044178971788063,  0.051823242678219,  0.002983800152995,
                                   0.029408597539190,  0.059988721256955, -0.027991681350926,
                                   0.015599395451521, -0.054715637008488, -0.021659370579305,
                                  -0.054715637008488, -0.054867712816987, -0.039708846062059,
                                  -0.018579478977668, -0.034062378125724,  0.038503569391985};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(5);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{ 0.007129352251451,  0.019605718691489,  0.016293028206027,
                                   0.013070479127660, -0.002970563438104};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{ 0.031482119893843,  0.086575829708068,  0.071947499481674,
                                   0.057717219805379, -0.013117549955768};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(5);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.010711063256865,  0.029455423956380,  0.024478472883023,
                                   0.019636949304253, -0.004462943023694};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.029700398508160,  0.081676095897440,  0.067875651750163,
                                   0.054450730598293, -0.012375166045067};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 1);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(15);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.026166985859894, -0.071959211114710,  0.010744833997060,
                                  -0.006932478908960,  0.035154013980002, -0.006932478908960,
                                  -0.052642767591242, -0.086927038856231, -0.035095178394161,
                                   0.002888532878733, -0.036099352013042, -0.006932478908960,
                                   0.007163222664707, -0.066182145357243, -0.021335258135468};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.010317693245241, -0.028373656424414,  0.087267806528214,
                                   0.045570475845650,  0.042404548322110,  0.045570475845650,
                                   0.072557333924872, -0.007549205123164,  0.048371555949914,
                                  -0.018987698269021, -0.036190392524976,  0.045570475845650,
                                   0.058178537685476, -0.066349052962455, -0.042078933986149};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePDSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(15);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.035644199512795,  0.098021548660186, -0.014706466132464,
                                   0.009402549007869, -0.047910247582724,  0.009402549007869,
                                   0.071630314398030,  0.118387909036415,  0.047753542932020,
                                  -0.003917728753279,  0.049192413356525,  0.009402549007869,
                                  -0.009804310754976,  0.090186091153628,  0.029090907780038};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.015451970508822, -0.042492918899260,  0.067033294236219,
                                   0.031212160982557,  0.041621505434743,  0.031212160982557,
                                   0.037118264316706, -0.031797094682732,  0.024745509544471,
                                  -0.013005067076065, -0.037365301664395,  0.031212160982557,
                                   0.044688862824146, -0.068503053051390, -0.037205900890604};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(15);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.001437379374511, -0.003952793279906,  0.009368058669662,
                                   0.004725762740340,  0.004948565226513,  0.004725762740340,
                                   0.006973245322980, -0.001949556321075,  0.004648830215320,
                                  -0.001969067808475, -0.004304143034648,  0.004725762740340,
                                   0.006245372446442, -0.007890928896855, -0.004731092799596};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.027661754243450, -0.076069824169487, -0.043347547948289,
                                  -0.039154218465543,  0.018356047726445, -0.039154218465543,
                                  -0.117131462491792, -0.109501661386788, -0.078087641661195,
                                   0.016314257693976, -0.023695259335382, -0.039154218465543,
                                  -0.028898365298859, -0.043441308781534, -0.000372450464435};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(15);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{ 0.001790472754338, -0.005806431999702,  0.005324996404098,
                                   0.004923800074429, -0.008531261777010, -0.005806431999702,
                                  -0.011531440463354,  0.002472850563674, -0.007687626975569,
                                  -0.005806431999702, -0.005687507851340,  0.009762493407514,
                                  -0.006116782585876,  0.002419346666543,  0.005837379905645};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.025753364783031, -0.038083906994492, -0.021319180358868,
                                  -0.070821753153335, -0.112201231141178, -0.038083906994492,
                                  -0.043160472049359, -0.102849514035386, -0.028773648032906,
                                  -0.038083906994492, -0.074800820760785, -0.039085163991258,
                                   0.016125920314164,  0.015868294581038,  0.000789964334342};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CMemBlock<double> fints(25);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 0.083187219132364,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.083187219132364,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.083187219132364,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.083187219132364,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.083187219132364};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.035447085056383,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.035447085056383,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.035447085056383,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.035447085056383,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.035447085056383};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(11));

    std::vector<double> r0123vals{ 0.001044915069104, -0.000670837481573,  0.003427651300041,
                                  -0.002720492103989, -0.000311679251205, -0.000670837481573,
                                  -0.000555947102831,  0.002338695614901,  0.002283362281660,
                                   0.003234762868478,  0.003427651300041,  0.002338695614901,
                                  -0.003649406880990,  0.001559130409934, -0.001428188041684,
                                  -0.002720492103989,  0.002283362281660,  0.001559130409934,
                                  -0.002458749004215, -0.003299332500071, -0.000311679251205,
                                   0.003234762868478, -0.001428188041684, -0.003299332500071,
                                   0.000426751220881};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(6));

    std::vector<double> r0312vals{ 0.021338339176323,  0.034215212887908,  0.017984236420477,
                                   0.028771074480771, -0.004100317245653,  0.034215212887908,
                                   0.102988279022468,  0.068040998640218,  0.053515539981455,
                                  -0.022005551025443,  0.017984236420477,  0.068040998640218,
                                   0.083486227948454,  0.045360665760145, -0.007493431841865,
                                   0.028771074480771,  0.053515539981455,  0.045360665760145,
                                   0.058391995704589, -0.000364129217099, -0.004100317245653,
                                  -0.022005551025443, -0.007493431841865, -0.000364129217099,
                                   0.013206043305778};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(25);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.385370186031906,  0.018261768434687, -0.005108780522769,
                                   0.023745770032733, -0.000663069567145,  0.018261768434687,
                                   0.428949406160136,  0.022026464935234,  0.015597055921471,
                                  -0.022651705246943, -0.005108780522769,  0.022026464935234,
                                   0.432063647748182,  0.014684309956822,  0.002128658551154,
                                   0.023745770032733,  0.015597055921471,  0.014684309956822,
                                   0.415951859558910,  0.012669881751764, -0.000663069567145,
                                  -0.022651705246943,  0.002128658551154,  0.012669881751764,
                                   0.384055098057068};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.015076639794480, -0.015698705577617,  0.063153907323239,
                                  -0.053933060493743, -0.005524545841330, -0.015698705577617,
                                  -0.022386180333923,  0.038155750177978,  0.038395709336246,
                                   0.063048561131938,  0.063153907323239,  0.038155750177978,
                                  -0.079093237014484,  0.025437166785318, -0.026314128051350,
                                  -0.053933060493743,  0.038395709336246,  0.025437166785318,
                                  -0.054382604780795, -0.062289042172837, -0.005524545841330,
                                   0.063048561131938, -0.026314128051350, -0.062289042172837,
                                   0.004119623875841};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(25);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.169084424724620};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.868177509086430};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.480186085555079,  0.022547697024725, -0.006394868186558,
                                   0.029368438901709, -0.000809655751887,  0.022547697024725,
                                   0.533993089818627,  0.027111335248510,  0.019180817722707,
                                  -0.028032506911096, -0.006394868186558,  0.027111335248510,
                                   0.537918299299586,  0.018074223499007,  0.002664528411066,
                                   0.029368438901709,  0.019180817722707,  0.018074223499007,
                                   0.518009075049705,  0.015719600183812, -0.000809655751887,
                                  -0.028032506911096,  0.002664528411066,  0.015719600183812,
                                   0.478580268313837};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.021662986735168, -0.015325507819405,  0.036560732032704,
                                  -0.038337645143337, -0.002790796984405, -0.015325507819405,
                                  -0.014909247833867,  0.012870534308492,  0.015362406465951,
                                   0.042942460167605,  0.036560732032704,  0.012870534308492,
                                  -0.047197089899741,  0.008580356205662, -0.015233638346960,
                                  -0.038337645143337,  0.015362406465951,  0.008580356205662,
                                  -0.027711253222159, -0.038861229054556, -0.002790796984405,
                                   0.042942460167605, -0.015233638346960, -0.038861229054556,
                                   0.016127906049431};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 1);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(45);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000, -0.109588304547766,
                                   0.000000000000000, -0.189812511392062,  0.000000000000000,
                                   0.189812511392062,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.189812511392062,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.189812511392062,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.219176609095532,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.189812511392062,  0.000000000000000,
                                   0.189812511392062,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.189812511392062,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.109588304547766,  0.000000000000000,  0.189812511392062};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000, -0.061950265531195,
                                   0.000000000000000, -0.107301007442413,  0.000000000000000,
                                   0.107301007442413,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.107301007442413,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.107301007442413,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.123900531062391,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.107301007442413,  0.000000000000000,
                                   0.107301007442413,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.107301007442413,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.061950265531195,  0.000000000000000,  0.107301007442413};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.002683710523490, -0.007380203939598, -0.007167623648171,
                                   0.000912583709132, -0.014592444906776,  0.000912583709132,
                                   0.013742566412044, -0.026830573950940, -0.010672789139108,
                                  -0.004553828377102,  0.011564809833433,  0.000912583709132,
                                  -0.000606374839587, -0.016487860693802, -0.006938806585022,
                                  -0.012625216264280, -0.012410001718412,  0.008341662114042,
                                  -0.003288534608156,  0.015752871583017, -0.003288534608156,
                                  -0.017792549675725, -0.003915904142357, -0.011861699783817,
                                   0.001370222753398, -0.016487860693802, -0.003288534608156,
                                   0.005561108076028, -0.009669556211615, -0.008868605253599,
                                   0.017718355275893, -0.012625216264280, -0.000606374839587,
                                   0.003818839266316,  0.007829702476883,  0.003818839266316,
                                  -0.010672789139108, -0.017887049300627,  0.022636557361301,
                                   0.004669194719813, -0.001792814631629,  0.003818839266316,
                                  -0.006662311281849, -0.003286826824654,  0.016457663721409};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{ 0.006312829792922,  0.017360281930535, -0.052516270103774,
                                  -0.024493271656511, -0.032510995227073, -0.039687968665287,
                                  -0.063994825996378,  0.013963665423634, -0.032534737638485,
                                   0.025897115775664,  0.030715099227601, -0.039687968665287,
                                  -0.037159672067074,  0.046362467968700,  0.033802889131605,
                                   0.013659677506177,  0.039063126834858, -0.124637585026139,
                                  -0.061633453732011, -0.054124353481531, -0.053345437181769,
                                  -0.092599762834564,  0.002508729127412, -0.061733175223043,
                                   0.022227265492404,  0.046362467968700, -0.061633453732011,
                                  -0.083091723350760,  0.090424338278200,  0.053331536793977,
                                   0.006466169149664,  0.013659677506177, -0.037159672067074,
                                  -0.033659001288495, -0.024394543055446, -0.033659001288495,
                                  -0.032534737638485,  0.009309110282423, -0.036882544630974,
                                  -0.000016109377486,  0.015721175869080, -0.010866955775331,
                                  -0.021549876714546,  0.028822155759980,  0.023330159534572};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePDSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(45);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.189812511392062,
                                   0.000000000000000,  0.189812511392062,  0.000000000000000,
                                  -0.109588304547766,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.189812511392062,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.189812511392062,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.219176609095532,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.189812511392062,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.189812511392062,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.109588304547766,
                                   0.000000000000000,  0.189812511392062,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.189812511392062};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.107301007442413,
                                   0.000000000000000,  0.107301007442413,  0.000000000000000,
                                  -0.061950265531195,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.107301007442413,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.107301007442413,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.123900531062391,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.107301007442413,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.107301007442413,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.061950265531195,
                                   0.000000000000000,  0.107301007442413,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.107301007442413};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.002970099659877, -0.014079475120361,  0.019766004348188,
                                  -0.008167774064661, -0.013832743098127, -0.014079475120361,
                                  -0.007965533717964,  0.009348894097912, -0.000679047548725,
                                   0.001041807263445, -0.003649981038987,  0.004292160013374,
                                  -0.016282840078501,  0.017570928826872,  0.008744235685386,
                                   0.001041807263445, -0.003649981038987,  0.004292160013374,
                                   0.015378825143315, -0.019767932817041, -0.011897953248335,
                                  -0.029937349948609, -0.004283408160093, -0.019958233299072,
                                  -0.011897953248335, -0.013178621878027,  0.025293786183594,
                                  -0.005110994748836,  0.001520825432911,  0.005226962578029,
                                   0.012892694173731, -0.018389763562334, -0.001977572707801,
                                   0.001041807263445, -0.003649981038987,  0.004292160013374,
                                  -0.000679047548725,  0.006232596065275, -0.007399660760693,
                                  -0.018389763562334, -0.010791092232305, -0.003625549964302,
                                  -0.007751708733311, -0.009894319805776,  0.018344370231701};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{ 0.011919274924097,  0.018790574747359,  0.009165925174276,
                                   0.032778006041266,  0.054007069457745,  0.018790574747359,
                                  -0.036914176376626, -0.097162163437753, -0.029462148373191,
                                  -0.012979042929044, -0.044186399359580, -0.024154654366061,
                                  -0.032964512923203, -0.052293543460966, -0.022404618770124,
                                  -0.028466016418363, -0.035738959274497, -0.024154654366061,
                                  -0.037284072566944, -0.046035567570035, -0.018410076905551,
                                   0.025614061899556,  0.037491034524145,  0.017076041266371,
                                  -0.018410076905551, -0.030690378380024, -0.021942341812318,
                                   0.018591343287616,  0.014891233031041, -0.000031315350756,
                                   0.030646193885967,  0.046730232382521,  0.018662682291909,
                                  -0.028466016418363, -0.044186399359580, -0.000924194132083,
                                  -0.029462148373191, -0.064774775625169, -0.012362386065634,
                                   0.046730232382521,  0.090829068924062,  0.034214917535166,
                                   0.029148026137936,  0.047225275814966,  0.020222030749867};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);

    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(45);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000, -0.046936997644915,
                                   0.000000000000000, -0.081297264675733,  0.000000000000000,
                                   0.081297264675733,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.081297264675733,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.081297264675733,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.093873995289830,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.081297264675733,  0.000000000000000,
                                   0.081297264675733,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.081297264675733,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.046936997644915,  0.000000000000000,  0.081297264675733};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000, -0.026492295871231,
                                   0.000000000000000, -0.045886002458120,  0.000000000000000,
                                   0.045886002458120,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.045886002458120,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.045886002458120,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052984591742463,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.045886002458120,  0.000000000000000,
                                   0.045886002458120,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.045886002458120,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.026492295871231,  0.000000000000000,  0.045886002458120};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{ 0.002624812181476,  0.007218233499058,  0.015414307779995,
                                   0.011786525031334, -0.001428570495409, -0.000826509981420,
                                  -0.000669936538451,  0.003060568611434,  0.002975802882928,
                                   0.004122162565127,  0.001302416826301, -0.000826509981420,
                                   0.004407469156370, -0.003456993889820, -0.000403132530489,
                                  -0.000826509981420, -0.000669936538451,  0.003060568611434,
                                   0.002975802882928,  0.004122162565127,  0.009855640162612,
                                   0.016641456399239,  0.005522419697482,  0.011094304266159,
                                  -0.004106516734422, -0.003456993889820,  0.002975802882928,
                                   0.002040379074289, -0.003149772274224, -0.004226260988544,
                                   0.001302416826301, -0.000826509981420,  0.004407469156370,
                                  -0.003456993889820, -0.000403132530489, -0.003456993889820,
                                   0.002975802882928,  0.002040379074289, -0.003149772274224,
                                  -0.004226260988544,  0.001818547120498,  0.015462558629311,
                                   0.011741416816354,  0.003334003054246, -0.000422829213747};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{ 0.011465618674111,  0.031530451353806,  0.017818843320120,
                                   0.017707847964169, -0.011394222497393, -0.012418933901961,
                                  -0.036112308977243, -0.078331202752149, -0.048087225528750,
                                  -0.007178992217622, -0.006476216711201, -0.012418933901961,
                                  -0.020840229163434, -0.028017705880633, -0.007023718341617,
                                  -0.031913031136096, -0.073771360452277, -0.027684052466334,
                                  -0.037454081582858,  0.022062153633580, -0.021021291474060,
                                  -0.052839872146112, -0.032142095191320, -0.035226581430742,
                                   0.008758871447525, -0.028017705880633, -0.037454081582858,
                                  -0.018456034977556, -0.042559625799895, -0.001473541873380,
                                   0.002384736577043, -0.031913031136096, -0.020840229163434,
                                   0.001223439970569,  0.014242569550166,  0.001223439970569,
                                  -0.048087225528750, -0.052220801834766,  0.003960378963382,
                                   0.018020555360755,  0.018684469882661,  0.046413612769764,
                                   0.035185700956314,  0.034254861451545, -0.001168314401262};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(45);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.081297264675733,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.081297264675733,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.081297264675733,  0.000000000000000,
                                   0.081297264675733,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.046936997644915,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.093873995289830,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.046936997644915,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.081297264675733,
                                   0.000000000000000,  0.081297264675733,  0.000000000000000,
                                  -0.081297264675733,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.081297264675733};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.045886002458120,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.045886002458120,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.045886002458120,  0.000000000000000,
                                   0.045886002458120,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.026492295871231,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052984591742463,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.026492295871231,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.045886002458120,
                                   0.000000000000000,  0.045886002458120,  0.000000000000000,
                                  -0.045886002458120,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.045886002458120};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.003511002740106, -0.001134249900888,  0.001776621900595,
                                  -0.001134249900888,  0.013226325977804, -0.004671360736810,
                                   0.001776621900595, -0.004671360736810,  0.002431007872660,
                                   0.009655257535292, -0.000930110817433, -0.001134249900888,
                                  -0.000930110817433,  0.022277697548576,  0.003971296422528,
                                  -0.001134249900888,  0.003971296422528,  0.020779970540200,
                                   0.020692508303827,  0.004076141115479,  0.005920912150271,
                                   0.004076141115479,  0.007355335500733,  0.002717427410319,
                                   0.005920912150271,  0.002717427410319,  0.015758414845268,
                                   0.015833304283785,  0.003971296422528, -0.004671360736810,
                                   0.003971296422528,  0.014851798365717, -0.004239524502873,
                                  -0.004671360736810, -0.004239524502873,  0.004456847766543,
                                  -0.001943545868036,  0.005562356502454, -0.000539997433723,
                                   0.005562356502454, -0.005510969157418, -0.005688228258621,
                                  -0.000539997433723, -0.005688228258621, -0.000532291887283};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.007245173595825, -0.025600330115785, -0.012790197649571,
                                  -0.037602809085266, -0.031322204357537, -0.030816315771425,
                                  -0.007334525390716, -0.012812597317203,  0.016549863725452,
                                   0.019924227388520, -0.073881894516784, -0.025600330115785,
                                  -0.097068501616919, -0.088824338501496, -0.053508749477707,
                                  -0.037602809085266, -0.060055556188333,  0.048200401763262,
                                   0.018389962307879, -0.086540319777981, -0.024952918718845,
                                  -0.055356964690103, -0.074776195710166, -0.036904643126735,
                                  -0.024952918718845, -0.057693546518654,  0.039184061240250,
                                   0.015075002604525, -0.060055556188333, -0.030816315771425,
                                  -0.053508749477707, -0.059216225667664, -0.052477877052163,
                                  -0.012812597317203, -0.023835597693173,  0.030341416829996,
                                  -0.004122088967003,  0.005136227960260, -0.001894461645813,
                                   0.023139946414482,  0.013050918482307,  0.001631967627995,
                                   0.011199151775440,  0.013634446597476, -0.005792509916863};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 1);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(75);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{ 0.001866729396007, -0.001424911591865, -0.000856077155325,
                                  -0.000645283005878, -0.002519880644459, -0.000645283005878,
                                   0.001766913719470, -0.007953574392973,  0.006319137377598,
                                   0.000388489410652,  0.000729802203233, -0.000645283005878,
                                  -0.000322486455228, -0.000647932770561,  0.005010970906799,
                                  -0.001424911591865, -0.001533627811852, -0.003787738621143,
                                   0.000200538148039, -0.006102024893930,  0.000200538148039,
                                   0.003203718971876, -0.005396853218936, -0.004727071582407,
                                  -0.002090442874976,  0.003662249815556,  0.000200538148039,
                                   0.000780591058590, -0.005714513561804, -0.004654183363955,
                                  -0.006634742176309,  0.000883046450490,  0.004481905141829,
                                  -0.002130797425264,  0.006981629281098, -0.002130797425264,
                                  -0.006885867445698,  0.006693690715109, -0.004590578297132,
                                   0.000887832260527, -0.007667125584952, -0.002130797425264,
                                   0.002987936761219,  0.002658710971544, -0.003131094401058,
                                   0.007864477127728, -0.004727437581161,  0.000780591058590,
                                   0.001677450031995,  0.005431162185259,  0.001677450031995,
                                  -0.004727071582407, -0.003597902145957,  0.007142945290549,
                                   0.002311390456608, -0.000647932770561,  0.001677450031995,
                                  -0.004438231169968, -0.000104656667944,  0.006965710044941,
                                   0.004127460614037,  0.000306621583590,  0.000033997671900,
                                   0.000388489410652,  0.002518612131242,  0.000388489410652,
                                  -0.007419767417117,  0.003313989330405,  0.007392356143674,
                                   0.000125220991915, -0.004960040980945,  0.000388489410652,
                                   0.000618421070573,  0.000557063567421,  0.002193758814933};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.003695048403619, -0.008683436598962,  0.023851032369390,
                                   0.016447021137961,  0.016956987396888,  0.030051416254592,
                                   0.028529278562117, -0.009736603891965,  0.024195427657649,
                                  -0.005986842467724, -0.013245669260045,  0.008611171799369,
                                   0.015383473220745, -0.020657121534409, -0.016337384888921,
                                  -0.008683436598962, -0.024416885742051,  0.073875760448810,
                                   0.033813922198894,  0.041848125046758,  0.045923871176719,
                                   0.081253883471484, -0.011569761691248,  0.042204724376688,
                                  -0.028477488850199, -0.039223957890785,  0.053720323705891,
                                   0.054123499768272, -0.063788959707467, -0.042856490506439,
                                  -0.007624428040697, -0.024174235713382,  0.088101098577040,
                                   0.038265634218152,  0.028577083078407,  0.021979127208275,
                                   0.047678550960900,  0.007570235582316,  0.031785700640600,
                                  -0.009157969670115, -0.024621596689785,  0.038265634218152,
                                   0.058734065718027, -0.056062264228508, -0.027841358471430,
                                  -0.008499251191856, -0.018954359063501,  0.054123499768272,
                                   0.045223064446620,  0.030882805570992,  0.037802485292094,
                                   0.042204724376688, -0.007713174460832,  0.046083279824244,
                                  -0.001737221748523, -0.020657121534409,  0.015363462186125,
                                   0.028772843975250, -0.039849537902844, -0.029622892491187,
                                   0.001311210141786,  0.005457424036772, -0.009265550620744,
                                   0.003416176738587, -0.005319686790949, -0.005986842467724,
                                  -0.018615880601991,  0.004056918288319,  0.000011593610976,
                                   0.018177512026940,  0.007131280993733, -0.015389861674034,
                                  -0.007418349808397,  0.006767808518799,  0.007235845797000};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePDSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(75);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.001659944960387,  0.001284486363008,  0.006020064409068,
                                  -0.007094623447449, -0.003730064387623,  0.001284486363008,
                                   0.001405306587701, -0.000738134233049,  0.004315951124284,
                                  -0.000273379098250,  0.000766990028467,  0.003416663877628,
                                  -0.004060838362723, -0.000703280422647, -0.000005232864383,
                                   0.000588510565793, -0.000182946320051,  0.001943306669870,
                                  -0.001529914240791, -0.000354305882832,  0.002274398832089,
                                   0.005499745202119, -0.006320952062433, -0.004938329193726,
                                  -0.002252421037292,  0.000588510565793, -0.000182946320051,
                                   0.001943306669870, -0.001529914240791, -0.000354305882832,
                                  -0.001587715260959, -0.002856513795386,  0.006233512963422,
                                   0.004296086367078,  0.006722457127713,  0.007236480854902,
                                   0.004929654754877, -0.006027203592032,  0.003286436503251,
                                  -0.003015200356209, -0.005720714611650,  0.004296086367078,
                                   0.004155675308948, -0.006436585767951, -0.006707732565283,
                                  -0.000354305882832,  0.001906562669005, -0.000809711112446,
                                  -0.002108038286479, -0.000114196101824, -0.000641782127074,
                                  -0.003274604021704,  0.006946139266549,  0.000575874597345,
                                   0.004479099677402,  0.000588510565793, -0.000182946320051,
                                   0.001943306669870, -0.001529914240791, -0.000354305882832,
                                   0.000269521828610, -0.000703280422647, -0.002707225575149,
                                   0.004002730896500, -0.000583820233138,  0.000575874597345,
                                   0.005218431925711, -0.002357556457940,  0.000099874551325,
                                  -0.000501771301897, -0.004527595152165,  0.004229717428062,
                                   0.002824663143587, -0.006274895602266, -0.001966461871379};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.007626832387046, -0.014882852715975, -0.007962693872822,
                                  -0.011397750084366,  0.003926627634868, -0.014882852715975,
                                  -0.043142730913806, -0.027171355791747, -0.024077538084213,
                                   0.013843141086391,  0.016546675499264,  0.049525249446383,
                                   0.076043600118660,  0.044027673906266, -0.003862580405843,
                                   0.008227715352978,  0.017084984622416,  0.031679442008879,
                                   0.037358697096156,  0.007751792432683,  0.016255351132035,
                                   0.042380175286200,  0.025604791167467,  0.029139557462600,
                                  -0.007575025612254,  0.025555963408998,  0.032847692637063,
                                   0.010480388648267,  0.027906925429237, -0.004487486731631,
                                   0.019230018832600,  0.048336442310112,  0.016694678699195,
                                   0.023260784502889, -0.013709255455039, -0.014216075146506,
                                  -0.029630099393380, -0.018561510448276, -0.019753399595587,
                                   0.005923364644378,  0.017202126099776,  0.023260784502889,
                                   0.011129785799463,  0.028952455224371,  0.001377568870610,
                                  -0.004487486731631, -0.021497541371329, -0.004366828603445,
                                   0.000088618563314,  0.016655781391262, -0.015987780506338,
                                  -0.039420286912892, -0.022452310315566, -0.023857826179238,
                                   0.005161503510050, -0.001971683950617,  0.042857746222378,
                                   0.031679442008879, -0.001300445303786, -0.016726765895946,
                                   0.008698911061932,  0.044027673906266,  0.050695733412440,
                                   0.012835521191162, -0.008172347854081, -0.023857826179238,
                                  -0.063044229363729, -0.053570890799146, -0.046713768767477,
                                   0.002298808453281, -0.013331581735701, -0.038114530925862,
                                  -0.024075373782534, -0.026238225746528,  0.005853208506804};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(75);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.020772183026680,  0.001929972313023,  0.002448131187524,
                                  -0.005000767798569,  0.001361982586540,  0.001929972313023,
                                   0.025377798773666, -0.001507220350795, -0.002399882822444,
                                   0.008502576488598,  0.002448131187524, -0.001507220350795,
                                   0.051160624687932,  0.007184109044555, -0.004228092331997,
                                  -0.005000767798569, -0.002399882822444,  0.007184109044555,
                                   0.046184774218459,  0.006522452732311,  0.001361982586540,
                                   0.008502576488598, -0.004228092331997,  0.006522452732311,
                                   0.020815466868072,  0.064863179899687,  0.006726460814620,
                                  -0.005089787393139,  0.008077166383046, -0.000563884679997,
                                   0.006726460814620,  0.063493735127235,  0.001646409320578,
                                   0.006722992093001, -0.007473408931381, -0.005089787393139,
                                   0.001646409320578,  0.051243588526133,  0.001097606213719,
                                   0.002120744747141,  0.008077166383046,  0.006722992093001,
                                   0.001097606213719,  0.057891241716401,  0.003640589394999,
                                  -0.000563884679997, -0.007473408931381,  0.002120744747141,
                                   0.003640589394999,  0.063744808617692,  0.012396334149650,
                                  -0.010906048975492,  0.004099808743371,  0.000802202953028,
                                  -0.000531701244875, -0.010906048975492,  0.034278907678322,
                                   0.007184109044555, -0.008833411532689, -0.007650222092306,
                                   0.004099808743371,  0.007184109044555,  0.034107083125288,
                                  -0.007493977887924,  0.003103802862722,  0.000802202953028,
                                  -0.008833411532689, -0.007493977887924,  0.013429474316428,
                                  -0.009640982849453, -0.000531701244875, -0.007650222092306,
                                   0.003103802862722, -0.009640982849453,  0.015328765780184};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.009743741952946,  0.009615943221590, -0.038235458263433,
                                   0.032374613257963, -0.000889172389662,  0.009615943221590,
                                   0.013203395280394, -0.025825347213221, -0.020257326666947,
                                  -0.048241896238800, -0.051251345809949, -0.025852542948527,
                                   0.057077767853155, -0.019834708411824,  0.028011298388668,
                                   0.045452318283654, -0.033571558988817, -0.013277085678464,
                                   0.037741830358228,  0.049766568792694,  0.009654443606391,
                                  -0.039896747776045,  0.012362425956787,  0.038408433041647,
                                  -0.002011080542707, -0.020074951346617,  0.021596975509768,
                                  -0.087472503469334,  0.073483977429285,  0.006929824358499,
                                   0.020283305728016,  0.030331810234924, -0.051395773720411,
                                  -0.049287563825652, -0.077128401334973, -0.074357258002615,
                                  -0.048996349598114,  0.101160807550546, -0.032664233065409,
                                   0.030982190834423,  0.066350683778908, -0.049287563825652,
                                  -0.034263849146941,  0.071404780089633,  0.075369417681237,
                                   0.006929824358499, -0.085810531678674,  0.036446876445556,
                                   0.084599363895537, -0.006330799702260, -0.007019942240648,
                                   0.006142439355340, -0.022744909365640,  0.023475591938588,
                                   0.010309722416199,  0.009755031255158,  0.010749709109123,
                                  -0.013277085678464, -0.028740071788682, -0.035906920075696,
                                  -0.039288003079342, -0.019834708411824,  0.038051845235436,
                                  -0.009323619272007,  0.006385144831276,  0.023475591938588,
                                  -0.008768723305878, -0.014761109147834,  0.023213774649958,
                                   0.034453756588725, -0.005505701577880, -0.024548784324649,
                                   0.014830568048482,  0.026108608125971, -0.000816606089788};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(75);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.021267483486391, -0.066199040817495, -0.012704909745513,
                                  -0.001963080758393, -0.006829393125022,  0.011080041329827,
                                  -0.002468538535437,  0.005206165714148, -0.004155227678955,
                                   0.005072110669815, -0.008221854879514, -0.000825882852183,
                                  -0.001384279582048,  0.000568598953105,  0.000538281085893,
                                  -0.001963080758393, -0.006829393125022,  0.011080041329827,
                                  -0.025952108023466, -0.064799374223798, -0.034930877542230,
                                   0.001553065879871, -0.001638787503607, -0.007281083667641,
                                   0.002455590622269, -0.006796171859789,  0.008982673828701,
                                  -0.008622763171277,  0.007615184437108,  0.007778115745119,
                                  -0.002468538535437,  0.005206165714148, -0.004155227678955,
                                   0.001553065879871, -0.001638787503607, -0.007281083667641,
                                  -0.052149451593998, -0.052393729359908, -0.034766301062666,
                                  -0.007281083667641, -0.001092525002404,  0.007620635602906,
                                   0.004290953641695, -0.002169235714228, -0.003162249011663,
                                   0.005072110669815, -0.008221854879514, -0.000825882852183,
                                   0.002455590622269, -0.006796171859789,  0.008982673828701,
                                  -0.007281083667641, -0.001092525002404,  0.007620635602906,
                                  -0.047097028418047, -0.059135897673974, -0.013768546752112,
                                  -0.006640917838909, -0.003728633086060,  0.009784831342350,
                                  -0.001384279582048,  0.000568598953105,  0.000538281085893,
                                  -0.008622763171277,  0.007615184437108,  0.007778115745119,
                                   0.004290953641695, -0.002169235714228, -0.003162249011663,
                                  -0.006640917838909, -0.003728633086060,  0.009784831342350,
                                  -0.021315398461195, -0.065071319560502, -0.015683678219544};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.016020629265359, -0.026136481249791, -0.010101575275081,
                                   0.012133618926262,  0.019691870733530,  0.002563267817061,
                                  -0.018560251274544, -0.052735735615173, -0.007270436720172,
                                   0.017074754073036,  0.052004362066092,  0.019589282576925,
                                  -0.004803870605645,  0.003727831825332,  0.010743158985667,
                                   0.012133618926262,  0.017196649808857,  0.009425125359910,
                                   0.012934597717766,  0.016948187150825,  0.006760377420911,
                                  -0.013163647764530, -0.014632987071823, -0.004301904279261,
                                  -0.005230192586741, -0.019552951589058, -0.014726205746766,
                                  -0.035078965276513, -0.048767700801434, -0.024235453734342,
                                  -0.030649482764357, -0.040554219915515, -0.022635757659515,
                                  -0.007115538138768, -0.015717134922301, -0.010392662129090,
                                   0.031006743516636,  0.059898719548207,  0.020671162344424,
                                  -0.010392662129090, -0.010478089948200,  0.001545013635474,
                                   0.015634217713408,  0.016897591631464,  0.005136165849075,
                                   0.032727863326305,  0.043466302473400,  0.019589282576925,
                                  -0.016321451578172, -0.019552951589058,  0.001910682740380,
                                  -0.004301904279261, -0.009755324714549, -0.009578727531813,
                                   0.019896537751858,  0.033242313475040,  0.015127046225505,
                                   0.031691117385006,  0.044292685707699,  0.021861089650185,
                                   0.007539254492415,  0.003727831825332, -0.007771528661423,
                                  -0.023821481878454, -0.058744328855832, -0.013685951681839,
                                   0.001099454662679,  0.021973223172989,  0.012980323685977,
                                   0.021141615332502,  0.054140589881086,  0.010603606252126,
                                  -0.012248272923775, -0.018742948129549, -0.008744359517675};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 1);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(135);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.021925786952995, -0.060295914120735,  0.010088471820883,
                                  -0.002015901298229,  0.022284419317089, -0.002015901298229,
                                  -0.023931040801607, -0.079551706310175, -0.031709798029153,
                                   0.000418428726388, -0.019120367249015, -0.002015901298229,
                                   0.005858973583816, -0.059459056667960, -0.019932051003531,
                                  -0.001790694869328, -0.001937361386024, -0.004799281856853,
                                   0.000234454387705, -0.007716938206989,  0.000234454387705,
                                   0.003959840522694, -0.007029773780628, -0.006117237227855,
                                  -0.002713128737115,  0.004555502331618,  0.000234454387705,
                                   0.000929952468580, -0.007363618860255, -0.005937267282732,
                                   0.002335887733962, -0.001790694869328, -0.001102604698330,
                                  -0.000800081812181, -0.003177158109322, -0.000800081812181,
                                   0.002181939530904, -0.010235318615707,  0.008019530865637,
                                   0.000509004983356,  0.000914548779755, -0.000800081812181,
                                  -0.000373955508565, -0.000772684902615,  0.006331353072280,
                                  -0.001790694869328, -0.001937361386024, -0.004799281856853,
                                   0.000234454387705, -0.007716938206989,  0.000234454387705,
                                   0.003959840522694, -0.007029773780628, -0.006117237227855,
                                  -0.002713128737115,  0.004555502331618,  0.000234454387705,
                                   0.000929952468580, -0.007363618860255, -0.005937267282732,
                                  -0.031419759446477, -0.058394846579243,  0.019922152505636,
                                  -0.006342443656757,  0.040825954603481, -0.006342443656757,
                                  -0.048600425309427, -0.060529164048083, -0.032400283539618,
                                   0.002642684856982, -0.042280651041843, -0.006342443656757,
                                   0.013281435003757, -0.053109476865279, -0.023984644317073,
                                   0.009927395886691, -0.006120077781482,  0.000929952468580,
                                   0.002168179393525,  0.006955277249444,  0.002168179393525,
                                  -0.006117237227855, -0.004686515853752,  0.009057538212573,
                                   0.003019751032722, -0.000772684902615,  0.002168179393525,
                                  -0.005574242247336, -0.000080575716038,  0.008785013111965,
                                   0.002335887733962, -0.001790694869328, -0.001102604698330,
                                  -0.000800081812181, -0.003177158109322, -0.000800081812181,
                                   0.002181939530904, -0.010235318615707,  0.008019530865637,
                                   0.000509004983356,  0.000914548779755, -0.000800081812181,
                                  -0.000373955508565, -0.000772684902615,  0.006331353072280,
                                   0.009927395886691, -0.006120077781482,  0.000929952468580,
                                   0.002168179393525,  0.006955277249444,  0.002168179393525,
                                  -0.006117237227855, -0.004686515853752,  0.009057538212573,
                                   0.003019751032722, -0.000772684902615,  0.002168179393525,
                                  -0.005574242247336, -0.000080575716038,  0.008785013111965,
                                  -0.011504449743082, -0.059646728692044,  0.010068411914655,
                                  -0.000997891331517,  0.028635935214872, -0.000997891331517,
                                  -0.042818075474449, -0.071022274130419, -0.012789612821551,
                                   0.000837318202673, -0.031621141247290, -0.000997891331517,
                                   0.007578948906543, -0.057972092286699, -0.014412287652569};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.007978088037634, -0.021939742103493,  0.028557277596778,
                                   0.012835210371382,  0.018833176428124,  0.017189542839427,
                                   0.025571087613032, -0.012316333308695,  0.014332057096611,
                                  -0.016050100299783, -0.015449412344690,  0.017189542839427,
                                   0.017485409508314, -0.025049737308736, -0.018517365395966,
                                   0.009213122817927,  0.026880161779638, -0.094506560141099,
                                  -0.051235430905311, -0.025564231407124, -0.017952504393853,
                                  -0.028785145128803, -0.014237973051064, -0.038092758207429,
                                  -0.007350448358516,  0.016576715615068, -0.019670373314159,
                                  -0.055178317489489,  0.052111039579911,  0.022640112038602,
                                   0.004894300509586,  0.009213122817927, -0.025543678461311,
                                  -0.030786205462557, -0.014389666489641, -0.009563522525746,
                                  -0.010164832172045,  0.000119546844491, -0.024076571105945,
                                  -0.009286124820233,  0.005725403257960, -0.004839382994905,
                                  -0.012276021515376,  0.016266529481465,  0.013512181208431,
                                   0.005912379445485,  0.015383253066296, -0.046506868374932,
                                  -0.022881909291136, -0.034651673072937, -0.050556983381627,
                                  -0.070621444692108,  0.024279290378175, -0.034792014834987,
                                   0.034960273120636,  0.033065800285323, -0.038297221590658,
                                  -0.030127320529666,  0.043116993396135,  0.037103311876825,
                                   0.004737842471394,  0.011800853968007, -0.042698617720851,
                                  -0.026815594813488, -0.024403122956878, -0.031169927281532,
                                  -0.053366014848344, -0.001227133018692, -0.035577343232229,
                                   0.012987469700638,  0.020411655437697, -0.026815594813488,
                                  -0.028465745147234,  0.034147182979247,  0.025138676458323,
                                   0.006460602902961,  0.012325442778743, -0.030127320529666,
                                  -0.032817003157967, -0.026749161840845, -0.044393012032591,
                                  -0.034792014834987,  0.016186193585450, -0.041628098996286,
                                  -0.002345206720525,  0.016266529481465, -0.009694034708684,
                                  -0.021400767933544,  0.026674721529544,  0.025410461020041,
                                   0.001274165753206,  0.005912379445485, -0.021753425251700,
                                  -0.001649768061369, -0.009805172638179,  0.002489144006106,
                                  -0.027178437489124,  0.000119546844491,  0.001443836869674,
                                   0.019640274856212,  0.013521839523021, -0.031225200919058,
                                  -0.017961401329794,  0.019567272853906,  0.012314402344413,
                                   0.002024184764557,  0.013938685738163, -0.055178317489489,
                                   0.001697451001427, -0.012285962002622, -0.000560137399014,
                                  -0.038092758207429, -0.009491982034043,  0.002958820044055,
                                   0.022479378366855,  0.019567272853906, -0.045650135385300,
                                  -0.048524628899859,  0.038722115168203,  0.013572399877194,
                                  -0.007298155966902, -0.018841716080654,  0.046582230624644,
                                   0.027543692875406,  0.029392231084508,  0.027543692875406,
                                   0.038189170772658, -0.012415955679104,  0.028174781827182,
                                  -0.002588747914730, -0.026934985589673,  0.023189360407362,
                                   0.032607595972634, -0.049380806914401, -0.026612927564162};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePDPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(135);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.026754743072893,  0.002193080941073, -0.002827371813216,
                                   0.002193080941073,  0.038338413001516, -0.012140289864041,
                                  -0.002827371813216, -0.012140289864041,  0.013988077942710,
                                   0.073575543450456,  0.002406116796164,  0.002193080941073,
                                   0.002406116796164,  0.071180341180274,  0.007481031626830,
                                   0.002193080941073,  0.007481031626830,  0.072717508916345,
                                  -0.012404344848634,  0.005857454440357,  0.001325911740109,
                                   0.005857454440357, -0.024444408588934, -0.001169265259405,
                                   0.001325911740109, -0.001169265259405, -0.012355636541845,
                                   0.002400414179179, -0.000301649057475,  0.000998133430620,
                                  -0.000301649057475,  0.007652351286095, -0.002649830106117,
                                   0.000998133430620, -0.002649830106117,  0.001164664271670,
                                  -0.027213216621845,  0.009419806486986,  0.003889615185342,
                                   0.009419806486986, -0.049863290673548, -0.008512913851240,
                                   0.003889615185342, -0.008512913851240, -0.034916862522091,
                                   0.002400414179179, -0.000301649057475,  0.000998133430620,
                                  -0.000301649057475,  0.007652351286095, -0.002649830106117,
                                   0.000998133430620, -0.002649830106117,  0.001164664271670,
                                   0.029018047094006, -0.004835488576072, -0.002686821498719,
                                  -0.004835488576072,  0.059151699593642,  0.007445018047989,
                                  -0.002686821498719,  0.007445018047989,  0.052122139057125,
                                   0.096964884734856,  0.008571460983712,  0.012505343238772,
                                   0.008571460983712,  0.073814531166216,  0.005714307322474,
                                   0.012505343238772,  0.005714307322474,  0.086543765369213,
                                   0.038605433080894,  0.007445018047989, -0.009816242812295,
                                   0.007445018047989,  0.039434466395761, -0.011039670282729,
                                  -0.009816242812295, -0.011039670282729,  0.015488024353193,
                                  -0.000515406116266,  0.003309037062088, -0.000617874953755,
                                   0.003309037062088, -0.003188479702540, -0.003670929054662,
                                  -0.000617874953755, -0.003670929054662, -0.000970043238254,
                                   0.023289785984601, -0.005561090782803, -0.001093752813599,
                                  -0.005561090782803,  0.051627366557100,  0.000957331033564,
                                  -0.001093752813599,  0.000957331033564,  0.038605866785366,
                                   0.002400414179179, -0.000301649057475,  0.000998133430620,
                                  -0.000301649057475,  0.007652351286095, -0.002649830106117,
                                   0.000998133430620, -0.002649830106117,  0.001164664271670,
                                  -0.007204669610688, -0.001169265259405,  0.000440235484350,
                                  -0.001169265259405, -0.016296272392623,  0.006831842156528,
                                   0.000440235484350,  0.006831842156528, -0.009301984649632,
                                   0.072544731217924,  0.009024190920340,  0.000957331033564,
                                   0.009024190920340,  0.064803381775195,  0.000139173517506,
                                   0.000957331033564,  0.000139173517506,  0.070777422439837,
                                   0.024329543825460,  0.007277163943730, -0.007751381225268,
                                   0.007277163943730,  0.029322025152250, -0.010732478143085,
                                  -0.007751381225268, -0.010732478143085,  0.017611456434872};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.007419476462815,  0.014013103213930,  0.006557291167604,
                                   0.018037456513800,  0.015148837119548,  0.013917300719600,
                                   0.004953898811605,  0.008294845867021, -0.009416249150621,
                                  -0.020403560272742,  0.039997651109936,  0.014013103213930,
                                   0.047997829128629,  0.042904176530000,  0.028034075476810,
                                   0.018037456513800,  0.030078604514111, -0.027139559615450,
                                   0.021247872875056, -0.092416730028507, -0.026919084451569,
                                  -0.053656148476532, -0.073218463362023, -0.039579118782556,
                                  -0.020540047984947, -0.053724606295066,  0.038473213648322,
                                   0.009693077980652, -0.045937120849830, -0.025917539001453,
                                  -0.019131814638900, -0.035049238745457, -0.031620540383199,
                                  -0.005539684584680, -0.012131344759735,  0.016654077493448,
                                   0.013671186554897, -0.038426854223489, -0.017834142321211,
                                  -0.049128998377333, -0.042758287761381, -0.033183530457618,
                                  -0.017271157819064, -0.025016629434452,  0.026475293832424,
                                   0.006805826781275, -0.023582469051294, -0.013579914244202,
                                  -0.045204896012391, -0.032161987546080, -0.038350153887154,
                                  -0.006774163391564, -0.011457992683253,  0.016654077493448,
                                   0.006171179874957, -0.026232495024118, -0.010525987511004,
                                  -0.055571011245057, -0.045091222170627, -0.027841962312803,
                                  -0.018999437064144, -0.023817609012932,  0.019578866871535,
                                  -0.012479954723146,  0.027985733095402,  0.012279789731618,
                                   0.040431846107582,  0.025021666631831,  0.026954564071721,
                                   0.012279789731618,  0.018657155396934, -0.022713112832828,
                                   0.005134513846114, -0.023817609012932, -0.016622059493853,
                                  -0.027841962312803, -0.030060814780418, -0.032369375984389,
                                  -0.003911885164144, -0.006384487513341,  0.012032183984880,
                                  -0.004476623195588,  0.004283279415147, -0.003242775667080,
                                   0.029512996847725,  0.013400828144200, -0.000037204477528,
                                   0.013091026379252,  0.013088287652361, -0.005298336918879,
                                  -0.010069997671269,  0.033245474859136,  0.013076161996466,
                                   0.045673680564353,  0.038134818819368,  0.027885707225973,
                                   0.015400407725706,  0.023861353926102, -0.026639661492860,
                                   0.006805826781275, -0.027553401730400, -0.002659849376660,
                                  -0.038071380320841, -0.035049238745457, -0.003211191860287,
                                  -0.026391331543326, -0.039706923438880,  0.019541328692824,
                                   0.011812050667725, -0.053724606295066, -0.010054934330686,
                                  -0.039579118782556, -0.048812308908015, -0.020673549491068,
                                  -0.019623489030620, -0.047646224782619,  0.028002007014527,
                                  -0.020659703379699,  0.074316664607032,  0.027885707225973,
                                   0.067537643313278,  0.072111875484547,  0.048989572817794,
                                   0.023861353926102,  0.046130939297201, -0.048839379403576,
                                  -0.011067664601591,  0.034864880146625,  0.017065813420877,
                                   0.043031781169790,  0.038779900600951,  0.029994346819094,
                                   0.016027313936236,  0.022645830415143, -0.023334036152606};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSDDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(125);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000, -0.060312414940408,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.030156207470204,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.030156207470204,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.060312414940408,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000, -0.033169580854008,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000, -0.028725699652453,  0.000000000000000,
                                   0.016584790427004,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.016584790427004,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.028725699652452,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.033169580854008,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652452,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.001318320086298, -0.001238651931187, -0.001698438054261,
                                  -0.000369291832801, -0.000658993963283, -0.001238651931187,
                                   0.003682356652316, -0.000097625158607,  0.001747588963362,
                                   0.002165336210181, -0.001698438054261, -0.000097625158607,
                                   0.011579702632372, -0.003475395063780, -0.000091523228930,
                                  -0.000369291832801,  0.001747588963362, -0.003475395063780,
                                   0.002051061046083, -0.004546128072753, -0.000658993963283,
                                   0.002165336210181, -0.000091523228930, -0.004546128072753,
                                  -0.001554875501769,  0.007042549020961,  0.000671010851535,
                                   0.001629128018547, -0.000631475768993,  0.001168394798527,
                                   0.000671010851535,  0.006709312010227,  0.001342677349037,
                                  -0.000709836954478,  0.002617992043414,  0.001629128018547,
                                   0.001342677349037,  0.010579583425925,  0.001484338551615,
                                  -0.003920806700328, -0.000631475768993, -0.000709836954478,
                                   0.001484338551615,  0.015580084762442,  0.003093447539126,
                                   0.001168394798527,  0.002617992043414, -0.003920806700328,
                                   0.003093447539126,  0.007049024297520,  0.013295455052255,
                                   0.002829561102694, -0.003318559278851,  0.002403784122053,
                                  -0.000632759517615,  0.002829561102694,  0.007686533462897,
                                   0.000262137782030,  0.004645272704577, -0.001851616862456,
                                  -0.003318559278851,  0.000262137782030,  0.002187895696563,
                                   0.000174758521354,  0.001382733032855,  0.002403784122053,
                                   0.004645272704577,  0.000174758521354,  0.003815472875750,
                                   0.000007372886979, -0.000632759517615, -0.001851616862456,
                                   0.001382733032855,  0.000007372886979,  0.012040482008987,
                                   0.003432852758646, -0.003498451574653,  0.003579934083596,
                                  -0.000022836003797, -0.000472731890294, -0.003498451574653,
                                   0.012115251864392,  0.001484338551615, -0.003657548465863,
                                  -0.003787294394457,  0.003579934083596,  0.001484338551615,
                                   0.007053055617283,  0.000105728556024,  0.003371365837402,
                                  -0.000022836003797, -0.003657548465863,  0.000105728556024,
                                   0.002744345984054, -0.002888064896638, -0.000472731890294,
                                  -0.003787294394457,  0.003371365837402, -0.002888064896638,
                                   0.005961529453675,  0.000562788172959, -0.003781455512760,
                                  -0.000091523228930,  0.004374059511658,  0.001377401625210,
                                  -0.003781455512760, -0.001324349507682,  0.004474082261541,
                                  -0.000815647803116, -0.002005571874122, -0.000091523228930,
                                   0.004474082261541, -0.004824876096822, -0.005202026391607,
                                  -0.001879959124973,  0.004374059511658, -0.000815647803116,
                                  -0.005202026391607, -0.001064574533318, -0.001793524037405,
                                   0.001377401625210, -0.002005571874122, -0.001879959124973,
                                  -0.001793524037405, -0.000464223416512};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.004755874418206, -0.004284752363439,  0.015312754723014,
                                  -0.014872833345269, -0.004290645825383, -0.005579959792007,
                                  -0.006636563626017,  0.010440066731281,  0.013757689361716,
                                   0.025147998718574,  0.031862077305011,  0.014278470979734,
                                  -0.028234265596005,  0.008827989931946, -0.009027992244198,
                                  -0.019640130004326,  0.006953780790440,  0.009832006872426,
                                  -0.016567502539779, -0.024308187796424,  0.001093226786622,
                                   0.014546543563305, -0.007798606187260, -0.014261742623539,
                                   0.002029833238983,  0.013321430040069, -0.013509796292227,
                                   0.052964536227694, -0.044831150652932,  0.000596484807674,
                                  -0.013038811772747, -0.018493325361549,  0.033544069980170,
                                   0.027164921770946,  0.059579755673781,  0.059613465818968,
                                   0.032592934257885, -0.070844633672512,  0.025962202335928,
                                  -0.030456856030001, -0.055507310227399,  0.042892706571949,
                                   0.019388705679125, -0.050250803188857, -0.060083290397878,
                                  -0.011627627910836,  0.055057024870016, -0.016088406489536,
                                  -0.052983279668713,  0.003472615895169,  0.012163794603356,
                                  -0.013850851438638,  0.060504611973329, -0.047669038993089,
                                  -0.003826884863297, -0.011040326095817, -0.018186518692894,
                                   0.032217040486581,  0.027832827254096,  0.038335914068736,
                                   0.034325265285145,  0.026961532362841, -0.059217183794802,
                                   0.017974354908561, -0.014302193868811, -0.033310816009118,
                                   0.027832827254096,  0.021478026991054, -0.041380541404641,
                                  -0.036724160622752, -0.003826884863297,  0.055736867543628,
                                  -0.025210254988887, -0.055086419585840,  0.004573806291150,
                                   0.009476530655851, -0.008569884742798,  0.030709575224614,
                                  -0.032449317445330, -0.012254518088263, -0.012442135745988,
                                  -0.014698236821537,  0.019388705679125,  0.036135611901372,
                                   0.043273128357508,  0.044063781268105,  0.025962202335928,
                                  -0.047229755781675,  0.010957765644611, -0.009933040953564,
                                  -0.029848973813117,  0.012543934699867,  0.017386815247566,
                                  -0.031131182212067, -0.041709845013005,  0.006081650989502,
                                   0.034043758515610, -0.021765881751594, -0.035907394146034,
                                   0.001719499967641, -0.002297400874089,  0.003429704184063,
                                  -0.007798606187260,  0.002413516445326, -0.006107673146725,
                                   0.001807046620821,  0.004326626464610, -0.008083578905077,
                                   0.003199220828651, -0.007996693589450, -0.009027992244198,
                                  -0.005051074970650,  0.011764277331669, -0.005025761044670,
                                   0.013956559354019,  0.010703240350206, -0.013130159742412,
                                   0.001503657120388,  0.005341734437805,  0.010589027750021,
                                   0.006813621122088, -0.008458207428887, -0.000154480881719,
                                   0.006377526817052, -0.000529977316406};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDSDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(125);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000, -0.060312414940408,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.030156207470204,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.060312414940408,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.030156207470204,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.030156207470204,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.052232083501982,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.060312414940408,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.052232083501982,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.060312414940408,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000, -0.033169580854008,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000, -0.028725699652452,  0.000000000000000,
                                   0.016584790427004,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.033169580854008,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.016584790427004,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.016584790427004,
                                   0.000000000000000,  0.028725699652452,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028725699652453,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.033169580854008,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028725699652453,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.033169580854008,
                                   0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{ 0.001124495835953,  0.005896954587924,  0.011081345963470,
                                   0.002882720098173,  0.000455445002010, -0.001007778744240,
                                   0.000574995488572,  0.002366592765555, -0.002899232135772,
                                  -0.003132012523433, -0.001387948519228,  0.001333208670836,
                                  -0.002754723153985,  0.002949822546561, -0.000074934487686,
                                  -0.000289170197979, -0.000502827091787,  0.002020874541595,
                                  -0.000001378136381,  0.003614424412979, -0.000547514859819,
                                   0.000970683823787, -0.000525666150785, -0.000392722219915,
                                   0.001144416250801, -0.001007778744240,  0.000574995488572,
                                   0.002366592765555, -0.002899232135772, -0.003132012523433,
                                   0.003085807609359,  0.005681379886685,  0.006458882784355,
                                   0.010113791532373, -0.001089280389109, -0.000071114312372,
                                   0.001128106490818,  0.000231718277944,  0.001243114227618,
                                   0.003705756152897,  0.001432003737767, -0.000572241371555,
                                   0.003857108400801, -0.003017412973224, -0.000678531883010,
                                   0.001781470076171,  0.002153096510151, -0.001562170159573,
                                  -0.003163872865258, -0.001646439163521, -0.001387948519228,
                                   0.001333208670836, -0.002754723153985,  0.002949822546561,
                                  -0.000074934487686, -0.000071114312372,  0.001128106490818,
                                   0.000231718277944,  0.001243114227618,  0.003705756152897,
                                   0.009644686720475,  0.008890543734429,  0.001904900285185,
                                   0.005927029156286, -0.004018619466865, -0.002875198148819,
                                   0.001243114227618,  0.000154478851963,  0.000092177967803,
                                  -0.004316188555438, -0.000074934487686, -0.003234825408653,
                                   0.001147801314160,  0.002789889965973, -0.001536568586471,
                                  -0.000289170197979, -0.000502827091787,  0.002020874541595,
                                  -0.000001378136381,  0.003614424412979,  0.001432003737767,
                                  -0.000572241371555,  0.003857108400801, -0.003017412973224,
                                  -0.000678531883010, -0.002875198148819,  0.001243114227618,
                                   0.000154478851963,  0.000092177967803, -0.004316188555438,
                                   0.001728743843339,  0.013011636382667,  0.003244625783687,
                                   0.002348219313862, -0.000916782716181, -0.003755799486427,
                                   0.002587499240306,  0.000022103701890, -0.002392103740230,
                                  -0.001480308407195, -0.000547514859819,  0.000970683823787,
                                  -0.000525666150785, -0.000392722219915,  0.001144416250801,
                                   0.001781470076171,  0.002153096510151, -0.001562170159573,
                                  -0.003163872865258, -0.001646439163521, -0.000074934487686,
                                  -0.003234825408653,  0.001147801314160,  0.002789889965973,
                                  -0.001536568586471, -0.003755799486427,  0.002587499240306,
                                   0.000022103701890, -0.002392103740230, -0.001480308407195,
                                  -0.001262687417588,  0.005902363535023,  0.010038774764413,
                                   0.004983491983792, -0.000397865176329};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{ 0.007802378710029,  0.019692948724320,  0.012156804530894,
                                   0.012297747344655, -0.004771825856929, -0.001263508628851,
                                  -0.016233131322256, -0.010753716189679, -0.002732813349670,
                                   0.007881798521074,  0.003454150962420,  0.024894698479164,
                                   0.044041955216727,  0.006800981720654, -0.004636992723075,
                                  -0.008254016717696, -0.021859575751596, -0.035273373329677,
                                  -0.027451630557633, -0.006709335485382, -0.004337487338856,
                                   0.006856484972423, -0.002564580549815, -0.015042553902761,
                                  -0.011221645542799, -0.008107079663225, -0.013744560037029,
                                  -0.005112506967926, -0.010505089135740,  0.004624844858340,
                                  -0.006124030090421, -0.015077490020398, -0.003727055287844,
                                  -0.007696126661651,  0.008391439574697,  0.006095788869125,
                                   0.014435660042533,  0.002955765109127,  0.005192646942669,
                                  -0.007571007505037,  0.004040954475720,  0.006301269867245,
                                   0.010026730940551,  0.019957061300607,  0.006609157950936,
                                   0.017710040496250,  0.039983513012179,  0.016650309769945,
                                   0.027766511244092, -0.008148958839924,  0.020647325705870,
                                   0.032914709393680,  0.013675688316887,  0.022908969743877,
                                  -0.006119891505420,  0.008474590561888,  0.006170831824811,
                                   0.005796851652040,  0.014349555162148, -0.000536483750415,
                                  -0.014440370819586, -0.034827036598169, -0.036538461395893,
                                  -0.023218024398779,  0.006016821174827,  0.003346191925588,
                                   0.014349555162148,  0.003864567768027, -0.005787130810313,
                                  -0.005886140111219, -0.006119891505420, -0.014970041439884,
                                  -0.005698203465370, -0.007662035271173,  0.008509540886787,
                                  -0.017525233468702, -0.037248033173264, -0.014577644542702,
                                  -0.022394603238903,  0.005239349100854, -0.003575227159878,
                                   0.021751442438712,  0.010026730940551, -0.003218197556592,
                                  -0.011669677974499,  0.007933932689341,  0.005192646942669,
                                   0.001970510072752,  0.010108454256976,  0.004240837927127,
                                  -0.011184550113984, -0.024214292349936, -0.012082664404969,
                                  -0.018498394918572, -0.001179531156195, -0.014388444356011,
                                  -0.036416554445966, -0.015706129240505, -0.024925528313875,
                                   0.006687874466603,  0.003842941652129, -0.010487745294319,
                                  -0.002564580549815,  0.010973791497351,  0.008411384035564,
                                   0.003622763892243,  0.032522198350693,  0.041016213043225,
                                   0.017342865264623, -0.003424364033191, -0.004636992723075,
                                   0.002361338078774, -0.018350814673636, -0.021934936051244,
                                  -0.005742551271679, -0.000320783542135, -0.028561364500001,
                                  -0.040106008058923, -0.014320214234073,  0.005128768276335,
                                   0.009453796478940,  0.014570736697873,  0.007070386440429,
                                   0.010544709603474, -0.002418247138475};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(225);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.062086251738290,  0.000000000000000,  0.058946090578727,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.028053710474980,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.058946090578727,  0.000000000000000,
                                   1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.029473045289363,
                                   0.000000000000000, -0.051048811894956,  0.000000000000000,
                                   0.029473045289363,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.058946090578727,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                  -0.058946090578727,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                  -0.051048811894956,  0.000000000000000,  0.029473045289363,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.028053710474980,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.130151334264890,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.164183875528200,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.130151334264890,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.028053710474980,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.029473045289363,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.058946090578727,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000, -0.058946090578727,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.029473045289363,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.029473045289363,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  1.130151334264890,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.028053710474980,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.062086251738290,
                                   0.000000000000000, -0.058946090578727,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.130151334264890,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.058946090578727,  0.000000000000000,  1.130151334264890};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.764629737666849,  0.000000000000000,  0.042611711291179,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.740027854682256,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.042611711291179,  0.000000000000000,
                                   0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.021305855645589,
                                   0.000000000000000, -0.036902824476889,  0.000000000000000,
                                   0.021305855645589,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.042611711291179,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                  -0.042611711291179,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                  -0.036902824476889,  0.000000000000000,  0.021305855645589,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.740027854682256,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.813833503636034,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.838435386620627,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.813833503636034,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.740027854682256,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.021305855645589,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.042611711291179,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000, -0.042611711291179,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.021305855645589,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.021305855645589,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.813833503636034,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.740027854682256,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.764629737666849,
                                   0.000000000000000, -0.042611711291179,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.813833503636034,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.042611711291179,  0.000000000000000,  0.813833503636034};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.355787127971241,  0.014096087696362, -0.003304835579206,
                                   0.009015464262110, -0.001424717166913,  0.014096087696362,
                                   0.389425519064833,  0.005271269251748,  0.005008379166556,
                                  -0.012208472136504, -0.003304835579206,  0.005271269251748,
                                   0.402703234758725,  0.014087063301909,  0.003833197675669,
                                   0.009015464262110,  0.005008379166556,  0.014087063301909,
                                   0.386048474718214,  0.012546349780913, -0.001424717166913,
                                  -0.012208472136504,  0.003833197675669,  0.012546349780913,
                                   0.357780576658355,  0.008822829523671,  0.000851880856359,
                                   0.002071773558931, -0.000759473237754,  0.001525624417867,
                                   0.000851880856359,  0.008425549201863,  0.001776839210863,
                                  -0.000858229112236,  0.003386080702174,  0.002071773558931,
                                   0.001776839210863,  0.013496681317875,  0.001982486770580,
                                  -0.004948021886145, -0.000759473237754, -0.000858229112236,
                                   0.001982486770580,  0.019854525601399,  0.004019867719567,
                                   0.001525624417867,  0.003386080702174, -0.004948021886145,
                                   0.004019867719567,  0.008821171119139,  0.001647311512457,
                                  -0.001552770187619, -0.002136574094252, -0.000465765846908,
                                  -0.000849626372227, -0.001552770187619,  0.004624826205918,
                                  -0.000054540434785,  0.002170093476987,  0.002753252374520,
                                  -0.002136574094252, -0.000054540434785,  0.014896668802178,
                                  -0.004441728540167, -0.000133170315312, -0.000465765846908,
                                   0.002170093476987, -0.004441728540167,  0.002484496242244,
                                  -0.005798719687028, -0.000849626372227,  0.002753252374520,
                                  -0.000133170315312, -0.005798719687028, -0.002045754848635,
                                   0.008822829523671,  0.000851880856359,  0.002071773558931,
                                  -0.000759473237754,  0.001525624417867,  0.000851880856359,
                                   0.008425549201863,  0.001776839210863, -0.000858229112236,
                                   0.003386080702174,  0.002071773558931,  0.001776839210863,
                                   0.013496681317875,  0.001982486770580, -0.004948021886145,
                                  -0.000759473237754, -0.000858229112236,  0.001982486770580,
                                   0.019854525601399,  0.004019867719567,  0.001525624417867,
                                   0.003386080702174, -0.004948021886145,  0.004019867719567,
                                   0.008821171119139,  0.385800384524770,  0.015744581663846,
                                  -0.010743514053092,  0.020047160339838, -0.001100022502298,
                                   0.015744581663846,  0.404928436642160,  0.011543404853431,
                                   0.014352313343316, -0.018976246693059, -0.010743514053092,
                                   0.011543404853431,  0.400955995449154,  0.007695603235620,
                                   0.004476464188788,  0.020047160339838,  0.014352313343316,
                                   0.007695603235620,  0.392968175522730,  0.010271023024752,
                                  -0.001100022502298, -0.018976246693059,  0.004476464188788,
                                   0.010271023024752,  0.383618673228545,  0.004228278611307,
                                  -0.004495934879233,  0.004523323067434, -0.000075891644055,
                                  -0.000622744728080, -0.004495934879233,  0.015506680929862,
                                   0.001982486770580, -0.004692839461750, -0.004947640219981,
                                   0.004523323067434,  0.001982486770580,  0.008997787545250,
                                   0.000124766902047,  0.004242456410122, -0.000075891644055,
                                  -0.004692839461750,  0.000124766902047,  0.003346702272312,
                                  -0.003709409540634, -0.000622744728080, -0.004947640219981,
                                   0.004242456410122, -0.003709409540634,  0.007534388483899,
                                   0.001647311512457, -0.001552770187619, -0.002136574094252,
                                  -0.000465765846908, -0.000849626372227, -0.001552770187619,
                                   0.004624826205918, -0.000054540434785,  0.002170093476987,
                                   0.002753252374520, -0.002136574094252, -0.000054540434785,
                                   0.014896668802178, -0.004441728540167, -0.000133170315312,
                                  -0.000465765846908,  0.002170093476987, -0.004441728540167,
                                   0.002484496242244, -0.005798719687028, -0.000849626372227,
                                   0.002753252374520, -0.000133170315312, -0.005798719687028,
                                  -0.002045754848635,  0.004228278611307, -0.004495934879233,
                                   0.004523323067434, -0.000075891644055, -0.000622744728080,
                                  -0.004495934879233,  0.015506680929862,  0.001982486770580,
                                  -0.004692839461750, -0.004947640219981,  0.004523323067434,
                                   0.001982486770580,  0.008997787545250,  0.000124766902047,
                                   0.004242456410122, -0.000075891644055, -0.004692839461750,
                                   0.000124766902047,  0.003346702272312, -0.003709409540634,
                                  -0.000622744728080, -0.004947640219981,  0.004242456410122,
                                  -0.003709409540634,  0.007534388483899,  0.357265929167363,
                                   0.004466179345290, -0.003571176209830,  0.020122644538223,
                                   0.002111102593197,  0.004466179345290,  0.386368102185412,
                                   0.016770677064874,  0.002868049202882, -0.017318049135186,
                                  -0.003571176209830,  0.016770677064874,  0.390289344090244,
                                   0.000607567575839, -0.000968192763571,  0.020122644538223,
                                   0.002868049202882,  0.000607567575839,  0.383181456224166,
                                   0.008001136558783,  0.002111102593197, -0.017318049135186,
                                  -0.000968192763571,  0.008001136558783,  0.356633811575714};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.009220754943649, -0.011432853295081,  0.021552408815691,
                                  -0.019186800553326,  0.000613963214196, -0.011432853295081,
                                  -0.018062190419611,  0.013374660091478,  0.010468208237208,
                                   0.028185580071319,  0.025589403177293,  0.013294917278293,
                                  -0.034120498690506,  0.007209897120979, -0.015493279304807,
                                  -0.023490828934847,  0.014037608324727,  0.004725889957902,
                                  -0.020025060926764, -0.027810313200183, -0.002223261619513,
                                   0.026838840712095, -0.009565534938373, -0.023507886587368,
                                   0.004421360205845, -0.012738189067186,  0.014861928347702,
                                  -0.064246361176305,  0.053223738495087,  0.012664600904099,
                                   0.013508769664702,  0.021600063247586, -0.034152544965954,
                                  -0.040673178999218, -0.036451417458714, -0.031770987928699,
                                  -0.031308672575964,  0.059160863439080, -0.014251477042897,
                                   0.004381959063446,  0.024671538715550, -0.015824985217469,
                                  -0.027386111636758,  0.041294274744880,  0.032700646244084,
                                  -0.006508375382681, -0.057829832286180,  0.036267033540861,
                                   0.059354706323047, -0.004992235558486, -0.004319346763401,
                                   0.003670849730992, -0.014399933327532,  0.016333817968106,
                                   0.010438624317984,  0.007392036109243,  0.007541846465364,
                                  -0.006278460867271, -0.026642561482567, -0.020803181943679,
                                  -0.016059997715934, -0.011032630944801,  0.019755255572317,
                                  -0.002021350263906, -0.003232721649334,  0.008749012724606,
                                  -0.003901099435732, -0.009384090085444,  0.013206153115487,
                                   0.019607026404353, -0.007329267897489, -0.018481671369885,
                                   0.012849537634112,  0.021623100393421,  0.001752374085096,
                                  -0.009197728189581,  0.008531070642492, -0.035079759666274,
                                   0.029722914374836, -0.001612442523678,  0.009298574902557,
                                   0.011856676582506, -0.027157644435611, -0.020252716319250,
                                  -0.055206709545562, -0.064501324484143, -0.029300658562923,
                                   0.064676020756078, -0.020773279620826,  0.038221210462059,
                                   0.052732047299384, -0.036406958370779, -0.010970593454627,
                                   0.037753813867274,  0.059150725718664,  0.011283929116419,
                                  -0.036954533911775,  0.012935467760772,  0.035577065601776,
                                  -0.000764006210404, -0.006668655323208,  0.006482667575392,
                                  -0.031635745007189,  0.027001173999614,  0.003421842981549,
                                   0.007559019369777,  0.010197886400051, -0.027113949892813,
                                  -0.026931887939032, -0.037295204379076, -0.039145085022275,
                                  -0.029207433876920,  0.039503674861469, -0.019471622584613,
                                   0.016310452092615,  0.031305202381135, -0.026931887939032,
                                  -0.018075966595209,  0.032641126349244,  0.038174585137030,
                                   0.003421842981549, -0.032184325790572,  0.013181560419662,
                                   0.032974332284732,  0.000117999923532, -0.006763863704456,
                                   0.005548964470981, -0.022093353546612,  0.021667722076793,
                                   0.012269329723846,  0.010317406563285,  0.010017318283240,
                                  -0.010970593454627, -0.031882601113755, -0.042644412047973,
                                  -0.051728312676864, -0.020773279620826,  0.043117347170719,
                                  -0.011989592212235,  0.004534975724860,  0.025804888573248,
                                  -0.007651238036461, -0.018015483223422,  0.023056342016613,
                                   0.040580266430643, -0.007075227736299, -0.022440414167475,
                                   0.011727212128018,  0.024371703612073,  0.000122707437799,
                                  -0.003485688011308,  0.004128762368427, -0.019090345278962,
                                   0.013039890574670, -0.007684733315710,  0.002018125653249,
                                   0.004452578336009, -0.012688779682854, -0.001126654345368,
                                  -0.014539880844525, -0.007583797226209, -0.006745097170279,
                                   0.019755255572317, -0.008452650925689,  0.017110159526006,
                                   0.019905395133589, -0.026842616465135,  0.000231388137931,
                                   0.016295421244843,  0.017641314895604,  0.012447512927854,
                                  -0.022061644270617,  0.001592548950679,  0.020736119495040,
                                   0.000918715333004, -0.007595796577357,  0.010609042776420,
                                  -0.050136843130922,  0.031296071268995, -0.011438796986073,
                                   0.003272033902060,  0.010849251156142, -0.027386111636758,
                                   0.001282920184464, -0.018658165308292, -0.014368387408739,
                                  -0.014251477042897,  0.039440575626053, -0.019432441706883,
                                   0.019270756947242,  0.027551250600494, -0.035989370488159,
                                  -0.011330785268656,  0.031080307505503,  0.019343075681549,
                                   0.017320667444097, -0.042920563401754,  0.006643776728450,
                                   0.038597332445890, -0.004224486506425,  0.010577286404001,
                                  -0.010054095490820,  0.035994495400482, -0.039710690832008,
                                  -0.005743059935737, -0.011130447285206, -0.014811981026530,
                                   0.016737914375689,  0.022791183104205,  0.043822910364462,
                                   0.039466841053966,  0.018911141172981, -0.050583211667437,
                                   0.014260808513204, -0.011613489124884, -0.039710690832008,
                                   0.019221783016686,  0.015349159686876, -0.040565269970555,
                                  -0.039766393021289, -0.002905835102028,  0.040058771135183,
                                  -0.014412341818366, -0.038868566781806,  0.005203972644081};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(225);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.028053710474980,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.130151334264890,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.058946090578727,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.058946090578727,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.130151334264890,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.028053710474980,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                   0.029473045289363,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.051048811894956,  0.000000000000000,
                                  -0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.058946090578727,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.058946090578727,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                   0.029473045289363,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.062086251738290,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.164183875528200,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.062086251738290,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.029473045289363,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                   0.058946090578727,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.058946090578727,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.029473045289363,
                                   0.000000000000000,  0.029473045289363,  0.000000000000000,
                                   1.028053710474980,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.130151334264890,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.130151334264890,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.051048811894956,  0.000000000000000,
                                  -0.051048811894956,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.058946090578727,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.058946090578727,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.051048811894956,
                                   0.000000000000000,  0.051048811894956,  0.000000000000000,
                                   1.130151334264890,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.028053710474980,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.130151334264890};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.740027854682256,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.813833503636034,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.042611711291179,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.042611711291179,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.813833503636034,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.740027854682256,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                   0.021305855645589,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.036902824476889,  0.000000000000000,
                                  -0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.042611711291179,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.042611711291179,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                   0.021305855645589,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.764629737666849,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.838435386620627,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.764629737666849,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.021305855645589,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                   0.042611711291179,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.042611711291179,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.021305855645589,
                                   0.000000000000000,  0.021305855645589,  0.000000000000000,
                                   0.740027854682256,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.813833503636034,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.813833503636034,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.036902824476889,  0.000000000000000,
                                  -0.036902824476889,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.042611711291179,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.042611711291179,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.036902824476889,
                                   0.000000000000000,  0.036902824476889,  0.000000000000000,
                                   0.813833503636034,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.740027854682256,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.813833503636034};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.396946505248092,  0.010052818861931,  0.001914666694555,
                                   0.010052818861931,  0.430704096024368,  0.004873321147611,
                                   0.001914666694555,  0.004873321147611,  0.398538452903950,
                                   0.015598694959847,  0.000987518134931, -0.001718117258579,
                                   0.000987518134931,  0.017387301064423, -0.005006588876445,
                                  -0.001718117258579, -0.005006588876445,  0.004831509652861,
                                  -0.003699760158805,  0.002289729070947, -0.002370719309518,
                                   0.002289729070947, -0.012032097539926,  0.005050336240100,
                                  -0.002370719309518,  0.005050336240100, -0.003975987070969,
                                   0.009884363501320, -0.000838581508292, -0.000494114807744,
                                  -0.000838581508292,  0.022205252327584, -0.000028124084560,
                                  -0.000494114807744, -0.000028124084560,  0.022291576962577,
                                  -0.001586580633917,  0.001689840021891, -0.000947058276240,
                                   0.001689840021891, -0.001203370256285, -0.000686760094969,
                                  -0.000947058276240, -0.000686760094969,  0.002363575303981,
                                   0.015598694959847,  0.000987518134931, -0.001718117258579,
                                   0.000987518134931,  0.017387301064423, -0.005006588876445,
                                  -0.001718117258579, -0.005006588876445,  0.004831509652861,
                                   0.434170663674999,  0.009699125539665,  0.005267858542390,
                                   0.009699125539665,  0.451660392273576,  0.017429112559173,
                                   0.005267858542390,  0.017429112559173,  0.430604318593715,
                                   0.005629146676860,  0.001977315841567, -0.000082776369306,
                                   0.001977315841567,  0.012687891246471,  0.002189120352974,
                                  -0.000082776369306,  0.002189120352974,  0.018427959465896,
                                   0.005408465728369, -0.000957561052355,  0.002430179826511,
                                  -0.000957561052355,  0.015803281683148, -0.005206302729096,
                                   0.002430179826511, -0.005206302729096,  0.003040175440704,
                                  -0.013487244282428,  0.003729889966203,  0.003059606676475,
                                   0.003729889966203, -0.021042542546869, -0.005509429951529,
                                   0.003059606676475, -0.005509429951529, -0.019147886244919,
                                  -0.003699760158805,  0.002289729070947, -0.002370719309518,
                                   0.002289729070947, -0.012032097539926,  0.005050336240100,
                                  -0.002370719309518,  0.005050336240100, -0.003975987070969,
                                   0.005629146676860,  0.001977315841567, -0.000082776369306,
                                   0.001977315841567,  0.012687891246471,  0.002189120352974,
                                  -0.000082776369306,  0.002189120352974,  0.018427959465896,
                                   0.448880778593399,  0.015295880495218,  0.016663590091395,
                                   0.015295880495218,  0.447402461005581,  0.010197253663478,
                                   0.016663590091395,  0.010197253663478,  0.434994453517236,
                                   0.015503379049345,  0.002189120352974, -0.004951273662081,
                                   0.002189120352974,  0.008458594164314,  0.000153048880756,
                                  -0.004951273662081,  0.000153048880756,  0.000534691712492,
                                   0.004243758336950, -0.005535059030204, -0.000138113456082,
                                  -0.005535059030204,  0.005013373974969,  0.004767201109256,
                                  -0.000138113456082,  0.004767201109256, -0.001045530324545,
                                   0.009884363501320, -0.000838581508292, -0.000494114807744,
                                  -0.000838581508292,  0.022205252327584, -0.000028124084560,
                                  -0.000494114807744, -0.000028124084560,  0.022291576962577,
                                   0.005408465728369, -0.000957561052355,  0.002430179826511,
                                  -0.000957561052355,  0.015803281683148, -0.005206302729096,
                                   0.002430179826511, -0.005206302729096,  0.003040175440704,
                                   0.015503379049345,  0.002189120352974, -0.004951273662081,
                                   0.002189120352974,  0.008458594164314,  0.000153048880756,
                                  -0.004951273662081,  0.000153048880756,  0.000534691712492,
                                   0.430487145938732,  0.022373707688196,  0.002899568254725,
                                   0.022373707688196,  0.438490990870954,  0.003952776259401,
                                   0.002899568254725,  0.003952776259401,  0.427247302022421,
                                   0.013926765376200,  0.004493787732038, -0.006461192264323,
                                   0.004493787732038,  0.011444562185213, -0.004129637487191,
                                  -0.006461192264323, -0.004129637487191,  0.008850410671922,
                                  -0.001586580633917,  0.001689840021891, -0.000947058276240,
                                   0.001689840021891, -0.001203370256285, -0.000686760094969,
                                  -0.000947058276240, -0.000686760094969,  0.002363575303981,
                                  -0.013487244282428,  0.003729889966203,  0.003059606676475,
                                   0.003729889966203, -0.021042542546869, -0.005509429951529,
                                   0.003059606676475, -0.005509429951529, -0.019147886244919,
                                   0.004243758336950, -0.005535059030204, -0.000138113456082,
                                  -0.005535059030204,  0.005013373974969,  0.004767201109256,
                                  -0.000138113456082,  0.004767201109256, -0.001045530324545,
                                   0.013926765376200,  0.004493787732038, -0.006461192264323,
                                   0.004493787732038,  0.011444562185213, -0.004129637487191,
                                  -0.006461192264323, -0.004129637487191,  0.008850410671922,
                                   0.399186667457387,  0.010056549795509, -0.002208199081057,
                                   0.010056549795509,  0.428317411682737,  0.008532924624016,
                                  -0.002208199081057,  0.008532924624016,  0.397839330123615};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.009302838556413, -0.020243533317730, -0.007682220260585,
                                  -0.024495620617241, -0.020996808366040, -0.015539840497723,
                                  -0.008443130768293, -0.014100118228845,  0.014672344633896,
                                  -0.008105142101738,  0.014223486401125,  0.002491038073235,
                                   0.017579619956947,  0.015708336006947,  0.004247237502934,
                                   0.009754905886976,  0.015326906124245, -0.010501003420599,
                                   0.011281395673340, -0.051623219781864, -0.007105971039211,
                                  -0.027030115525063, -0.040470224510258, -0.012712711159763,
                                  -0.013814123045937, -0.041406392676628,  0.022696971500912,
                                  -0.013064943767515,  0.052796423630378,  0.019868242492231,
                                   0.025584717435425,  0.041414888002866,  0.027575557017382,
                                   0.009081026557140,  0.025940131933016, -0.027952340893446,
                                  -0.000978468366576,  0.009815424730956,  0.011025901955978,
                                  -0.006662165958021,  0.003133017902914,  0.015155774321186,
                                  -0.007306248551252, -0.007468107816084, -0.001291946333613,
                                  -0.008105142101738,  0.012747267455687,  0.006550640173189,
                                   0.019383454098874,  0.014381548955384,  0.015192153876312,
                                   0.004794361996677,  0.007548186123915, -0.009174216369036,
                                  -0.010038977822735,  0.011437841713294,  0.004480150586245,
                                   0.020330983272123,  0.015209188437933,  0.010782434147139,
                                   0.006928717188249,  0.009342711878261, -0.009107063043564,
                                   0.002266364335899, -0.008996812485281, -0.001310240999962,
                                  -0.021485369444892, -0.017853822840749, -0.006741107036398,
                                  -0.007191283877909, -0.009197314339560,  0.004997439345396,
                                   0.003692522777709, -0.021264798064282, -0.016523172814250,
                                  -0.008166431324920, -0.019633022329581, -0.021427388328489,
                                  -0.001748583772498, -0.003072886244852,  0.007027748956676,
                                   0.013912347211806, -0.035584051730379, -0.019232787112507,
                                  -0.056016822728420, -0.042803399655918, -0.038642143165835,
                                  -0.018116120749913, -0.022973074169051,  0.027070153817711,
                                   0.008886606514505, -0.029912705070341, -0.015135279034747,
                                  -0.050467416128913, -0.035999928638009, -0.038299539028036,
                                  -0.012239663788986, -0.019287770371759,  0.020621464787499,
                                   0.003869132339824, -0.016437524096685, -0.007678006071879,
                                  -0.012925112978098, -0.017919574326121, -0.014886630630576,
                                  -0.000252242114188, -0.003614337453868,  0.003460422826843,
                                  -0.015453139101848,  0.053922076917315,  0.016905597996394,
                                   0.049908773826857,  0.047277689454512,  0.033272515884571,
                                   0.016905597996394,  0.035948051278210, -0.029541137432176,
                                  -0.000847944215972, -0.003614337453868,  0.004353975979518,
                                  -0.014886630630576, -0.011946382884080, -0.000519587452618,
                                  -0.006784669957018, -0.013425576218462,  0.005734314327083,
                                  -0.004097111230103,  0.011613384201634,  0.002392690841584,
                                   0.027079063478435,  0.014999970265837,  0.006881681124599,
                                   0.009342167431410,  0.009311935354745, -0.008197918479065,
                                  -0.009553976292220,  0.030583654222233,  0.013216008701613,
                                   0.052319187379936,  0.037903920527570,  0.031204048175901,
                                   0.016146490419495,  0.022086242553550, -0.027952340893446,
                                   0.002065917781727, -0.010354708070203, -0.000050595606671,
                                  -0.023986045390091, -0.019633022329581,  0.002302032769268,
                                  -0.016865656816759, -0.019438021235970,  0.008654353952658,
                                   0.000762084613005, -0.009197314339560, -0.006932613390975,
                                  -0.006741107036398, -0.011902548560499, -0.015867780247894,
                                   0.001888950925945, -0.001332383868981,  0.004080451174525,
                                  -0.006752718370623,  0.026473243677844,  0.011890548363199,
                                   0.030725330977356,  0.031570040379250,  0.023255108685847,
                                   0.009441981761194,  0.015931345049164, -0.021326882274330,
                                  -0.011009613534377,  0.032312049266914,  0.018492915323938,
                                   0.050462737242862,  0.039423331166941,  0.034938696780597,
                                   0.017528470467519,  0.021070089226530, -0.024036906804971,
                                   0.001215207937472, -0.002565189297613, -0.007465710757795,
                                   0.010703895417198,  0.003133017902914, -0.010893317741643,
                                   0.009357300575815,  0.011102813226769, -0.003485622637661,
                                   0.015217360054281, -0.054048092182716, -0.015402398749725,
                                  -0.034235943590430, -0.046770603200599, -0.020607007388772,
                                  -0.018287848817607, -0.039105098241414,  0.029732344519917,
                                  -0.007939225481070,  0.030597861577778,  0.013757570321857,
                                   0.004362972521414,  0.016862593545941,  0.015646326404278,
                                  -0.002341994494285,  0.003620383612260, -0.006218427508201,
                                  -0.013798265267815,  0.050821743773305,  0.013732700629852,
                                   0.030602944449207,  0.043081991462029,  0.015421286517793,
                                   0.018321481616924,  0.037313251210158, -0.024906915366621,
                                   0.007257135885661, -0.014160325915279, -0.004710693748727,
                                  -0.019040164383649, -0.014782989525260, -0.013484016169538,
                                  -0.003949783241020, -0.008835787926494,  0.012215058149273};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePDDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoBlock bgtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(375);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.011352502854703, -0.008705968077814,  0.002359317836368,
                                   0.000456002356125,  0.000619172935808, -0.008705968077814,
                                   0.027790974982564,  0.001150583882502, -0.002972146048470,
                                  -0.003050171049804,  0.002359317836368,  0.001150583882502,
                                   0.031919290895117, -0.007096153569423,  0.003186778465974,
                                   0.000456002356125, -0.002972146048470, -0.007096153569423,
                                   0.012037810461540, -0.009239059382840,  0.000619172935808,
                                  -0.003050171049804,  0.003186778465974, -0.009239059382840,
                                   0.012921208602682,  0.055123316015345,  0.004580326221416,
                                  -0.002171388861818,  0.001604760658535, -0.001132229191631,
                                   0.004580326221416,  0.052521248037142, -0.000919047399897,
                                   0.001663159432445, -0.003181641995016, -0.002171388861818,
                                  -0.000919047399897,  0.046903645149815,  0.001570141275567,
                                   0.002043908498621,  0.001604760658535,  0.001663159432445,
                                   0.001570141275567,  0.052594303492351,  0.004191369342343,
                                  -0.001132229191631, -0.003181641995016,  0.002043908498621,
                                   0.004191369342343,  0.056917403746014, -0.009840648111012,
                                  -0.001091483051821,  0.000063393299938,  0.004300701203259,
                                   0.000895151034671, -0.001091483051821, -0.014074784660083,
                                   0.003864235684785,  0.001977448038747, -0.003732991673281,
                                   0.000063393299938,  0.003864235684785, -0.024024032285376,
                                  -0.001880035533911, -0.000891504974952,  0.004300701203259,
                                   0.001977448038747, -0.001880035533911, -0.018798224354986,
                                  -0.001263756218049,  0.000895151034671, -0.003732991673281,
                                  -0.000891504974952, -0.001263756218049, -0.010001083533555,
                                   0.000636985404623, -0.000474029529208, -0.001322635547852,
                                  -0.000071056093417, -0.000728885255743, -0.000474029529208,
                                   0.001750797404373,  0.000521238604441,  0.000641125657445,
                                   0.001194309607697, -0.001322635547852,  0.000521238604441,
                                   0.003229482845524, -0.000562024072941,  0.000076446836801,
                                  -0.000071056093417,  0.000641125657445, -0.000562024072941,
                                   0.000608600336121, -0.002694864463782, -0.000728885255743,
                                   0.001194309607697,  0.000076446836801, -0.002694864463782,
                                  -0.002491835197352, -0.020446120905705, -0.003004316101104,
                                  -0.002234360315519,  0.006629929160592, -0.000840687760221,
                                  -0.003004316101104, -0.026674741805752,  0.005759812726975,
                                   0.004582389441525, -0.006838325364904, -0.002234360315519,
                                   0.005759812726975, -0.039927622883654, -0.007798613590978,
                                   0.002298488490113,  0.006629929160592,  0.004582389441525,
                                  -0.007798613590978, -0.035624054952129, -0.005375985102382,
                                  -0.000840687760221, -0.006838325364904,  0.002298488490113,
                                  -0.005375985102382, -0.019153668406180,  0.000636985404623,
                                  -0.000474029529208, -0.001322635547852, -0.000071056093417,
                                  -0.000728885255743, -0.000474029529208,  0.001750797404373,
                                   0.000521238604441,  0.000641125657445,  0.001194309607697,
                                  -0.001322635547852,  0.000521238604441,  0.003229482845524,
                                  -0.000562024072941,  0.000076446836801, -0.000071056093417,
                                   0.000641125657445, -0.000562024072941,  0.000608600336121,
                                  -0.002694864463782, -0.000728885255743,  0.001194309607697,
                                   0.000076446836801, -0.002694864463782, -0.002491835197352,
                                   0.022042225095994,  0.001587594370460,  0.003367939237311,
                                  -0.003616478215175,  0.002243729770461,  0.001587594370460,
                                   0.023949262979416, -0.000371906802382, -0.002175157990811,
                                   0.007039261359563,  0.003367939237311, -0.000371906802382,
                                   0.040079319736906,  0.003994069920070, -0.008135525989601,
                                  -0.003616478215175, -0.002175157990811,  0.003994069920070,
                                   0.043120850588024,  0.005881551377006,  0.002243729770461,
                                   0.007039261359563, -0.008135525989601,  0.005881551377006,
                                   0.022054491872514,  0.070843376512221,  0.005329533543956,
                                  -0.008669368315634,  0.005510491000814, -0.001168287900900,
                                   0.005329533543956,  0.061040077915571, -0.000394793145856,
                                   0.007720453154382, -0.004765348206278, -0.008669368315634,
                                  -0.000394793145856,  0.051960050985119, -0.000263195430571,
                                   0.003612236798181,  0.005510491000814,  0.007720453154382,
                                  -0.000263195430571,  0.054606366953586,  0.001521025927439,
                                  -0.001168287900900, -0.004765348206278,  0.003612236798181,
                                   0.001521025927439,  0.068526272175436,  0.012270894690666,
                                  -0.009001094813861,  0.007423922035685,  0.000427617287623,
                                  -0.000907902842690, -0.009001094813861,  0.031989826965394,
                                   0.003994069920070, -0.008126626846284, -0.007041528459842,
                                   0.007423922035685,  0.003994069920070,  0.026719546491271,
                                  -0.003700298402440,  0.007005026112880,  0.000427617287623,
                                  -0.008126626846284, -0.003700298402440,  0.012723582079566,
                                  -0.007878932883765, -0.000907902842690, -0.007041528459842,
                                   0.007005026112880, -0.007878932883765,  0.017126916621673,
                                   0.000929759718706, -0.002252453252329,  0.000076446836801,
                                   0.002475279826258,  0.001488912741326, -0.002252453252329,
                                   0.000000011905968,  0.000965188666822, -0.000571098534126,
                                  -0.001074867799126,  0.000076446836801,  0.000965188666822,
                                  -0.001345617852302, -0.001539380430950, -0.001171015988197,
                                   0.002475279826258, -0.000571098534126, -0.001539380430950,
                                  -0.000983094297841, -0.000961542125472,  0.001488912741326,
                                  -0.001074867799126, -0.001171015988197, -0.000961542125472,
                                  -0.000156905638402,  0.018640606697969,  0.000473782370710,
                                   0.001023778505440, -0.004435582046799, -0.000147405062259,
                                   0.000473782370710,  0.021652763310520, -0.001305146613875,
                                  -0.001017809021724,  0.006666667237014,  0.001023778505440,
                                  -0.001305146613875,  0.041762541374198, -0.000548226606138,
                                  -0.004292905587671, -0.004435582046799, -0.001017809021724,
                                  -0.000548226606138,  0.040056275310809,  0.002694456252664,
                                  -0.000147405062259,  0.006666667237014, -0.004292905587671,
                                   0.002694456252664,  0.020223208277191,  0.000636985404623,
                                  -0.000474029529208, -0.001322635547852, -0.000071056093417,
                                  -0.000728885255743, -0.000474029529208,  0.001750797404373,
                                   0.000521238604441,  0.000641125657445,  0.001194309607697,
                                  -0.001322635547852,  0.000521238604441,  0.003229482845524,
                                  -0.000562024072941,  0.000076446836801, -0.000071056093417,
                                   0.000641125657445, -0.000562024072941,  0.000608600336121,
                                  -0.002694864463782, -0.000728885255743,  0.001194309607697,
                                   0.000076446836801, -0.002694864463782, -0.002491835197352,
                                  -0.007617774800257,  0.005164633155408,  0.000707716892249,
                                  -0.001545558802760, -0.000451764180416,  0.005164633155408,
                                  -0.012222174077340, -0.001880035533911,  0.002501208819144,
                                   0.000809680467109,  0.000707716892249, -0.001880035533911,
                                  -0.016016021523584,  0.005430931963044,  0.001002754611529,
                                  -0.001545558802760,  0.002501208819144,  0.005430931963044,
                                  -0.009693165266040,  0.004831759112689, -0.000451764180416,
                                   0.000809680467109,  0.001002754611529,  0.004831759112689,
                                  -0.005610046296121,  0.056982835452757,  0.000075419716757,
                                  -0.002018495188217,  0.006555320311052,  0.001845596291020,
                                   0.000075419716757,  0.052521271849078,  0.001011329933748,
                                   0.000520962364193, -0.005331377593269, -0.002018495188217,
                                   0.001011329933748,  0.044212409445212, -0.001508619586333,
                                  -0.000298123477773,  0.006555320311052,  0.000520962364193,
                                  -0.001508619586333,  0.050628114896670,  0.002268285091398,
                                   0.001845596291020, -0.005331377593269, -0.000298123477773,
                                   0.002268285091398,  0.056603592469210,  0.014112764830626,
                                  -0.009240533763668,  0.003838802584567,  0.002648601303964,
                                  -0.000744597264977, -0.009240533763668,  0.025368278463733,
                                   0.006099803102338, -0.009486701552911, -0.005393765116967,
                                   0.003838802584567,  0.006099803102338,  0.022540818610118,
                                  -0.006010794925459,  0.002081425712213,  0.002648601303964,
                                  -0.009486701552911, -0.006010794925459,  0.017474826344969,
                                  -0.007008725774037, -0.000744597264977, -0.005393765116967,
                                   0.002081425712213, -0.007008725774037,  0.013922593236686};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.003996053996195,  0.002388695986808, -0.009101002341234,
                                   0.008691261101780,  0.001998708860675,  0.008474154752748,
                                   0.007906043876367, -0.005590229503081, -0.006507336498756,
                                  -0.014256118677437, -0.014971544209827, -0.006682869356015,
                                   0.016441131386674, -0.002151358171286,  0.006369207919866,
                                   0.010665225293603, -0.004803958402176, -0.004587443716754,
                                   0.008989785633939,  0.012777068913283, -0.001586915125126,
                                  -0.011322338683772,  0.004249056288919,  0.011047365880799,
                                  -0.001996363096157, -0.011669547868134,  0.015231026647597,
                                  -0.029188597395098,  0.026005732268600, -0.000631594622545,
                                   0.013018132550891,  0.022422020038606, -0.017254222024964,
                                  -0.013984229855874, -0.035016348344630, -0.031919001834450,
                                  -0.017906283369783,  0.043256470228449, -0.010027418958297,
                                   0.017061934726234,  0.029965569492634, -0.018727681047048,
                                  -0.007336830743804,  0.026827609853686,  0.034270396665437,
                                   0.003247749195958, -0.036119201455436,  0.013526135611367,
                                   0.031773480073443, -0.005639056789882,  0.013237292096994,
                                  -0.015909707840064,  0.063428183050629, -0.048930841791870,
                                  -0.005800053644913, -0.013204703733523, -0.022107647775539,
                                   0.035027190259408,  0.030040714728293,  0.039836740137900,
                                   0.037169054620037,  0.029604257437140, -0.064787673168850,
                                   0.015004163890621, -0.015662142195169, -0.032657765706736,
                                   0.021808718950580,  0.021666005505607, -0.040196052640080,
                                  -0.037173599305643, -0.000214245657526,  0.055943081932375,
                                  -0.033394788857554, -0.055942354790562,  0.005906585687836,
                                   0.005723949268035, -0.005106469892608,  0.018779945110729,
                                  -0.022399392778659, -0.012243601753701, -0.008833337146072,
                                  -0.010073636559447,  0.010927855019324,  0.032334717772335,
                                   0.023863304136042,  0.019364612096039,  0.014897509345164,
                                  -0.025790669464135,  0.001437226688902,  0.000661071750564,
                                  -0.011445581295398,  0.005418785553338,  0.010918721601390,
                                  -0.017488014496999, -0.022708490275189,  0.008433485708870,
                                   0.025206603833362, -0.019436922664212, -0.029169512928472,
                                  -0.001394256144293,  0.008568312493728, -0.008635366786944,
                                   0.025255441644745, -0.025841760679319, -0.002250595029858,
                                  -0.010197101636015, -0.013552283628396,  0.015396771736633,
                                   0.019265797608037,  0.039278128884069,  0.033301171651804,
                                   0.015104931045284, -0.040114518642773,  0.009994380227227,
                                  -0.015116786368516, -0.029290564533204,  0.018129405923382,
                                   0.008665141297370, -0.026749700249304, -0.035693887924560,
                                  -0.002969106772268,  0.027551465028541, -0.010885260716922,
                                  -0.029217274146063,  0.001863458130762,  0.004484338980619,
                                  -0.003838966348651,  0.014157698276879, -0.013514928497937,
                                  -0.004871931559521, -0.005361536106524, -0.005984959175834,
                                   0.010348354065899,  0.014129087058103,  0.028698840654497,
                                   0.043283366361356,  0.015484385480444, -0.033195859672191,
                                   0.010328316024390, -0.010567530085060, -0.021926513370518,
                                   0.005954607264636,  0.011789519528758, -0.016248404209584,
                                  -0.027682276119940,  0.001605324697903,  0.012925063065420,
                                  -0.006739340672722, -0.012779334530114,  0.001922256680150,
                                   0.004125499420773, -0.003633268009917,  0.018054476350907,
                                  -0.014775900060120,  0.001254620484559, -0.004327724509740,
                                  -0.005235827005065,  0.016378937421944,  0.011899360951363,
                                   0.028823528548101,  0.032879522767156,  0.021847855187000,
                                  -0.029966456541840,  0.014655624833888, -0.019010801698168,
                                  -0.027567416561897,  0.024685490182417,  0.008515130409831,
                                  -0.022028404234026, -0.031575439130262, -0.006643999282079,
                                   0.019339193618528, -0.005628074587039, -0.018619010858193,
                                  -0.000580798858614, -0.008029381803460,  0.008253099465944,
                                  -0.021553213660867,  0.021277898080675,  0.001991848013566,
                                   0.008982606218654,  0.011934279683065, -0.004985138654774,
                                  -0.009746620305576, -0.030667797675451, -0.031845948824176,
                                  -0.010024737001792,  0.039824184594368, -0.006683158001195,
                                   0.013269145343407,  0.026699977287573, -0.009746620305576,
                                  -0.003323425769849,  0.020056463271045,  0.029262577090029,
                                   0.001991848013566, -0.023947372745203,  0.008980505692028,
                                   0.021897081084642, -0.004078883243220,  0.003098900463697,
                                  -0.002457499417500,  0.010578914573036, -0.011307128184016,
                                  -0.007374761636941, -0.004985120564220, -0.004920295265835,
                                   0.008515130409831,  0.022380273938353,  0.022622071264164,
                                   0.026005066879530,  0.014655624833888, -0.019977637694560,
                                   0.009634834492094, -0.002868943715358, -0.013281092375838,
                                   0.003201080091771,  0.009282995413752, -0.013255858893558,
                                  -0.021486521512443,  0.004473168013017,  0.010945150684095,
                                  -0.007249816910857, -0.013114028537722, -0.000735766755592,
                                  -0.002349399666833,  0.003203473107068, -0.006739340672722,
                                   0.001803050407414, -0.007569081795217,  0.002284426602355,
                                   0.004768019869061, -0.010669616057731,  0.004677653235285,
                                  -0.008236201624873, -0.010567530085060, -0.006458837365507,
                                   0.013831608196746, -0.004292949887179,  0.022324431692654,
                                   0.011939344958274, -0.014941098269035,  0.004624402991760,
                                   0.004495881541530,  0.012861354882920,  0.007976333222601,
                                  -0.007902744531470,  0.000791339275981,  0.005687121244524,
                                  -0.000320015191821, -0.006769368699679,  0.006419142196808,
                                  -0.023022939169024,  0.022692789464581,  0.001435451136205,
                                   0.007113598696631,  0.009240066581239, -0.013249101953622,
                                  -0.013264916970436, -0.031121337141276, -0.031547730556838,
                                  -0.012256552408269,  0.036548482897520, -0.010946575070320,
                                   0.012035454081587,  0.033642298174269, -0.016148177046907,
                                  -0.009047926676826,  0.027610015847483,  0.030455554112090,
                                   0.005111947317603, -0.026975373636405,  0.008775314088846,
                                   0.024846218682842, -0.003329368466513,  0.004761907530714,
                                  -0.005858959228359,  0.026621498144251, -0.017906054182672,
                                   0.009095854243736, -0.002679163732272, -0.006262527725929,
                                   0.015884142616218,  0.001854385543049,  0.017806837458681,
                                   0.012052915372388,  0.007791137187213, -0.025790669464135,
                                   0.012096784925828, -0.016887000386197, -0.023363170458363,
                                   0.032723193754692,  0.003484290206050, -0.021299123330518,
                                  -0.020070979960583, -0.014814019734254,  0.029894127026660,
                                  -0.000617195383760, -0.028125746696754, -0.000432214406972,
                                   0.008078465621442, -0.009321831055124,  0.047644208638721,
                                  -0.032732791652996,  0.001444413744645, -0.006758468727047,
                                  -0.011493331079600,  0.021666005505607,  0.009757022405190,
                                   0.024984713445610,  0.024914012903220,  0.015004163890621,
                                  -0.043191782112567,  0.017100787528290, -0.010118284554444,
                                  -0.025393589593556,  0.022105016071759,  0.016972185671402,
                                  -0.030042469197479, -0.023754117918284, -0.006934298236436,
                                   0.039119270977630, -0.009402184719779, -0.035675373810267,
                                   0.004684119568444, -0.014005437238648,  0.013738946130486,
                                  -0.049242715443070,  0.053599689709074,  0.007249081997880,
                                   0.014394186708268,  0.020004382096376, -0.023912738438643,
                                  -0.029953167818119, -0.053518747865917, -0.048144930470084,
                                  -0.025290607529363,  0.064748694781895, -0.018770508307801,
                                   0.016298037067322,  0.047627230412880, -0.025209716626944,
                                  -0.020107809565267,  0.052213290342956,  0.049069214819444,
                                   0.003369738179377, -0.054331359365403,  0.019153578071204,
                                   0.052627202347483, -0.006911245022486, -0.007061576489534,
                                   0.007674609676356, -0.022734263944893,  0.022159246786014,
                                   0.004394341724786,  0.007563782104027,  0.010732464178929,
                                  -0.012122838471115, -0.018507083942834, -0.032142261327588,
                                  -0.028362491516321, -0.014258085941532,  0.034667536306855,
                                  -0.007759329200978,  0.011488843504610,  0.026396537070668,
                                  -0.017590159451691, -0.009297539995005,  0.023806417181713,
                                   0.034967204321679,  0.001387453824875, -0.026945155241091,
                                   0.009698763796920,  0.024352382421003, -0.000633011611902};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 1);
    
    CGtoBlock dgtos(mlih, mbas, 2);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(375);

    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.010416260504712, -0.050418752082308,  0.008925950778969,
                                  -0.000625184469979,  0.018712128428764, -0.000625184469979,
                                  -0.020324003601835, -0.064901083487216, -0.011344503645188,
                                  -0.000828249033647, -0.017072075428104, -0.000625184469979,
                                   0.006915571106267, -0.052075250149601, -0.012906363659885,
                                   0.007930414758108, -0.004171742914141,  0.000970454305506,
                                   0.000435999174610,  0.002740436341174,  0.000435999174610,
                                  -0.001482553078744, -0.004860969159267,  0.008226951135949,
                                   0.002073089341584, -0.000422992901130,  0.000435999174610,
                                  -0.004692521037848, -0.000025564230973,  0.008424221230978,
                                  -0.002165783599953,  0.002002057374270, -0.000046722382957,
                                   0.001188190457811,  0.002017312034038,  0.001188190457811,
                                  -0.003035554068153,  0.007951378031797, -0.006753386857995,
                                  -0.000074591586149, -0.000937316492982,  0.001188190457811,
                                  -0.000636459335653,  0.001852874201973, -0.003485223351739,
                                  -0.000424077779178, -0.001439948126617, -0.003925216615958,
                                   0.000053318830580, -0.006049588867321,  0.000053318830580,
                                   0.003291621264544, -0.005053988314107, -0.000434666500868,
                                  -0.002265430278345,  0.004040707006777,  0.000053318830580,
                                   0.001389081846024, -0.005970808683307, -0.002416549054876,
                                  -0.000566376325463,  0.001027071036270, -0.000816103346933,
                                   0.000663968994773,  0.000764623662203,  0.000663968994773,
                                  -0.002042612512529,  0.001059150775585,  0.000824716813605,
                                  -0.001356323455936,  0.000135543610790,  0.000663968994773,
                                   0.000412827211885, -0.001685575875603,  0.000676055433336,
                                   0.007930414758108, -0.004171742914141,  0.000970454305506,
                                   0.000435999174610,  0.002740436341174,  0.000435999174610,
                                  -0.001482553078744, -0.004860969159267,  0.008226951135949,
                                   0.002073089341584, -0.000422992901130,  0.000435999174610,
                                  -0.004692521037848, -0.000025564230973,  0.008424221230978,
                                  -0.025389621765202, -0.048047424159955,  0.012766021611601,
                                  -0.001684744647592,  0.024371800716835, -0.001684744647592,
                                  -0.022101903669374, -0.055871188021721, -0.029329935180268,
                                   0.000012866837096, -0.019841500994355, -0.001684744647592,
                                   0.011083508855226, -0.048021690485763, -0.023137513878459,
                                  -0.001024232566153,  0.000882697304514, -0.003538599321361,
                                  -0.000481367511744, -0.005277132541976, -0.000481367511744,
                                   0.000361633656377,  0.000377592298666, -0.003635496974705,
                                  -0.000885250252602,  0.001213534926627, -0.000481367511744,
                                   0.001685049972017, -0.000887803200690, -0.005529856573377,
                                   0.002729576367329, -0.001477942347041, -0.001827225021724,
                                  -0.000569655647708, -0.004186971340120, -0.000569655647708,
                                   0.001981578122841, -0.007012715839017,  0.007402440721069,
                                   0.000524485727962,  0.000959453357208, -0.000569655647708,
                                  -0.002290161589827, -0.000428970891117,  0.008639540584862,
                                   0.002774300348019,  0.002882917336109,  0.003402159817653,
                                  -0.001076560302375,  0.006229110271091, -0.001076560302375,
                                  -0.006399019644222,  0.004382748686669,  0.006446339134624,
                                   0.000989680482692, -0.006062623032776, -0.001076560302375,
                                  -0.000711451536008,  0.004862278301492,  0.004904308994093,
                                  -0.002165783599953,  0.002002057374270, -0.000046722382957,
                                   0.001188190457811,  0.002017312034038,  0.001188190457811,
                                  -0.003035554068153,  0.007951378031797, -0.006753386857995,
                                  -0.000074591586149, -0.000937316492982,  0.001188190457811,
                                  -0.000636459335653,  0.001852874201973, -0.003485223351739,
                                  -0.001024232566153,  0.000882697304514, -0.003538599321361,
                                  -0.000481367511744, -0.005277132541976, -0.000481367511744,
                                   0.000361633656377,  0.000377592298666, -0.003635496974705,
                                  -0.000885250252602,  0.001213534926627, -0.000481367511744,
                                   0.001685049972017, -0.000887803200690, -0.005529856573377,
                                  -0.029171101217619, -0.043013451867195,  0.021819260645557,
                                  -0.003040784128222,  0.036449332901114, -0.003040784128222,
                                  -0.036792902147915, -0.047613161404248, -0.024528601431944,
                                   0.001266993386759, -0.038135637268494, -0.003040784128222,
                                   0.014546173763705, -0.040479465093677, -0.020552212228786,
                                   0.006468569530380, -0.001416125990315,  0.001685049972017,
                                   0.000514334517167,  0.007082226475630,  0.000514334517167,
                                  -0.003635496974705,  0.000251728199111,  0.003391214468631,
                                   0.001414424024924,  0.000528137336100,  0.000514334517167,
                                  -0.004942807631375,  0.001412722059532,  0.005483838887832,
                                  -0.002904918231576, -0.001843361223666,  0.000806372064019,
                                  -0.000074591586149, -0.002118604153359, -0.000074591586149,
                                   0.007413403584725, -0.003313074179915, -0.006408972893661,
                                   0.001040250478615,  0.003898122046682, -0.000074591586149,
                                  -0.000915165216824,  0.000237139733565, -0.001918942009226,
                                  -0.000424077779178, -0.001439948126617, -0.003925216615958,
                                   0.000053318830580, -0.006049588867321,  0.000053318830580,
                                   0.003291621264544, -0.005053988314107, -0.000434666500868,
                                  -0.002265430278345,  0.004040707006777,  0.000053318830580,
                                   0.001389081846024, -0.005970808683307, -0.002416549054876,
                                   0.002729576367329, -0.001477942347041, -0.001827225021724,
                                  -0.000569655647708, -0.004186971340120, -0.000569655647708,
                                   0.001981578122841, -0.007012715839017,  0.007402440721069,
                                   0.000524485727962,  0.000959453357208, -0.000569655647708,
                                  -0.002290161589827, -0.000428970891117,  0.008639540584862,
                                   0.006468569530380, -0.001416125990315,  0.001685049972017,
                                   0.000514334517167,  0.007082226475630,  0.000514334517167,
                                  -0.003635496974705,  0.000251728199111,  0.003391214468631,
                                   0.001414424024924,  0.000528137336100,  0.000514334517167,
                                  -0.004942807631375,  0.001412722059532,  0.005483838887832,
                                  -0.011085502604005, -0.048194025736223,  0.017075939225599,
                                  -0.000635773191668,  0.032561710972506, -0.000635773191668,
                                  -0.039564829233596, -0.050027258155873, -0.011781220088378,
                                   0.000954015595929, -0.036579512586786, -0.000635773191668,
                                   0.008811131702907, -0.046285994544364, -0.015957045265989,
                                   0.008409035390729, -0.003850464850906,  0.001130079076526,
                                   0.002459044187595,  0.004905393872141,  0.002459044187595,
                                  -0.005398452556747, -0.001430188841250,  0.007197861135727,
                                   0.000897134274288, -0.002450413061720,  0.002459044187595,
                                  -0.004385299880791, -0.002056196302331,  0.006384387944953,
                                  -0.000566376325463,  0.001027071036270, -0.000816103346933,
                                   0.000663968994773,  0.000764623662203,  0.000663968994773,
                                  -0.002042612512529,  0.001059150775585,  0.000824716813605,
                                  -0.001356323455936,  0.000135543610790,  0.000663968994773,
                                   0.000412827211885, -0.001685575875603,  0.000676055433336,
                                   0.002774300348019,  0.002882917336109,  0.003402159817653,
                                  -0.001076560302375,  0.006229110271091, -0.001076560302375,
                                  -0.006399019644222,  0.004382748686669,  0.006446339134624,
                                   0.000989680482692, -0.006062623032776, -0.001076560302375,
                                  -0.000711451536008,  0.004862278301492,  0.004904308994093,
                                  -0.002904918231576, -0.001843361223666,  0.000806372064019,
                                  -0.000074591586149, -0.002118604153359, -0.000074591586149,
                                   0.007413403584725, -0.003313074179915, -0.006408972893661,
                                   0.001040250478615,  0.003898122046682, -0.000074591586149,
                                  -0.000915165216824,  0.000237139733565, -0.001918942009226,
                                   0.008409035390729, -0.003850464850906,  0.001130079076526,
                                   0.002459044187595,  0.004905393872141,  0.002459044187595,
                                  -0.005398452556747, -0.001430188841250,  0.007197861135727,
                                   0.000897134274288, -0.002450413061720,  0.002459044187595,
                                  -0.004385299880791, -0.002056196302331,  0.006384387944953,
                                  -0.011841270842401, -0.052061677722652,  0.009073985096056,
                                   0.002224999812271,  0.017537965739533,  0.002224999812271,
                                  -0.020338505001398, -0.062800434448972, -0.015763835423634,
                                   0.000161659307691, -0.018517142475356,  0.002224999812271,
                                   0.005084386143750, -0.051738359107269, -0.012728979707506};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.004603402610068, -0.011731355374054,  0.018615232944613,
                                   0.009149974058749,  0.012152823741784,  0.011787736168944,
                                   0.015857911569309, -0.007836906070043,  0.009887894893405,
                                  -0.006937528669964, -0.010901903897188,  0.010242200937712,
                                   0.012092137441551, -0.017716257658288, -0.011562808498822,
                                  -0.000305898428028,  0.012104859564151, -0.014816766167360,
                                  -0.002873600726256, -0.007982213900809, -0.002326253866845,
                                  -0.012837808248364,  0.003411008759543, -0.002104746649913,
                                   0.010183156909333,  0.008012212387008, -0.013507227811928,
                                  -0.009081670582550,  0.013425319220867,  0.009440433537613,
                                  -0.002918638662648, -0.014892906944737,  0.053219792866015,
                                   0.006407645300153,  0.011073217080683,  0.005791954812595,
                                   0.019849124424374,  0.006592485888275,  0.006777286484544,
                                  -0.006102136566010, -0.011853717495151,  0.019284143327279,
                                   0.040843124573897, -0.031511916802556, -0.010108736340985,
                                   0.006230922176175,  0.016800038915879, -0.043799309327034,
                                  -0.028349626897213, -0.020697358072299, -0.012155530428850,
                                  -0.017351897071047,  0.000415430739956, -0.021712622563739,
                                  -0.007228149676125,  0.011768866348513, -0.011420293763727,
                                  -0.027908889155827,  0.039360213099508,  0.016282759359225,
                                   0.002995166314292,  0.000958106289274, -0.003144527426200,
                                  -0.015416916224463, -0.002913358398313, -0.006096220093283,
                                   0.005928563007647, -0.000241855527738, -0.012931864660606,
                                  -0.015463749804650, -0.002951387319499,  0.010365756778954,
                                  -0.001082960864745,  0.001922755250748,  0.002772262895085,
                                   0.007298589687640,  0.009339591158453, -0.011693424659566,
                                  -0.007176839733781, -0.009785475961105, -0.012905435232963,
                                  -0.011796632511884,  0.008855936483983, -0.009070371870322,
                                   0.007213751091626,  0.006971036650527, -0.005412670362411,
                                  -0.007910438327307,  0.010337027547047,  0.008285104438044,
                                   0.005184879541717,  0.013330416936087, -0.011504632171363,
                                  -0.004924679725585, -0.009919990679771, -0.009566643961179,
                                  -0.014009550302041,  0.008527211142273, -0.006955290836598,
                                   0.012179310229450,  0.007449931486835, -0.008021108729947,
                                  -0.005697349850771,  0.010317974701510,  0.009921805026879,
                                  -0.000779891529669, -0.001770143871970,  0.008499503910799,
                                   0.002003638874395,  0.006092991754845,  0.009868451907009,
                                   0.014961254692874,  0.001862603047277,  0.006770781220576,
                                  -0.012575309939655, -0.006432733854762,  0.007828505822046,
                                   0.003853738700906, -0.004854381876840, -0.005773609874837,
                                  -0.002644417700974, -0.004226149115845,  0.012095264033803,
                                   0.022623950022570,  0.009775800643446,  0.006128328097409,
                                   0.005933479998085,  0.002022158028137,  0.021069750426582,
                                   0.009182240179275, -0.002743648295798,  0.001911867884064,
                                   0.006238765216755, -0.008511277944660, -0.008864337964534,
                                  -0.007190046959693, -0.015900493420225,  0.030748811832354,
                                   0.021998612939960,  0.029257766885479,  0.028179416182252,
                                   0.036087700368428, -0.018867277863633,  0.025526473977846,
                                  -0.012941027460187, -0.020129997584807,  0.020538381224211,
                                   0.019132584537819, -0.030731246394288, -0.025926347236360,
                                  -0.004541731583007, -0.010314515624337,  0.025544155878478,
                                   0.016743262943485,  0.019145397512841,  0.036388266730017,
                                   0.030789155112747, -0.016949301846590,  0.022380792439422,
                                  -0.009675891865794, -0.016356420950230,  0.015473646242299,
                                   0.018917569029300, -0.023359770182371, -0.018517162834135,
                                  -0.002109659056999, -0.005642407802896,  0.012492155492627,
                                   0.010959016817823,  0.003372445338127,  0.012519624123899,
                                   0.012389263216973, -0.001031016281444,  0.014450213729790,
                                  -0.002163212024490,  0.000583572369289, -0.002496117245681,
                                   0.003710215607684, -0.002919451445400, -0.004979549437591,
                                   0.006726291575890,  0.020128277071980, -0.052522815586900,
                                  -0.020953471719689, -0.022156443750574, -0.025233580111445,
                                  -0.029619975953517,  0.009716977270701, -0.019746650635678,
                                   0.010513991713102,  0.019371749969019, -0.020953471719689,
                                  -0.035015210391267,  0.037589503505054,  0.020959170903838,
                                   0.003208994129610,  0.002006854348120,  0.003710215607684,
                                  -0.009673908836102,  0.001321462787180,  0.005997725292496,
                                   0.014450213729790, -0.000687344187629,  0.000347418442148,
                                  -0.007078999245909, -0.005767745707410,  0.010508792259154,
                                   0.009400309152890, -0.007714760513651, -0.000090193352956,
                                   0.002471644023369,  0.003843565720222, -0.013097969592775,
                                  -0.004999087237594, -0.007507277033830, -0.009675891865795,
                                  -0.015239910370512,  0.007062209102746, -0.005708686572791,
                                   0.017197747862858,  0.006493467400073, -0.008046167320440,
                                  -0.004200463463761,  0.010187386699239,  0.008323748539006,
                                   0.004456494818322,  0.010806673370964, -0.026811123570205,
                                  -0.014996713040851, -0.021542156151156, -0.027793691640022,
                                  -0.033734751307694,  0.017578807139558, -0.019938195205885,
                                   0.008491820325810,  0.025491696601831, -0.018200783362786,
                                  -0.018789900166882,  0.031875032961335,  0.021547610106730,
                                  -0.000497712890866, -0.002151005461491,  0.011169502109290,
                                  -0.000242954057067,  0.005575685231064, -0.004628746786208,
                                   0.022648188610248,  0.002022158028137, -0.004002312491663,
                                  -0.016634739541407, -0.008506517106989,  0.023049485585074,
                                   0.007627408103525, -0.010586421599014, -0.008666912676429,
                                  -0.002136473321166, -0.000245490728203,  0.003853738700906,
                                   0.007256612948534,  0.004537415405793,  0.013089317921392,
                                   0.006770781220576,  0.001241735364851,  0.009318937009060,
                                   0.007241300000356, -0.002016085998713, -0.001480687472942,
                                   0.005288054993377, -0.004170859771004, -0.004928001118020,
                                   0.003287327281120,  0.007581043918928, -0.019061512723336,
                                  -0.014381381213116, -0.016482336645932, -0.017019143323312,
                                  -0.023335468743687,  0.006842079452159, -0.017941388527221,
                                  -0.001101898860912,  0.015792238048534, -0.011284952208754,
                                  -0.014680080079028,  0.026681870269089,  0.014702319390999,
                                   0.004710088370519,  0.012474178851120, -0.026229060045135,
                                  -0.021046282805453, -0.023411805404154, -0.022956677310129,
                                  -0.033668036671847,  0.015441453518144, -0.022864855486024,
                                   0.010562716511256,  0.018952181105318, -0.020194694267544,
                                  -0.018047426488001,  0.027391499874588,  0.026860752126148,
                                  -0.002654652035872, -0.001783612541372, -0.001917981156548,
                                   0.010425014636751, -0.000366660655912,  0.005261000742259,
                                  -0.009030723974669, -0.000241855527738,  0.009507065812869,
                                   0.011793580200651,  0.005671499254336, -0.013191408390055,
                                  -0.002922780269223,  0.004664474081393, -0.000692202001902,
                                  -0.006520688689187, -0.020486530777040,  0.048324690110686,
                                   0.021018837483379,  0.015747949564156,  0.006951834860988,
                                   0.025727185809081,  0.000994893979851,  0.013859173347674,
                                  -0.004749481907286, -0.019539880784732,  0.026722075692214,
                                   0.035115104999148, -0.041033216389704, -0.018217710860495,
                                   0.003997381085658,  0.013302295409307, -0.029147155156359,
                                  -0.023761302561461, -0.006113449294045, -0.006102136566010,
                                   0.000121633561061, -0.002746869120114, -0.015412022475385,
                                  -0.006310616043325,  0.000505017651552,  0.007142292703642,
                                  -0.006559606379177,  0.006033047818732,  0.005832325524022,
                                   0.007065018573695,  0.018885788054681, -0.048207227987615,
                                  -0.020432495396691, -0.016436420649662, -0.002381166603608,
                                  -0.022733987663049, -0.003797317586139, -0.011331172970705,
                                   0.007438373246096,  0.014821234147603, -0.026638010697996,
                                  -0.030522333726025,  0.034037365948388,  0.012177766848854,
                                  -0.003213554511612, -0.009009144041065,  0.013012674244118,
                                   0.004890648157334,  0.009840642884834,  0.013812720784685,
                                   0.011529227445593, -0.008316586200058,  0.008370197783196,
                                  -0.003729328394048, -0.008722017126946,  0.003798421278370,
                                   0.008993134017604, -0.014724760270234, -0.008365159479998};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeDDDDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 2);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CMemBlock<double> fints(625);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 1.101230025915860,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.976188030537569,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.036667069492300,  0.000000000000000,  0.054961117819803,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.022679639608025,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.062520997689144,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.026188192065634,  0.000000000000000,
                                   0.062520997689144,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.026188192065634,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.032281478211778,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.032281478211778,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.054961117819803,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.022679639608025,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.101230025915860,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.021547309753620,  0.000000000000000,  0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.026188192065634,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.039841358081119,
                                   0.000000000000000, -0.013094096032817,  0.000000000000000,
                                   0.039841358081119,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.013094096032817,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.026188192065634,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.054961117819803,  0.000000000000000,
                                  -0.026188192065634,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.013094096032817,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000, -0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.062520997689144,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.026188192065634,  0.000000000000000,  0.062520997689144,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.039841358081119,  0.000000000000000, -0.013094096032817,
                                   0.000000000000000,  0.039841358081119,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.976188030537569,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.021547309753620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.101230025915860,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.021547309753620,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.976188030537569,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.039841358081119,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.039841358081119,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.062520997689144,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.062520997689144,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.026188192065634,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000, -0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.054961117819803,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.039841358081119,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.039841358081119,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.021547309753620,  0.000000000000000, -0.026188192065634,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.101230025915860,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.026188192065634,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.022679639608025,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.054961117819803,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.032281478211778,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.032281478211778,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.013094096032817,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000, -0.013094096032817,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.026188192065634,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.062520997689144,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.026188192065634,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.062520997689144,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.022679639608025,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.022679639608025,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.013094096032817,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.013094096032817,  0.000000000000000,
                                   0.054961117819803,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.054961117819803,  0.000000000000000,
                                   1.036667069492300,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.976188030537569,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.991307790276252,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.101230025915860};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.768988866552595,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.681005874064739,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.723431229256007,  0.000000000000000,  0.038688326845020,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.015909508196726,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.043991496243929,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.018370717680108,  0.000000000000000,
                                   0.043991496243929,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.018370717680108,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.022778818648294,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.022778818648294,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.038688326845020,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.015909508196726,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.768988866552595,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.712824890458190,  0.000000000000000,  0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.018370717680108,  0.000000000000000,
                                   0.691612212862555,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028081988047203,
                                   0.000000000000000, -0.009185358840054,  0.000000000000000,
                                   0.028081988047203,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.009185358840054,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.018370717680108,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.038688326845020,  0.000000000000000,
                                  -0.018370717680108,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.009185358840054,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000, -0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.043991496243929,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.018370717680108,  0.000000000000000,  0.043991496243929,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.028081988047203,  0.000000000000000, -0.009185358840054,
                                   0.000000000000000,  0.028081988047203,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.681005874064739,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.712824890458190,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.768988866552596,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.712824890458190,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.681005874064739,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.028081988047203,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.028081988047203,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.043991496243929,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000, -0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.043991496243929,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.018370717680108,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000, -0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.038688326845020,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.028081988047203,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.028081988047203,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862556,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.712824890458190,  0.000000000000000, -0.018370717680108,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.768988866552595,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.018370717680108,  0.000000000000000,
                                   0.691612212862555,  0.000000000000000,  0.015909508196726,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.038688326845020,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.022778818648294,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.022778818648294,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.009185358840054,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000, -0.009185358840054,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  -0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.018370717680108,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.043991496243929,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000, -0.018370717680108,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.043991496243929,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.015909508196726,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.015909508196726,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.009185358840054,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.009185358840054,  0.000000000000000,
                                   0.038688326845020,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.038688326845020,  0.000000000000000,
                                   0.723431229256007,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862555,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.681005874064739,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.691612212862555,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.768988866552595};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));

    std::vector<double> r0123vals{ 0.331015487530059,  0.004002296508697, -0.001777762433551,
                                   0.007797087388452,  0.000329884591335,  0.004002296508697,
                                   0.353297184403464,  0.004228159293716,  0.001850122684312,
                                  -0.009638045610313, -0.001777762433551,  0.004228159293716,
                                   0.364429362026196,  0.000313686264262,  0.000675456756960,
                                   0.007797087388452,  0.001850122684312,  0.000313686264262,
                                   0.356334913762142,  0.007113503657013,  0.000329884591335,
                                  -0.009638045610313,  0.000675456756960,  0.007113503657013,
                                   0.333188602619478,  0.003894644180792, -0.003132212204246,
                                   0.001844662282100,  0.000009332088705,  0.000601738236393,
                                  -0.003132212204246,  0.011610160652775, -0.000091926132477,
                                  -0.001247479371853, -0.001665168015158,  0.001844662282100,
                                  -0.000091926132477,  0.008728377660179, -0.000198564939942,
                                   0.003811428970873,  0.000009332088705, -0.001247479371853,
                                  -0.000198564939942,  0.002915111360809, -0.003937203897498,
                                   0.000601738236393, -0.001665168015158,  0.003811428970873,
                                  -0.003937203897498,  0.005362508476496, -0.001755383569002,
                                   0.001862317239002,  0.001490281545567,  0.000548557588061,
                                  -0.000078254242790,  0.001862317239002, -0.004767748159621,
                                   0.001655640332556, -0.002792024067967, -0.002334927706517,
                                   0.001490281545567,  0.001655640332556, -0.012581878918335,
                                   0.004715614802553,  0.000074726912567,  0.000548557588061,
                                  -0.002792024067967,  0.004715614802553, -0.004332914606589,
                                   0.004433623505234, -0.000078254242790, -0.002334927706517,
                                   0.000074726912567,  0.004433623505234, -0.001823151227687,
                                   0.007663705969549,  0.000025901607010,  0.000541724374004,
                                  -0.000667432202681, -0.000155384999179,  0.000025901607010,
                                   0.006855331571019,  0.001387803039313, -0.000324786533962,
                                   0.002471548419976,  0.000541724374004,  0.001387803039313,
                                   0.011114926170552, -0.000957677616326, -0.004807193233753,
                                  -0.000667432202681, -0.000324786533962, -0.000957677616326,
                                   0.016642061739330,  0.001451694081180, -0.000155384999179,
                                   0.002471548419976, -0.004807193233753,  0.001451694081180,
                                   0.009219475026384,  0.000327457649350,  0.000608544315643,
                                  -0.000078986961466, -0.000165725818679, -0.000745043019839,
                                   0.000608544315643, -0.002174356530556,  0.002207801575602,
                                  -0.001395716400426, -0.000802719715964, -0.000078986961466,
                                   0.002207801575602, -0.001542917494167, -0.001012343159495,
                                  -0.000031822376489, -0.000165725818679, -0.001395716400426,
                                  -0.001012343159495,  0.003530076246552, -0.000687273368975,
                                  -0.000745043019839, -0.000802719715964, -0.000031822376489,
                                  -0.000687273368975,  0.000355975943413,  0.003894644180792,
                                  -0.003132212204246,  0.001844662282100,  0.000009332088705,
                                   0.000601738236393, -0.003132212204246,  0.011610160652775,
                                  -0.000091926132477, -0.001247479371853, -0.001665168015158,
                                   0.001844662282100, -0.000091926132477,  0.008728377660179,
                                  -0.000198564939942,  0.003811428970873,  0.000009332088705,
                                  -0.001247479371853, -0.000198564939942,  0.002915111360809,
                                  -0.003937203897498,  0.000601738236393, -0.001665168015158,
                                   0.003811428970873, -0.003937203897498,  0.005362508476496,
                                   0.353040286802783,  0.011717812980679, -0.004774265248194,
                                   0.007067221897021, -0.002184971431875,  0.011717812980679,
                                   0.368272487574481,  0.002766856797790,  0.004564497622183,
                                  -0.009904758725521, -0.004774265248194,  0.002766856797790,
                                   0.372772904424830,  0.007702676398920,  0.004324112452608,
                                   0.007067221897021,  0.004564497622183,  0.007702676398920,
                                   0.365325700613905,  0.010614396482074, -0.002184971431875,
                                  -0.009904758725521,  0.004324112452608,  0.010614396482074,
                                   0.356436767397167,  0.004010138743175, -0.000129567347360,
                                   0.001687955096984,  0.001379377862221,  0.002213771850524,
                                  -0.000129567347360,  0.002582711869159,  0.001555432721499,
                                   0.000292518647247,  0.000510207666792,  0.001687955096984,
                                   0.001555432721499,  0.004432447330962,  0.000395528552290,
                                  -0.005363180517371,  0.001379377862221,  0.000292518647247,
                                   0.000395528552290,  0.008681602539714,  0.001446863421274,
                                   0.002213771850524,  0.000510207666792, -0.005363180517371,
                                   0.001446863421274,  0.003799787211659,  0.001741345125770,
                                  -0.001234874113907, -0.002812172181348, -0.000346760478890,
                                  -0.001396350953413, -0.001234874113907,  0.004387351238532,
                                   0.000311233998168,  0.001987170967824,  0.002391089479923,
                                  -0.002812172181348,  0.000311233998168,  0.009261076421628,
                                  -0.002233387834719,  0.000198890381453, -0.000346760478890,
                                   0.001987170967824, -0.002233387834719,  0.002374308324550,
                                  -0.005181007406299, -0.001396350953413,  0.002391089479923,
                                   0.000198890381453, -0.005181007406299, -0.004248920626667,
                                  -0.009506876525216, -0.001624802192064, -0.002341644094614,
                                   0.002438297888035, -0.000813060535465, -0.001624802192064,
                                  -0.009735342087429,  0.000482585286133,  0.002380849834672,
                                  -0.002204477396364, -0.002341644094614,  0.000482585286133,
                                  -0.010521634074637, -0.002336517317031,  0.002042493195922,
                                   0.002438297888035,  0.002380849834672, -0.002336517317031,
                                  -0.014274896809913, -0.003190172590676, -0.000813060535465,
                                  -0.002204477396364,  0.002042493195922, -0.003190172590676,
                                  -0.008275540073717, -0.001755383569002,  0.001862317239002,
                                   0.001490281545567,  0.000548557588061, -0.000078254242790,
                                   0.001862317239002, -0.004767748159621,  0.001655640332556,
                                  -0.002792024067967, -0.002334927706517,  0.001490281545567,
                                   0.001655640332556, -0.012581878918335,  0.004715614802553,
                                   0.000074726912567,  0.000548557588061, -0.002792024067967,
                                   0.004715614802553, -0.004332914606589,  0.004433623505234,
                                  -0.000078254242790, -0.002334927706517,  0.000074726912567,
                                   0.004433623505234, -0.001823151227687,  0.004010138743175,
                                  -0.000129567347360,  0.001687955096984,  0.001379377862221,
                                   0.002213771850524, -0.000129567347360,  0.002582711869159,
                                   0.001555432721499,  0.000292518647247,  0.000510207666792,
                                   0.001687955096984,  0.001555432721499,  0.004432447330962,
                                   0.000395528552290, -0.005363180517371,  0.001379377862221,
                                   0.000292518647247,  0.000395528552290,  0.008681602539714,
                                   0.001446863421274,  0.002213771850524,  0.000510207666792,
                                  -0.005363180517371,  0.001446863421274,  0.003799787211659,
                                   0.364058805621932,  0.008882395827818, -0.012522269472199,
                                   0.011266217208873, -0.001561182737816,  0.008882395827818,
                                   0.372845498163956,  0.004577441282160,  0.009386931616688,
                                  -0.010649004249017, -0.012522269472199,  0.004577441282160,
                                   0.373366703566746,  0.003051627521440,  0.005217612280083,
                                   0.011266217208873,  0.009386931616688,  0.003051627521440,
                                   0.365023055150050,  0.005727751810776, -0.001561182737816,
                                  -0.010649004249017,  0.005217612280083,  0.005727751810776,
                                   0.360962459858596,  0.000160363087743, -0.000191057447051,
                                   0.004709815623857, -0.001004885247063, -0.001016272651183,
                                  -0.000191057447051,  0.007576567491167,  0.000395528552290,
                                  -0.002244470726166, -0.002322181320977,  0.004709815623857,
                                   0.000395528552290,  0.002954964887308,  0.001225825594591,
                                   0.005027375663834, -0.001004885247063, -0.002244470726166,
                                   0.001225825594591, -0.000067024551918, -0.000037519085783,
                                  -0.001016272651183, -0.002322181320977,  0.005027375663834,
                                  -0.000037519085783,  0.005046254215479,  0.000669324668099,
                                   0.003828876281671,  0.000074726912567, -0.004812082477047,
                                  -0.000028961841554,  0.003828876281671,  0.004256785537504,
                                  -0.005385261093668,  0.000217416776516,  0.002075072107900,
                                   0.000074726912567, -0.005385261093668,  0.005242449549306,
                                   0.005078276931591,  0.001638489922159, -0.004812082477047,
                                   0.000217416776516,  0.005078276931591, -0.000464842718250,
                                   0.002076368643207, -0.000028961841554,  0.002075072107900,
                                   0.001638489922159,  0.002076368643207,  0.000821731497188,
                                   0.007663705969549,  0.000025901607010,  0.000541724374004,
                                  -0.000667432202681, -0.000155384999179,  0.000025901607010,
                                   0.006855331571019,  0.001387803039313, -0.000324786533962,
                                   0.002471548419976,  0.000541724374004,  0.001387803039313,
                                   0.011114926170552, -0.000957677616326, -0.004807193233753,
                                  -0.000667432202681, -0.000324786533962, -0.000957677616326,
                                   0.016642061739330,  0.001451694081180, -0.000155384999179,
                                   0.002471548419976, -0.004807193233753,  0.001451694081180,
                                   0.009219475026384,  0.001741345125770, -0.001234874113907,
                                  -0.002812172181348, -0.000346760478890, -0.001396350953413,
                                  -0.001234874113907,  0.004387351238532,  0.000311233998168,
                                   0.001987170967824,  0.002391089479923, -0.002812172181348,
                                   0.000311233998168,  0.009261076421628, -0.002233387834719,
                                   0.000198890381453, -0.000346760478890,  0.001987170967824,
                                  -0.002233387834719,  0.002374308324550, -0.005181007406299,
                                  -0.001396350953413,  0.002391089479923,  0.000198890381453,
                                  -0.005181007406299, -0.004248920626667,  0.000160363087743,
                                  -0.000191057447051,  0.004709815623857, -0.001004885247063,
                                  -0.001016272651183, -0.000191057447051,  0.007576567491167,
                                   0.000395528552290, -0.002244470726166, -0.002322181320977,
                                   0.004709815623857,  0.000395528552290,  0.002954964887308,
                                   0.001225825594591,  0.005027375663834, -0.001004885247063,
                                  -0.002244470726166,  0.001225825594591, -0.000067024551918,
                                  -0.000037519085783, -0.001016272651183, -0.002322181320977,
                                   0.005027375663834, -0.000037519085783,  0.005046254215479,
                                   0.356163134662040,  0.003119834557986, -0.004376484485287,
                                   0.016775443158234,  0.003514506753389,  0.003119834557986,
                                   0.365473322600281,  0.008853775766409,  0.002551454708201,
                                  -0.014413764453911, -0.004376484485287,  0.008853775766409,
                                   0.365055340740140,  0.000044411977213, -0.000511300063657,
                                   0.016775443158234,  0.002551454708201,  0.000044411977213,
                                   0.362490149285537,  0.006655414802385,  0.003514506753389,
                                  -0.014413764453911, -0.000511300063657,  0.006655414802385,
                                   0.355403565788658,  0.007017158813006, -0.003937591551468,
                                   0.004422416364479,  0.001411328258086, -0.000694079448225,
                                  -0.003937591551468,  0.010593292683873,  0.001470765833181,
                                  -0.005199196636827, -0.003206742108981,  0.004422416364479,
                                   0.001470765833181,  0.005695995836613, -0.000045250390798,
                                   0.002049712611344,  0.001411328258086, -0.005199196636827,
                                  -0.000045250390798,  0.006553190107670, -0.002018556979619,
                                  -0.000694079448225, -0.003206742108981,  0.002049712611344,
                                  -0.002018556979619,  0.006836086926283,  0.000327457649350,
                                   0.000608544315643, -0.000078986961466, -0.000165725818679,
                                  -0.000745043019839,  0.000608544315643, -0.002174356530556,
                                   0.002207801575602, -0.001395716400426, -0.000802719715964,
                                  -0.000078986961466,  0.002207801575602, -0.001542917494167,
                                  -0.001012343159495, -0.000031822376489, -0.000165725818679,
                                  -0.001395716400426, -0.001012343159495,  0.003530076246552,
                                  -0.000687273368975, -0.000745043019839, -0.000802719715964,
                                  -0.000031822376489, -0.000687273368975,  0.000355975943413,
                                  -0.009506876525216, -0.001624802192064, -0.002341644094614,
                                   0.002438297888035, -0.000813060535465, -0.001624802192064,
                                  -0.009735342087429,  0.000482585286133,  0.002380849834672,
                                  -0.002204477396364, -0.002341644094614,  0.000482585286133,
                                  -0.010521634074637, -0.002336517317031,  0.002042493195922,
                                   0.002438297888035,  0.002380849834672, -0.002336517317031,
                                  -0.014274896809913, -0.003190172590676, -0.000813060535465,
                                  -0.002204477396364,  0.002042493195922, -0.003190172590676,
                                  -0.008275540073717,  0.000669324668099,  0.003828876281671,
                                   0.000074726912567, -0.004812082477047, -0.000028961841554,
                                   0.003828876281671,  0.004256785537504, -0.005385261093668,
                                   0.000217416776516,  0.002075072107900,  0.000074726912567,
                                  -0.005385261093668,  0.005242449549306,  0.005078276931591,
                                   0.001638489922159, -0.004812082477047,  0.000217416776516,
                                   0.005078276931591, -0.000464842718250,  0.002076368643207,
                                  -0.000028961841554,  0.002075072107900,  0.001638489922159,
                                   0.002076368643207,  0.000821731497188,  0.007017158813006,
                                  -0.003937591551468,  0.004422416364479,  0.001411328258086,
                                  -0.000694079448225, -0.003937591551468,  0.010593292683873,
                                   0.001470765833181, -0.005199196636827, -0.003206742108981,
                                   0.004422416364479,  0.001470765833181,  0.005695995836613,
                                  -0.000045250390798,  0.002049712611344,  0.001411328258086,
                                  -0.005199196636827, -0.000045250390798,  0.006553190107670,
                                  -0.002018556979619, -0.000694079448225, -0.003206742108981,
                                   0.002049712611344, -0.002018556979619,  0.006836086926283,
                                   0.333183789184543,  0.005445719852102, -0.001842487300712,
                                   0.009357640270707,  0.000358402885398,  0.005445719852102,
                                   0.356724062681100,  0.004020569705094, -0.004142778202882,
                                  -0.008411492984234, -0.001842487300712,  0.004020569705094,
                                   0.361369242329431,  0.005185466401611,  0.000832980632316,
                                   0.009357640270707, -0.004142778202882,  0.005185466401611,
                                   0.355596879587897,  0.006907990817992,  0.000358402885398,
                                  -0.008411492984234,  0.000832980632316,  0.006907990817992,
                                   0.332375777590643};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));

    std::vector<double> r0312vals{ 0.007288959625488, -0.005948029413270,  0.014226645671457,
                                  -0.013575467913362, -0.000652322799705, -0.006761289447279,
                                  -0.009419788637837,  0.007800737156122,  0.007402679112683,
                                   0.018532300970240,  0.017293403242787,  0.007757832486206,
                                  -0.021911144222168,  0.004817828524136, -0.008217957830366,
                                  -0.016560204634388,  0.008301881132516,  0.004257497721682,
                                  -0.014122743298611, -0.017670574424808, -0.001357439037125,
                                   0.017096419541614, -0.006140140812297, -0.015453581395548,
                                   0.002350338281913, -0.004779190769011,  0.011653837867980,
                                  -0.012639524783521,  0.010099142419064, -0.004493319643623,
                                   0.000308461580951,  0.009296465858980, -0.007447795077824,
                                  -0.002194654673217, -0.011251363815515, -0.007073658807722,
                                  -0.005631986112670,  0.015410281220362, -0.006129107569007,
                                   0.007570934719781,  0.011760367613021, -0.010678384816528,
                                   0.000959559288956,  0.010675586922995,  0.012683336501353,
                                   0.005989182789889, -0.015972259143399,  0.005058781173549,
                                   0.011885430551506, -0.002066343445260,  0.007416712698617,
                                  -0.010687168096580,  0.048212540538411, -0.030370629453054,
                                   0.001657484047478, -0.004669536991898, -0.011775191943855,
                                   0.020779836977740,  0.006864032451797,  0.018019797376781,
                                   0.013856426314767,  0.014621465905123, -0.037070853640842,
                                   0.014320395888424, -0.008298875542685, -0.017514352296394,
                                   0.015349959115104,  0.013476585958411, -0.024582961363605,
                                  -0.016691826038346, -0.004993568723184,  0.033091514179371,
                                  -0.013262369139371, -0.029549922785729,  0.005454265791536,
                                  -0.008849809515038,  0.008916132713614, -0.032660934552567,
                                   0.040808806805092,  0.009044435521719,  0.009209709435489,
                                   0.013355256923896, -0.012281793047719, -0.022496408347628,
                                  -0.028633204966405, -0.021004417501954, -0.016727434382956,
                                   0.036005939163869, -0.008457630476336,  0.003620446159557,
                                   0.017431186737861, -0.009160865467928, -0.014307431374410,
                                   0.030888035263938,  0.023322119827483, -0.003789676376871,
                                  -0.033091591781535,  0.014787839091683,  0.034828656268416,
                                  -0.004302491448097, -0.000541555415727, -0.002014090326409,
                                  -0.000213285296109,  0.005359849116552,  0.012590045993803,
                                   0.004364556066989,  0.000861181670710,  0.004089191060044,
                                  -0.013309695260823, -0.004200260816798, -0.005313451809982,
                                  -0.003664237872508,  0.002066415262625,  0.003809556689509,
                                  -0.009992515537928, -0.004280460507225,  0.008885214612306,
                                  -0.007939626649260,  0.001497399132839,  0.004035530911964,
                                  -0.010034093055202, -0.001285930116757,  0.004641768479258,
                                   0.001936725895133,  0.000533237659821, -0.004100474186636,
                                   0.002088660952773, -0.008342500257686,  0.007888656308824,
                                   0.002145656988207,  0.011493303067735,  0.008617749276605,
                                  -0.005628889150868, -0.006049708281252, -0.015617824150631,
                                  -0.017221952677687, -0.006277644062959,  0.017665276261505,
                                  -0.001754302009397,  0.008246168272490,  0.011100334412488,
                                  -0.004464445580218, -0.004892254770523,  0.008430031149170,
                                   0.013337821897719, -0.002218545989872, -0.010484270542031,
                                   0.003097974917396,  0.010267424830398, -0.002075782712579,
                                  -0.004850315228578,  0.006635650598336, -0.008661174713003,
                                   0.007388247377474, -0.001202838229792,  0.007448910632345,
                                   0.013499347226382, -0.007103468368793, -0.004291413342342,
                                  -0.015200190194335, -0.013586933067140, -0.007743913944832,
                                   0.016734881466909, -0.003177621730730,  0.009956621411298,
                                   0.011779097716903, -0.007425800851320, -0.000987255162193,
                                   0.009045743424952,  0.015309687227554,  0.000888550247666,
                                  -0.012509087813690,  0.003988946315076,  0.010269389433283,
                                  -0.002144853828453,  0.000982713135641, -0.001128672167041,
                                   0.005414741626845, -0.001558829542113,  0.003641863651161,
                                  -0.001456918207219, -0.002610408768380,  0.016866093334196,
                                   0.003854207623769,  0.009925095178931,  0.011565432899076,
                                   0.004309606733090, -0.008636920296810,  0.004794303277005,
                                  -0.012224927004425, -0.010680333976250,  0.008065171262496,
                                   0.003364983350120, -0.006112562592319, -0.009668982460833,
                                  -0.004516152962530,  0.003903424034508, -0.001836480804379,
                                  -0.005119769122422, -0.000282497843220,  0.002012338132186,
                                  -0.001527761752900,  0.007448400092507, -0.008908586825074,
                                  -0.007992749843818, -0.004197142025169, -0.003893443957967,
                                   0.004846621510298,  0.024114595089837,  0.014374240051762,
                                   0.010008773886409,  0.006966110678354, -0.009398426207724,
                                   0.003675778351103,  0.004001874916044, -0.004332356879172,
                                   0.001617530217825,  0.006079020956423, -0.008088690754804,
                                  -0.013141947782577,  0.005717815125867,  0.009937711417551,
                                  -0.006835614273054, -0.013486923217925, -0.002608015691836,
                                   0.008496433135613, -0.007974754916107,  0.024259838152979,
                                  -0.024014984490526, -0.002422986988457, -0.011368124642475,
                                  -0.013597725627803,  0.014893104917139,  0.020306839980515,
                                   0.046510991740498,  0.038777848512514,  0.015174630028548,
                                  -0.043670777163777,  0.010384549214461, -0.017682980712954,
                                  -0.032130588544347,  0.018377094577992,  0.008006818693598,
                                  -0.026631641183246, -0.041338891841014, -0.003058405699466,
                                   0.025288664493010, -0.009698520358350, -0.027248103007039,
                                   0.001019261238706,  0.005667071760240, -0.005069144500330,
                                   0.019498876733775, -0.017607789850813, -0.003202309430528,
                                  -0.006975434784673, -0.008147764126428,  0.011716792270917,
                                   0.015063298140211,  0.031377864712531,  0.046625258335438,
                                   0.014785047518243, -0.039863847956236,  0.010424182207873,
                                  -0.013311550770558, -0.026319209633723,  0.012625246093665,
                                   0.011197956811486, -0.019861253576424, -0.030288886774177,
                                  -0.001142239097282,  0.021487595364083, -0.007771878926346,
                                  -0.021593354747415,  0.001647765440699,  0.001981002193967,
                                  -0.002726889808932,  0.007256385320776, -0.006434793645782,
                                  -0.004955034149873, -0.002866365915587, -0.004481647885506,
                                   0.004732357590438,  0.009771926400809,  0.006179951626858,
                                   0.011943892361582,  0.018122327331497, -0.010217354440867,
                                   0.004563278399365, -0.002903747071319, -0.001190677753910,
                                   0.001634033447113,  0.004764706938122, -0.004344720809341,
                                  -0.008644523555905,  0.005369700800903,  0.007388113261441,
                                  -0.008643736098526, -0.005957916674566,  0.000403044576659,
                                  -0.009210265404079,  0.011348337133741, -0.040669542556001,
                                   0.032910409351035,  0.002527705620285,  0.009919022047758,
                                   0.015630674772295, -0.020351657453802, -0.016589077879408,
                                  -0.034774651146129, -0.037813487236726, -0.017711340986976,
                                   0.053813762037499, -0.011807560657984,  0.015755619681969,
                                   0.030183239282784, -0.016589077879408, -0.013567771635868,
                                   0.029454906338468,  0.033386238238184,  0.002527705620285,
                                  -0.037676780446162,  0.016945642731667,  0.035709789397724,
                                  -0.004196982590515,  0.000234185445882, -0.002150212435618,
                                   0.009160853456058, -0.004359728253670,  0.006804344769024,
                                   0.003824677806066,  0.001525536600964,  0.004764706938122,
                                  -0.004181975540118,  0.003392329788038,  0.006368075872331,
                                   0.004563278399365, -0.006811569627244,  0.014319595332026,
                                  -0.005762677065815, -0.008118559683454,  0.008024863890426,
                                   0.000761768475336, -0.007409782397529, -0.001602899293923,
                                  -0.008682757657140,  0.004325078229828,  0.004613341048947,
                                  -0.004435355735356,  0.001355179067869, -0.002566513725334,
                                   0.001356554310781, -0.007771878926346,  0.008603206832537,
                                  -0.001770493666625,  0.003944832894008,  0.003606816942402,
                                  -0.009284787333019, -0.002931082269142, -0.011619273617081,
                                  -0.013311550770558, -0.006898165487025,  0.016609936648432,
                                  -0.003236815721645,  0.020224015973831,  0.010860660953609,
                                  -0.008782407180854,  0.001938370492085,  0.008063607100453,
                                   0.013922737365887,  0.003173675133166, -0.008465118842829,
                                   0.004084650196521,  0.008584947622041, -0.000481335108391,
                                  -0.006356866493908,  0.005766077985780, -0.021120620718540,
                                   0.020958480511219,  0.001351171129055,  0.007035643996151,
                                   0.008451729297691, -0.012842504631979, -0.012941449354871,
                                  -0.034390698412149, -0.037203565606916, -0.012324586360765,
                                   0.039707808831888, -0.010877228077321,  0.013065520818410,
                                   0.040492121468037, -0.016295665927261, -0.008695851431168,
                                   0.028395092242808,  0.034392549395351,  0.006009908966525,
                                  -0.025074243150076,  0.009561372519063,  0.022985063926281,
                                  -0.003100833489912,  0.001454085781107, -0.001474551177624,
                                   0.009985959454883, -0.006307453477849,  0.006035626743321,
                                  -0.000788982178453, -0.001646903148571,  0.008917177777967,
                                  -0.000624523248164,  0.008472847274140,  0.002935201997384,
                                   0.006835756264163, -0.009398426207724,  0.003871309972388,
                                  -0.012974697617616, -0.012891770629881,  0.024484531214663,
                                  -0.000026813445081, -0.010335231564200, -0.012563076169665,
                                  -0.009417761957579,  0.013400774317701, -0.000745471803352,
                                  -0.012247239463814, -0.002049763340757,  0.000437641455879,
                                   0.001521978585407,  0.003287011181892, -0.003974270905837,
                                  -0.006625630359985, -0.003411740498047, -0.000388653393479,
                                   0.003364983350120,  0.006568174391644,  0.007775641886855,
                                   0.013407209083853,  0.004794303277005, -0.005757946864540,
                                   0.000314354002253,  0.005522657826442, -0.003350258781197,
                                   0.000251728933554,  0.014061940542429, -0.005426660846987,
                                  -0.007922104767090,  0.005611394560550,  0.002274170383626,
                                  -0.001999080302664, -0.001777197754142,  0.000029168739069,
                                  -0.004509374315290,  0.003953728048440, -0.016242260789408,
                                   0.020684697123457,  0.003837513325615,  0.005232424571699,
                                   0.006472048805318, -0.008777628453818, -0.013867588457552,
                                  -0.025176299993649, -0.022559755768712, -0.011344288641238,
                                   0.024566903306679, -0.009547846659983,  0.005104498936974,
                                   0.023669433844484, -0.010733200948575, -0.009600142719548,
                                   0.026058154106659,  0.022558992067982,  0.001746124848157,
                                  -0.019220519646471,  0.006387485144262,  0.019009599894003,
                                  -0.001989396775366, -0.006895351752671,  0.007116903149244,
                                  -0.022181422334311,  0.020475548300356,  0.005021272049580,
                                   0.007536295773120,  0.010398971402201, -0.011073780973897,
                                  -0.019457780068782, -0.037274201424677, -0.032731999541626,
                                  -0.014984133228823,  0.037920447735605, -0.008026766636765,
                                   0.013933610013155,  0.030142849887195, -0.018906336015014,
                                  -0.008924463007957,  0.023937373739724,  0.042710971197301,
                                   0.001444240204939, -0.025125524683995,  0.009093273641183,
                                   0.022551717183077,  0.000422679079793,  0.000195158555757,
                                   0.001054934981531, -0.003595653312830,  0.000634167344258,
                                  -0.008604278098038, -0.003779641632501, -0.001722800183700,
                                  -0.004892918405667,  0.009070165048642, -0.000501453758953,
                                   0.002640998194111,  0.002431292502429,  0.002066415262625,
                                  -0.005333738872898,  0.009098164471896,  0.008170752464875,
                                  -0.011341327485152,  0.005533537549308,  0.004081380987248,
                                  -0.000948880259928,  0.012621380413417, -0.002723477615230,
                                  -0.003475914760872,  0.003799901997923, -0.000203476311664,
                                   0.009453337379310, -0.011544899907389,  0.037779964066223,
                                  -0.033763173505301, -0.002675042685834, -0.009306128313174,
                                  -0.016149231507831,  0.016995928384379,  0.016754720642697,
                                   0.023626116656837,  0.021164223275300,  0.017165705092660,
                                  -0.039873094363867,  0.011305169994906, -0.008399754616557,
                                  -0.026914317066327,  0.021082540166564,  0.012459804655718,
                                  -0.032379624584465, -0.025620935115961, -0.003286217755949,
                                   0.047027808885595, -0.014527254569963, -0.040051762583396,
                                   0.004564163476344, -0.004045920907812,  0.006729547773840,
                                  -0.013262369139371,  0.009993026822568,  0.007728620548429,
                                   0.004705605637060,  0.009168657496855, -0.008168635105952,
                                  -0.016586996705843, -0.008499943015154, -0.008298875542685,
                                  -0.012036854997646,  0.015446189017017,  0.002950034685591,
                                  -0.002603010178224,  0.003279461571235,  0.003779227286093,
                                  -0.006349689601498,  0.005980573047920,  0.008734698958360,
                                  -0.008233906101162, -0.013757226689842,  0.021908841745325,
                                   0.014074366262106, -0.001316986796418, -0.008343832487911,
                                   0.008867546245264, -0.034293935247327,  0.030401297032405,
                                  -0.001848460731818,  0.009343419015992,  0.014810530316737,
                                  -0.019319440952304, -0.016063288193694, -0.023070276938429,
                                  -0.016769486673722, -0.014130814358572,  0.035175741131462,
                                  -0.011295031513907,  0.010688808552434,  0.021376884325143,
                                  -0.021091143084132, -0.007931250695948,  0.027065515754708,
                                   0.018053864419044,  0.006003442077382, -0.038868752877946,
                                   0.014549748316911,  0.042833356752094, -0.002923834510344,
                                   0.003999841883965, -0.004968799747137,  0.010748538477981,
                                  -0.011677717100063, -0.001019359933061, -0.004818181184412,
                                  -0.007163592533933,  0.005277243210492,  0.003685860410098,
                                   0.015129019852817,  0.016068978920817,  0.005886937177041,
                                  -0.017812753951296,  0.004278684584695, -0.005683034737802,
                                  -0.013341397371241,  0.002786658390265,  0.004461155856061,
                                  -0.011701087475477, -0.013125236829291, -0.000314243695642,
                                   0.010820147853654, -0.004266185916635, -0.009076091049841,
                                   0.005622959140887};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeMaxQValues)
{
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(kgtos, 1.0e-13);
    
    CGtoPairsContainer bcont({bpairs.pick(0), bpairs.pick(1)});
    
    CGtoPairsContainer kcont({kpairs.pick(0), kpairs.pick(1)});
    
    CVecMemBlock<double> bvec({CMemBlock<double>(1), CMemBlock<double>(1)});
    
    CVecMemBlock<double> kvec({CMemBlock<double>(1), CMemBlock<double>(1)});
    
    eridrv.computeMaxQValues(&bvec, &kvec, &bcont, &kcont);
    
    vlxtest::compare({1.27641389857182689641}, bvec[0].data());
    
    vlxtest::compare({1.13863555260159077820}, bvec[1].data());
    
    vlxtest::compare({1.15070414827136605747}, kvec[0].data());
    
    vlxtest::compare({0.99288328066163043835}, kvec[1].data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestJKForHeAtom)
{
    auto mhe = vlxmol::getMoleculeHeAtom();

    auto mbas = vlxbas::getMolecularBasisForHeAtom();

    const int nrows = 5, ncols = 5;

    std::vector<double> denvals{  0.351213608529393,  0.304029605741231, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.304029605741231,
                                  0.263184566094148, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000};

    const CAODensityMatrix dmat({CDenseMatrix(denvals, nrows, ncols)}, denmat::rest);

    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);

    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mhe, mbas);

    CAOFockMatrix fock(dmat);

    eridrv.compute(fock, dmat, mhe, mbas, qqdata);

    ASSERT_EQ(fock.getNumberOfRows(0), nrows);

    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);

    std::vector<double> jkvals {  1.329815389795996,  0.545891434518611,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.545891434518611,
                                  0.865839035335250,  0.000000000000000, -0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  1.614984494992040,  0.000000000000000, -0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  1.614984494992040,  0.000000000000000, -0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  1.614984494992040};

    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals.size()));

    CAOFockMatrix jk ({CDenseMatrix(jkvals, nrows, ncols)}, {fockmat::restjk}, {1.0}, {0});

    ASSERT_EQ(jk, fock);
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestJKForHeAnionWithSPD)
{
    auto mhe = vlxmol::getMoleculeHeAtom();
    
    auto mbas = vlxbas::getMolecularBasisSPDForHeAtom();
    
    const int nrows = 9, ncols = 9;
    
    std::vector<double> denvals{  1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00,
                                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00};
    
    const CAODensityMatrix dmat({CDenseMatrix(denvals, nrows, ncols)}, denmat::rest);
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mhe, mbas);
    
    CAOFockMatrix fock(dmat);
    
    eridrv.compute(fock, dmat, mhe, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> jkvals { 15.5086818237907,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000, 13.6451976996879,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000, 13.6451976996879,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                 13.6451976996879,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000, 11.1814020170352,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000, 11.1814020170352,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                 11.1814020170352,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000, 11.1814020170352,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000,  0.0000000000000,
                                  0.0000000000000,  0.0000000000000, 11.1814020170352};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals.size()));
    
    CAOFockMatrix jk ({CDenseMatrix(jkvals, nrows, ncols)}, {fockmat::restjk}, {1.0}, {0});
    
    ASSERT_EQ(jk, fock);
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestJKForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getRestDensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> jkvals {
    12.5167971574448949,  -3.6658315880217520,  -1.6447438883350372,  -0.5140361149152390,  -0.6527772822764955,
    -0.5140361149152397,  -0.6527772822764958,   0.0000000000000005,  -0.0000000000000001,   0.7933779079230698,
    -0.7933779079230704,  -0.0040937126161196,   0.0012546716139458,   0.6242803742030969,   0.6242803742030973,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000003,  -0.0002408453277244,  -0.0000000000000000,   0.0045540032341421,  -3.6658315880217511,
     7.7833897590559324,   5.3597706735459720,   2.1591397431695505,   2.5095629959661094,   2.1591397431695505,
     2.5095629959661081,  -0.0000000000000000,   0.0000000000000000,  -2.6086424092976292,   2.6086424092976310,
     0.1239397763729299,   0.0556812082323447,  -2.0184051429381134,  -2.0184051429381147,  -0.0000000000000001,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0023221868883481,  -0.0000000000000000,   0.0020012845931223,  -1.6447438883350372,   5.3597706735459720,
     5.5119107894364481,   2.6937998841993211,   3.3986261138240828,   2.6937998841993220,   3.3986261138240828,
     0.0000000000000002,   0.0000000000000001,  -1.8679699639745320,   1.8679699639745322,   0.0756373318876606,
     0.0674099981127682,  -1.4340459501635281,  -1.4340459501635281,  -0.0000000000000000,  -0.0000000000000001,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0003364275824354,
    -0.0000000000000000,   0.0069795700710596,  -0.5140361149152392,   2.1591397431695505,   2.6937998841993207,
     4.6984758636102502,   2.9807907662378388,   0.7015391222557474,   1.5161116635199221,   1.3284500274338313,
     2.2520408051656227,  -0.7285847892739213,   0.8767886167460910,   1.0688009913637599,   1.7971807795394827,
    -0.5358628411743045,  -0.1101918318052870,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   1.1648256211196386,   0.1116222785865695,  -0.0000000000000001,
    -0.7389928339991191,  -0.6527772822764955,   2.5095629959661090,   3.3986261138240828,   2.9807907662378388,
     3.7436963584841090,   1.5161116635199221,   2.5461241562725809,   0.4566133209227201,   1.2253691290857067,
    -0.5599691986665823,   1.0302870613699144,   0.3823749166752148,   1.0044623015714094,  -0.4112125656361736,
    -0.2314624600919833,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0984088291831277,   0.0133987389938632,  -0.0000000000000000,  -0.0570207593897959,
    -0.5140361149152397,   2.1591397431695505,   2.6937998841993220,   0.7015391222557474,   1.5161116635199219,
     4.6984758636102502,   2.9807907662378392,  -1.3284500274338322,  -2.2520408051656235,  -0.8767886167460905,
     0.7285847892739212,   1.0688009913637599,   1.7971807795394830,  -0.1101918318052870,  -0.5358628411743045,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,
    -1.1648256211196382,   0.1116222785865690,   0.0000000000000000,  -0.7389928339991202,  -0.6527772822764957,
     2.5095629959661081,   3.3986261138240828,   1.5161116635199221,   2.5461241562725809,   2.9807907662378392,
     3.7436963584841081,  -0.4566133209227199,  -1.2253691290857069,  -1.0302870613699142,   0.5599691986665823,
     0.3823749166752150,   1.0044623015714096,  -0.2314624600919833,  -0.4112125656361735,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0984088291831277,
     0.0133987389938635,   0.0000000000000000,  -0.0570207593897960,   0.0000000000000005,  -0.0000000000000000,
     0.0000000000000001,   1.3284500274338322,   0.4566133209227202,  -1.3284500274338327,  -0.4566133209227200,
     8.4874200961780133,   3.1586526036021225,  -1.3363761800341809,  -1.3363761800341816,  -0.0000000000000003,
    -0.0000000000000002,  -1.9981360354136575,   1.9981360354136577,  -0.0000000000000004,  -0.0000000000000002,
    -0.0000000000000001,  -0.0000000000000001,   0.0000000000000001,  -0.0162208529550673,   0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,  -0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
     2.2520408051656227,   1.2253691290857067,  -2.2520408051656231,  -1.2253691290857069,   3.1586526036021225,
     4.6468099021321869,  -0.0188551414255343,  -0.0188551414255340,  -0.0000000000000002,  -0.0000000000000001,
    -1.4246081676669780,   1.4246081676669777,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,  -0.0290807180424306,   0.0000000000000001,  -0.0000000000000001,
     0.0000000000000001,   0.7933779079230701,  -2.6086424092976292,  -1.8679699639745317,  -0.7285847892739215,
    -0.5599691986665823,  -0.8767886167460907,  -1.0302870613699144,  -1.3363761800341809,  -0.0188551414255343,
     5.1625290629115046,  -1.4483843395850213,  -2.0267384197998890,  -1.4593371224984597,   0.2495270972564686,
     0.1116711153347607,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -1.2690027214023281,  -1.1668467253894739,   0.0000000000000000,  -0.0123557686110778,
    -0.7933779079230706,   2.6086424092976310,   1.8679699639745322,   0.8767886167460910,   1.0302870613699144,
     0.7285847892739212,   0.5599691986665825,  -1.3363761800341816,  -0.0188551414255340,  -1.4483843395850213,
     5.1625290629115046,   2.0267384197998899,   1.4593371224984597,  -0.1116711153347607,  -0.2495270972564684,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -1.2690027214023272,   1.1668467253894748,  -0.0000000000000000,   0.0123557686110778,  -0.0040937126161196,
     0.1239397763729298,   0.0756373318876606,   1.0688009913637602,   0.3823749166752148,   1.0688009913637602,
     0.3823749166752151,  -0.0000000000000003,  -0.0000000000000002,  -2.0267384197998890,   2.0267384197998899,
     8.4141882312829530,   3.1289263564853074,  -0.3900592393287933,  -0.3900592393287938,  -0.0000000000000001,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0570205597096819,  -0.0000000000000001,   0.0337887649434693,   0.0012546716139458,   0.0556812082323448,
     0.0674099981127681,   1.7971807795394830,   1.0044623015714094,   1.7971807795394830,   1.0044623015714096,
    -0.0000000000000002,  -0.0000000000000001,  -1.4593371224984597,   1.4593371224984597,   3.1289263564853069,
     4.6368894032500183,   0.6831145527113701,   0.6831145527113699,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0448380769025145,
    -0.0000000000000000,   0.0341073697712625,   0.6242803742030971,  -2.0184051429381134,  -1.4340459501635276,
    -0.5358628411743045,  -0.4112125656361736,  -0.1101918318052871,  -0.2314624600919833,  -1.9981360354136575,
    -1.4246081676669780,   0.2495270972564686,  -0.1116711153347608,  -0.3900592393287933,   0.6831145527113701,
     5.0140198152668827,   0.2697082126478946,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.3666260309293115,   1.2811551100351930,   0.0000000000000001,
     1.2571297739771972,   0.6242803742030976,  -2.0184051429381147,  -1.4340459501635281,  -0.1101918318052870,
    -0.2314624600919833,  -0.5358628411743045,  -0.4112125656361736,   1.9981360354136575,   1.4246081676669777,
     0.1116711153347607,  -0.2495270972564684,  -0.3900592393287938,   0.6831145527113699,   0.2697082126478946,
     5.0140198152668845,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.3666260309293113,   1.2811551100351928,  -0.0000000000000000,   1.2571297739771972,
    -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000004,  -0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,   0.0000000000000001,  -0.0000000000000001,
     8.3211161111573198,   3.1016068336523843,   1.1660593068764644,   1.1660593068764649,  -0.0000000000000002,
     0.0000000000000000,   0.0000000000000000,   0.0501425346888325,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000002,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,   3.1016068336523843,
     4.6339027900460161,   1.7953727912815483,   1.7953727912815480,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0414598716562782,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   1.1660593068764644,   1.7953727912815483,
     4.8073849654371248,   0.2039462944955927,   1.6011125976176674,   0.0000000000000000,   0.0000000000000000,
     1.2671672733227579,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   1.1660593068764646,   1.7953727912815480,   0.2039462944955927,
     4.8073849654371248,  -1.6011125976176674,   0.0000000000000000,   0.0000000000000000,   1.2671672733227575,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000002,  -0.0000000000000001,   1.6011125976176674,  -1.6011125976176674,
     7.2828583222221344,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000003,  -0.0000000000000000,  -0.0000000000000000,   1.1648256211196386,   0.0984088291831277,
    -1.1648256211196382,  -0.0984088291831277,  -0.0162208529550673,  -0.0290807180424306,  -1.2690027214023281,
    -1.2690027214023272,  -0.0000000000000000,   0.0000000000000000,  -0.3666260309293115,   0.3666260309293113,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     7.2760112552851215,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0002408453277242,
     0.0023221868883481,   0.0003364275824351,   0.1116222785865695,   0.0133987389938630,   0.1116222785865693,
     0.0133987389938635,   0.0000000000000001,   0.0000000000000001,  -1.1668467253894739,   1.1668467253894748,
     0.0570205597096820,   0.0448380769025145,   1.2811551100351930,   1.2811551100351930,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     7.2881208080016986,   0.0000000000000000,   0.0189989538120016,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,   0.0501425346888325,   0.0414598716562782,
     1.2671672733227579,   1.2671672733227575,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     7.2725545859718439,   0.0000000000000000,   0.0045540032341425,   0.0020012845931223,   0.0069795700710595,
    -0.7389928339991190,  -0.0570207593897962,  -0.7389928339991196,  -0.0570207593897961,   0.0000000000000002,
     0.0000000000000001,  -0.0123557686110778,   0.0123557686110778,   0.0337887649434693,   0.0341073697712625,
     1.2571297739771972,   1.2571297739771974,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0189989538120020,   0.0000000000000000,
     7.2757950881870475};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals.size()));
    
    vlxtest::compare(jkvals, fock.getFock(0));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestJKForH2Se)
{
    auto mh2o = vlxmol::getMoleculeH2Se();
    
    auto mbas = vlxbas::getMolecularBasisForH2Se();
    
    const int nrows = 42, ncols = 42;

    std::vector<double> denvals{   1.130133663073508,  -0.421201198203008,   0.163219785957157, 
                                  -0.028083939914295,  -0.017723587480734,  -0.003180613872843, 
                                   0.001730560940715,  -0.003180614694659,   0.001730558116180, 
                                   0.000618302770613,  -0.001888574599778,   0.003299547500286, 
                                   0.003188876216579,   0.000094069158093,   0.000801073671636, 
                                   0.000618302233159,  -0.001888572937492,   0.003299544722008, 
                                   0.003188872438378,   0.000801073689138,   0.000094069372157, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,   0.000000000000001,  -0.000000000000000, 
                                  -0.000022026671931,  -0.000020407271410,   0.000140759814639, 
                                   0.000043315659442,   0.000017944764312,   0.000120833156143, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000075025796807,  -0.000031081149809,  -0.000209287285852, 
                                  -0.421201198203008,   1.370865976440662,  -0.609537008437353, 
                                   0.081773490902150,   0.046444457025718,   0.008893306069080, 
                                  -0.004106816737625,   0.008893308470274,  -0.004106809039648, 
                                  -0.001620076972551,   0.004956467896001,  -0.008557392644885, 
                                  -0.008890443640922,  -0.000261172497288,  -0.001766469999201, 
                                  -0.001620075523043,   0.004956463408190,  -0.008557385077944, 
                                  -0.008890433506053,  -0.001766470017660,  -0.000261173100361, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000051703907782,   0.000072086125936,  -0.000387605492669, 
                                  -0.000131433999634,  -0.000067538181553,  -0.000218656402485, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000227652847569,   0.000116979093501,   0.000378718933310, 
                                   0.163219785957157,  -0.609537008437353,   1.529824812536579, 
                                  -0.433982810818172,  -0.206620554369395,  -0.038174427636490, 
                                   0.017696268411754,  -0.038174438011355,   0.017696229732907, 
                                   0.008375558988431,  -0.025397087113290,   0.042456812037313, 
                                   0.047747679045610,   0.001376245534378,   0.006661430551956, 
                                   0.008375551600561,  -0.025397064288295,   0.042456774296006, 
                                   0.047747626510755,   0.006661431029182,   0.001376248389117, 
                                   0.000000000000000,  -0.000000000000001,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,   0.000000000000000,  -0.000000000000001, 
                                  -0.000365831757144,  -0.000137604540654,   0.001809362463855, 
                                   0.000647745034643,  -0.000022674016300,   0.001214772996106, 
                                   0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                  -0.001121938574427,   0.000039272423106,  -0.002104022178383, 
                                  -0.028083939914295,   0.081773490902150,  -0.433982810818172, 
                                   0.629762989103894,   0.281259843491581,   0.049369169984606, 
                                  -0.024674216952522,   0.049369173711990,  -0.024674165869942, 
                                  -0.012958734091383,   0.038113164769331,  -0.061784160514848, 
                                  -0.070728161134257,  -0.002016918226736,  -0.008842128221279, 
                                  -0.012958721174751,   0.038113125272843,  -0.061784098346394, 
                                  -0.070728077927939,  -0.008842130026143,  -0.002016921966566, 
                                  -0.000000000000001,   0.000000000000001,  -0.000000000000001, 
                                  -0.000000000000002,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000276944176721,  -0.000115765146289,  -0.002500714568656, 
                                   0.000411726857927,   0.000764492171729,  -0.002086741990648, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000713113849490,  -0.001324135241077,   0.003614306588939, 
                                  -0.017723587480734,   0.046444457025718,  -0.206620554369395, 
                                   0.281259843491581,   0.148938075796371,  -0.019318652316124, 
                                  -0.030267550224105,  -0.019318624576403,  -0.030267488829977, 
                                  -0.012982746798948,   0.042661263584894,  -0.070248336349736, 
                                  -0.068920799098955,  -0.001870719226099,  -0.000351844741179, 
                                  -0.012982741044365,   0.042661244940319,  -0.070248301028742, 
                                  -0.068920735383946,  -0.000351843644392,  -0.001870720241717, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.001372255545177,   0.000355726748506,  -0.002556769037454, 
                                  -0.001963764727687,  -0.001149489843282,  -0.003707197825470, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.003401343051464,   0.001990972007775,   0.006421026457229, 
                                  -0.003180613872843,   0.008893306069080,  -0.038174427636490, 
                                   0.049369169984606,  -0.019318652316124,   0.175148764323550, 
                                   0.079864218604927,  -0.018869553056406,  -0.014768899346370, 
                                  -0.002015853565957,   0.005647168744329,  -0.009159681745087, 
                                   0.001608030732764,  -0.000312179418082,   0.000812364173577, 
                                   0.033140038159731,  -0.099328447577652,   0.152417610004031, 
                                   0.121402415127438,  -0.015095947766034,   0.003477979347193, 
                                   0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000001, 
                                  -0.000845794217369,  -0.000298311925869,   0.002416564651907, 
                                  -0.011045894924994,   0.003844966431369,   0.023365635072761, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001610311920115,   0.000264239260440,   0.002645162337113, 
                                   0.001730560940715,  -0.004106816737625,   0.017696268411754, 
                                  -0.024674216952522,  -0.030267550224105,   0.079864218604927, 
                                   0.040077046605490,  -0.014768852459195,  -0.006108897962559, 
                                  -0.000620636764477,   0.001298885377507,  -0.001496797477882, 
                                   0.004661027551196,  -0.000036599197126,   0.001221952149875, 
                                   0.015277162328067,  -0.049469461167416,   0.077391241198689, 
                                   0.063149430359726,  -0.006514551892468,   0.001812163570314, 
                                   0.000000000000001,  -0.000000000000002,   0.000000000000001, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                  -0.000284267234222,  -0.000088529538746,   0.001314099433636, 
                                  -0.008010315174323,   0.000956115588644,   0.011406655322777, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.002471327320525,  -0.000265593918285,   0.001182945165001, 
                                  -0.003180614694659,   0.008893308470274,  -0.038174438011355, 
                                   0.049369173711990,  -0.019318624576403,  -0.018869553056406, 
                                  -0.014768852459195,   0.175148743539541,   0.079864170087819, 
                                   0.033140037399829,  -0.099328445632606,   0.152417610847087, 
                                   0.121402429086300,   0.003477977517841,  -0.015095945731129, 
                                  -0.002015869486221,   0.005647217820743,  -0.009159757826493, 
                                   0.001607918207531,   0.000812365148129,  -0.000312185586182, 
                                   0.000000000000000,  -0.000000000000001,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                  -0.000845785724925,  -0.000298311471946,   0.002416525409529, 
                                   0.006917514460130,  -0.002151320197124,  -0.013973585958644, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.008760866137960,  -0.003197718772714,  -0.018912647100344, 
                                   0.001730558116180,  -0.004106809039648,   0.017696229732907, 
                                  -0.024674165869942,  -0.030267488829977,  -0.014768899346370, 
                                  -0.006108897962559,   0.079864170087819,   0.040076998781513, 
                                   0.015277152398327,  -0.049469430805724,   0.077391196252631, 
                                   0.063149397129303,   0.001812161180957,  -0.006514548050792, 
                                  -0.000620656899561,   0.001298947554462,  -0.001496894751448, 
                                   0.004660917662890,   0.001221956462592,  -0.000036603203484, 
                                   0.000000000000000,   0.000000000000001,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000284261910439,  -0.000088528954174,   0.001314078461734, 
                                   0.006145391326434,  -0.000248046219549,  -0.006727785473157, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.005701467356404,  -0.000960817046125,  -0.009286971810792, 
                                   0.000618302770613,  -0.001620076972551,   0.008375558988431, 
                                  -0.012958734091383,  -0.012982746798948,  -0.002015853565957, 
                                  -0.000620636764477,   0.033140037399829,   0.015277152398327, 
                                   1.178465527160369,  -0.466818953948649,   0.020186854801810, 
                                   0.029313308397886,   0.000836342819716,  -0.003675777996171, 
                                   0.000778917635252,  -0.002399430136192,   0.003717991032045, 
                                   0.005611275726618,   0.000249739852281,   0.000066378539561, 
                                   0.000000000000001,  -0.000000000000002,   0.000000000000005, 
                                  -0.000000000000004,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                   0.000010948336592,  -0.000173030779244,   0.000701338604507, 
                                   0.002009686499321,  -0.000206108338126,  -0.001783794803275, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.003108331878095,  -0.000215283455036,  -0.003403637205848, 
                                  -0.001888574599778,   0.004956467896001,  -0.025397087113290, 
                                   0.038113164769331,   0.042661263584894,   0.005647168744329, 
                                   0.001298885377507,  -0.099328445632606,  -0.049469430805724, 
                                  -0.466818953948649,   1.224915685134512,  -0.082893886908445, 
                                  -0.083407480522353,  -0.002248637021290,   0.008080590821895, 
                                  -0.002399430215637,   0.007385477222431,  -0.011446716358901, 
                                  -0.016737611900872,  -0.000508162073976,  -0.000243094071009, 
                                  -0.000000000000001,   0.000000000000004,  -0.000000000000013, 
                                   0.000000000000011,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000001,   0.000000000000000, 
                                   0.000030065658153,   0.000411111870806,  -0.002181316144919, 
                                  -0.006190517072232,   0.000425559520459,   0.006732861641002, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.009629717197440,   0.000574107583058,   0.011828390455359, 
                                   0.003299547500286,  -0.008557392644885,   0.042456812037313, 
                                  -0.061784160514848,  -0.070248336349736,  -0.009159681745087, 
                                  -0.001496797477882,   0.152417610847087,   0.077391196252631, 
                                   0.020186854801810,  -0.082893886908445,   0.152458847233981, 
                                   0.125935425485599,   0.003571222408570,  -0.012356800234144, 
                                   0.003717991698091,  -0.011446718015693,   0.017477220899588, 
                                   0.026059940498965,   0.000832557450151,   0.000414384296847, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000071582126691,   0.000035733111393,   0.002954940054824, 
                                   0.009900302018246,  -0.000280204317557,  -0.010284249217750, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                   0.015006588843529,  -0.000528137570271,  -0.017880579610272, 
                                   0.003188876216579,  -0.008890443640922,   0.047747679045610, 
                                  -0.070728161134257,  -0.068920799098955,   0.001608030732764, 
                                   0.004661027551196,   0.121402429086300,   0.063149397129303, 
                                   0.029313308397886,  -0.083407480522353,   0.125935425485599, 
                                   0.105573765431549,   0.002967704502851,  -0.009708171710812, 
                                   0.005611286888442,  -0.016737645914008,   0.026059990767888, 
                                   0.031468266379135,   0.000067014072075,   0.000627024806316, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.001253572369597,  -0.000387396864904,   0.002631608855800, 
                                   0.007091127614498,  -0.000108324513905,  -0.006868772939759, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.014478377160544,   0.000362026227885,  -0.014518238937695, 
                                   0.000094069158093,  -0.000261172497288,   0.001376245534378, 
                                  -0.002016918226736,  -0.001870719226099,  -0.000312179418082, 
                                  -0.000036599197126,   0.003477977517841,   0.001812161180957, 
                                   0.000836342819716,  -0.002248637021290,   0.003571222408570, 
                                   0.002967704502851,   0.000086754778785,  -0.000277980560876, 
                                   0.000066379126352,  -0.000243095914634,   0.000414387242345, 
                                   0.000627025718371,   0.000032740984094,   0.000012700542550, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001302607049501,   0.000415235514624,   0.000097860712113, 
                                   0.000166400227129,  -0.000029639128056,  -0.000243101843803, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000160889875395,  -0.000069906107095,  -0.000420405481079, 
                                   0.000801073671636,  -0.001766469999201,   0.006661430551956, 
                                  -0.008842128221279,  -0.000351844741179,   0.000812364173577, 
                                   0.001221952149875,  -0.015095945731129,  -0.006514548050792, 
                                  -0.003675777996171,   0.008080590821895,  -0.012356800234144, 
                                  -0.009708171710812,  -0.000277980560876,   0.001347985613512, 
                                   0.000249740956378,  -0.000508165400901,   0.000832562155507, 
                                   0.000067021931589,   0.000025058884343,   0.000032741496528, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000241430639577,   0.000081230344740,  -0.000187788595619, 
                                   0.001033888715980,   0.000676720916674,   0.001157766485042, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.001517969562850,   0.000993351268582,   0.001601900755406, 
                                   0.000618302233159,  -0.001620075523043,   0.008375551600561, 
                                  -0.012958721174751,  -0.012982741044365,   0.033140038159731, 
                                   0.015277162328067,  -0.002015869486221,  -0.000620656899561, 
                                   0.000778917635252,  -0.002399430215637,   0.003717991698091, 
                                   0.005611286888442,   0.000066379126352,   0.000249740956378, 
                                   1.178465524742984,  -0.466818946502035,   0.020186842555999, 
                                   0.029313293317148,  -0.003675778194075,   0.000836343003130, 
                                  -0.000000000000000,   0.000000000000002,  -0.000000000000005, 
                                   0.000000000000004,   0.000000000000000,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000001,   0.000000000000000, 
                                   0.000010946357627,  -0.000173030789848,   0.000701345321657, 
                                  -0.003696740354257,   0.000289496032902,   0.003839536617017, 
                                   0.000000000000000,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000186276328918,   0.000070853682100,  -0.000157003335327, 
                                  -0.001888572937492,   0.004956463408190,  -0.025397064288295, 
                                   0.038113125272843,   0.042661244940319,  -0.099328447577652, 
                                  -0.049469461167416,   0.005647217820743,   0.001298947554462, 
                                  -0.002399430136192,   0.007385477222431,  -0.011446718015693, 
                                  -0.016737645914008,  -0.000243095914634,  -0.000508165400901, 
                                  -0.466818946502035,   1.224915662215716,  -0.082893849540772, 
                                  -0.083407433589175,   0.008080591325009,  -0.002248637590973, 
                                   0.000000000000001,  -0.000000000000004,   0.000000000000012, 
                                  -0.000000000000011,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000030071829446,   0.000411111750795,  -0.002181336595023, 
                                   0.011434846261419,  -0.000709973469283,  -0.013610126939056, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000546294688455,  -0.000081492393346,   0.000083354665401, 
                                   0.003299544722008,  -0.008557385077944,   0.042456774296006, 
                                  -0.061784098346394,  -0.070248301028742,   0.152417610004031, 
                                   0.077391241198689,  -0.009159757826493,  -0.001496894751448, 
                                   0.003717991032045,  -0.011446716358901,   0.017477220899588, 
                                   0.026059990767888,   0.000414387242345,   0.000832562155507, 
                                   0.020186842555999,  -0.082893849540772,   0.152458789641470, 
                                   0.125935348202675,  -0.012356800232682,   0.003571223251788, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                   0.000071572867280,   0.000035732026009,   0.002954970599919, 
                                  -0.017946250700459,   0.000597481893273,   0.020627177398946, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001070631975989,  -0.000021404956360,  -0.000033850739019, 
                                   0.003188872438378,  -0.008890433506053,   0.047747626510755, 
                                  -0.070728077927939,  -0.068920735383946,   0.121402415127438, 
                                   0.063149430359726,   0.001607918207531,   0.004660917662890, 
                                   0.005611275726618,  -0.016737611900872,   0.026059940498965, 
                                   0.031468266379135,   0.000627025718371,   0.000067021931589, 
                                   0.029313293317148,  -0.083407433589175,   0.125935348202675, 
                                   0.105573669087531,  -0.009708170597064,   0.002967704402686, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000001, 
                                  -0.001253577092852,  -0.000387396982805,   0.002631630078765, 
                                  -0.016084211002539,  -0.000259359911650,   0.016007566734574, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.001098079482870,   0.000274826054815,  -0.001310565926614, 
                                   0.000801073689138,  -0.001766470017660,   0.006661431029182, 
                                  -0.008842130026143,  -0.000351843644392,  -0.015095947766034, 
                                  -0.006514551892468,   0.000812365148129,   0.001221956462592, 
                                   0.000249739852281,  -0.000508162073976,   0.000832557450151, 
                                   0.000067014072075,   0.000032740984094,   0.000025058884343, 
                                  -0.003675778194075,   0.008080591325009,  -0.012356800232682, 
                                  -0.009708170597064,   0.001347985830443,  -0.000277980666316, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000241431357552,   0.000081230393604,  -0.000187791776830, 
                                  -0.001831543613629,  -0.001198627805073,  -0.001966170631561, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000136388941926,  -0.000089381954102,  -0.000201705293681, 
                                   0.000094069372157,  -0.000261173100361,   0.001376248389117, 
                                  -0.002016921966566,  -0.001870720241717,   0.003477979347193, 
                                   0.001812163570314,  -0.000312185586182,  -0.000036603203484, 
                                   0.000066378539561,  -0.000243094071009,   0.000414384296847, 
                                   0.000627024806316,   0.000012700542550,   0.000032741496528, 
                                   0.000836343003130,  -0.002248637590973,   0.003571223251788, 
                                   0.002967704402686,  -0.000277980666316,   0.000086754861197, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001302607031553,   0.000415235554460,   0.000097861422749, 
                                  -0.000222535588470,   0.000075360027722,   0.000485633714586, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000063662164137,  -0.000009284664369,   0.000000330484401, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000001,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   1.182678290332052,  -0.479757793827607,   0.041285187963634, 
                                   0.067876283335214,   0.001692541819219,   0.001692542953229, 
                                  -0.000225139144518,  -0.000202692415904,   0.001470546131387, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000002, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000225136772336,  -0.000202692352925,   0.001470537087093, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000001, 
                                  -0.000000000000002,  -0.000000000000001,   0.000000000000001, 
                                  -0.000000000000002,   0.000000000000004,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000001,   0.000000000000000, 
                                   0.000000000000002,  -0.000000000000004,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000001, 
                                  -0.479757793827607,   1.264522211729289,  -0.146197983684989, 
                                  -0.204331974218631,  -0.005048126414011,  -0.005048129944513, 
                                   0.000762798615986,   0.000502588845533,  -0.004536012728482, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000004, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000001, 
                                   0.000762791185729,   0.000502588840465,  -0.004535984945703, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                  -0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000001,   0.000000000000001, 
                                   0.000000000000005,  -0.000000000000013,   0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000005,   0.000000000000012,   0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000001, 
                                   0.041285187963634,  -0.146197983684989,   0.250998435354553, 
                                   0.319446049748368,   0.008031247147271,   0.008031252741148, 
                                  -0.001132245966276,  -0.000225787746325,   0.006702318663475, 
                                   0.000000000000002,   0.000000000000001,  -0.000000000000007, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001132234512109,  -0.000225787174411,   0.006702275982904, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                  -0.000000000000002,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000001,   0.000000000000001, 
                                  -0.000000000000004,   0.000000000000011,   0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000001,  -0.000000000000000, 
                                   0.000000000000004,  -0.000000000000011,   0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000001, 
                                   0.067876283335214,  -0.204331974218631,   0.319446049748368, 
                                   0.406915693597961,   0.010226527456445,   0.010226534583722, 
                                  -0.002769992432593,  -0.000713293227986,   0.008511062503724, 
                                   0.000000000000002,   0.000000000000001,  -0.000000000000009, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.002769979285124,  -0.000713292968842,   0.008511008099177, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001692541819219,  -0.005048126414011,   0.008031247147271, 
                                   0.010226527456445,   0.000259830151464,   0.000256966097269, 
                                  -0.000082669309901,  -0.000022450870524,   0.000213810076855, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.001352621977624,   0.000432533769123,   0.000245694632101, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001692542953229,  -0.005048129944513,   0.008031252741148, 
                                   0.010226534583722,   0.000256966097269,   0.000259830510952, 
                                   0.001352621913991,   0.000432533842825,   0.000245696155911, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000082669257932,  -0.000022450949387,   0.000213808853560, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000001,   0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000001,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000225139144518,   0.000762798615986,  -0.001132245966276, 
                                  -0.002769992432593,  -0.000082669309901,   0.001352621913991, 
                                   0.719439404527231,   0.227849549759394,   0.016037592810119, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000001, 
                                   0.000203139720555,  -0.000147039491035,   0.000059291330793, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000001,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000202692415904,   0.000502588845533,  -0.000225787746325, 
                                  -0.000713293227986,  -0.000022450870524,   0.000432533842825, 
                                   0.227849549759394,   0.072161160746783,   0.005082524244797, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000147039537845,  -0.000113343843586,   0.000017431681913, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.001470546131387,  -0.004536012728482,   0.006702318663475, 
                                   0.008511062503724,   0.000213810076855,   0.000245696155911, 
                                   0.016037592810119,   0.005082524244797,   0.000538194691525, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000059291248450,   0.000017431596934,   0.000183225176131, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000022026671931,   0.000051703907782,  -0.000365831757144, 
                                   0.000276944176721,   0.001372255545177,  -0.000845794217369, 
                                  -0.000284267234222,  -0.000845785724925,  -0.000284261910439, 
                                   0.000010948336592,   0.000030065658153,   0.000071582126691, 
                                  -0.001253572369597,   0.001302607049501,   0.000241430639577, 
                                   0.000010946357627,   0.000030071829446,   0.000071572867280, 
                                  -0.001253577092852,   0.000241431357552,   0.001302607031553, 
                                   0.000000000000000,  -0.000000000000001,   0.000000000000002, 
                                   0.000000000000002,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.718179558358057,   0.228740396024795,   0.015495336685894, 
                                   0.000155172934398,  -0.000097329117055,   0.000063234182785, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000268764597324,   0.000168578264546,  -0.000109527260864, 
                                  -0.000020407271410,   0.000072086125936,  -0.000137604540654, 
                                  -0.000115765146289,   0.000355726748506,  -0.000298311925869, 
                                  -0.000088529538746,  -0.000298311471946,  -0.000088528954174, 
                                  -0.000173030779244,   0.000411111870806,   0.000035733111393, 
                                  -0.000387396864904,   0.000415235514624,   0.000081230344740, 
                                  -0.000173030789848,   0.000411111750795,   0.000035732026009, 
                                  -0.000387396982805,   0.000081230393604,   0.000415235554460, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.228740396024795,   0.072854325900428,   0.004935442169249, 
                                  -0.000062348827983,  -0.000066882198390,   0.000017631420949, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000107990865396,   0.000115842858866,  -0.000030538799582, 
                                   0.000140759814639,  -0.000387605492669,   0.001809362463855, 
                                  -0.002500714568656,  -0.002556769037454,   0.002416564651907, 
                                   0.001314099433636,   0.002416525409529,   0.001314078461734, 
                                   0.000701338604507,  -0.002181316144919,   0.002954940054824, 
                                   0.002631608855800,   0.000097860712113,  -0.000187788595619, 
                                   0.000701345321657,  -0.002181336595023,   0.002954970599919, 
                                   0.002631630078765,  -0.000187791776830,   0.000097861422749, 
                                  -0.000000000000002,   0.000000000000004,  -0.000000000000007, 
                                  -0.000000000000009,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.015495336685894,   0.004935442169249,   0.000437716615756, 
                                  -0.000127045357179,   0.000005571250277,   0.000180794944795, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000220041635011,  -0.000009649248029,  -0.000313136824430, 
                                   0.000043315659442,  -0.000131433999634,   0.000647745034643, 
                                   0.000411726857927,  -0.001963764727687,  -0.011045894924994, 
                                  -0.008010315174323,   0.006917514460130,   0.006145391326434, 
                                   0.002009686499321,  -0.006190517072232,   0.009900302018246, 
                                   0.007091127614498,   0.000166400227129,   0.001033888715980, 
                                  -0.003696740354257,   0.011434846261419,  -0.017946250700459, 
                                  -0.016084211002539,  -0.001831543613629,  -0.000222535588470, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.000155172934398,  -0.000062348827983,  -0.000127045357179, 
                                   0.720330357209242,   0.228694828282184,   0.009715316186519, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000026403869134,   0.000565114444114,  -0.002553941034585, 
                                   0.000017944764312,  -0.000067538181553,  -0.000022674016300, 
                                   0.000764492171729,  -0.001149489843282,   0.003844966431369, 
                                   0.000956115588644,  -0.002151320197124,  -0.000248046219549, 
                                  -0.000206108338126,   0.000425559520459,  -0.000280204317557, 
                                  -0.000108324513905,  -0.000029639128056,   0.000676720916674, 
                                   0.000289496032902,  -0.000709973469283,   0.000597481893273, 
                                  -0.000259359911650,  -0.001198627805073,   0.000075360027722, 
                                  -0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000097329117055,  -0.000066882198390,   0.000005571250277, 
                                   0.228694828282184,   0.072990727605863,   0.004310613587103, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000565114574565,   0.000509665422617,  -0.000308111982823, 
                                   0.000120833156143,  -0.000218656402485,   0.001214772996106, 
                                  -0.002086741990648,  -0.003707197825470,   0.023365635072761, 
                                   0.011406655322777,  -0.013973585958644,  -0.006727785473157, 
                                  -0.001783794803275,   0.006732861641002,  -0.010284249217750, 
                                  -0.006868772939759,  -0.000243101843803,   0.001157766485042, 
                                   0.003839536617017,  -0.013610126939056,   0.020627177398946, 
                                   0.016007566734574,  -0.001966170631561,   0.000485633714586, 
                                  -0.000000000000000,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000063234182785,   0.000017631420949,   0.000180794944795, 
                                   0.009715316186519,   0.004310613587103,   0.004109639117620, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.002553940640306,  -0.000308111921880,   0.001449136219520, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                  -0.000225136772336,   0.000762791185729,  -0.001132234512109, 
                                  -0.002769979285124,   0.001352621977624,  -0.000082669257932, 
                                   0.000203139720555,  -0.000147039537845,   0.000059291248450, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.719439400444076,   0.227849552140461,   0.016037593124298, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000202692352925,   0.000502588840465,  -0.000225787174411, 
                                  -0.000713292968842,   0.000432533769123,  -0.000022450949387, 
                                  -0.000147039491035,  -0.000113343843586,   0.000017431596934, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.227849552140461,   0.072161162661593,   0.005082524322250, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001470537087093,  -0.004535984945703,   0.006702275982904, 
                                   0.008511008099177,   0.000245694632101,   0.000213808853560, 
                                   0.000059291330793,   0.000017431681913,   0.000183225176131, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.016037593124298,   0.005082524322250,   0.000538192402583, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000075025796807,   0.000227652847569,  -0.001121938574427, 
                                  -0.000713113849490,   0.003401343051464,  -0.001610311920115, 
                                  -0.002471327320525,   0.008760866137960,   0.005701467356404, 
                                   0.003108331878095,  -0.009629717197440,   0.015006588843529, 
                                   0.014478377160544,   0.000160889875395,   0.001517969562850, 
                                  -0.000186276328918,   0.000546294688455,  -0.001070631975989, 
                                   0.001098079482870,  -0.000136388941926,  -0.000063662164137, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                  -0.000268764597324,   0.000107990865396,   0.000220041635011, 
                                   0.000026403869134,   0.000565114574565,  -0.002553940640306, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.720299866650672,   0.228042290066672,   0.012664356596310, 
                                  -0.000031081149809,   0.000116979093501,   0.000039272423106, 
                                  -0.001324135241077,   0.001990972007775,   0.000264239260440, 
                                  -0.000265593918285,  -0.003197718772714,  -0.000960817046125, 
                                  -0.000215283455036,   0.000574107583058,  -0.000528137570271, 
                                   0.000362026227885,  -0.000069906107095,   0.000993351268582, 
                                   0.000070853682100,  -0.000081492393346,  -0.000021404956360, 
                                   0.000274826054815,  -0.000089381954102,  -0.000009284664369, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000168578264546,   0.000115842858866,  -0.000009649248029, 
                                   0.000565114444114,   0.000509665422617,  -0.000308111921880, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.228042290066672,   0.072402216659709,   0.004666391389492, 
                                  -0.000209287285852,   0.000378718933310,  -0.002104022178383, 
                                   0.003614306588939,   0.006421026457229,   0.002645162337113, 
                                   0.001182945165001,  -0.018912647100344,  -0.009286971810792, 
                                  -0.003403637205848,   0.011828390455359,  -0.017880579610272, 
                                  -0.014518238937695,  -0.000420405481079,   0.001601900755406, 
                                  -0.000157003335327,   0.000083354665401,  -0.000033850739019, 
                                  -0.001310565926614,  -0.000201705293681,   0.000000330484401, 
                                   0.000000000000000,  -0.000000000000001,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000109527260864,  -0.000030538799582,  -0.000313136824430, 
                                  -0.002553941034585,  -0.000308111982823,   0.001449136219520, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.012664356596310,   0.004666391389492,   0.002436318893508};

    const CAODensityMatrix dmat({CDenseMatrix(denvals, nrows, ncols)}, denmat::rest);
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-12, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> jkvals { 117.857768727795246,  32.609856553679066,   2.524727088717546, 
                                   2.073474526407391,   0.988378747689814,   0.070615494990968, 
                                   0.410428569563903,   0.070615494985720,   0.410428569555598, 
                                   0.000842591302766,   0.000217983428433,  -0.000028126782322, 
                                  -0.000044965264832,  -0.000009489558653,  -0.042181454619972, 
                                   0.000842593295200,   0.000217983610989,  -0.000028126753648, 
                                  -0.000044965243435,  -0.042181454636970,  -0.000009489543004, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000061471353724,  -0.000017054701545,  -0.000023562087915, 
                                  -0.000086370708979,   0.000020215809996,   0.000025296519345, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000149598299673,  -0.000035014772696,  -0.000043814821989, 
                                  32.609856553679066,  83.061592388366392,  28.053818148794846, 
                                   9.630662717298723,   4.651141062240152,   0.349133848186671, 
                                   1.940759666699717,   0.349133848183611,   1.940759666684476, 
                                   0.009012834364007,   0.006168812356577,   0.000809943079483, 
                                   0.000338494350981,   0.000024184292712,  -0.228542649817976, 
                                   0.009012841280139,   0.006168816971486,   0.000809943562283, 
                                   0.000338494442472,  -0.228542649784721,   0.000024184211660, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000609586843736,   0.000176563881884,   0.000126609750030, 
                                  -0.000413474804192,  -0.000172795677023,  -0.000119281271301, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000716156331716,   0.000299290276865,   0.000206600896872, 
                                   2.524727088717546,  28.053818148794846,  49.335835244927836, 
                                  21.923928029598066,  11.472376364930529,   1.129502461261824, 
                                   4.932417792128861,   1.129502460912963,   4.932417791467648, 
                                   0.004746168075247,   0.031549521077770,   0.016293496866849, 
                                   0.004647004869250,   0.000210118889103,  -1.119730324395056, 
                                   0.004746172060157,   0.031549533914077,   0.016293503649725, 
                                   0.004647006975591,  -1.119730324996343,   0.000210119210717, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.002444639735731,   0.002268996701883,   0.000432752350941, 
                                  -0.000506775192747,   0.000219162289004,   0.000417938472960, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000877753343313,  -0.000379608509687,  -0.000723892735289, 
                                   2.073474526407391,   9.630662717298723,  21.923928029598066, 
                                  25.350579738053760,  18.181216066061204,   3.781290991056523, 
                                   8.980440496103746,   3.781290983500770,   8.980440493545085, 
                                   0.001941730663084,   0.024336446439559,   0.093284159932147, 
                                   0.060013420525597,   0.011879140551586,  -4.213124948300206, 
                                   0.001941731749484,   0.024336455607841,   0.093284155524529, 
                                   0.060013415382703,  -4.213124950236165,   0.011879145582388, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.001564157786959,   0.002334207404951,   0.000710644790127, 
                                   0.000668016737652,   0.006146421749064,   0.004491677471710, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001157049994426,  -0.010645912861443,  -0.007779790243670, 
                                   0.988378747689814,   4.651141062240152,  11.472376364930529, 
                                  18.181216066061204,  16.764404431581550,   4.614373743957179, 
                                   9.537339146235455,   4.614373732026622,   9.537339138512788, 
                                   0.001073246856103,   0.014341178341603,   0.070445711693932, 
                                   0.071260269217909,   0.014493044396626,  -3.229902960511601, 
                                   0.001073247341955,   0.014341183236479,   0.070445712447604, 
                                   0.071260274651458,  -3.229902963778099,   0.014493051852879, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000769701874973,   0.000277029891903,   0.000224777045591, 
                                   0.001024630991709,   0.004368422822374,   0.003695911256960, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001774720064801,  -0.007566322787821,  -0.006401483381471, 
                                   0.070615494990968,   0.349133848186671,   1.129502461261824, 
                                   3.781290991056523,   4.614373743957179,  12.494899365146384, 
                                   8.360529848383344,   0.336023170678084,   2.102349246882248, 
                                   0.000127645919132,   0.001253310409837,   0.016194663627884, 
                                   0.017777625663503,   0.028427560417647,  -0.274763220927076, 
                                   0.108393420466371,   1.049467095692840,   4.951046387476951, 
                                   6.489587471862320,  -2.411542820879911,   0.210945451244569, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   0.001095892727442,   0.006392394292714,   0.015940873488793, 
                                   0.339359226095313,   2.055207811276272,   5.992674202150579, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000478812130913,  -0.002793953214152,  -0.006267772523248, 
                                   0.410428569563903,   1.940759666699717,   4.932417792128861, 
                                   8.980440496103746,   9.537339146235455,   8.360529848383344, 
                                  11.526892444130739,   2.102349253166844,   5.328926206342184, 
                                   0.000459702338025,   0.006410775058633,   0.031599181184735, 
                                   0.037503924249580,   0.024001211170335,  -1.460364184852949, 
                                   0.170382809860090,   1.407733150725250,   4.969353402478769, 
                                   7.086329590701152,  -2.071428213083649,   0.908140072927748, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.001584984901569,   0.008799788644018,   0.026994160232387, 
                                   0.146399499111725,   0.783040291434189,   2.391307763904538, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000457563103277,  -0.004705709143368,  -0.006567026245229, 
                                   0.070615494985720,   0.349133848183611,   1.129502460912963, 
                                   3.781290983500770,   4.614373732026622,   0.336023170678084, 
                                   2.102349253166844,  12.494899360057213,   8.360529846113939, 
                                   0.108393420356982,   1.049467095120744,   4.951046383113243, 
                                   6.489587469881652,   0.210945457856467,  -2.411542824056427, 
                                   0.000127646014754,   0.001253309998293,   0.016194711579447, 
                                   0.017777670836737,  -0.274763234814444,   0.028427584819176, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001095888831268,   0.006392412933216,   0.015940911897444, 
                                  -0.170094273652441,  -1.025184276423810,  -2.990909067453096, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.293654303033034,  -1.781259152766736,  -5.192941986059282, 
                                   0.410428569555598,   1.940759666684476,   4.932417791467648, 
                                   8.980440493545085,   9.537339138512788,   2.102349246882248, 
                                   5.328926206342184,   8.360529846113939,  11.526892434938365, 
                                   0.170382809624717,   1.407733149096074,   4.969353391310212, 
                                   7.086329579693339,   0.908140074992716,  -2.071428212879054, 
                                   0.000459702590491,   0.006410776739697,   0.031599208092944, 
                                   0.037503949583236,  -1.460364194555921,   0.024001231433518, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.001584982114382,   0.008799801492873,   0.026994187973493, 
                                  -0.072803487374764,  -0.387444885118817,  -1.189966682029321, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.127014466868015,  -0.680485639467545,  -2.074216786397519, 
                                   0.000842591302766,   0.009012834364007,   0.004746168075247, 
                                   0.001941730663084,   0.001073246856103,   0.000127645919132, 
                                   0.000459702338025,   0.108393420356982,   0.170382809624716, 
                                  88.570583017168104,  27.486402643454905,   2.241362565695363, 
                                   0.504020528274561,   0.009894993596382,  -0.112475532457305, 
                                   0.000238502220001,   0.000495328366326,  -0.000027952211524, 
                                  -0.000051989418100,  -0.000095817604301,  -0.000048611052021, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.004784482472229,   0.000290215691857,  -0.000129138994621, 
                                  -0.003123188017088,  -0.000363437337951,  -0.000043508421716, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.005333639289730,  -0.000611069646923,  -0.000068501257130, 
                                   0.000217983428433,   0.006168812356577,   0.031549521077770, 
                                   0.024336446439559,   0.014341178341603,   0.001253310409837, 
                                   0.006410775058633,   1.049467095120744,   1.407733149096074, 
                                  27.486402643454905,  49.453801234465701,  14.183699395194147, 
                                   3.968568018808944,   0.173318590010071,  -1.468925416182689, 
                                   0.000495328373318,   0.001517434749548,   0.001988673745953, 
                                   0.000791320001368,  -0.001138011804650,   0.000719210037112, 
                                   0.000000000000000,  -0.000000000000001,   0.000000000000003, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.021881117305557,   0.018646731840302,   0.006679933703819, 
                                  -0.012325870260963,  -0.014132209170082,  -0.005488656716481, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.021590764066772,  -0.024271575806446,  -0.009333029527092, 
                                  -0.000028126782322,   0.000809943079483,   0.016293496866849, 
                                   0.093284159932147,   0.070445711693932,   0.016194663627884, 
                                   0.031599181184735,   4.951046383113244,   4.969353391310211, 
                                   2.241362565695363,  14.183699395194147,  22.571344785914619, 
                                  12.023038675794472,   1.546998890448782,  -6.497783334802272, 
                                  -0.000027952187760,   0.001988673917288,  -0.013203312292400, 
                                  -0.014117042311481,  -0.012890121010359,   0.009008658575075, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.009910061948506,   0.049662846692087,   0.059880307234960, 
                                  -0.013056803863772,  -0.023581836836335,  -0.002843018329669, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                  -0.021420903118261,  -0.042008814364454,  -0.009965949366088, 
                                  -0.000044965264832,   0.000338494350981,   0.004647004869250, 
                                   0.060013420525597,   0.071260269217909,   0.017777625663503, 
                                   0.037503924249580,   6.489587469881652,   7.086329579693339, 
                                   0.504020528274561,   3.968568018808944,  12.023038675794472, 
                                  13.155513636060956,   2.202737447851170,  -3.289732595608442, 
                                  -0.000051989368161,   0.000791319919334,  -0.014117030728298, 
                                  -0.012955670024300,  -0.009835416240766,   0.016555725864068, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.003256414288573,   0.023710601474494,   0.061694137210839, 
                                  -0.007763757262002,  -0.001169677615586,   0.003402275017596, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                  -0.012344088175682,  -0.003741457558729,  -0.000345604330298, 
                                  -0.000009489558653,   0.000024184292712,   0.000210118889103, 
                                   0.011879140551586,   0.014493044396626,   0.028427560417647, 
                                   0.024001211170335,   0.210945457856467,   0.908140074992716, 
                                   0.009894993596382,   0.173318590010071,   1.546998890448782, 
                                   2.202737447851170,  12.419032458969346,  -0.159831231850929, 
                                  -0.000048610991375,   0.000719210076361,   0.009008644164129, 
                                   0.016555707808128,  -0.008618413704740,   0.172520373971464, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                   0.000000000000002,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.124757737865327,   1.205683686360233,   4.105360213977148, 
                                  -0.000751132885379,   0.003498861594793,   0.019258015062522, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.001470753029586,   0.000782032192792,   0.003669577771402, 
                                  -0.042181454619972,  -0.228542649817976,  -1.119730324395056, 
                                  -4.213124948300206,  -3.229902960511601,  -0.274763220927076, 
                                  -1.460364184852949,  -2.411542824056427,  -2.071428212879054, 
                                  -0.112475532457305,  -1.468925416182689,  -6.497783334802272, 
                                  -3.289732595608442,  -0.159831231850929,  13.645886474024312, 
                                  -0.000095817559292,  -0.001138011875843,  -0.012890119827290, 
                                  -0.009835410969445,   0.217730366498301,  -0.008618419585509, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.001014707749876,  -0.008285933389673,  -0.013892290436904, 
                                   0.324269190108425,   2.033317926464628,   2.703082478375381, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.561514481324687,   3.524892957702942,   4.685834594993277, 
                                   0.000842593295200,   0.009012841280139,   0.004746172060157, 
                                   0.001941731749484,   0.001073247341955,   0.108393420466371, 
                                   0.170382809860090,   0.000127646014754,   0.000459702590491, 
                                   0.000238502220001,   0.000495328373318,  -0.000027952187760, 
                                  -0.000051989368161,  -0.000048610991375,  -0.000095817559292, 
                                  88.570583015645894,  27.486402641710413,   2.241362565557243, 
                                   0.504020528307776,  -0.112475532563608,   0.009894993647857, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000002,  -0.000000000000000,   0.000000000000000, 
                                   0.004784479194438,   0.000290215608672,  -0.000129138861109, 
                                   0.006180665466861,   0.000710921112566,   0.000081078326823, 
                                  -0.000000000000000,   0.000000000000001,   0.000000000000000, 
                                   0.000037941204147,   0.000009211114295,   0.000003428627161, 
                                   0.000217983610989,   0.006168816971486,   0.031549533914077, 
                                   0.024336455607841,   0.014341183236479,   1.049467095692840, 
                                   1.407733150725250,   0.001253309998293,   0.006410776739697, 
                                   0.000495328366326,   0.001517434749548,   0.001988673917288, 
                                   0.000791319919334,   0.000719210076361,  -0.001138011875843, 
                                  27.486402641710413,  49.453801229074941,  14.183699388630654, 
                                   3.968568016745553,  -1.468925416580938,   0.173318589943629, 
                                  -0.000000000000001,   0.000000000000000,  -0.000000000000002, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.021881105541992,   0.018646727894000,   0.006679933880673, 
                                   0.024861096119497,   0.028085917116554,   0.010826972476841, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000120867204383,   0.000103065148431,   0.000086802502056, 
                                  -0.000028126753648,   0.000809943562283,   0.016293503649725, 
                                   0.093284155524529,   0.070445712447604,   4.951046387476951, 
                                   4.969353402478769,   0.016194711579447,   0.031599208092944, 
                                  -0.000027952211524,   0.001988673745953,  -0.013203312292400, 
                                  -0.014117030728298,   0.009008644164129,  -0.012890119827290, 
                                   2.241362565557243,  14.183699388630654,  22.571344827658255, 
                                  12.023038715484967,  -6.497783347643201,   1.546998888320734, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   0.009910064222828,   0.049662827726772,   0.059880275410369, 
                                   0.025079460498785,   0.048171621264256,   0.010052250701834, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000597079896750,  -0.000581946335793,  -0.002520881803036, 
                                  -0.000044965243435,   0.000338494442472,   0.004647006975591, 
                                   0.060013415382703,   0.071260274651458,   6.489587471862320, 
                                   7.086329590701152,   0.017777670836737,   0.037503949583236, 
                                  -0.000051989418100,   0.000791320001368,  -0.014117042311481, 
                                  -0.012955670024300,   0.016555707808128,  -0.009835410969445, 
                                   0.504020528307776,   3.968568016745553,  12.023038715484967, 
                                  13.155513667730688,  -3.289732603339671,   2.202737441067176, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000001,   0.000000000000001, 
                                   0.003256418300098,   0.023710583382813,   0.061694098934713, 
                                   0.014572179179418,   0.003825026629936,  -0.001401863181513, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000551572807294,  -0.000857768112767,  -0.003119292180612, 
                                  -0.042181454636970,  -0.228542649784721,  -1.119730324996342, 
                                  -4.213124950236165,  -3.229902963778099,  -2.411542820879911, 
                                  -2.071428213083649,  -0.274763234814444,  -1.460364194555921, 
                                  -0.000095817604301,  -0.001138011804650,  -0.012890121010359, 
                                  -0.009835416240766,  -0.008618413704740,   0.217730366498301, 
                                  -0.112475532563608,  -1.468925416580938,  -6.497783347643201, 
                                  -3.289732603339671,  13.645886510153678,  -0.159831232577100, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001014709891907,  -0.008285927478409,  -0.013892281298798, 
                                  -0.648420402384140,  -4.069305812839722,  -5.409593038449139, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000068117147048,   0.001541502696662,   0.001979210436194, 
                                  -0.000009489543004,   0.000024184211660,   0.000210119210717, 
                                   0.011879145582389,   0.014493051852879,   0.210945451244569, 
                                   0.908140072927748,   0.028427584819176,   0.024001231433518, 
                                  -0.000048611052021,   0.000719210037112,   0.009008658575075, 
                                   0.016555725864068,   0.172520373971464,  -0.008618419585509, 
                                   0.009894993647857,   0.173318589943629,   1.546998888320734, 
                                   2.202737441067176,  -0.159831232577100,  12.419032424006367, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                   0.000000000000002,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.124757738386229,   1.205683683845024,   4.105360204141695, 
                                   0.001649277355717,  -0.002426692805988,  -0.012806963687668, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000084875133086,  -0.002639091649222,  -0.014843153328253, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000001,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  88.570636594293262,  27.488096834465640,   2.240994069699742, 
                                   0.503672360570630,   0.009823025800347,   0.009823025787744, 
                                   0.004813133819611,   0.000315443166109,  -0.000110612144291, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.004813136742540,   0.000315443460869,  -0.000110612016108, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  27.488096834465640,  49.452186794645854,  14.185297337875763, 
                                   3.969650058022597,   0.173625708022963,   0.173625708180023, 
                                   0.021738255216554,   0.018901307197351,   0.006825916264432, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.021738267231155,   0.018901312621533,   0.006825916618409, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000003,   0.000000000000001, 
                                   0.000000000000000,   0.000000000000001,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000002,  -0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   2.240994069699742,  14.185297337875763,  22.481504719083535, 
                                  11.957976858731014,   1.527285433133446,   1.527285421051051, 
                                   0.010622706313969,   0.048562027988382,   0.056846919921712, 
                                  -0.000000000000001,   0.000000000000002,   0.000000000000004, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.010622703701966,   0.048562049080608,   0.056846953338863, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000002,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000002, 
                                   0.503672360570630,   3.969650058022597,  11.957976858731014, 
                                  13.097683403715873,   2.184087252519409,   2.184087237826530, 
                                   0.003775782898387,   0.022357369885220,   0.058559648636651, 
                                  -0.000000000000001,   0.000000000000002,   0.000000000000004, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.003775778506984,   0.022357385909400,   0.058559683001073, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.009823025800347,   0.173625708022963,   1.527285433133446, 
                                   2.184087252519409,  12.409642162563140,   0.020968418657434, 
                                   0.000279085602975,   0.001186888496035,   0.003422727159813, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.124780349145123,   1.202175660491695,   4.097692967674701, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.009823025787744,   0.173625708180023,   1.527285421051051, 
                                   2.184087237826530,   0.020968418657434,  12.409642122944087, 
                                   0.124780350011933,   1.202175655931682,   4.097692954089526, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000279084846218,   0.001186892583928,   0.003422734979779, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000002,  -0.000000000000001,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.004813133819611,   0.021738255216553,   0.010622706313969, 
                                   0.003775782898387,   0.000279085602975,   0.124780350011933, 
                                  53.043866154994106,  16.848265588064670,   3.279900339772641, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.001221765011844,   0.001163243320619,   0.000288396835838, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                   0.000315443166109,   0.018901307197351,   0.048562027988382, 
                                   0.022357369885220,   0.001186888496035,   1.202175655931682, 
                                  16.848265588064670,  28.162873729098990,  12.681472415182409, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.001163243454352,   0.002488135672623,   0.000336228529170, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000001, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                  -0.000110612144291,   0.006825916264432,   0.056846919921712, 
                                   0.058559648636651,   0.003422727159813,   4.097692954089526, 
                                   3.279900339772641,  12.681472415182409,  16.762919645600487, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000288396972560,   0.000336228965980,   0.000926729460314, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000061471353724,   0.000609586843736,   0.002444639735731, 
                                   0.001564157786959,   0.000769701874973,   0.001095892727442, 
                                   0.001584984901569,   0.001095888831268,   0.001584982114382, 
                                   0.004784482472229,   0.021881117305557,   0.009910061948506, 
                                   0.003256414288573,   0.124757737865327,  -0.001014707749876, 
                                   0.004784479194438,   0.021881105541992,   0.009910064222828, 
                                   0.003256418300098,  -0.001014709891907,   0.124757738386229, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000001, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  53.043413445429721,  16.847375938825863,   3.279720453439758, 
                                   0.000656711597808,   0.000686804378049,   0.000233613176149, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.001137460399120,  -0.001189578793481,  -0.000404626931465, 
                                  -0.000017054701545,   0.000176563881884,   0.002268996701883, 
                                   0.002334207404951,   0.000277029891903,   0.006392394292714, 
                                   0.008799788644018,   0.006392412933216,   0.008799801492873, 
                                   0.000290215691857,   0.018646731840302,   0.049662846692087, 
                                   0.023710601474494,   1.205683686360233,  -0.008285933389673, 
                                   0.000290215608672,   0.018646727894000,   0.049662827726772, 
                                   0.023710583382813,  -0.008285927478409,   1.205683683845024, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000002, 
                                   0.000000000000002,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  16.847375938825859,  28.173165809543644,  12.695560224273802, 
                                   0.000691454549861,   0.001437198889057,   0.000202557761922, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.001197632518145,  -0.002489303174187,  -0.000350852108812, 
                                  -0.000023562087915,   0.000126609750030,   0.000432752350941, 
                                   0.000710644790127,   0.000224777045591,   0.015940873488793, 
                                   0.026994160232387,   0.015940911897444,   0.026994187973493, 
                                  -0.000129138994621,   0.006679933703819,   0.059880307234960, 
                                   0.061694137210839,   4.105360213977148,  -0.013892290436904, 
                                  -0.000129138861109,   0.006679933880673,   0.059880275410369, 
                                   0.061694098934713,  -0.013892281298798,   4.105360204141695, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000004, 
                                   0.000000000000004,   0.000000000000001,   0.000000000000001, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   3.279720453439758,  12.695560224273802,  16.796462514357145, 
                                   0.000221879644186,   0.000163274172431,   0.000427098405110, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000384303924943,  -0.000282804371388,  -0.000739779253320, 
                                  -0.000086370708979,  -0.000413474804194,  -0.000506775192761, 
                                   0.000668016737651,   0.001024630991710,   0.339359226095313, 
                                   0.146399499111725,  -0.170094273652441,  -0.072803487374763, 
                                  -0.003123188017088,  -0.012325870260963,  -0.013056803863772, 
                                  -0.007763757262002,  -0.000751132885379,   0.324269190108425, 
                                   0.006180665466861,   0.024861096119497,   0.025079460498785, 
                                   0.014572179179418,  -0.648420402384139,   0.001649277355717, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                   0.000656711597808,   0.000691454549861,   0.000221879644186, 
                                  53.042595732647001,  16.849978923186139,   3.283545628467981, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000741060316429,   0.000372774864051,   0.001443598188326, 
                                   0.000020215809996,  -0.000172795677023,   0.000219162289005, 
                                   0.006146421749063,   0.004368422822374,   2.055207811276272, 
                                   0.783040291434189,  -1.025184276423810,  -0.387444885118817, 
                                  -0.000363437337951,  -0.014132209170082,  -0.023581836836335, 
                                  -0.001169677615586,   0.003498861594793,   2.033317926464628, 
                                   0.000710921112566,   0.028085917116554,   0.048171621264256, 
                                   0.003825026629936,  -4.069305812839723,  -0.002426692805989, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000686804378049,   0.001437198889057,   0.000163274172431, 
                                  16.849978923186139,  28.169621007839503,  12.674358532440564, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000372774639199,   0.005760660965885,   0.001843374508220, 
                                   0.000025296519345,  -0.000119281271301,   0.000417938472961, 
                                   0.004491677471710,   0.003695911256960,   5.992674202150579, 
                                   2.391307763904538,  -2.990909067453096,  -1.189966682029321, 
                                  -0.000043508421716,  -0.005488656716481,  -0.002843018329669, 
                                   0.003402275017597,   0.019258015062522,   2.703082478375381, 
                                   0.000081078326823,   0.010826972476841,   0.010052250701834, 
                                  -0.001401863181513,  -5.409593038449140,  -0.012806963687668, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000233613176149,   0.000202557761922,   0.000427098405110, 
                                   3.283545628467981,  12.674358532440564,  16.733952547051310, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.001443597879339,   0.001843374195670,  -0.000298626118140, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.004813136742540,   0.021738267231155,   0.010622703701966, 
                                   0.003775778506984,   0.124780349145123,   0.000279084846218, 
                                   0.001221765011844,   0.001163243454352,   0.000288396972560, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  53.043866153475314,  16.848265582998117,   3.279900337902300, 
                                  -0.000000000000001,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000001,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   0.000315443460869,   0.018901312621533,   0.048562049080608, 
                                   0.022357385909400,   1.202175660491695,   0.001186892583928, 
                                   0.001163243320619,   0.002488135672623,   0.000336228965980, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  16.848265582998117,  28.162873719441901,  12.681472414232271, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000000000000000,   0.000000000000000,   0.000000000000001, 
                                   0.000000000000001,   0.000000000000000,  -0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                  -0.000110612016108,   0.006825916618409,   0.056846953338863, 
                                   0.058559683001073,   4.097692967674701,   0.003422734979779, 
                                   0.000288396835838,   0.000336228529170,   0.000926729460314, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,   0.000000000000000,   0.000000000000000, 
                                   3.279900337902300,  12.681472414232271,  16.762919646034860, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                   0.000149598299673,   0.000716156331708,   0.000877753343310, 
                                  -0.001157049994425,  -0.001774720064802,   0.000478812130912, 
                                  -0.000457563103276,  -0.293654303033033,  -0.127014466868014, 
                                  -0.005333639289730,  -0.021590764066772,  -0.021420903118261, 
                                  -0.012344088175682,  -0.001470753029586,   0.561514481324687, 
                                   0.000037941204147,  -0.000120867204383,   0.000597079896750, 
                                   0.000551572807294,  -0.000068117147047,  -0.000084875133086, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001137460399120,  -0.001197632518145,  -0.000384303924943, 
                                  -0.000741060316426,   0.000372774639199,   0.001443597879338, 
                                  -0.000000000000001,  -0.000000000000000,  -0.000000000000000, 
                                  53.043451444792510,  16.849548484368132,   3.281878703713686, 
                                  -0.000035014772696,   0.000299290276865,  -0.000379608509689, 
                                  -0.010645912861443,  -0.007566322787821,  -0.002793953214152, 
                                  -0.004705709143368,  -1.781259152766736,  -0.680485639467545, 
                                  -0.000611069646923,  -0.024271575806446,  -0.042008814364454, 
                                  -0.003741457558729,   0.000782032192792,   3.524892957702942, 
                                   0.000009211114295,   0.000103065148431,  -0.000581946335793, 
                                  -0.000857768112767,   0.001541502696662,  -0.002639091649222, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.001189578793481,  -0.002489303174187,  -0.000282804371388, 
                                   0.000372774864050,   0.005760660965883,   0.001843374195669, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  16.849548484368132,  28.162969178956097,  12.672229989186372, 
                                  -0.000043814821989,   0.000206600896872,  -0.000723892735288, 
                                  -0.007779790243669,  -0.006401483381472,  -0.006267772523248, 
                                  -0.006567026245229,  -5.192941986059282,  -2.074216786397519, 
                                  -0.000068501257130,  -0.009333029527092,  -0.009965949366088, 
                                  -0.000345604330298,   0.003669577771402,   4.685834594993277, 
                                   0.000003428627161,   0.000086802502056,  -0.002520881803036, 
                                  -0.003119292180612,   0.001979210436194,  -0.014843153328253, 
                                  -0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000000000000000,  -0.000000000000000,   0.000000000000000, 
                                   0.000000000000000,  -0.000000000000000,  -0.000000000000000, 
                                  -0.000404626931465,  -0.000350852108812,  -0.000739779253320, 
                                   0.001443598188326,   0.001843374508219,  -0.000298626118140, 
                                   0.000000000000000,   0.000000000000000,  -0.000000000000000, 
                                   3.281878703713686,  12.672229989186372,  16.734297379699314};

    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals.size()));
    
    // Maximum element is 117.9 and we compare up to its 14-th significant digit.

    vlxtest::compare(jkvals, fock.getFock(0), 1.0e-11);

    // Small deviations:
    //Row 0 Col 0 Diff= 5.40012e-13
    //Row 0 Col 5 Diff= 8.86653e-12
    //Row 0 Col 7 Diff= 8.86623e-12
    //Row 1 Col 5 Diff= 2.13718e-13
    //Row 1 Col 7 Diff= 2.1394e-13
    //Row 4 Col 6 Diff= 1.36779e-13
    //Row 4 Col 8 Diff= 1.40332e-13
    //Row 5 Col 0 Diff= 8.86653e-12
    //Row 5 Col 1 Diff= 2.13718e-13
    //Row 5 Col 7 Diff= 3.07365e-13
    //Row 6 Col 4 Diff= 1.36779e-13
    //Row 7 Col 0 Diff= 8.86623e-12
    //Row 7 Col 1 Diff= 2.1394e-13
    //Row 7 Col 5 Diff= 3.07365e-13
    //Row 8 Col 4 Diff= 1.40332e-13
    //Row 8 Col 8 Diff= 1.1724e-13
    //Row 9 Col 9 Diff= 2.70006e-13
    //Row 10 Col 10 Diff= 1.27898e-13
    //Row 15 Col 15 Diff= 3.12639e-13
    //Row 16 Col 16 Diff= 1.06581e-13
    //Row 21 Col 21 Diff= 3.55271e-13
    //Row 27 Col 27 Diff= 1.63425e-13
    //Row 30 Col 30 Diff= 1.7053e-13
    //Row 31 Col 31 Diff= 1.03029e-13
    //Row 33 Col 33 Diff= 1.13687e-13
    //Row 36 Col 36 Diff= 1.7053e-13
    //Row 37 Col 37 Diff= 1.13687e-13
    //Row 39 Col 39 Diff= 1.42109e-13
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestJForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getRestDensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    fock.setFockType(fockmat::restj, 0);
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> jvals {
     8.7050379187043561,  -2.4657043786446717,  -1.1137995569755967,  -0.3511318732750601,  -0.4446320177607397,
    -0.3511318732750600,  -0.4446320177607397,  -0.0000000000000001,  -0.0000000000000000,   0.5361288159064489,
    -0.5361288159064490,  -0.0044770319812429,  -0.0006299498865047,   0.4210136824140648,   0.4210136824140649,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000235185375759,  -0.0000000000000000,  -0.0006317274157580,  -2.4657043786446717,
     4.6107988401885738,   3.2399595033624133,   1.3380164196748638,   1.5519940396258918,   1.3380164196748638,
     1.5519940396258911,  -0.0000000000000000,  -0.0000000000000001,  -1.5278603190511189,   1.5278603190511197,
     0.0426591930930834,   0.0323882948027741,  -1.1900843008517898,  -1.1900843008517907,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0022396741511996,   0.0000000000000000,  -0.0067024628386738,  -1.1137995569755967,   3.2399595033624133,
     3.2510234849410455,   1.6016184595970253,   1.9832829274635904,   1.6016184595970262,   1.9832829274635904,
    -0.0000000000000000,  -0.0000000000000001,  -1.1203065357796762,   1.1203065357796764,   0.0394685894113138,
     0.0593780468191573,  -0.8628237714505190,  -0.8628237714505193,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0027616119843075,
     0.0000000000000000,  -0.0103712485702743,  -0.3511318732750601,   1.3380164196748638,   1.6016184595970253,
     2.6736148307623582,   1.7151234216931157,   0.4305685278263138,   0.8916189352925326,   0.8358379607934615,
     1.2946077512285303,  -0.5092833429303614,   0.5432278109275681,   0.6688777249219887,   1.0343712373243215,
    -0.3771699317197293,  -0.0674918927665985,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.6556417187043254,   0.0649272942106411,   0.0000000000000000,
    -0.4126231745655132,  -0.4446320177607397,   1.5519940396258918,   1.9832829274635904,   1.7151234216931157,
     2.0730952363039421,   0.8916189352925326,   1.4369218507590551,   0.2959165689551334,   0.6810628720504702,
    -0.4022917927272958,   0.6269996274430015,   0.2521880185372851,   0.5753797135624102,  -0.2969946963794896,
    -0.1674320398271054,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0686013019904289,   0.0111720312088013,   0.0000000000000000,  -0.0457391853379266,
    -0.3511318732750600,   1.3380164196748638,   1.6016184595970262,   0.4305685278263138,   0.8916189352925326,
     2.6736148307623586,   1.7151234216931157,  -0.8358379607934618,  -1.2946077512285310,  -0.5432278109275680,
     0.5092833429303614,   0.6688777249219888,   1.0343712373243217,  -0.0674918927665984,  -0.3771699317197293,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.6556417187043254,   0.0649272942106411,  -0.0000000000000000,  -0.4126231745655135,  -0.4446320177607397,
     1.5519940396258911,   1.9832829274635904,   0.8916189352925326,   1.4369218507590551,   1.7151234216931157,
     2.0730952363039417,  -0.2959165689551335,  -0.6810628720504704,  -0.6269996274430014,   0.4022917927272958,
     0.2521880185372851,   0.5753797135624102,  -0.1674320398271053,  -0.2969946963794897,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0686013019904289,
     0.0111720312088014,   0.0000000000000000,  -0.0457391853379268,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,   0.8358379607934615,   0.2959165689551334,  -0.8358379607934618,  -0.2959165689551335,
     4.7076345675032147,   1.8801553141599465,  -0.7297834861388904,  -0.7297834861388908,   0.0000000000000000,
     0.0000000000000000,  -1.1263972132955591,   1.1263972132955593,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0487109928656118,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,
     1.2946077512285303,   0.6810628720504702,  -1.2946077512285310,  -0.6810628720504704,   1.8801553141599465,
     2.5960648818262251,  -0.0412100329723808,  -0.0412100329723807,   0.0000000000000000,  -0.0000000000000000,
    -0.8070479620307962,   0.8070479620307962,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0452456368688820,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.5361288159064489,  -1.5278603190511189,  -1.1203065357796762,  -0.5092833429303614,
    -0.4022917927272958,  -0.5432278109275680,  -0.6269996274430014,  -0.7297834861388904,  -0.0412100329723808,
     2.7404754067902717,  -0.8331893089142363,  -1.1460059830560330,  -0.8324963086962188,   0.1813765630566836,
     0.0848158427995801,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.6653025800873301,  -0.6151565796374144,   0.0000000000000000,  -0.0012575959103357,
    -0.5361288159064490,   1.5278603190511197,   1.1203065357796764,   0.5432278109275681,   0.6269996274430015,
     0.5092833429303614,   0.4022917927272958,  -0.7297834861388908,  -0.0412100329723807,  -0.8331893089142363,
     2.7404754067902721,   1.1460059830560334,   0.8324963086962188,  -0.0848158427995801,  -0.1813765630566836,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.6653025800873297,   0.6151565796374148,  -0.0000000000000000,   0.0012575959103357,  -0.0044770319812429,
     0.0426591930930834,   0.0394685894113138,   0.6688777249219887,   0.2521880185372851,   0.6688777249219888,
     0.2521880185372851,   0.0000000000000000,   0.0000000000000000,  -1.1460059830560330,   1.1460059830560334,
     4.7086781283286561,   1.8777970967376274,  -0.1827852235097767,  -0.1827852235097768,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
     0.0362978447333053,  -0.0000000000000000,  -0.0090481046098741,  -0.0006299498865047,   0.0323882948027741,
     0.0593780468191573,   1.0343712373243215,   0.5753797135624102,   1.0343712373243217,   0.5753797135624102,
     0.0000000000000000,  -0.0000000000000000,  -0.8324963086962188,   0.8324963086962188,   1.8777970967376274,
     2.5873491028090645,   0.3624137967490834,   0.3624137967490834,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0329115572610546,
    -0.0000000000000000,  -0.0088431046618546,   0.4210136824140648,  -1.1900843008517898,  -0.8628237714505190,
    -0.3771699317197293,  -0.2969946963794896,  -0.0674918927665984,  -0.1674320398271053,  -1.1263972132955591,
    -0.8070479620307962,   0.1813765630566836,  -0.0848158427995801,  -0.1827852235097767,   0.3624137967490834,
     2.6383614315022412,   0.1565872541926126,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.2015306815839267,   0.6742742998009658,  -0.0000000000000000,
     0.6568305949296774,   0.4210136824140649,  -1.1900843008517907,  -0.8628237714505193,  -0.0674918927665985,
    -0.1674320398271054,  -0.3771699317197293,  -0.2969946963794897,   1.1263972132955593,   0.8070479620307962,
     0.0848158427995801,  -0.1813765630566836,  -0.1827852235097768,   0.3624137967490834,   0.1565872541926126,
     2.6383614315022421,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.2015306815839266,   0.6742742998009656,   0.0000000000000000,   0.6568305949296774,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     4.7060131273126666,   1.8712403651288989,   0.7119430700412535,   0.7119430700412537,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0306147836458636,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   1.8712403651288989,
     2.5680728042510141,   1.0030967709139023,   1.0030967709139023,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0275594275451729,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.7119430700412535,   1.0030967709139023,
     2.5002601160998239,   0.1361000572995744,   0.8432718608386539,   0.0000000000000000,   0.0000000000000000,
     0.6696522336226355,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.7119430700412537,   1.0030967709139023,   0.1361000572995744,
     2.5002601160998239,  -0.8432718608386539,   0.0000000000000000,   0.0000000000000000,   0.6696522336226353,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.8432718608386539,  -0.8432718608386539,
     3.7888295812519379,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.6556417187043254,   0.0686013019904289,
    -0.6556417187043254,  -0.0686013019904289,   0.0487109928656118,   0.0452456368688820,  -0.6653025800873301,
    -0.6653025800873297,  -0.0000000000000001,  -0.0000000000000001,  -0.2015306815839267,   0.2015306815839266,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     3.8055595971784735,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000235185375759,
     0.0022396741511996,   0.0027616119843075,   0.0649272942106411,   0.0111720312088013,   0.0649272942106411,
     0.0111720312088014,   0.0000000000000000,   0.0000000000000000,  -0.6151565796374144,   0.6151565796374148,
     0.0362978447333053,   0.0329115572610546,   0.6742742998009658,   0.6742742998009656,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     3.7920055763533203,   0.0000000000000000,   0.0041449416324419,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0306147836458636,   0.0275594275451729,
     0.6696522336226355,   0.6696522336226353,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     3.7869595774205220,   0.0000000000000000,  -0.0006317274157580,  -0.0067024628386738,  -0.0103712485702743,
    -0.4126231745655132,  -0.0457391853379266,  -0.4126231745655135,  -0.0457391853379268,   0.0000000000000000,
     0.0000000000000001,  -0.0012575959103357,   0.0012575959103357,  -0.0090481046098741,  -0.0088431046618546,
     0.6568305949296774,   0.6568305949296774,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0041449416324422,   0.0000000000000000,
     3.7925725062365911};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jvals.size()));
    
    vlxtest::compare(jvals, fock.getFock(0));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestKForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getRestDensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    fock.setFockType(fockmat::restk, 0);
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> kvals {
     4.8932786799638182,  -1.2655771692675917,  -0.5828552256161563,  -0.1882276316348811,  -0.2364867532449839,
    -0.1882276316348804,  -0.2364867532449837,  -0.0000000000000006,   0.0000000000000000,   0.2788797238898280,
    -0.2788797238898276,  -0.0048603513463662,  -0.0025145713869552,   0.2177469906250327,   0.2177469906250324,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,
     0.0000000000000003,   0.0002878824028762,   0.0000000000000000,  -0.0058174580656580,  -1.2655771692675921,
     1.4382079213212147,   1.1201483331788544,   0.5168930961801773,   0.5944250832856743,   0.5168930961801773,
     0.5944250832856740,  -0.0000000000000000,  -0.0000000000000002,  -0.4470782288046084,   0.4470782288046086,
    -0.0386213901867632,   0.0090953813732034,  -0.3617634587654666,  -0.3617634587654666,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0021571614140512,   0.0000000000000001,  -0.0154062102704700,  -0.5828552256161562,   1.1201483331788546,
     0.9901361804456427,   0.5094370349947297,   0.5679397411030980,   0.5094370349947306,   0.5679397411030981,
    -0.0000000000000003,  -0.0000000000000003,  -0.3726431075848204,   0.3726431075848207,   0.0032998469349670,
     0.0513460955255465,  -0.2916015927375100,  -0.2916015927375104,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0051867963861795,
     0.0000000000000000,  -0.0277220672116082,  -0.1882276316348810,   0.5168930961801772,   0.5094370349947298,
     0.6487537979144666,   0.4494560771483926,   0.1595979333968801,   0.2671262070651429,   0.3432258941530916,
     0.3371746972914381,  -0.2899818965868014,   0.2096670051090453,   0.2689544584802174,   0.2715616951091601,
    -0.2184770222651541,  -0.0247919537279099,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.1464578162890123,   0.0182323098347128,   0.0000000000000001,
    -0.0862535151319074,  -0.2364867532449839,   0.5944250832856744,   0.5679397411030980,   0.4494560771483926,
     0.4024941141237751,   0.2671262070651431,   0.3277195452455293,   0.1352198169875467,   0.1367566150152336,
    -0.2446143867880092,   0.2237121935160887,   0.1220011203993555,   0.1462971255534110,  -0.1827768271228056,
    -0.1034016195622274,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,
     0.0000000000000000,   0.0387937747977301,   0.0089453234237394,   0.0000000000000000,  -0.0344576112860573,
    -0.1882276316348803,   0.5168930961801772,   0.5094370349947305,   0.1595979333968801,   0.2671262070651432,
     0.6487537979144675,   0.4494560771483924,  -0.3432258941530914,  -0.3371746972914385,  -0.2096670051090454,
     0.2899818965868015,   0.2689544584802178,   0.2715616951091605,  -0.0247919537279099,  -0.2184770222651542,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
    -0.1464578162890126,   0.0182323098347133,  -0.0000000000000000,  -0.0862535151319067,  -0.2364867532449838,
     0.5944250832856740,   0.5679397411030980,   0.2671262070651430,   0.3277195452455293,   0.4494560771483924,
     0.4024941141237752,  -0.1352198169875470,  -0.1367566150152340,  -0.2237121935160885,   0.2446143867880093,
     0.1220011203993552,   0.1462971255534109,  -0.1034016195622273,  -0.1827768271228059,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,  -0.0387937747977301,
     0.0089453234237393,  -0.0000000000000000,  -0.0344576112860575,  -0.0000000000000006,   0.0000000000000000,
    -0.0000000000000002,   0.3432258941530909,   0.1352198169875466,  -0.3432258941530909,  -0.1352198169875469,
     0.9278490388284155,   0.6016580247177705,  -0.1231907922436001,  -0.1231907922435999,   0.0000000000000004,
     0.0000000000000002,  -0.2546583911774608,   0.2546583911774610,   0.0000000000000005,   0.0000000000000002,
     0.0000000000000001,   0.0000000000000001,  -0.0000000000000002,   0.1136428386862909,  -0.0000000000000000,
     0.0000000000000002,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000002,  -0.0000000000000003,
     0.3371746972914381,   0.1367566150152336,  -0.3371746972914387,  -0.1367566150152340,   0.6016580247177704,
     0.5453198615202637,  -0.0635649245192273,  -0.0635649245192274,   0.0000000000000002,   0.0000000000000001,
    -0.1894877563946145,   0.1894877563946148,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,   0.1195719917801947,  -0.0000000000000000,   0.0000000000000001,
     0.0000000000000000,   0.2788797238898277,  -0.4470782288046084,  -0.3726431075848207,  -0.2899818965868013,
    -0.2446143867880091,  -0.2096670051090453,  -0.2237121935160883,  -0.1231907922436001,  -0.0635649245192273,
     0.3184217506690389,  -0.2179942782434513,  -0.2652735463121769,  -0.2056554948939779,   0.1132260288568987,
     0.0579605702643994,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0616024387723321,  -0.0634664338853550,  -0.0000000000000000,   0.0098405767904064,
    -0.2788797238898274,   0.4470782288046084,   0.3726431075848207,   0.2096670051090453,   0.2237121935160885,
     0.2899818965868015,   0.2446143867880092,  -0.1231907922435999,  -0.0635649245192274,  -0.2179942782434513,
     0.3184217506690394,   0.2652735463121770,   0.2056554948939779,  -0.0579605702643994,  -0.1132260288568988,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000001,
    -0.0616024387723321,   0.0634664338853548,  -0.0000000000000000,  -0.0098405767904064,  -0.0048603513463661,
    -0.0386213901867630,   0.0032998469349670,   0.2689544584802173,   0.1220011203993554,   0.2689544584802175,
     0.1220011203993552,   0.0000000000000004,   0.0000000000000002,  -0.2652735463121769,   0.2652735463121770,
     1.0031680253743585,   0.6266678369899475,   0.0244887923092399,   0.0244887923092403,   0.0000000000000001,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
     0.0155751297569286,  -0.0000000000000000,  -0.0518849741632175,  -0.0025145713869552,   0.0090953813732034,
     0.0513460955255466,   0.2715616951091600,   0.1462971255534111,   0.2715616951091604,   0.1462971255534109,
     0.0000000000000002,   0.0000000000000001,  -0.2056554948939779,   0.2056554948939779,   0.6266678369899477,
     0.5378088023681107,   0.0417130407867966,   0.0417130407867968,   0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000002,   0.0209850376195948,
    -0.0000000000000000,  -0.0517935790949717,   0.2177469906250324,  -0.3617634587654666,  -0.2916015927375104,
    -0.2184770222651541,  -0.1827768271228056,  -0.0247919537279098,  -0.1034016195622273,  -0.2546583911774608,
    -0.1894877563946145,   0.1132260288568987,  -0.0579605702643994,   0.0244887923092399,   0.0417130407867966,
     0.2627030477375997,   0.0434662957373307,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0364353322385420,   0.0673934895667385,  -0.0000000000000001,
     0.0565314158821576,   0.2177469906250322,  -0.3617634587654667,  -0.2916015927375105,  -0.0247919537279099,
    -0.1034016195622274,  -0.2184770222651542,  -0.1827768271228059,   0.2546583911774610,   0.1894877563946147,
     0.0579605702643994,  -0.1132260288568988,   0.0244887923092403,   0.0417130407867968,   0.0434662957373307,
     0.2627030477375998,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,   0.0364353322385419,   0.0673934895667385,   0.0000000000000000,   0.0565314158821577,
     0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000005,   0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,   0.0000000000000001,  -0.0000000000000001,   0.0000000000000001,
     1.0909101434680140,   0.6408738966054136,   0.2578268332060425,   0.2578268332060426,   0.0000000000000001,
    -0.0000000000000000,   0.0000000000000000,   0.0110870326028947,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000002,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000001,   0.6408738966054134,
     0.5022428184560123,   0.2108207505462563,   0.2108207505462565,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0136589834340676,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.2578268332060425,   0.2108207505462563,
     0.1931352667625228,   0.0682538201035560,   0.0854311240596404,  -0.0000000000000000,   0.0000000000000000,
     0.0721371939225131,  -0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.2578268332060427,   0.2108207505462565,   0.0682538201035560,
     0.1931352667625230,  -0.0854311240596404,  -0.0000000000000000,   0.0000000000000000,   0.0721371939225131,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000001,  -0.0000000000000002,  -0.0000000000000001,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000001,  -0.0000000000000000,   0.0854311240596403,  -0.0854311240596404,
     0.2948008402817411,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000003,   0.0000000000000001,   0.0000000000000001,   0.1464578162890123,   0.0387937747977301,
    -0.1464578162890127,  -0.0387937747977301,   0.1136428386862909,   0.1195719917801946,  -0.0616024387723321,
    -0.0616024387723321,  -0.0000000000000001,  -0.0000000000000002,  -0.0364353322385420,   0.0364353322385419,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.3351079390718259,  -0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,   0.0002878824028760,
     0.0021571614140511,   0.0051867963861799,   0.0182323098347127,   0.0089453234237396,   0.0182323098347130,
     0.0089453234237393,  -0.0000000000000000,  -0.0000000000000000,  -0.0634664338853550,   0.0634664338853548,
     0.0155751297569286,   0.0209850376195948,   0.0673934895667385,   0.0673934895667383,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
     0.2958903447049419,   0.0000000000000000,  -0.0107090705471177,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000002,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0110870326028947,   0.0136589834340676,
     0.0721371939225131,   0.0721371939225131,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
     0.3013645688692004,   0.0000000000000000,  -0.0058174580656584,  -0.0154062102704700,  -0.0277220672116081,
    -0.0862535151319075,  -0.0344576112860571,  -0.0862535151319073,  -0.0344576112860574,  -0.0000000000000001,
     0.0000000000000000,   0.0098405767904064,  -0.0098405767904064,  -0.0518849741632175,  -0.0517935790949716,
     0.0565314158821576,   0.0565314158821576,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0107090705471177,   0.0000000000000000,
     0.3093499242861350};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(kvals.size()));
    
    vlxtest::compare(kvals, fock.getFock(0));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestGenJForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getGenOODensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    fock.setFockType(fockmat::rgenj, 0);

    fock.setFockType(fockmat::rgenj, 1);
    
    ASSERT_TRUE(fockmat::rgenj == fock.getFockType(0));

    ASSERT_TRUE(fockmat::rgenj == fock.getFockType(1));
    
    ASSERT_TRUE(fockmat::rgenj == fock.getFockType(0));
    
    ASSERT_TRUE(fockmat::rgenj == fock.getFockType(1));
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    ASSERT_EQ(fock.getNumberOfRows(1), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(1), ncols);
    
    std::vector<double> jvals00{
     0.9903111380024709,  -0.3311052090308504,  -0.1513736205733688,  -0.0472898591165016,  -0.0604013564553225,
    -0.0472898591165016,  -0.0604013564553225,   0.0000000000000000,   0.0000000000000000,   0.0714903991346079,
    -0.0714903991346079,   0.0199519230731715,   0.0025181179725149,   0.0566414756881583,   0.0566414756881583,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0036484587011672,   0.0000000000000000,   0.0000197107354208,  -0.3311052090308504,
     0.8039658202674774,   0.5840669618046817,   0.2321585004461648,   0.2783392199672013,   0.2321585004461650,
     0.2783392199672012,  -0.0000000000000002,  -0.0000000000000001,  -0.2580718053356358,   0.2580718053356359,
    -0.0745336932827087,  -0.0319906745306454,  -0.2023576943597737,  -0.2023576943597739,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,
     0.0494758194960815,  -0.0000000000000000,  -0.0014776383358171,  -0.1513736205733688,   0.5840669618046817,
     0.6117261883414851,   0.2885975889562702,   0.3716854883680290,   0.2885975889562703,   0.3716854883680291,
    -0.0000000000000001,  -0.0000000000000001,  -0.1943019618051428,   0.1943019618051428,  -0.0516714434949631,
    -0.0360872961116262,  -0.1477523656041203,  -0.1477523656041205,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0433971985646316,
    -0.0000000000000000,  -0.0018137035136946,  -0.0472898591165016,   0.2321585004461648,   0.2885975889562702,
     0.4859786204035437,   0.3148626028014539,   0.0770921955224147,   0.1624727855709966,   0.1457158565949765,
     0.2363689882942751,  -0.0859068922008617,   0.0965635393722776,   0.1073099953234435,   0.1854425537275453,
    -0.0519079977248084,  -0.0088120361147318,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.1099952791007642,   0.0197068229881318,  -0.0000000000000000,
    -0.0730248269109102,  -0.0604013564553225,   0.2783392199672013,   0.3716854883680290,   0.3148626028014539,
     0.3917876099946511,   0.1624727855709965,   0.2703940998275061,   0.0531800387011594,   0.1289501655751335,
    -0.0677351047235564,   0.1110484164078130,   0.0220124651070852,   0.0875042431262234,  -0.0444651617770854,
    -0.0244588618627793,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0076307264258489,   0.0191316520276135,  -0.0000000000000000,  -0.0082843545210679,
    -0.0472898591165016,   0.2321585004461650,   0.2885975889562703,   0.0770921955224147,   0.1624727855709965,
     0.4859786204035438,   0.3148626028014540,  -0.1457158565949766,  -0.2363689882942752,  -0.0965635393722777,
     0.0859068922008617,   0.1073099953234434,   0.1854425537275452,  -0.0088120361147318,  -0.0519079977248085,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.1099952791007642,   0.0197068229881318,  -0.0000000000000000,  -0.0730248269109103,  -0.0604013564553225,
     0.2783392199672012,   0.3716854883680291,   0.1624727855709966,   0.2703940998275061,   0.3148626028014540,
     0.3917876099946513,  -0.0531800387011596,  -0.1289501655751336,  -0.1110484164078131,   0.0677351047235564,
     0.0220124651070852,   0.0875042431262233,  -0.0244588618627793,  -0.0444651617770854,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0076307264258488,
     0.0191316520276135,  -0.0000000000000000,  -0.0082843545210679,   0.0000000000000000,  -0.0000000000000002,
    -0.0000000000000001,   0.1457158565949765,   0.0531800387011594,  -0.1457158565949766,  -0.0531800387011596,
     0.7933856712510483,   0.3387554940427114,  -0.1229869632075393,  -0.1229869632075394,   0.0000000000000001,
     0.0000000000000001,  -0.1942880835962469,   0.1942880835962472,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0476071515016019,   0.0000000000000001,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,
     0.2363689882942751,   0.1289501655751335,  -0.2363689882942752,  -0.1289501655751336,   0.3387554940427114,
     0.4923615292063780,  -0.0038319613847527,  -0.0038319613847527,   0.0000000000000001,   0.0000000000000001,
    -0.1426600259742565,   0.1426600259742567,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0270379041540919,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,   0.0714903991346079,  -0.2580718053356358,  -0.1943019618051428,  -0.0859068922008617,
    -0.0677351047235564,  -0.0965635393722777,  -0.1110484164078131,  -0.1229869632075393,  -0.0038319613847527,
     0.4959282365923101,  -0.1479844046183797,  -0.1904754793127637,  -0.1496120101584177,   0.0214161932231745,
     0.0103780370500320,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.1129010818380865,  -0.1182191309093775,   0.0000000000000000,  -0.0008734935429353,
    -0.0714903991346079,   0.2580718053356359,   0.1943019618051428,   0.0965635393722776,   0.1110484164078130,
     0.0859068922008617,   0.0677351047235564,  -0.1229869632075394,  -0.0038319613847527,  -0.1479844046183797,
     0.4959282365923099,   0.1904754793127639,   0.1496120101584177,  -0.0103780370500320,  -0.0214161932231745,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.1129010818380865,   0.1182191309093776,  -0.0000000000000000,   0.0008734935429353,   0.0199519230731715,
    -0.0745336932827087,  -0.0516714434949631,   0.1073099953234435,   0.0220124651070852,   0.1073099953234434,
     0.0220124651070852,   0.0000000000000001,   0.0000000000000001,  -0.1904754793127637,   0.1904754793127639,
     0.8652295015575346,   0.3719983965318763,  -0.0170848209137793,  -0.0170848209137793,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0545135027927105,   0.0000000000000000,  -0.0017285237294440,   0.0025181179725149,  -0.0319906745306454,
    -0.0360872961116262,   0.1854425537275453,   0.0875042431262234,   0.1854425537275452,   0.0875042431262233,
     0.0000000000000001,   0.0000000000000001,  -0.1496120101584177,   0.1496120101584177,   0.3719983965318763,
     0.5339004594878535,   0.0787930040441641,   0.0787930040441642,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0315536172720997,
     0.0000000000000000,  -0.0017060816959773,   0.0566414756881583,  -0.2023576943597737,  -0.1477523656041203,
    -0.0519079977248084,  -0.0444651617770854,  -0.0088120361147318,  -0.0244588618627793,  -0.1942880835962469,
    -0.1426600259742565,   0.0214161932231745,  -0.0103780370500320,  -0.0170848209137793,   0.0787930040441641,
     0.4817022886806471,   0.0272500177030089,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0222636570994396,   0.1167283449647423,   0.0000000000000000,
     0.1161280577269636,   0.0566414756881583,  -0.2023576943597739,  -0.1477523656041205,  -0.0088120361147318,
    -0.0244588618627793,  -0.0519079977248085,  -0.0444651617770854,   0.1942880835962472,   0.1426600259742567,
     0.0103780370500320,  -0.0214161932231745,  -0.0170848209137793,   0.0787930040441642,   0.0272500177030089,
     0.4817022886806471,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0222636570994395,   0.1167283449647424,   0.0000000000000000,   0.1161280577269637,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.7917198217367811,   0.3371333557843082,   0.1225523345839089,   0.1225523345839090,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0510641989604900,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.3371333557843082,
     0.4879822974974724,   0.1811891771835378,   0.1811891771835378,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,  -0.0304500675460464,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.1225523345839089,   0.1811891771835378,
     0.4550116507506123,   0.0243494013844303,   0.1480805991296236,  -0.0000000000000000,  -0.0000000000000000,
     0.1146355498247648,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.1225523345839090,   0.1811891771835378,   0.0243494013844303,
     0.4550116507506121,  -0.1480805991296238,   0.0000000000000000,  -0.0000000000000000,   0.1146355498247649,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.1480805991296236,  -0.1480805991296238,
     0.6786544431561625,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.1099952791007642,   0.0076307264258489,
    -0.1099952791007642,  -0.0076307264258488,  -0.0476071515016019,  -0.0270379041540919,  -0.1129010818380865,
    -0.1129010818380865,  -0.0000000000000001,  -0.0000000000000001,  -0.0222636570994396,   0.0222636570994395,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.7332004704140767,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0036484587011672,
     0.0494758194960815,   0.0433971985646316,   0.0197068229881318,   0.0191316520276135,   0.0197068229881318,
     0.0191316520276135,   0.0000000000000001,   0.0000000000000000,  -0.1182191309093775,   0.1182191309093776,
    -0.0545135027927105,  -0.0315536172720997,   0.1167283449647423,   0.1167283449647424,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.7474281785169019,  -0.0000000000000000,   0.0003853747547328,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0510641989604900,  -0.0304500675460464,
     0.1146355498247648,   0.1146355498247649,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.7290323148863876,  -0.0000000000000000,   0.0000197107354208,  -0.0014776383358171,  -0.0018137035136946,
    -0.0730248269109102,  -0.0082843545210679,  -0.0730248269109103,  -0.0082843545210679,   0.0000000000000001,
     0.0000000000000001,  -0.0008734935429353,   0.0008734935429353,  -0.0017285237294440,  -0.0017060816959773,
     0.1161280577269636,   0.1161280577269637,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0003853747547327,  -0.0000000000000000,
     0.6790668631458123};

    std::vector<double> jvals01{
     0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0131336386046957,  -0.0016651372155552,  -0.0003444217856378,  -0.0003444217856378,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0035434925070312,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0521029999545313,
     0.0240875907371629,   0.0020144833160349,   0.0020144833160349,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,  -0.0467770759652298,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0373287282629393,   0.0301279689042306,
     0.0007177209071113,   0.0007177209071114,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
    -0.0406463373234668,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0067733153081290,   0.0034385602962241,  -0.0043132752731360,
    -0.0008128311775260,   0.0040750797682010,  -0.0000000000000000,  -0.0000000000000000,  -0.0091989800056414,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0151938118383802,   0.0145628945002607,  -0.0012482753574777,  -0.0008231662774004,
     0.0029937890223691,  -0.0000000000000000,  -0.0000000000000000,  -0.0177145792739417,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0067733153081291,   0.0034385602962241,  -0.0008128311775260,  -0.0043132752731359,  -0.0040750797682011,
    -0.0000000000000000,  -0.0000000000000000,  -0.0091989800056414,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0151938118383802,
     0.0145628945002608,  -0.0008231662774004,  -0.0012482753574776,  -0.0029937890223692,  -0.0000000000000000,
    -0.0000000000000000,  -0.0177145792739417,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
     0.0019344226233444,  -0.0019344226233444,   0.0343668992671324,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0008348282358765,
    -0.0008348282358765,   0.0212444685965512,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0047514582955382,  -0.0003589333227862,   0.0036588784722432,   0.0014640764791231,
    -0.0023508829431383,   0.0000000000000000,   0.0000000000000000,   0.0100745299607064,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0047514582955382,   0.0003589333227862,  -0.0014640764791231,  -0.0036588784722432,  -0.0023508829431384,
    -0.0000000000000000,  -0.0000000000000000,  -0.0100745299607064,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0400800102207108,
    -0.0184343231458730,  -0.0041051314180200,  -0.0041051314180199,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0332450227500307,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0184343231458730,  -0.0229357807827095,
    -0.0055444763226152,  -0.0055444763226151,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0205129237239060,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,  -0.0093583235685706,  -0.0064824318329930,  -0.0057552297930153,
    -0.0004384205182883,  -0.0087750045743895,   0.0000000000000000,   0.0000000000000000,   0.0059203292338797,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0093583235685706,  -0.0064824318329930,  -0.0004384205182883,  -0.0057552297930154,
     0.0087750045743896,   0.0000000000000000,   0.0000000000000000,   0.0059203292338797,  -0.0000000000000000,
    -0.0131336386046957,   0.0521029999545313,   0.0373287282629393,   0.0067733153081290,   0.0151938118383802,
     0.0067733153081291,   0.0151938118383802,  -0.0000000000000001,  -0.0000000000000000,  -0.0047514582955382,
     0.0047514582955382,  -0.0400800102207108,  -0.0184343231458730,  -0.0093583235685706,  -0.0093583235685706,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0216060081523868,   0.0000000000000000,   0.0363008273515796,  -0.0016651372155552,
     0.0240875907371629,   0.0301279689042306,   0.0034385602962241,   0.0145628945002607,   0.0034385602962241,
     0.0145628945002608,  -0.0000000000000000,  -0.0000000000000001,  -0.0003589333227862,   0.0003589333227862,
    -0.0184343231458730,  -0.0229357807827095,  -0.0064824318329930,  -0.0064824318329930,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0136849862314935,   0.0000000000000000,   0.0229715465811821,  -0.0003444217856378,   0.0020144833160349,
     0.0007177209071113,  -0.0043132752731360,  -0.0012482753574777,  -0.0008128311775260,  -0.0008231662774004,
     0.0019344226233444,   0.0008348282358765,   0.0036588784722432,  -0.0014640764791231,  -0.0041051314180200,
    -0.0055444763226152,  -0.0057552297930153,  -0.0004384205182883,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0036028505506500,  -0.0061679110806137,
     0.0000000000000000,   0.0016250928015258,  -0.0003444217856378,   0.0020144833160349,   0.0007177209071114,
    -0.0008128311775260,  -0.0008231662774004,  -0.0043132752731359,  -0.0012482753574776,  -0.0019344226233444,
    -0.0008348282358765,   0.0014640764791231,  -0.0036588784722432,  -0.0041051314180199,  -0.0055444763226151,
    -0.0004384205182883,  -0.0057552297930154,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,  -0.0000000000000000,   0.0036028505506500,  -0.0061679110806137,   0.0000000000000000,
     0.0016250928015258,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.0040750797682010,
     0.0029937890223691,  -0.0040750797682011,  -0.0029937890223692,   0.0343668992671324,   0.0212444685965512,
    -0.0023508829431383,  -0.0023508829431384,  -0.0000000000000000,  -0.0000000000000000,  -0.0087750045743895,
     0.0087750045743896,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,  -0.0284803375911141,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0036028505506500,   0.0036028505506500,  -0.0284803375911141,
     0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0216060081523868,
    -0.0136849862314935,  -0.0061679110806137,  -0.0061679110806137,   0.0000000000000001,  -0.0000000000000000,
     0.0000000000000001,  -0.0163223349276599,  -0.0000000000000000,   0.0035434925070312,  -0.0467770759652298,
    -0.0406463373234668,  -0.0091989800056414,  -0.0177145792739417,  -0.0091989800056414,  -0.0177145792739417,
    -0.0000000000000000,  -0.0000000000000000,   0.0100745299607064,  -0.0100745299607064,   0.0332450227500307,
     0.0205129237239060,   0.0059203292338797,   0.0059203292338797,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0163223349276599,
     0.0000000000000001,  -0.0276572108545747,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0363008273515796,   0.0229715465811821,   0.0016250928015258,
     0.0016250928015258,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0276572108545747,
     0.0000000000000001};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jvals00.size()));

    ASSERT_EQ(fock.getNumberOfElements(1), static_cast<int32_t>(jvals01.size()));
    
    vlxtest::compare(jvals00, fock.getFock(0));

    vlxtest::compare(jvals01, fock.getFock(1));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestGenKForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getGenOODensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    fock.setFockType(fockmat::rgenk, 0);
    
    fock.setFockType(fockmat::rgenk, 1);
    
    ASSERT_TRUE(fockmat::rgenk == fock.getFockType(0));
    
    ASSERT_TRUE(fockmat::rgenk == fock.getFockType(1));
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    ASSERT_EQ(fock.getNumberOfRows(1), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(1), ncols);
    
    std::vector<double> kvals00{
     0.0299990177031620,  -0.0321081056615648,  -0.0207249277815864,  -0.0089955719965280,  -0.0102458268553737,
    -0.0089955719965279,  -0.0102458268553737,  -0.0000000000000000,   0.0000000000000000,   0.0113270928218322,
    -0.0113270928218321,   0.0001957110628680,  -0.0003730807868611,   0.0040568089603356,   0.0040568089603356,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0192916846409663,   0.0000000000000000,   0.0004892502662002,  -0.0321081056615648,
     0.1922919725138338,   0.1576675476035700,   0.0172974518937634,   0.0589840727101024,   0.0172974518937636,
     0.0589840727101025,  -0.0000000000000004,  -0.0000000000000003,  -0.0223443065373079,   0.0223443065373079,
    -0.2213813176066475,  -0.1405274188301533,  -0.0322921270349251,  -0.0322921270349254,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000003,
     0.1041034233989530,  -0.0000000000000000,   0.0063801035333442,  -0.0207249277815864,   0.1576675476035701,
     0.1454218874902447,   0.0297397874327186,   0.0639830366188770,   0.0297397874327185,   0.0639830366188771,
    -0.0000000000000003,  -0.0000000000000003,  -0.0272980237660741,   0.0272980237660742,  -0.1654192760955300,
    -0.1087439608826968,  -0.0186896909315823,  -0.0186896909315825,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000003,   0.0970881443053711,
    -0.0000000000000000,   0.0039901080506775,  -0.0089955719965280,   0.0172974518937634,   0.0297397874327186,
     0.0899990520596242,   0.0514576165140274,   0.0516445324874090,   0.0404326027975860,   0.0080908682328547,
     0.0145386941812365,  -0.0574416876473621,   0.0497933080547688,   0.1160081613849310,   0.1126751709142935,
     0.0310717391643595,   0.0191556241405285,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0171616958693438,   0.0500658006279216,  -0.0000000000000000,
    -0.0107109931168154,  -0.0102458268553737,   0.0589840727101024,   0.0639830366188770,   0.0514576165140274,
     0.0486908343378983,   0.0404326027975861,   0.0446058617502175,   0.0040683989831157,   0.0060862629734933,
    -0.0358750183532276,   0.0346269757824696,   0.0017774314209434,   0.0208490162195430,   0.0106270446498704,
     0.0073367587431949,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0007987891769662,   0.0562340055503058,  -0.0000000000000000,  -0.0032195770592885,
    -0.0089955719965279,   0.0172974518937636,   0.0297397874327185,   0.0516445324874090,   0.0404326027975861,
     0.0899990520596239,   0.0514576165140272,  -0.0080908682328544,  -0.0145386941812362,  -0.0497933080547686,
     0.0574416876473620,   0.1160081613849306,   0.1126751709142930,   0.0191556241405283,   0.0310717391643596,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0171616958693436,   0.0500658006279215,  -0.0000000000000000,  -0.0107109931168153,  -0.0102458268553737,
     0.0589840727101025,   0.0639830366188771,   0.0404326027975860,   0.0446058617502175,   0.0514576165140272,
     0.0486908343378983,  -0.0040683989831157,  -0.0060862629734933,  -0.0346269757824695,   0.0358750183532275,
     0.0017774314209432,   0.0208490162195427,   0.0073367587431948,   0.0106270446498703,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0007987891769659,
     0.0562340055503058,  -0.0000000000000000,  -0.0032195770592885,  -0.0000000000000000,  -0.0000000000000004,
    -0.0000000000000003,   0.0080908682328547,   0.0040683989831157,  -0.0080908682328544,  -0.0040683989831157,
     0.0497637082660318,   0.0273218195849261,  -0.0066073415563007,  -0.0066073415563003,   0.0000000000000011,
     0.0000000000000008,  -0.0077122054432206,   0.0077122054432209,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0403964960673728,   0.0000000000000001,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000003,  -0.0000000000000003,
     0.0145386941812365,   0.0060862629734933,  -0.0145386941812362,  -0.0060862629734933,   0.0273218195849261,
     0.0274235768213140,  -0.0044195878504926,  -0.0044195878504923,   0.0000000000000008,   0.0000000000000007,
    -0.0008116145934039,   0.0008116145934042,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0184095319173107,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0113270928218322,  -0.0223443065373079,  -0.0272980237660740,  -0.0574416876473622,
    -0.0358750183532276,  -0.0497933080547686,  -0.0346269757824695,  -0.0066073415563007,  -0.0044195878504926,
     0.0660981710031994,  -0.0524999459166119,  -0.1148087634019676,  -0.0943910171591704,  -0.0154985092199650,
    -0.0155605921605301,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0074893165914186,  -0.0535995930172967,   0.0000000000000000,   0.0072672624267722,
    -0.0113270928218322,   0.0223443065373079,   0.0272980237660741,   0.0497933080547687,   0.0346269757824696,
     0.0574416876473620,   0.0358750183532275,  -0.0066073415563003,  -0.0044195878504923,  -0.0524999459166119,
     0.0660981710031993,   0.1148087634019674,   0.0943910171591705,   0.0155605921605300,   0.0154985092199652,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0074893165914185,   0.0535995930172967,  -0.0000000000000000,  -0.0072672624267722,   0.0001957110628681,
    -0.2213813176066477,  -0.1654192760955298,   0.1160081613849310,   0.0017774314209435,   0.1160081613849305,
     0.0017774314209433,   0.0000000000000011,   0.0000000000000008,  -0.1148087634019676,   0.1148087634019674,
     0.6827905892207583,   0.4575021591931676,   0.0772583967587062,   0.0772583967587069,  -0.0000000000000002,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,
    -0.0206478537719072,  -0.0000000000000000,  -0.0286493828934799,  -0.0003730807868611,  -0.1405274188301533,
    -0.1087439608826967,   0.1126751709142935,   0.0208490162195430,   0.1126751709142931,   0.0208490162195426,
     0.0000000000000008,   0.0000000000000007,  -0.0943910171591704,   0.0943910171591704,   0.4575021591931677,
     0.3810810344541768,   0.0724754470842340,   0.0724754470842345,  -0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,  -0.0034436784384035,
    -0.0000000000000000,  -0.0235007685753581,   0.0040568089603356,  -0.0322921270349251,  -0.0186896909315823,
     0.0310717391643596,   0.0106270446498704,   0.0191556241405283,   0.0073367587431948,  -0.0077122054432206,
    -0.0008116145934039,  -0.0154985092199650,   0.0155605921605300,   0.0772583967587062,   0.0724754470842340,
     0.0308382718476638,   0.0149369651799740,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0119350582453950,   0.0007822186080120,   0.0000000000000000,
    -0.0023739785590416,   0.0040568089603356,  -0.0322921270349253,  -0.0186896909315825,   0.0191556241405285,
     0.0073367587431949,   0.0310717391643596,   0.0106270446498703,   0.0077122054432209,   0.0008116145934042,
    -0.0155605921605301,   0.0154985092199652,   0.0772583967587069,   0.0724754470842346,   0.0149369651799740,
     0.0308382718476640,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0119350582453952,   0.0007822186080120,   0.0000000000000000,  -0.0023739785590416,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000002,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0521467655241788,   0.0287757634295119,   0.0041464148099868,   0.0041464148099868,  -0.0000000000000001,
     0.0000000000000000,  -0.0000000000000000,  -0.0475492433468698,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,   0.0287757634295119,
     0.0258045349351589,   0.0058323009208816,   0.0058323009208816,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,  -0.0295570100384178,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,   0.0041464148099868,   0.0058323009208816,
     0.0166768358774699,   0.0040186341215168,   0.0030153233567606,  -0.0000000000000000,  -0.0000000000000000,
     0.0101700418110080,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,   0.0041464148099868,   0.0058323009208816,   0.0040186341215168,
     0.0166768358774698,  -0.0030153233567605,   0.0000000000000000,  -0.0000000000000000,   0.0101700418110079,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.0030153233567606,  -0.0030153233567605,
     0.0212953362638715,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000003,   0.0000000000000003,   0.0171616958693438,   0.0007987891769662,
    -0.0171616958693436,  -0.0007987891769659,  -0.0403964960673727,  -0.0184095319173107,  -0.0074893165914186,
    -0.0074893165914185,  -0.0000000000000001,  -0.0000000000000001,   0.0119350582453950,  -0.0119350582453952,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.1100524241949620,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0192916846409663,
     0.1041034233989533,   0.0970881443053710,   0.0500658006279216,   0.0562340055503059,   0.0500658006279215,
     0.0562340055503059,   0.0000000000000001,   0.0000000000000000,  -0.0535995930172967,   0.0535995930172967,
    -0.0206478537719073,  -0.0034436784384034,   0.0007822186080120,   0.0007822186080120,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.1305002002036307,  -0.0000000000000000,  -0.0043689973404840,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0475492433468698,  -0.0295570100384178,
     0.0101700418110080,   0.0101700418110079,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.1009943317751939,  -0.0000000000000000,   0.0004892502662002,   0.0063801035333442,   0.0039901080506775,
    -0.0107109931168154,  -0.0032195770592885,  -0.0107109931168153,  -0.0032195770592885,   0.0000000000000000,
     0.0000000000000000,   0.0072672624267722,  -0.0072672624267722,  -0.0286493828934799,  -0.0235007685753581,
    -0.0023739785590416,  -0.0023739785590416,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0043689973404840,  -0.0000000000000000,
     0.0232370717506740};
    
    std::vector<double> kvals01{
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0115766226069109,  -0.0044082240172522,   0.0012433609677109,   0.0012433609677109,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0186750249197919,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0518951774952845,
     0.0282072858836438,  -0.0121991192010193,  -0.0121991192010192,  -0.0000000000000002,  -0.0000000000000000,
    -0.0000000000000000,  -0.1039794169557690,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0424969707483207,   0.0274196992926053,
    -0.0151684096646303,  -0.0151684096646303,  -0.0000000000000002,  -0.0000000000000000,  -0.0000000000000000,
    -0.0930499441849744,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
    -0.0000000000000002,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,   0.0082412131111433,   0.0032059178635208,  -0.0134096192146708,
    -0.0057015965053197,   0.0013551039344826,  -0.0000000000000000,  -0.0000000000000000,  -0.0304509224521330,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0191394421037141,   0.0119478164322832,  -0.0129466745792242,  -0.0098145026648313,
     0.0014468260596332,  -0.0000000000000000,  -0.0000000000000000,  -0.0488505966683656,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0082412131111433,   0.0032059178635208,  -0.0057015965053197,  -0.0134096192146707,  -0.0013551039344828,
    -0.0000000000000000,  -0.0000000000000000,  -0.0304509224521329,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0191394421037141,
     0.0119478164322832,  -0.0098145026648313,  -0.0129466745792241,  -0.0014468260596334,  -0.0000000000000000,
    -0.0000000000000000,  -0.0488505966683656,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000002,
    -0.0000000000000002,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
    -0.0029911134851925,   0.0029911134851925,   0.0156976651924537,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000002,  -0.0000000000000002,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0051505436978847,
     0.0051505436978846,   0.0099058659914311,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0052460724924203,  -0.0011876453187414,   0.0097246449180823,   0.0064725883940889,
    -0.0004351976998665,   0.0000000000000000,   0.0000000000000000,   0.0312145268256564,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
     0.0052460724924203,   0.0011876453187415,  -0.0064725883940890,  -0.0097246449180822,  -0.0004351976998667,
     0.0000000000000000,  -0.0000000000000000,  -0.0312145268256564,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0373940234877797,
    -0.0205165396938501,  -0.0087995814066925,  -0.0087995814066925,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0126132788749459,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000002,
    -0.0000000000000002,  -0.0000000000000000,  -0.0000000000000000,  -0.0200807169189198,  -0.0198085661642814,
    -0.0100010461017224,  -0.0100010461017223,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0073121907127538,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000000,  -0.0102598519895891,  -0.0066249394618745,   0.0033561589691749,
     0.0017712078187861,  -0.0038304165810353,  -0.0000000000000000,  -0.0000000000000000,   0.0222116581753981,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0000000000000000,
    -0.0000000000000000,  -0.0102598519895892,  -0.0066249394618745,   0.0017712078187861,   0.0033561589691749,
     0.0038304165810354,   0.0000000000000000,   0.0000000000000000,   0.0222116581753980,  -0.0000000000000000,
     0.0097271149234484,   0.2020057108491713,   0.1490909010856937,  -0.1243546094304262,  -0.0115278640224865,
    -0.1243546094304257,  -0.0115278640224862,  -0.0000000000000011,  -0.0000000000000009,   0.1174425907410178,
    -0.1174425907410177,  -0.6573262353954621,  -0.4490108445650662,  -0.0776671299197574,  -0.0776671299197582,
     0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000001,  -0.0376968571270681,   0.0000000000000000,   0.0486431343537278,   0.0033775938489869,
     0.1297686137954843,   0.1006845249530000,  -0.1022520995507086,  -0.0181899899834886,  -0.1022520995507082,
    -0.0181899899834883,  -0.0000000000000007,  -0.0000000000000007,   0.0843441320080655,  -0.0843441320080655,
    -0.4233410961932380,  -0.3428912062158293,  -0.0632940317355481,  -0.0632940317355487,   0.0000000000000001,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000001,
    -0.0270741276682119,   0.0000000000000000,   0.0347662217014413,   0.0038851212677551,   0.0308894296058742,
     0.0157185756480167,  -0.0686240362732267,  -0.0270304576254541,  -0.0454825170291763,  -0.0211666345119982,
     0.0019347800315930,  -0.0041652040662049,   0.0470222663823368,  -0.0424530997636404,  -0.1682135455909137,
    -0.1407056032274311,  -0.0314997781898025,  -0.0230476575412153,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0223786014928474,  -0.0271408554656654,
     0.0000000000000000,   0.0110559411822781,   0.0038851212677551,   0.0308894296058742,   0.0157185756480168,
    -0.0454825170291764,  -0.0211666345119982,  -0.0686240362732264,  -0.0270304576254539,  -0.0019347800315936,
     0.0041652040662043,   0.0424530997636404,  -0.0470222663823367,  -0.1682135455909134,  -0.1407056032274309,
    -0.0230476575412151,  -0.0314997781898028,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0223786014928473,  -0.0271408554656655,   0.0000000000000000,
     0.0110559411822780,   0.0000000000000000,  -0.0000000000000002,  -0.0000000000000002,  -0.0132429735323381,
    -0.0008620625107151,   0.0132429735323379,   0.0008620625107149,   0.0316688540888159,   0.0139190699725128,
     0.0041679614456731,   0.0041679614456730,  -0.0000000000000001,  -0.0000000000000001,  -0.0090904954020718,
     0.0090904954020720,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0783014259428255,   0.0000000000000001,  -0.0000000000000000,   0.0000000000000001,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0046780080111981,   0.0046780080111981,  -0.0150088579506998,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0229966836505982,
    -0.0144537826790626,   0.0019527917439096,   0.0019527917439095,   0.0000000000000001,   0.0000000000000000,
    -0.0000000000000000,   0.0253168268874070,  -0.0000000000000000,   0.0184364135801722,  -0.0997565742598086,
    -0.0923022217692810,  -0.0462902022376409,  -0.0533659076774428,  -0.0462902022376408,  -0.0533659076774427,
    -0.0000000000000000,  -0.0000000000000000,   0.0473044896210186,  -0.0473044896210186,   0.0108687583420571,
    -0.0009174527373442,   0.0048090347381838,   0.0048090347381838,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0784856559121120,
     0.0000000000000000,  -0.0111169913717673,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,   0.0379103422202137,   0.0237269450608379,  -0.0077458522071919,
    -0.0077458522071918,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,  -0.0723702242993719,
     0.0000000000000000};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(kvals00.size()));
    
    ASSERT_EQ(fock.getNumberOfElements(1), static_cast<int32_t>(kvals01.size()));
    
    vlxtest::compare(kvals00, fock.getFock(0));
    
    vlxtest::compare(kvals01, fock.getFock(1));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeRestGenJKForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;
    
    auto dmat = vlxden::getGenOODensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    fock.setFockType(fockmat::rgenjk, 0);
    
    fock.setFockType(fockmat::rgenjk, 1);
    
    ASSERT_TRUE(fockmat::rgenjk == fock.getFockType(0));
    
    ASSERT_TRUE(fockmat::rgenjk == fock.getFockType(1));
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    ASSERT_EQ(fock.getNumberOfRows(1), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(1), ncols);
    
    std::vector<double> jkvals00 {
     1.9506232583017800,  -0.6301023124001359,  -0.2820223133651512,  -0.0855841462364752,  -0.1105568860552713,
    -0.0855841462364753,  -0.1105568860552713,   0.0000000000000000,   0.0000000000000000,   0.1316537054473836,
    -0.1316537054473837,   0.0397081350834750,   0.0054093167318909,   0.1092261424159810,   0.1092261424159810,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0119947672386319,   0.0000000000000000,  -0.0004498287953586,  -0.6301023124001359,
     1.4156396680211212,   1.0104663760057933,   0.4470195489985662,   0.4976943672243002,   0.4470195489985664,
     0.4976943672242998,   0.0000000000000000,   0.0000000000000001,  -0.4937993041339637,   0.4937993041339639,
     0.0723139310412301,   0.0765460697688625,  -0.3724232616846223,  -0.3724232616846224,  -0.0000000000000001,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0051517844067900,   0.0000000000000000,  -0.0093353802049784,  -0.2820223133651512,   1.0104663760057933,
     1.0780304891927255,   0.5474553904798218,   0.6793879401171810,   0.5474553904798221,   0.6793879401171812,
     0.0000000000000001,   0.0000000000000001,  -0.3613058998442115,   0.3613058998442114,   0.0620763891056038,
     0.0365693686594444,  -0.2768150402766583,  -0.2768150402766585,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,  -0.0102937471761079,
     0.0000000000000000,  -0.0076175150780667,  -0.0855841462364752,   0.4470195489985662,   0.5474553904798218,
     0.8819581887474632,   0.5782675890888804,   0.1025398585574204,   0.2845129683444072,   0.2833408449570983,
     0.4581992824073137,  -0.1143720967543613,   0.1433337706897864,   0.0986118292619560,   0.2582099365407972,
    -0.1348877346139763,  -0.0367796963699921,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.2028288623321846,  -0.0106521546516580,   0.0000000000000000,
    -0.1353386607050050,  -0.1105568860552713,   0.4976943672243002,   0.6793879401171810,   0.5782675890888804,
     0.7348843856514039,   0.2845129683444069,   0.4961823379047947,   0.1022916784192031,   0.2518140681767737,
    -0.0995951910938852,   0.1874698570331564,   0.0422474987932270,   0.1541594700329038,  -0.0995573682040412,
    -0.0562544824687535,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0144626636747316,  -0.0179707014950788,   0.0000000000000000,  -0.0133491319828473,
    -0.0855841462364753,   0.4470195489985664,   0.5474553904798221,   0.1025398585574204,   0.2845129683444069,
     0.8819581887474637,   0.5782675890888808,  -0.2833408449570988,  -0.4581992824073142,  -0.1433337706897868,
     0.1143720967543614,   0.0986118292619562,   0.2582099365407974,  -0.0367796963699919,  -0.1348877346139766,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.2028288623321848,  -0.0106521546516579,   0.0000000000000000,  -0.1353386607050053,  -0.1105568860552713,
     0.4976943672242998,   0.6793879401171812,   0.2845129683444072,   0.4961823379047947,   0.5782675890888808,
     0.7348843856514042,  -0.1022916784192035,  -0.2518140681767739,  -0.1874698570331567,   0.0995951910938853,
     0.0422474987932272,   0.1541594700329039,  -0.0562544824687534,  -0.0995573682040411,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0144626636747317,
    -0.0179707014950788,   0.0000000000000000,  -0.0133491319828473,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.2833408449570983,   0.1022916784192031,  -0.2833408449570988,  -0.1022916784192035,
     1.5370076342360648,   0.6501891685004967,  -0.2393665848587779,  -0.2393665848587785,  -0.0000000000000009,
    -0.0000000000000006,  -0.3808639617492732,   0.3808639617492735,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0548178069358310,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000002,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.4581992824073137,   0.2518140681767737,  -0.4581992824073142,  -0.2518140681767739,   0.6501891685004967,
     0.9572994815914420,  -0.0032443349190128,  -0.0032443349190131,  -0.0000000000000006,  -0.0000000000000005,
    -0.2845084373551091,   0.2845084373551092,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0356662763908731,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000002,   0.1316537054473836,  -0.4937993041339637,  -0.3613058998442116,  -0.1143720967543612,
    -0.0995951910938852,  -0.1433337706897868,  -0.1874698570331567,  -0.2393665848587779,  -0.0032443349190128,
     0.9257583021814209,  -0.2434688633201475,  -0.2661421952235598,  -0.2048330031576650,   0.0583308956663140,
     0.0363166662605941,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.2183128470847544,  -0.1828386688014583,   0.0000000000000000,  -0.0090142495126428,
    -0.1316537054473836,   0.4937993041339639,   0.3613058998442115,   0.1433337706897865,   0.1874698570331564,
     0.1143720967543614,   0.0995951910938853,  -0.2393665848587785,  -0.0032443349190131,  -0.2434688633201475,
     0.9257583021814205,   0.2661421952235604,   0.2048330031576649,  -0.0363166662605940,  -0.0583308956663142,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.2183128470847545,   0.1828386688014585,   0.0000000000000000,   0.0090142495126428,   0.0397081350834749,
     0.0723139310412303,   0.0620763891056036,   0.0986118292619560,   0.0422474987932269,   0.0986118292619563,
     0.0422474987932271,  -0.0000000000000009,  -0.0000000000000006,  -0.2661421952235598,   0.2661421952235604,
     1.0476684138943111,   0.2864946338705849,  -0.1114280385862648,  -0.1114280385862655,   0.0000000000000002,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000001,  -0.0000000000000000,  -0.0000000000000001,
    -0.0883791518135138,   0.0000000000000000,   0.0251923354345919,   0.0054093167318909,   0.0765460697688625,
     0.0365693686594443,   0.2582099365407972,   0.1541594700329038,   0.2582099365407973,   0.1541594700329040,
    -0.0000000000000006,  -0.0000000000000005,  -0.2048330031576650,   0.2048330031576650,   0.2864946338705848,
     0.6867198845215301,   0.0851105610040942,   0.0851105610040939,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,  -0.0596635561057959,
     0.0000000000000000,   0.0200886051834035,   0.1092261424159810,  -0.3724232616846223,  -0.2768150402766583,
    -0.1348877346139764,  -0.0995573682040412,  -0.0367796963699919,  -0.0562544824687534,  -0.3808639617492732,
    -0.2845084373551091,   0.0583308956663140,  -0.0363166662605940,  -0.1114280385862648,   0.0851105610040942,
     0.9325663055136304,   0.0395630702260438,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0564623724442742,   0.2326744713214726,   0.0000000000000000,
     0.2346300940129688,   0.1092261424159810,  -0.3724232616846225,  -0.2768150402766585,  -0.0367796963699921,
    -0.0562544824687535,  -0.1348877346139766,  -0.0995573682040411,   0.3808639617492735,   0.2845084373551092,
     0.0363166662605941,  -0.0583308956663142,  -0.1114280385862655,   0.0851105610040938,   0.0395630702260438,
     0.9325663055136302,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0564623724442742,   0.2326744713214728,   0.0000000000000000,   0.2346300940129690,
     0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000002,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     1.5312928779493835,   0.6454909481391045,   0.2409582543578310,   0.2409582543578312,  -0.0000000000000001,
    -0.0000000000000000,   0.0000000000000000,  -0.0545791545741102,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.6454909481391045,
     0.9501600600597859,   0.3565460534461940,   0.3565460534461940,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,  -0.0313431250536750,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.2409582543578310,   0.3565460534461940,
     0.8933464656237546,   0.0446801686473438,   0.2931458749024866,   0.0000000000000000,   0.0000000000000000,
     0.2191010578385216,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.2409582543578312,   0.3565460534461940,   0.0446801686473438,
     0.8933464656237544,  -0.2931458749024871,   0.0000000000000000,   0.0000000000000000,   0.2191010578385219,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.2931458749024866,  -0.2931458749024871,
     1.3360135500484536,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,   0.2028288623321846,   0.0144626636747316,
    -0.2028288623321848,  -0.0144626636747317,  -0.0548178069358311,  -0.0356662763908731,  -0.2183128470847544,
    -0.2183128470847545,  -0.0000000000000001,  -0.0000000000000001,  -0.0564623724442742,   0.0564623724442742,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     1.3563485166331912,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0119947672386319,
    -0.0051517844067903,  -0.0102937471761078,  -0.0106521546516580,  -0.0179707014950789,  -0.0106521546516579,
    -0.0179707014950789,   0.0000000000000001,   0.0000000000000000,  -0.1828386688014583,   0.1828386688014585,
    -0.0883791518135137,  -0.0596635561057960,   0.2326744713214726,   0.2326744713214728,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     1.3643561568301732,   0.0000000000000000,   0.0051397468499496,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0545791545741102,  -0.0313431250536750,
     0.2191010578385216,   0.2191010578385219,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     1.3570702979975815,   0.0000000000000000,  -0.0004498287953586,  -0.0093353802049784,  -0.0076175150780667,
    -0.1353386607050050,  -0.0133491319828473,  -0.1353386607050053,  -0.0133491319828473,   0.0000000000000002,
     0.0000000000000002,  -0.0090142495126428,   0.0090142495126428,   0.0251923354345919,   0.0200886051834035,
     0.2346300940129688,   0.2346300940129690,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0051397468499494,   0.0000000000000000,
     1.3348966545409506};

    std::vector<double> jkvals01 {
     0.0000000000000002,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0146906546024805,   0.0010779495861418,  -0.0019322045389865,  -0.0019322045389865,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0115880399057295,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0523108224137781,
     0.0199678955906820,   0.0162280858330891,   0.0162280858330890,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0104252650253094,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0321604857775579,   0.0328362385158559,
     0.0166038514788529,   0.0166038514788531,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0117572695380408,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000001,   0.0000000000000001,   0.0053054175051147,   0.0036712027289274,   0.0047830686683988,
     0.0040759341502677,   0.0067950556019194,   0.0000000000000000,   0.0000000000000000,   0.0120529624408502,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0112481815730463,   0.0171779725682382,   0.0104501238642688,   0.0081681701100305,
     0.0045407519851050,   0.0000000000000000,   0.0000000000000000,   0.0134214381204822,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0053054175051149,   0.0036712027289274,   0.0040759341502677,   0.0047830686683989,  -0.0067950556019194,
     0.0000000000000000,   0.0000000000000000,   0.0120529624408501,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0112481815730463,
     0.0171779725682384,   0.0081681701100305,   0.0104501238642689,  -0.0045407519851050,   0.0000000000000000,
     0.0000000000000000,   0.0134214381204822,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000002,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000002,
     0.0000000000000002,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0068599587318813,  -0.0068599587318813,   0.0530361333418111,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000002,   0.0000000000000002,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000002,   0.0068202001696377,
    -0.0068202001696376,   0.0325830712016713,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000002,  -0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0042568440986561,   0.0004697786731690,  -0.0024068879735959,  -0.0035444354358427,
    -0.0042665681864101,   0.0000000000000000,   0.0000000000000000,  -0.0110654669042436,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0042568440986561,  -0.0004697786731691,   0.0035444354358428,   0.0024068879735958,  -0.0042665681864101,
    -0.0000000000000000,   0.0000000000000000,   0.0110654669042436,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000003,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0427659969536419,
    -0.0163521065978959,   0.0005893185706525,   0.0005893185706527,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0538767666251155,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000002,
     0.0000000000000002,   0.0000000000000000,   0.0000000000000000,  -0.0167879293728262,  -0.0260629954011376,
    -0.0010879065435080,  -0.0010879065435079,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0337136567350582,  -0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000002,   0.0000000000000000,  -0.0084567951475521,  -0.0063399242041115,  -0.0148666185552055,
    -0.0026480488553627,  -0.0137195925677437,   0.0000000000000000,   0.0000000000000000,  -0.0103709997076387,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,  -0.0084567951475520,  -0.0063399242041115,  -0.0026480488553627,  -0.0148666185552057,
     0.0137195925677438,   0.0000000000000000,   0.0000000000000000,  -0.0103709997076386,   0.0000000000000000,
    -0.0359943921328398,  -0.0977997109401087,  -0.0744334445598151,   0.1379012400466842,   0.0419154876992469,
     0.1379012400466839,   0.0419154876992466,   0.0000000000000009,   0.0000000000000009,  -0.1269455073320942,
     0.1269455073320941,   0.5771662149540405,   0.4121421982733202,   0.0589504827826162,   0.0589504827826170,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,
    -0.0000000000000001,  -0.0055151591777055,   0.0000000000000000,   0.0239585203494314,  -0.0067078682800973,
    -0.0815934323211585,  -0.0404285871445388,   0.1091292201431568,   0.0473157789840100,   0.1091292201431564,
     0.0473157789840099,   0.0000000000000007,   0.0000000000000005,  -0.0850619986536379,   0.0850619986536379,
     0.3864724499014920,   0.2970196446504103,   0.0503291680695621,   0.0503291680695627,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,
    -0.0002958447947751,   0.0000000000000000,   0.0111768714609229,  -0.0045739648390307,  -0.0268604629738044,
    -0.0142831338337941,   0.0599974857269547,   0.0245339069104987,   0.0438568546741243,   0.0195203019571974,
     0.0019340652150958,   0.0058348605379579,  -0.0397045094378504,   0.0395249468053942,   0.1600032827548737,
     0.1296166505822007,   0.0199893186037719,   0.0221708165046387,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000002,   0.0000000000000000,   0.0000000000000000,   0.0151729003915474,   0.0148050333044380,
     0.0000000000000000,  -0.0078057555792265,  -0.0045739648390307,  -0.0268604629738044,  -0.0142831338337940,
     0.0438568546741244,   0.0195203019571974,   0.0599974857269546,   0.0245339069104987,  -0.0019340652150952,
    -0.0058348605379573,  -0.0395249468053942,   0.0397045094378503,   0.1600032827548736,   0.1296166505822007,
     0.0221708165046385,   0.0199893186037720,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000002,   0.0000000000000000,  -0.0151729003915473,   0.0148050333044381,   0.0000000000000000,
    -0.0078057555792264,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0213931330687401,
     0.0068496405554533,  -0.0213931330687401,  -0.0068496405554533,   0.0370649444454489,   0.0285698672205896,
    -0.0088697273319497,  -0.0088697273319498,   0.0000000000000001,   0.0000000000000001,  -0.0084595137467072,
     0.0084595137467072,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000002,   0.0213407507605973,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0025276930901019,   0.0025276930901019,  -0.0419518172315284,
     0.0000000000000002,   0.0000000000000000,  -0.0000000000000002,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0202153326541754,
    -0.0129161897839244,  -0.0142886139051370,  -0.0142886139051369,   0.0000000000000001,  -0.0000000000000000,
     0.0000000000000002,  -0.0579614967427268,   0.0000000000000000,  -0.0113494285661098,   0.0062024223293490,
     0.0110095471223474,   0.0278922422263581,   0.0179367491295594,   0.0278922422263580,   0.0179367491295593,
     0.0000000000000000,   0.0000000000000000,  -0.0271554296996058,   0.0271554296996058,   0.0556212871580043,
     0.0419433001851562,   0.0070316237295756,   0.0070316237295756,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0458409860567922,
     0.0000000000000002,  -0.0441974303373821,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0346913124829455,   0.0222161481015263,   0.0109960378102435,
     0.0109960378102434,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0170558025902225,
     0.0000000000000002};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals00.size()));
    
    ASSERT_EQ(fock.getNumberOfElements(1), static_cast<int32_t>(jkvals01.size()));
    
    vlxtest::compare(jkvals00, fock.getFock(0));
    
    vlxtest::compare(jkvals01, fock.getFock(1));
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeUnrestJKForH2O)
{
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    const int nrows = 24, ncols = 24;

    auto dmat = vlxden::getUnrestDensityMatrixForH2O();
    
    CElectronRepulsionIntegralsDriver eridrv(MPI_COMM_WORLD);
    
    auto qqdata = eridrv.compute(ericut::qq, 1.0e-16, mh2o, mbas);
    
    CAOFockMatrix fock(dmat);
    
    eridrv.compute(fock, dmat, mh2o, mbas, qqdata);
    
    ASSERT_EQ(fock.getNumberOfRows(0), nrows);
    
    ASSERT_EQ(fock.getNumberOfColumns(0), ncols);
    
    std::vector<double> jkvals_alpha {
    11.9265743865906479,  -3.4468889676290653,  -1.5420674441930364,  -0.4809335070297658,  -0.6112914875112113,
    -0.4809335070297655,  -0.6112914875112112,  -0.0000000000000002,  -0.0000000000000002,   0.7436130796047931,
    -0.7436130796047922,  -0.0027946806253360,   0.0014671912568964,   0.5851688823225205,   0.5851688823225198,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
    -0.0000000000000002,  -0.0031369781416341,  -0.0000000000000000,   0.0091133647207835,  -3.4468889676290653,
     7.1305548456278380,   4.8855154023547716,   1.9655823921043849,   2.2850315060244624,   1.9655823921043849,
     2.2850315060244619,  -0.0000000000000004,  -0.0000000000000001,  -2.3779392738920935,   2.3779392738920944,
     0.1016906647191691,   0.0400014312491542,  -1.8435301793553340,  -1.8435301793553345,   0.0000000000000002,
     0.0000000000000002,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0414454207450744,   0.0000000000000000,  -0.0622783952085829,  -1.5420674441930367,   4.8855154023547716,
     5.0154680032449086,   2.4575748225156047,   3.0941849698205099,   2.4575748225156056,   3.0941849698205104,
    -0.0000000000000001,   0.0000000000000000,  -1.6966548176957570,   1.6966548176957568,   0.0596181730937638,
     0.0537388339496559,  -1.3056156736456648,  -1.3056156736456646,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0335407091892086,
     0.0000000000000000,  -0.0481539515152068,  -0.4809335070297657,   1.9655823921043849,   2.4575748225156047,
     4.3089562407405397,   2.7297654861551020,   0.6391925554059588,   1.3843264595298159,   1.2053412286132730,
     2.0627499118823258,  -0.6609355756708808,   0.7979993345502649,   0.9699888473315637,   1.6445969356046997,
    -0.4881367102257520,  -0.0996038324629373,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,   1.0763182591993061,   0.1121368783786796,   0.0000000000000000,
    -0.7008881086890134,  -0.6112914875112113,   2.2850315060244624,   3.0941849698205099,   2.7297654861551020,
     3.4077322625669244,   1.3843264595298159,   2.3164854448889329,   0.4166376838975686,   1.1214608088905633,
    -0.5083876882884004,   0.9379476741600675,   0.3472838508622426,   0.9175218393256526,  -0.3747509162321329,
    -0.2088005452654261,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0924904261741575,   0.0282062685716005,   0.0000000000000000,  -0.0809166566269998,
    -0.4809335070297655,   1.9655823921043851,   2.4575748225156056,   0.6391925554059589,   1.3843264595298159,
     4.3089562407405371,   2.7297654861551020,  -1.2053412286132728,  -2.0627499118823254,  -0.7979993345502646,
     0.6609355756708803,   0.9699888473315632,   1.6445969356046988,  -0.0996038324629373,  -0.4881367102257515,
     0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
    -1.0763182591993057,   0.1121368783786799,   0.0000000000000000,  -0.7008881086890143,  -0.6112914875112112,
     2.2850315060244619,   3.0941849698205104,   1.3843264595298161,   2.3164854448889329,   2.7297654861551024,
     3.4077322625669240,  -0.4166376838975684,  -1.1214608088905638,  -0.9379476741600674,   0.5083876882884006,
     0.3472838508622422,   0.9175218393256527,  -0.2088005452654261,  -0.3747509162321327,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0924904261741574,
     0.0282062685716008,   0.0000000000000000,  -0.0809166566269997,  -0.0000000000000002,  -0.0000000000000004,
    -0.0000000000000001,   1.2053412286132732,   0.4166376838975686,  -1.2053412286132730,  -0.4166376838975686,
     7.8298865572612302,   2.8718654854332204,  -1.2189892156860487,  -1.2189892156860489,   0.0000000000000004,
     0.0000000000000001,  -1.8201380083148926,   1.8201380083148937,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000001,  -0.0266127528329694,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,  -0.0000000000000002,  -0.0000000000000001,   0.0000000000000000,
     2.0627499118823263,   1.1214608088905633,  -2.0627499118823254,  -1.1214608088905638,   2.8718654854332204,
     4.2419447867381166,  -0.0113435039577359,  -0.0113435039577359,   0.0000000000000001,  -0.0000000000000001,
    -1.2967691932050367,   1.2967691932050369,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,  -0.0283478418602922,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000001,   0.7436130796047928,  -2.3779392738920939,  -1.6966548176957574,  -0.6609355756708803,
    -0.5083876882884004,  -0.7979993345502646,  -0.9379476741600674,  -1.2189892156860487,  -0.0113435039577359,
     4.7018083764533261,  -1.3144033652914253,  -1.8427303994134177,  -1.3252729320206491,   0.2184328842889648,
     0.1000827389901412,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -1.1701829171033504,  -1.0854092101493136,   0.0000000000000000,   0.0078490577256688,
    -0.7436130796047921,   2.3779392738920939,   1.6966548176957565,   0.7979993345502648,   0.9379476741600675,
     0.6609355756708801,   0.5083876882884006,  -1.2189892156860489,  -0.0113435039577359,  -1.3144033652914253,
     4.7018083764533278,   1.8427303994134181,   1.3252729320206496,  -0.1000827389901412,  -0.2184328842889646,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -1.1701829171033498,   1.0854092101493131,   0.0000000000000000,  -0.0078490577256693,  -0.0027946806253360,
     0.1016906647191692,   0.0596181730937639,   0.9699888473315637,   0.3472838508622425,   0.9699888473315632,
     0.3472838508622422,   0.0000000000000004,   0.0000000000000001,  -1.8427303994134177,   1.8427303994134181,
     7.7694871894246527,   2.8525126070449023,  -0.3560970252022504,  -0.3560970252022503,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000002,
     0.0424736935693398,   0.0000000000000001,   0.0338582498847020,   0.0014671912568964,   0.0400014312491542,
     0.0537388339496559,   1.6445969356046997,   0.9175218393256526,   1.6445969356046988,   0.9175218393256527,
     0.0000000000000001,  -0.0000000000000001,  -1.3252729320206491,   1.3252729320206496,   2.8525126070449023,
     4.2377761455213161,   0.6286241521174318,   0.6286241521174318,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0350649743986403,
     0.0000000000000001,   0.0303831798694353,   0.5851688823225207,  -1.8435301793553340,  -1.3056156736456646,
    -0.4881367102257519,  -0.3747509162321327,  -0.0996038324629373,  -0.2088005452654262,  -1.8201380083148926,
    -1.2967691932050367,   0.2184328842889648,  -0.1000827389901412,  -0.3560970252022504,   0.6286241521174318,
     4.5745373094756498,   0.2521689482626353,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.3370948921074133,   1.1736546255219396,   0.0000000000000000,
     1.1759600479724523,   0.5851688823225201,  -1.8435301793553345,  -1.3056156736456641,  -0.0996038324629373,
    -0.2088005452654261,  -0.4881367102257514,  -0.3747509162321327,   1.8201380083148937,   1.2967691932050369,
     0.1000827389901412,  -0.2184328842889646,  -0.3560970252022504,   0.6286241521174318,   0.2521689482626353,
     4.5745373094756498,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.3370948921074132,   1.1736546255219391,  -0.0000000000000000,   1.1759600479724528,
     0.0000000000000000,   0.0000000000000002,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,  -0.0000000000000001,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     7.5717174486177568,   2.7900611165563864,   1.0518168406967665,   1.0518168406967674,   0.0000000000000001,
     0.0000000000000002,  -0.0000000000000002,   0.0378955335341061,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   2.7900611165563864,
     4.1956960546057864,   1.6342890240331325,   1.6342890240331329,   0.0000000000000001,   0.0000000000000001,
    -0.0000000000000001,   0.0333371775355829,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   1.0518168406967665,   1.6342890240331325,
     4.3790771785365523,   0.1868185409249087,   1.4669126876077387,   0.0000000000000001,   0.0000000000000000,
     1.1598092295923990,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   1.0518168406967672,   1.6342890240331329,   0.1868185409249087,
     4.3790771785365532,  -1.4669126876077381,   0.0000000000000001,   0.0000000000000000,   1.1598092295923981,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   1.4669126876077387,  -1.4669126876077381,
     6.6947154893414744,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   1.0763182591993061,   0.0924904261741575,
    -1.0763182591993057,  -0.0924904261741574,  -0.0266127528329694,  -0.0283478418602922,  -1.1701829171033504,
    -1.1701829171033498,   0.0000000000000002,   0.0000000000000000,  -0.3370948921074133,   0.3370948921074132,
     0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,
     6.7556748903898658,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0031369781416341,
     0.0414454207450741,   0.0335407091892090,   0.1121368783786805,   0.0282062685716003,   0.1121368783786807,
     0.0282062685716005,   0.0000000000000000,   0.0000000000000001,  -1.0854092101493134,   1.0854092101493131,
     0.0424736935693397,   0.0350649743986403,   1.1736546255219396,   1.1736546255219389,  -0.0000000000000002,
    -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     6.7443931301791276,   0.0000000000000000,   0.0553561297037818,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,   0.0378955335341061,   0.0333371775355829,
     1.1598092295923990,   1.1598092295923981,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     6.6862475984754148,   0.0000000000000000,   0.0091133647207836,  -0.0622783952085830,  -0.0481539515152070,
    -0.7008881086890135,  -0.0809166566269997,  -0.7008881086890141,  -0.0809166566269995,   0.0000000000000001,
     0.0000000000000001,   0.0078490577256686,  -0.0078490577256694,   0.0338582498847019,   0.0303831798694353,
     1.1759600479724523,   1.1759600479724523,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0553561297037813,   0.0000000000000000,
     6.6892130929551508};

    std::vector<double> jkvals_beta {
    11.9705288176907381,  -3.5018617599240747,  -1.5757385650936362,  -0.4921203536122575,  -0.6267033400337321,
    -0.4921203536122572,  -0.6267033400337320,  -0.0000000000000002,  -0.0000000000000002,   0.7580259999222546,
    -0.7580259999222541,  -0.0042936323747411,   0.0007834900182936,   0.5963921199459382,   0.5963921199459374,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
    -0.0000000000000002,   0.0110442395854444,  -0.0000000000000000,  -0.0144952468343042,  -3.5018617599240747,
     7.3496776787467324,   5.0462243450985982,   2.0138504662252883,   2.3617560199086687,   2.0138504662252892,
     2.3617560199086678,   0.0000000000000000,   0.0000000000000001,  -2.4486058259877721,   2.4486058259877717,
     0.1209211689368835,   0.0483965088165928,  -1.8943992109784198,  -1.8943992109784198,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
    -0.0460621522700096,   0.0000000000000000,   0.0886169833222800,  -1.5757385650936362,   5.0462243450985982,
     5.1415342365497061,   2.4901106184669000,   3.1550950443004782,   2.4901106184669000,   3.1550950443004782,
     0.0000000000000001,   0.0000000000000002,  -1.7513883333851434,   1.7513883333851430,   0.0743984288909749,
     0.0587478755092083,  -1.3445354127984332,  -1.3445354127984328,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,  -0.0418466330211103,
     0.0000000000000000,   0.0833732964318586,  -0.4921203536122575,   2.0138504662252883,   2.4901106184668995,
     4.2962228464302568,   2.7372323687596967,   0.6491328196255739,   1.3994317819830269,   1.2124177501432041,
     2.0560330969088252,  -0.6866843264579396,   0.8130608359564387,   0.9811318414849883,   1.6427856977133461,
    -0.5060443987754097,  -0.1100866202505442,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,   1.0716184602564285,   0.0867344674446972,   0.0000000000000000,
    -0.6518756374163015,  -0.6267033400337321,   2.3617560199086687,   3.1550950443004782,   2.7372323687596967,
     3.4347132104746532,   1.3994317819830266,   2.3460826263209822,   0.4198126964106246,   1.1197759975178736,
    -0.5388492343172994,   0.9639009369485960,   0.3565571196754161,   0.9183955192305364,  -0.3964753000943080,
    -0.2274089970286871,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0905267320481700,  -0.0108627709432746,   0.0000000000000000,  -0.0105542981914609,
    -0.4921203536122572,   2.0138504662252892,   2.4901106184669000,   0.6491328196255739,   1.3994317819830266,
     4.2962228464302541,   2.7372323687596962,  -1.2124177501432032,  -2.0560330969088239,  -0.8130608359564382,
     0.6866843264579388,   0.9811318414849882,   1.6427856977133448,  -0.1100866202505443,  -0.5060443987754093,
    -0.0000000000000002,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000002,
    -1.0716184602564278,   0.0867344674446977,   0.0000000000000000,  -0.6518756374163016,  -0.6267033400337317,
     2.3617560199086678,   3.1550950443004782,   1.3994317819830269,   2.3460826263209822,   2.7372323687596962,
     3.4347132104746523,  -0.4198126964106240,  -1.1197759975178736,  -0.9639009369485961,   0.5388492343172995,
     0.3565571196754155,   0.9183955192305362,  -0.2274089970286871,  -0.3964753000943078,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,  -0.0905267320481697,
    -0.0108627709432741,   0.0000000000000000,  -0.0105542981914610,  -0.0000000000000002,   0.0000000000000000,
     0.0000000000000001,   1.2124177501432039,   0.4198126964106246,  -1.2124177501432032,  -0.4198126964106240,
     7.9183201262142235,   2.9031283844359748,  -1.2360203803872225,  -1.2360203803872227,   0.0000000000000001,
     0.0000000000000000,  -1.8417760670987875,   1.8417760670987882,   0.0000000000000001,   0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,  -0.0000000000000004,  -0.0273482662397237,   0.0000000000000001,
     0.0000000000000000,  -0.0000000000000003,  -0.0000000000000002,   0.0000000000000001,   0.0000000000000001,
     2.0560330969088252,   1.1197759975178736,  -2.0560330969088243,  -1.1197759975178736,   2.9031283844359748,
     4.2489891538814826,  -0.0217859906222834,  -0.0217859906222833,   0.0000000000000000,   0.0000000000000000,
    -1.3068374562777001,   1.3068374562777001,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
    -0.0000000000000000,  -0.0000000000000002,  -0.0339789760010130,   0.0000000000000002,   0.0000000000000000,
    -0.0000000000000001,   0.7580259999222544,  -2.4486058259877721,  -1.7513883333851434,  -0.6866843264579396,
    -0.5388492343172994,  -0.8130608359564381,  -0.9639009369485960,  -1.2360203803872225,  -0.0217859906222834,
     4.7312048537838178,  -1.3395828136771246,  -1.8677578840692088,  -1.3369933325497021,   0.2423713565353374,
     0.1137691619977393,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,  -1.1744981093565121,  -1.0626563801774036,   0.0000000000000000,  -0.0370922668403683,
    -0.7580259999222534,   2.4486058259877717,   1.7513883333851428,   0.8130608359564386,   0.9639009369485962,
     0.6866843264579391,   0.5388492343172994,  -1.2360203803872227,  -0.0217859906222833,  -1.3395828136771246,
     4.7312048537838187,   1.8677578840692091,   1.3369933325497025,  -0.1137691619977393,  -0.2423713565353371,
    -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000001,
    -1.1744981093565112,   1.0626563801774029,   0.0000000000000000,   0.0370922668403681,  -0.0042936323747411,
     0.1209211689368835,   0.0743984288909750,   0.9811318414849880,   0.3565571196754161,   0.9811318414849878,
     0.3565571196754155,   0.0000000000000001,   0.0000000000000000,  -1.8677578840692088,   1.8677578840692091,
     7.8440476946749458,   2.8778278347571655,  -0.3664781959483035,  -0.3664781959483036,  -0.0000000000000001,
    -0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000003,
     0.0498949124104971,   0.0000000000000000,   0.0354116214742682,   0.0007834900182936,   0.0483965088165929,
     0.0587478755092083,   1.6427856977133461,   0.9183955192305364,   1.6427856977133453,   0.9183955192305362,
     0.0000000000000000,   0.0000000000000000,  -1.3369933325497021,   1.3369933325497025,   2.8778278347571655,
     4.2442031690058855,   0.6215738874012864,   0.6215738874012864,  -0.0000000000000001,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,   0.0387199876004132,
     0.0000000000000001,   0.0345063681032741,   0.5963921199459382,  -1.8943992109784198,  -1.3445354127984332,
    -0.5060443987754097,  -0.3964753000943080,  -0.1100866202505444,  -0.2274089970286871,  -1.8417760670987875,
    -1.3068374562777001,   0.2423713565353374,  -0.1137691619977393,  -0.3664781959483035,   0.6215738874012864,
     4.5889794031862081,   0.2582677751180372,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,  -0.3389715467414841,   1.1960806793267982,   0.0000000000000000,
     1.1422614696910087,   0.5963921199459379,  -1.8943992109784198,  -1.3445354127984330,  -0.1100866202505443,
    -0.2274089970286871,  -0.5060443987754094,  -0.3964753000943078,   1.8417760670987882,   1.3068374562777001,
     0.1137691619977393,  -0.2423713565353371,  -0.3664781959483036,   0.6215738874012864,   0.2582677751180372,
     4.5889794031862090,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000001,   0.3389715467414840,   1.1960806793267977,   0.0000000000000000,   1.1422614696910085,
     0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000002,   0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000001,
     8.4992622715432589,   3.3155613426521322,   1.2554012757055508,   1.2554012757055513,   0.0000000000000001,
     0.0000000000000000,  -0.0000000000000001,   0.0678629055999586,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   0.0000000000000000,  -0.0000000000000001,
    -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   3.3155613426521322,
     4.5557896302417094,   1.7726737208178982,   1.7726737208178982,   0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0531481322859097,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   1.2554012757055508,   1.7726737208178982,
     4.4476219585760921,   0.2373767209853494,   1.4949991779410523,   0.0000000000000000,   0.0000000000000000,
     1.1888970861937065,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,   0.0000000000000000,
    -0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
    -0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,  -0.0000000000000001,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   1.2554012757055513,   1.7726737208178982,   0.2373767209853494,
     4.4476219585760921,  -1.4949991779410516,   0.0000000000000000,   0.0000000000000000,   1.1888970861937056,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000001,   0.0000000000000001,  -0.0000000000000001,
    -0.0000000000000000,   0.0000000000000002,   0.0000000000000001,  -0.0000000000000004,  -0.0000000000000002,
     0.0000000000000000,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
    -0.0000000000000001,   0.0000000000000001,   0.0000000000000001,   1.4949991779410523,  -1.4949991779410516,
     6.8236087991130585,  -0.0000000000000001,   0.0000000000000000,  -0.0000000000000000,   0.0000000000000000,
    -0.0000000000000002,   0.0000000000000001,   0.0000000000000001,   1.0716184602564285,   0.0905267320481700,
    -1.0716184602564278,  -0.0905267320481697,  -0.0273482662397237,  -0.0339789760010130,  -1.1744981093565121,
    -1.1744981093565112,   0.0000000000000003,   0.0000000000000001,  -0.3389715467414841,   0.3389715467414840,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0000000000000001,
     6.7717763984195214,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0110442395854444,
    -0.0460621522700094,  -0.0418466330211099,   0.0867344674446967,  -0.0108627709432743,   0.0867344674446971,
    -0.0108627709432742,   0.0000000000000001,   0.0000000000000002,  -1.0626563801774036,   1.0626563801774032,
     0.0498949124104972,   0.0387199876004132,   1.1960806793267984,   1.1960806793267977,  -0.0000000000000001,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000001,
     6.7993385576021996,   0.0000000000000000,  -0.0043866054320444,  -0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000001,   0.0000000000000000,   0.0000000000000000,   0.0678629055999586,   0.0531481322859097,
     1.1888970861937065,   1.1888970861937056,  -0.0000000000000000,   0.0000000000000000,   0.0000000000000000,
     6.8135903486522080,   0.0000000000000000,  -0.0144952468343043,   0.0886169833222798,   0.0833732964318581,
    -0.6518756374163017,  -0.0105542981914609,  -0.6518756374163017,  -0.0105542981914609,  -0.0000000000000003,
    -0.0000000000000001,  -0.0370922668403686,   0.0370922668403683,   0.0354116214742681,   0.0345063681032741,
     1.1422614696910089,   1.1422614696910085,   0.0000000000000001,   0.0000000000000000,   0.0000000000000000,
     0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  -0.0043866054320447,   0.0000000000000000,
     6.8118753176166233};
    
    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals_alpha.size()));

    ASSERT_EQ(fock.getNumberOfElements(0), static_cast<int32_t>(jkvals_beta.size()));
    
    vlxtest::compare(jkvals_alpha, fock.getFock(0), 1.0e-12);

    vlxtest::compare(jkvals_beta, fock.getFock(0, true), 1.0e-12);
}
