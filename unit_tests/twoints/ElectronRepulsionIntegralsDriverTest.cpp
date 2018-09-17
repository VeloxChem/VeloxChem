//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronRepulsionIntegralsDriverTest.hpp"

#include "mpi.h"

#include "ElectronRepulsionIntegralsDriver.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPPDForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
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
