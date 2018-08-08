//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "BoysFunction.hpp"

#include <cmath> 
#include <array>

#include "MathConst.hpp"

CBoysFunction::CBoysFunction()

    : _order(-1)
{

}

CBoysFunction::CBoysFunction(const int32_t order)

    : _order(order)
{
    _table = CMemBlock2D<double>(7, 121 * (_order + 1)); 
    
    _setTable(); 
}

void CBoysFunction::compute(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments,
                            const int32_t              iOrder) const
{
    if (iOrder < 0) return;

    if (iOrder > 28) return;

    if (iOrder == 0)
    {
        _computeBF00(values, arguments);

        return;
    }

    if (iOrder == 1)
    {
        _computeBF01(values, arguments);

        return;
    }

    if (iOrder == 2)
    {
        _computeBF02(values, arguments);

        return;
    }

    if (iOrder == 3)
    {
        _computeBF03(values, arguments);

        return;
    }

    if (iOrder == 4)
    {
        _computeBF04(values, arguments);

        return;
    }

    if (iOrder == 5)
    {
        _computeBF05(values, arguments);

        return;
    }

    if (iOrder == 6)
    {
        _computeBF06(values, arguments);

        return;
    }

    if (iOrder == 7)
    {
        _computeBF07(values, arguments);

        return;
    }

    if (iOrder == 8)
    {
        _computeBF08(values, arguments);

        return;
    }

    if (iOrder == 9)
    {
        _computeBF09(values, arguments);

        return;
    }

    if (iOrder == 10)
    {
        _computeBF10(values, arguments);

        return;
    }

    if (iOrder == 11)
    {
        _computeBF11(values, arguments);

        return;
    }

    if (iOrder == 12)
    {
        _computeBF12(values, arguments);

        return;
    }

    if (iOrder == 13)
    {
        _computeBF13(values, arguments);

        return;
    }

    if (iOrder == 14)
    {
        _computeBF14(values, arguments);

        return;
    }

    if (iOrder == 15)
    {
        _computeBF15(values, arguments);

        return;
    }

    if (iOrder == 16)
    {
        _computeBF16(values, arguments);

        return;
    }

    if (iOrder == 17)
    {
        _computeBF17(values, arguments);

        return;
    }

    if (iOrder == 18)
    {
        _computeBF18(values, arguments);

        return;
    }

    if (iOrder == 19)
    {
        _computeBF19(values, arguments);

        return;
    }

    if (iOrder == 20)
    {
        _computeBF20(values, arguments);

        return;
    }

    if (iOrder == 21)
    {
        _computeBF21(values, arguments);

        return;
    }

    if (iOrder == 22)
    {
        _computeBF22(values, arguments);

        return;
    }

    if (iOrder == 23)
    {
        _computeBF23(values, arguments);

        return;
    }

    if (iOrder == 24)
    {
        _computeBF23(values, arguments);

        return;
    }

    if (iOrder == 25)
    {
        _computeBF25(values, arguments);

        return;
    }

    if (iOrder == 26)
    {
        _computeBF26(values, arguments);

        return;
    }

    if (iOrder == 27)
    {
        _computeBF27(values, arguments);

        return;
    }

    if (iOrder == 28)
    {
        _computeBF28(values, arguments);

        return;
    }
}

CBoysFunction::~CBoysFunction()
{

}

void
CBoysFunction::_setTable()
{
    if (_order < 0 ) return;

    if (_order > 28) return;

    _generateTable00();

    if (_order == 0) return;

    _generateTable01();

    if (_order == 1) return;

    _generateTable02();

    if (_order == 2) return;

    _generateTable03();

    if (_order == 3) return;

    _generateTable04();

    if (_order == 4) return;

    _generateTable05();

    if (_order == 5) return;

    _generateTable06();

    if (_order == 6) return;

    _generateTable07();

    if (_order == 7) return;

    _generateTable08();

    if (_order == 8) return;

    _generateTable09();

    if (_order == 9) return;

    _generateTable10();

    if (_order == 10) return;

    _generateTable11();

    if (_order == 11) return;

    _generateTable12();

    if (_order == 12) return;

    _generateTable13();

    if (_order == 13) return;

    _generateTable14();

    if (_order == 14) return;

    _generateTable15();

    if (_order == 15) return;

    _generateTable16();

    if (_order == 16) return;

    _generateTable17();

    if (_order == 17) return;

    _generateTable18();

    if (_order == 18) return;

    _generateTable19();

    if (_order == 19) return;

    _generateTable20();

    if (_order == 20) return;

    _generateTable21();

    if (_order == 21) return;

    _generateTable22();

    if (_order == 22) return;

    _generateTable23();

    if (_order == 23) return;

    _generateTable24();

    if (_order == 24) return;

    _generateTable25();
    
    if (_order == 25) return;
    
    _generateTable26();
    
    if (_order == 26) return;
    
    _generateTable27();
    
    if (_order == 27) return;
    
    _generateTable28();
}

void
CBoysFunction::_computeBF00(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double w  = argv[i] - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(pnt);

            val00[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;
        }
        else
        {
            double fia = 1.0 / argv[i];

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 361)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                val00[i] -= f * std::exp(-argv[i]);
            }
        }
    }
}

void
CBoysFunction::_computeBF01(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(121 + pnt);

            val01[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            val00[i] = 2.0 * fa * val01[i] + std::exp(-fa);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 381)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                val01[i] = pf * (val00[i] - fx);
            }
            else
            {
                val01[i] = pf * val00[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF02(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 2> ft{1.0, 1.0 / 3.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(242 + pnt);

            val02[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 401)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;
            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF03(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 3> ft{1.0, 1.0 / 3.0, 1.0 / 5.0 };

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(363 + pnt);

            val03[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 421)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;
                
                pf += fia;
                
                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];
                
                pf += fia;
                
                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF04(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 4> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(484 + pnt);

            val04[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 441)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;
                
            }
            else
            {
                val01[i] = pf * val00[i];
                
                pf += fia;
                
                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF05(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 5> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(605 + pnt);

            val05[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 461)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];
                
                pf += fia;
                
                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF06(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 6> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                             1.0 / 11.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(726 + pnt);

            val06[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 481)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];
                
                pf += fia;
                
                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF07(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 7> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                             1.0 / 11.0, 1.0 / 13.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(847 + pnt);

            val07[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 501)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;
                
                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF08(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 8> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                             1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(968 + pnt);

            val08[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 521)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;
                
                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF09(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 9> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                             1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1089 + pnt);

            val09[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 541)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];
                
                pf += fia;
                
                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];

                pf += fia;

                val09[i] = pf * val08[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF10(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 10> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1210 + pnt);

            val10[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 561)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];

                pf += fia;

                val10[i] = pf * val09[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF11(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 11> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1331 + pnt);

            val11[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 581)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];
                
                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];

                pf += fia;

                val11[i] = pf * val10[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF12(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 12> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1452 + pnt);

            val12[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 601)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;
                
                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];

                pf += fia;

                val12[i] = pf * val11[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF13(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 13> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1573 + pnt);

            val13[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 621)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];

                pf += fia;

                val13[i] = pf * val12[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF14(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 14> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1694 + pnt);

            val14[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 641)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];
                
                pf += fia;
                
                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];

                pf += fia;

                val14[i] = pf * val13[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF15(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 15> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1815 + pnt);

            val15[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 661)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];

                pf += fia;

                val15[i] = pf * val14[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF16(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 16> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(1936 + pnt);

            val16[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 681)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];
                
                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];

                pf += fia;

                val16[i] = pf * val15[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF17(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 17> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2057 + pnt);

            val17[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 701)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;
                
                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];

                pf += fia;

                val17[i] = pf * val16[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF18(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 18> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2178 + pnt);

            val18[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 721)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];
                
                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];

                pf += fia;

                val18[i] = pf * val17[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF19(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 19> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2299 + pnt);

            val19[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 741)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;
                
                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];

                pf += fia;

                val19[i] = pf * val18[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF20(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 20> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2420 + pnt);

            val20[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 761)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];

                pf += fia;

                val20[i] = pf * val19[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF21(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 21> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2541 + pnt);

            val21[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 781)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];
                
                pf += fia;
                
                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];

                pf += fia;

                val21[i] = pf * val20[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF22(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 22> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2662 + pnt);

            val22[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 801)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];

                pf += fia;

                val22[i] = pf * val21[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF23(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 23> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2783 + pnt);

            val23[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 821)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];
                
                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];

                pf += fia;

                val23[i] = pf * val22[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF24(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 24> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    auto val24 = values.data(24);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(2904 + pnt);

            val24[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val23[i] = ft[23] * (f2a * val24[i] + fx);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 841)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

                pf += fia;

                val24[i] = pf * val23[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];

                pf += fia;
                
                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];
                
                pf += fia;
                
                val23[i] = pf * val22[i];

                pf += fia;

                val24[i] = pf * val23[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF25(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 25> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    auto val24 = values.data(24);

    auto val25 = values.data(25);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(3025 + pnt);

            val25[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val24[i] = ft[24] * (f2a * val25[i] + fx);

            val23[i] = ft[23] * (f2a * val24[i] + fx);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 861)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

                pf += fia;

                val24[i] = pf * val23[i] - rterm;

                pf += fia;

                val25[i] = pf * val24[i] - rterm;

            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];

                pf += fia;

                val09[i] = pf * val08[i];
                
                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];
                
                pf += fia;
                
                val23[i] = pf * val22[i];
                
                pf += fia;
                
                val24[i] = pf * val23[i];

                pf += fia;

                val25[i] = pf * val24[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF26(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 26> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0,
                              1.0 / 51.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    auto val24 = values.data(24);

    auto val25 = values.data(25);

    auto val26 = values.data(26);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(3146 + pnt);

            val26[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val25[i] = ft[25] * (f2a * val26[i] + fx);

            val24[i] = ft[24] * (f2a * val25[i] + fx);

            val23[i] = ft[23] * (f2a * val24[i] + fx);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 881)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

                pf += fia;

                val24[i] = pf * val23[i] - rterm;

                pf += fia;

                val25[i] = pf * val24[i] - rterm;

                pf += fia;

                val26[i] = pf * val25[i] - rterm;
            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];

                pf += fia;

                val09[i] = pf * val08[i];

                pf += fia;
                
                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];
                
                pf += fia;
                
                val23[i] = pf * val22[i];
                
                pf += fia;
                
                val24[i] = pf * val23[i];
                
                pf += fia;
                
                val25[i] = pf * val24[i];

                pf += fia;

                val26[i] = pf * val25[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF27(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 27> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0,
                              1.0 / 51.0, 1.0 / 53.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    auto val24 = values.data(24);

    auto val25 = values.data(25);

    auto val26 = values.data(26);

    auto val27 = values.data(27);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(3267 + pnt);

            val27[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val26[i] = ft[26] * (f2a * val27[i] + fx);

            val25[i] = ft[25] * (f2a * val26[i] + fx);

            val24[i] = ft[24] * (f2a * val25[i] + fx);

            val23[i] = ft[23] * (f2a * val24[i] + fx);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 901)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

                pf += fia;

                val24[i] = pf * val23[i] - rterm;

                pf += fia;

                val25[i] = pf * val24[i] - rterm;

                pf += fia;

                val26[i] = pf * val25[i] - rterm;

                pf += fia;

                val27[i] = pf * val26[i] - rterm;
            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];

                pf += fia;

                val09[i] = pf * val08[i];

                pf += fia;

                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];
                
                pf += fia;
                
                val23[i] = pf * val22[i];
                
                pf += fia;
                
                val24[i] = pf * val23[i];
                
                pf += fia;
                
                val25[i] = pf * val24[i];
                
                pf += fia;
                
                val26[i] = pf * val25[i];

                pf += fia;

                val27[i] = pf * val26[i];
            }
        }
    }
}

void
CBoysFunction::_computeBF28(      CMemBlock2D<double>& values,
                            const CMemBlock<double>&   arguments) const
{
    auto fpi = 0.5 * std::sqrt(mathconst::getPiValue());

    std::array<double, 28> ft{1.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0,
                              1.0 / 11.0, 1.0 / 13.0, 1.0 / 15.0, 1.0 / 17.0,
                              1.0 / 19.0, 1.0 / 21.0, 1.0 / 23.0, 1.0 / 25.0,
                              1.0 / 27.0, 1.0 / 29.0, 1.0 / 31.0, 1.0 / 33.0,
                              1.0 / 35.0, 1.0 / 37.0, 1.0 / 39.0, 1.0 / 41.0,
                              1.0 / 43.0, 1.0 / 45.0, 1.0 / 47.0, 1.0 / 49.0,
                              1.0 / 51.0, 1.0 / 53.0, 1.0 / 55.0};

    auto ndim = arguments.size();

    auto argv = arguments.data();

    auto val00 = values.data(0);

    auto val01 = values.data(1);

    auto val02 = values.data(2);

    auto val03 = values.data(3);

    auto val04 = values.data(4);

    auto val05 = values.data(5);

    auto val06 = values.data(6);

    auto val07 = values.data(7);

    auto val08 = values.data(8);

    auto val09 = values.data(9);

    auto val10 = values.data(10);

    auto val11 = values.data(11);

    auto val12 = values.data(12);

    auto val13 = values.data(13);

    auto val14 = values.data(14);

    auto val15 = values.data(15);

    auto val16 = values.data(16);

    auto val17 = values.data(17);

    auto val18 = values.data(18);

    auto val19 = values.data(19);

    auto val20 = values.data(20);

    auto val21 = values.data(21);

    auto val22 = values.data(22);

    auto val23 = values.data(23);

    auto val24 = values.data(24);

    auto val25 = values.data(25);

    auto val26 = values.data(26);

    auto val27 = values.data(27);

    auto val28 = values.data(28);

    for (int32_t i = 0; i < ndim; i++)
    {
        int32_t pnt = (argv[i] > 1.0e5) ? 1000000 : static_cast<int32_t>(10.0 * argv[i] + 0.5);

        if (pnt < 121)
        {
            double fa = argv[i];

            double w  = fa - 0.1 * pnt;

            double w2 = w *  w;

            double w4 = w2 * w2;

            auto tbvec = _table.data(3388 + pnt);

            val28[i]  = tbvec[0] + tbvec[1] * w + tbvec[2] * w2 + tbvec[3] * w2 * w

                      + tbvec[4] * w4 + tbvec[5] * w4 * w + tbvec[6] * w4 * w2;

            double f2a  = fa + fa;

            double fx = std::exp(-fa);

            val27[i] = ft[27] * (f2a * val28[i] + fx);

            val26[i] = ft[26] * (f2a * val27[i] + fx);

            val25[i] = ft[25] * (f2a * val26[i] + fx);

            val24[i] = ft[24] * (f2a * val25[i] + fx);

            val23[i] = ft[23] * (f2a * val24[i] + fx);

            val22[i] = ft[22] * (f2a * val23[i] + fx);

            val21[i] = ft[21] * (f2a * val22[i] + fx);

            val20[i] = ft[20] * (f2a * val21[i] + fx);

            val19[i] = ft[19] * (f2a * val20[i] + fx);

            val18[i] = ft[18] * (f2a * val19[i] + fx);

            val17[i] = ft[17] * (f2a * val18[i] + fx);

            val16[i] = ft[16] * (f2a * val17[i] + fx);

            val15[i] = ft[15] * (f2a * val16[i] + fx);

            val14[i] = ft[14] * (f2a * val15[i] + fx);

            val13[i] = ft[13] * (f2a * val14[i] + fx);

            val12[i] = ft[12] * (f2a * val13[i] + fx);

            val11[i] = ft[11] * (f2a * val12[i] + fx);

            val10[i] = ft[10] * (f2a * val11[i] + fx);

            val09[i] = ft[9] * (f2a * val10[i] + fx);

            val08[i] = ft[8] * (f2a * val09[i] + fx);

            val07[i] = ft[7] * (f2a * val08[i] + fx);

            val06[i] = ft[6] * (f2a * val07[i] + fx);

            val05[i] = ft[5] * (f2a * val06[i] + fx);

            val04[i] = ft[4] * (f2a * val05[i] + fx);

            val03[i] = ft[3] * (f2a * val04[i] + fx);

            val02[i] = ft[2] * (f2a * val03[i] + fx);

            val01[i] = ft[1] * (f2a * val02[i] + fx);

            val00[i] = ft[0] * (f2a * val01[i] + fx);
        }
        else
        {
            double fia = 1.0 / argv[i];

            double pf = 0.5 * fia;

            val00[i] = fpi * std::sqrt(fia);

            if (pnt < 921)
            {
                double fia2 = fia * fia;

                double f = 0.4999489092 * fia - 0.2473631686 * fia2

                         + 0.3211809090 * fia2 * fia - 0.3811559346 * fia2 * fia2;

                double fx = std::exp(-argv[i]);

                val00[i] -= f * fx;

                double rterm = pf * fx;

                val01[i] = pf * val00[i] - rterm;

                pf += fia;

                val02[i] = pf * val01[i] - rterm;

                pf += fia;

                val03[i] = pf * val02[i] - rterm;

                pf += fia;

                val04[i] = pf * val03[i] - rterm;

                pf += fia;

                val05[i] = pf * val04[i] - rterm;

                pf += fia;

                val06[i] = pf * val05[i] - rterm;

                pf += fia;

                val07[i] = pf * val06[i] - rterm;

                pf += fia;

                val08[i] = pf * val07[i] - rterm;

                pf += fia;

                val09[i] = pf * val08[i] - rterm;

                pf += fia;

                val10[i] = pf * val09[i] - rterm;

                pf += fia;

                val11[i] = pf * val10[i] - rterm;

                pf += fia;

                val12[i] = pf * val11[i] - rterm;

                pf += fia;

                val13[i] = pf * val12[i] - rterm;

                pf += fia;

                val14[i] = pf * val13[i] - rterm;

                pf += fia;

                val15[i] = pf * val14[i] - rterm;

                pf += fia;

                val16[i] = pf * val15[i] - rterm;

                pf += fia;

                val17[i] = pf * val16[i] - rterm;

                pf += fia;

                val18[i] = pf * val17[i] - rterm;

                pf += fia;

                val19[i] = pf * val18[i] - rterm;

                pf += fia;

                val20[i] = pf * val19[i] - rterm;

                pf += fia;

                val21[i] = pf * val20[i] - rterm;

                pf += fia;

                val22[i] = pf * val21[i] - rterm;

                pf += fia;

                val23[i] = pf * val22[i] - rterm;

                pf += fia;

                val24[i] = pf * val23[i] - rterm;

                pf += fia;

                val25[i] = pf * val24[i] - rterm;

                pf += fia;

                val26[i] = pf * val25[i] - rterm;

                pf += fia;

                val27[i] = pf * val26[i] - rterm;

                pf += fia;

                val28[i] = pf * val27[i] - rterm;
            }
            else
            {
                val01[i] = pf * val00[i];

                pf += fia;

                val02[i] = pf * val01[i];

                pf += fia;

                val03[i] = pf * val02[i];

                pf += fia;

                val04[i] = pf * val03[i];

                pf += fia;

                val05[i] = pf * val04[i];

                pf += fia;

                val06[i] = pf * val05[i];

                pf += fia;

                val07[i] = pf * val06[i];

                pf += fia;

                val08[i] = pf * val07[i];

                pf += fia;

                val09[i] = pf * val08[i];

                pf += fia;

                val10[i] = pf * val09[i];
                
                pf += fia;
                
                val11[i] = pf * val10[i];
                
                pf += fia;
                
                val12[i] = pf * val11[i];
                
                pf += fia;
                
                val13[i] = pf * val12[i];
                
                pf += fia;
                
                val14[i] = pf * val13[i];
                
                pf += fia;
                
                val15[i] = pf * val14[i];
                
                pf += fia;
                
                val16[i] = pf * val15[i];
                
                pf += fia;
                
                val17[i] = pf * val16[i];
                
                pf += fia;
                
                val18[i] = pf * val17[i];
                
                pf += fia;
                
                val19[i] = pf * val18[i];
                
                pf += fia;
                
                val20[i] = pf * val19[i];
                
                pf += fia;
                
                val21[i] = pf * val20[i];
                
                pf += fia;
                
                val22[i] = pf * val21[i];
                
                pf += fia;
                
                val23[i] = pf * val22[i];
                
                pf += fia;
                
                val24[i] = pf * val23[i];
                
                pf += fia;
                
                val25[i] = pf * val24[i];
                
                pf += fia;
                
                val26[i] = pf * val25[i];
                
                pf += fia;
                
                val27[i] = pf * val26[i];

                pf += fia;

                val28[i] = pf * val27[i];
            }
        }
    }
}

void
CBoysFunction::_loadTable(const CBFTable& bfData,
                          const int32_t   identifier)
{
    for (int32_t i = 0; i < 121; i++)
    {
        auto mat = _table.data(121 * identifier + i);
        
        for (int32_t j = 0; j < 7; j++) mat[j] = bfData[i][j];
    }
}

void
CBoysFunction::_generateTable00()
{
#include "Table_BFunc_00.tab"

    _loadTable(BFTable, 0);
}

void
CBoysFunction::_generateTable01()
{
#include "Table_BFunc_01.tab"

    _loadTable(BFTable, 1);
}

void
CBoysFunction::_generateTable02()
{
#include "Table_BFunc_02.tab"

    _loadTable(BFTable, 2);
}

void
CBoysFunction::_generateTable03()
{
#include "Table_BFunc_03.tab"

    _loadTable(BFTable, 3);
}

void
CBoysFunction::_generateTable04()
{
#include "Table_BFunc_04.tab"

    _loadTable(BFTable, 4);
}

void
CBoysFunction::_generateTable05()
{
#include "Table_BFunc_05.tab"

    _loadTable(BFTable, 5);
}

void
CBoysFunction::_generateTable06()
{
#include "Table_BFunc_06.tab"

    _loadTable(BFTable, 6);
}

void
CBoysFunction::_generateTable07()
{
#include "Table_BFunc_07.tab"

    _loadTable(BFTable, 7);
}

void
CBoysFunction::_generateTable08()
{
#include "Table_BFunc_08.tab"

    _loadTable(BFTable, 8);
}

void
CBoysFunction::_generateTable09()
{
#include "Table_BFunc_09.tab"

    _loadTable(BFTable, 9);
}

void
CBoysFunction::_generateTable10()
{
#include "Table_BFunc_10.tab"

    _loadTable(BFTable, 10);
}

void
CBoysFunction::_generateTable11()
{
#include "Table_BFunc_11.tab"

    _loadTable(BFTable, 11);
}

void
CBoysFunction::_generateTable12()
{
#include "Table_BFunc_12.tab"

    _loadTable(BFTable, 12);
}

void
CBoysFunction::_generateTable13()
{
#include "Table_BFunc_13.tab"

    _loadTable(BFTable, 13);
}

void
CBoysFunction::_generateTable14()
{
#include "Table_BFunc_14.tab"

    _loadTable(BFTable, 14);
}

void
CBoysFunction::_generateTable15()
{
#include "Table_BFunc_15.tab"

    _loadTable(BFTable, 15);
}

void
CBoysFunction::_generateTable16()
{
#include "Table_BFunc_16.tab"

    _loadTable(BFTable, 16);
}

void
CBoysFunction::_generateTable17()
{
#include "Table_BFunc_17.tab"

    _loadTable(BFTable, 17);
}

void
CBoysFunction::_generateTable18()
{
#include "Table_BFunc_18.tab"

    _loadTable(BFTable, 18);
}

void
CBoysFunction::_generateTable19()
{
#include "Table_BFunc_19.tab"

    _loadTable(BFTable, 19);
}

void
CBoysFunction::_generateTable20()
{
#include "Table_BFunc_20.tab"

    _loadTable(BFTable, 20);
}

void
CBoysFunction::_generateTable21()
{
#include "Table_BFunc_21.tab"

    _loadTable(BFTable, 21);
}

void
CBoysFunction::_generateTable22()
{
#include "Table_BFunc_22.tab"

    _loadTable(BFTable, 22);
}

void
CBoysFunction::_generateTable23()
{
#include "Table_BFunc_23.tab"

    _loadTable(BFTable, 23);
}

void
CBoysFunction::_generateTable24()
{
#include "Table_BFunc_24.tab"

    _loadTable(BFTable, 24);
}

void
CBoysFunction::_generateTable25()
{
#include "Table_BFunc_25.tab"

    _loadTable(BFTable, 25);
}

void
CBoysFunction::_generateTable26()
{
#include "Table_BFunc_26.tab"

    _loadTable(BFTable, 26);
}

void
CBoysFunction::_generateTable27()
{
#include "Table_BFunc_27.tab"

    _loadTable(BFTable, 27);
}

void
CBoysFunction::_generateTable28()
{
#include "Table_BFunc_28.tab"

    _loadTable(BFTable, 28);
}













