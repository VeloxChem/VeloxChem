#include "BasisFunction.hpp"

#include <cmath>

#include "MathConst.hpp"

CBasisFunction::CBasisFunction(const std::vector<double>& exponents, const std::vector<double>& norms, const int64_t angmom)

    : _exponents(exponents)

    , _norms(norms)

    , _angmom(angmom)

{
}

auto
CBasisFunction::setExponents(const std::vector<double>& exponents) -> void
{
    _exponents = exponents;
}

auto
CBasisFunction::setNormalizationFactors(const std::vector<double>& norms) -> void
{
    _norms = norms;
}

auto
CBasisFunction::setAngularMomentum(const int64_t angmom) -> void
{
    _angmom = angmom;
}

auto
CBasisFunction::add(const double exponent, const double norm) -> void
{
    _exponents.push_back(exponent);

    _norms.push_back(norm);
}

auto
CBasisFunction::normalize() -> void
{
    // NOTE: Currently implemented for l = 0..6

    if (_angmom > 6) return;

    if (const auto npgtos = _exponents.size(); npgtos > 0)
    {
        // uncontracted GTO

        if (npgtos == 1) _norms[0] = 1.0;

        // normalize primitive GTOs

        _rescale();

        // compute contracted GTO self-overlap

        double fovl = 0.0;

        for (size_t i = 0; i < npgtos; i++)
        {
            fovl += _overlap(i, i);

            for (size_t j = i + 1; j < npgtos; j++)
            {
                fovl += 2.0 * _overlap(i, j);
            }
        }

        // renormaliza primitive GTOs

        fovl = 1.0 / std::sqrt(fovl);

        for (size_t i = 0; i < npgtos; i++)
        {
            _norms[i] *= fovl;
        }
    }
}

auto
CBasisFunction::getExponents() const -> std::vector<double>
{
    return _exponents;
}

auto
CBasisFunction::getNormalizationFactors() const -> std::vector<double>
{
    return _norms;
}

auto
CBasisFunction::getAngularMomentum() const -> int64_t
{
    return _angmom;
}

auto
CBasisFunction::getNumberOfPrimitiveFunctions() const -> int64_t
{
    return static_cast<int64_t>(_exponents.size());
}

auto
CBasisFunction::_rescale() -> void
{
    if (const auto npgtos = _exponents.size(); npgtos > 0)
    {
        const auto fpi = 2.0 / mathconst::getPiValue();

        for (size_t i = 0; i < npgtos; i++)
        {
            _norms[i] *= std::pow(_exponents[i] * fpi, 0.75);
        }

        if (_angmom == 1)
        {
            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= 2.0 * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 2)
        {
            const double fact = 2.0 / std::sqrt(3.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i];
            }
        }
        else if (_angmom == 3)
        {
            const double fact = 4.0 / std::sqrt(15.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 4)
        {
            const double fact = 2.0 / std::sqrt(105.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i];
            }
        }
        else if (_angmom == 5)
        {
            const double fact = 4.0 / std::sqrt(945.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i]

                             * std::sqrt(_exponents[i]);
            }
        }
        else if (_angmom == 6)
        {
            const double fact = 4.0 / std::sqrt(10395.0);

            for (size_t i = 0; i < npgtos; i++)
            {
                _norms[i] *= fact * _exponents[i] * _exponents[i]

                             * _exponents[i];
            }
        }
        else
        {
            // implement l > 6
        }
    }
}

auto
CBasisFunction::_overlap(const size_t i, const size_t j) const -> double
{
    const auto fab = 1.0 / (_exponents[i] + _exponents[j]);

    const auto fovl = _norms[i] * _norms[j]

                      * std::pow(mathconst::getPiValue() * fab, 1.5);

    if (_angmom == 0)
    {
        return fovl;
    }
    else if (_angmom == 1)
    {
        return 0.5 * fab * fovl;
    }
    else if (_angmom == 2)
    {
        return 3.0 * fab * fab * fovl;
    }
    else if (_angmom == 3)
    {
        return 7.5 * fab * fab * fab * fovl;
    }
    else if (_angmom == 4)
    {
        return 420.0 * fab * fab * fab * fab * fovl;
    }
    else if (_angmom == 5)
    {
        return 1890.0 * fab * fab * fab * fab * fab * fovl;
    }
    else if (_angmom == 6)
    {
        return 41580.0 * fab * fab * fab * fab * fab * fab * fovl;
    }
    else
    {
        return 0.0;
    }
}
