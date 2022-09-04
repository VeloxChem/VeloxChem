//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "GridPartitioner.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>

CGridPartitioner::CGridPartitioner(const CGridBox& box)

    : _numberOfPointsThreshold(1024)
{
    _boxes.clear();

    _boxes.push_back(box);
}

CGridPartitioner::~CGridPartitioner()
{
}

void
CGridPartitioner::partitionGridPoints()
{
    auto box = _boxes.begin();

    while (box != _boxes.end())
    {
        auto npoints = box->getNumberOfGridPoints();

        if (npoints > _numberOfPointsThreshold)
        {
            auto newboxes = divideBoxIntoTwo(*box);

            _boxes.insert(_boxes.end(), newboxes.begin(), newboxes.end());

            box = _boxes.erase(box);
        }
        else
        {
            ++box;
        }
    }
}

std::list<CGridBox>
CGridPartitioner::divideBoxIntoEight(const CGridBox& box) const
{
    // box info

    auto boxdim = box.getBoxDimension();

    auto xmin = boxdim[0];

    auto ymin = boxdim[1];

    auto zmin = boxdim[2];

    auto xmax = boxdim[3];

    auto ymax = boxdim[4];

    auto zmax = boxdim[5];

    // grid points info

    auto npoints = box.getNumberOfGridPoints();

    auto xcoords = box.getCoordinatesX();

    auto ycoords = box.getCoordinatesY();

    auto zcoords = box.getCoordinatesZ();

    auto weights = box.getWeights();

    // find the center for dividing the box

    auto xhalf = findCenter(xcoords, npoints);

    auto yhalf = findCenter(ycoords, npoints);

    auto zhalf = findCenter(zcoords, npoints);

    // sub boxes and grid points

    std::vector<std::array<double, 6>> subboxdims;

    std::vector<CMemBlock2D<double>> subgridpoints;

    std::vector<int32_t> subnumpoints;

    std::vector<std::array<double, 2>> xrange({std::array<double, 2>({xmin, xhalf}), std::array<double, 2>({xhalf, xmax})});

    std::vector<std::array<double, 2>> yrange({std::array<double, 2>({ymin, yhalf}), std::array<double, 2>({yhalf, ymax})});

    std::vector<std::array<double, 2>> zrange({std::array<double, 2>({zmin, zhalf}), std::array<double, 2>({zhalf, zmax})});

    for (const auto& x : xrange)
    {
        for (const auto& y : yrange)
        {
            for (const auto& z : zrange)
            {
                subboxdims.push_back(std::array<double, 6>({x[0], y[0], z[0], x[1], y[1], z[1]}));

                subgridpoints.push_back(CMemBlock2D<double>(npoints, 4));

                subnumpoints.push_back(0);
            }
        }
    }

    for (int32_t g = 0; g < npoints; g++)
    {
        int32_t xval = (xcoords[g] < xhalf) ? 0 : 1;

        int32_t yval = (ycoords[g] < yhalf) ? 0 : 1;

        int32_t zval = (zcoords[g] < zhalf) ? 0 : 1;

        // x y z  id
        // ---------
        // 0 0 0   0
        // 0 0 1   1
        // 0 1 0   2
        // 0 1 1   3
        // 1 0 0   4
        // 1 0 1   5
        // 1 1 0   6
        // 1 1 1   7

        int32_t box_id = (xval << 2) | (yval << 1) | (zval << 0);

        auto count = subnumpoints[box_id];

        subgridpoints[box_id].data(0)[count] = xcoords[g];

        subgridpoints[box_id].data(1)[count] = ycoords[g];

        subgridpoints[box_id].data(2)[count] = zcoords[g];

        subgridpoints[box_id].data(3)[count] = weights[g];

        subnumpoints[box_id]++;
    }

    std::list<CGridBox> newboxes;

    for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxdims.size()); box_id++)
    {
        auto count = subnumpoints[box_id];

        if (count > 0) newboxes.push_back(CGridBox(subboxdims[box_id], subgridpoints[box_id].slice(0, count)));
    }

    return newboxes;
}

std::list<CGridBox>
CGridPartitioner::divideBoxIntoTwo(const CGridBox& box) const
{
    // box info

    auto boxdim = box.getBoxDimension();

    auto xmin = boxdim[0];

    auto ymin = boxdim[1];

    auto zmin = boxdim[2];

    auto xmax = boxdim[3];

    auto ymax = boxdim[4];

    auto zmax = boxdim[5];

    auto xlen = xmax - xmin;

    auto ylen = ymax - ymin;

    auto zlen = zmax - zmin;

    // the dimension that will be divided

    int32_t icart = 0;

    if ((ylen >= xlen) && (ylen >= zlen)) icart = 1;

    else if ((zlen >= xlen) && (zlen >= ylen)) icart = 2;

    // grid points info

    auto npoints = box.getNumberOfGridPoints();

    auto xcoords = box.getCoordinatesX();

    auto ycoords = box.getCoordinatesY();

    auto zcoords = box.getCoordinatesZ();

    auto weights = box.getWeights();

    std::vector<const double*> coords({xcoords, ycoords, zcoords});

    // find the center for dividing the box

    auto center = findMedian(coords[icart], npoints);

    // sub boxes and grid points

    std::vector<std::array<double, 6>> subboxdims;

    std::vector<CMemBlock2D<double>> subgridpoints;

    std::vector<int32_t> subnumpoints;

    for (int32_t box_id = 0; box_id < 2; box_id++)
    {
        subboxdims.push_back(std::array<double, 6>(boxdim));

        subgridpoints.push_back(CMemBlock2D<double>(npoints, 4));

        subnumpoints.push_back(0);
    }

    // update subboxdims

    // icart == 0
    // std::array<double, 6>({xmin, ymin, zmin, center, ymax, zmax})
    // std::array<double, 6>({center, ymin, zmin, xmax, ymax, zmax})

    // icart == 1
    // std::array<double, 6>({xmin, ymin, zmin, xmax, center, zmax})
    // std::array<double, 6>({xmin, center, zmin, xmax, ymax, zmax})

    // icart == 2
    // std::array<double, 6>({xmin, ymin, zmin, xmax, ymax, center})
    // std::array<double, 6>({xmin, ymin, center, xmax, ymax, zmax})

    subboxdims[0][3 + icart] = center;

    subboxdims[1][icart] = center;

    // update subgridpoints

    for (int32_t g = 0; g < npoints; g++)
    {
        int32_t box_id = (coords[icart][g] < center) ? 0 : 1;

        auto count = subnumpoints[box_id];

        subgridpoints[box_id].data(0)[count] = xcoords[g];

        subgridpoints[box_id].data(1)[count] = ycoords[g];

        subgridpoints[box_id].data(2)[count] = zcoords[g];

        subgridpoints[box_id].data(3)[count] = weights[g];

        ++subnumpoints[box_id];
    }

    std::list<CGridBox> newboxes;

    for (int32_t box_id = 0; box_id < static_cast<int32_t>(subboxdims.size()); box_id++)
    {
        auto count = subnumpoints[box_id];

        if (count > 0) newboxes.push_back(CGridBox(subboxdims[box_id], subgridpoints[box_id].slice(0, count)));
    }

    return newboxes;
}

double
CGridPartitioner::findCenter(const double* ptr, const int32_t num) const
{
    double center = 0.0;

    for (int32_t i = 0; i < num; i++)
    {
        center += ptr[i];
    }

    center /= num;

    return center;
}

double
CGridPartitioner::findMedian(const double* ptr, const int32_t num) const
{
    // Note: ptr and num are not checked
    // Assuming ptr is a valid pointer and num > 0

    // use approximate median for large number (2^20 or ~1 million)

    if (num > 1048576) return estimateMedianP2(ptr, num);

    // sort data

    std::vector<double> vec(ptr, ptr + num);

    std::sort(vec.begin(), vec.end());

    // divide at threshold for special cases

    if (1024 * 2 <= num && num <= 1024 * 3) return 0.5 * (vec[1024 - 1] + vec[1024]);

    if (1024 * 4 <= num && num <= 1024 * 5) return 0.5 * (vec[1024 - 1] + vec[1024]);

    // divide evenly otherwise

    if (num % 2 == 1)
    {
        return vec[num / 2];
    }
    else
    {
        return 0.5 * (vec[num / 2 - 1] + vec[num / 2]);
    }
}

double
CGridPartitioner::estimateMedianP2(const double* ptr, const int32_t num) const
{
    // Note: ptr and num are not checked
    // Assuming ptr is a valid pointer and num is reasonably large (e.g. > 100)

    // Reference: The P2 algorithm (https://doi.org/10.1145/4372.4378)

    const double p = 0.5;

    std::vector<double> q(ptr, ptr + 5);

    std::sort(q.begin(), q.end());

    std::vector<int32_t> n({0, 1, 2, 3, 4});

    std::vector<double> nprime({0.0, 2.0 * p, 4.0 * p, 2.0 + 2.0 * p, 4.0});

    std::vector<double> dnprime({0.0, p / 2.0, p, (1.0 + p) / 2.0, 1.0});

    for (int32_t j = 5; j < num; j++)
    {
        int32_t k = -1;

        if (ptr[j] < q[0])
        {
            q[0] = ptr[j];

            k = 0;
        }
        else if (ptr[j] < q[1])
        {
            k = 0;
        }
        else if (ptr[j] < q[2])
        {
            k = 1;
        }
        else if (ptr[j] < q[3])
        {
            k = 2;
        }
        else if (ptr[j] < q[4])
        {
            k = 3;
        }
        else
        {
            q[4] = ptr[j];

            k = 3;
        }

        for (int32_t i = k + 1; i < 5; i++)
        {
            ++n[i];
        }

        for (int32_t i = 0; i < 5; i++)
        {
            nprime[i] += dnprime[i];
        }

        for (int32_t i = 1; i < 4; i++)
        {
            double di = nprime[i] - n[i];

            if (((di >= 1.0) && (n[i + 1] - n[i] > 1.0)) || ((di <= -1.0) && (n[i - 1] - n[i] < -1.0)))
            {
                int32_t di_sign = 0;

                if (di > 0.0) di_sign = 1;

                else if (di < 0.0) di_sign = -1;

                double qiprime = q[i] + static_cast<double>(di_sign) / (n[i + 1] - n[i - 1]) * (

                        (n[i] - n[i - 1] + di_sign) * (q[i + 1] - q[i]) / (n[i + 1] - n[i]) +

                        (n[i + 1] - n[i] - di_sign) * (q[i] - q[i - 1]) / (n[i] - n[i - 1]));

                if ((q[i - 1] < qiprime) && (qiprime < q[i + 1]))
                {
                    q[i] = qiprime;
                }
                else
                {
                    q[i] += static_cast<double>(di_sign) * (q[i + di_sign] - q[i]) / (n[i + di_sign] - n[i]);

                }

                n[i] += di_sign;
            }
        }
    }

    return q[2];
}

int32_t
CGridPartitioner::getNumberOfBoxes() const
{
    return static_cast<int32_t>(_boxes.size());
}

std::string
CGridPartitioner::getGridInformation() const
{
    std::stringstream ss;

    ss << "Number of grid boxes: " << getNumberOfBoxes() << "\n";

    int32_t boxind = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        auto boxdim = box->getBoxDimension();

        auto npoints = box->getNumberOfGridPoints();

        ss << "  Grid box " << boxind << ", number of points: " << npoints << ", xyz: ";

        ss << "(" << boxdim[0] << ", " << boxdim[1] << ", " << boxdim[2] << "), ";

        ss << "(" << boxdim[3] << ", " << boxdim[4] << ", " << boxdim[5] << ")\n";

        ++boxind;
    }

    return ss.str();
}

std::string
CGridPartitioner::getGridStatistics() const
{
    std::stringstream ss;

    ss << "Threshold for number of points per box: " << _numberOfPointsThreshold << "\n";

    ss << "Total number of boxes: " << getNumberOfBoxes() << "\n";

    int32_t npoints_max = -1, npoints_min = -1, npoints_sum = 0;

    int32_t nboxes_dense = 0, nboxes_sparse = 0;

    int32_t npoints_bin = 100;

    std::vector<int32_t> nboxes_count(_numberOfPointsThreshold / npoints_bin + 1, 0);

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        auto npoints = box->getNumberOfGridPoints();

        if ((npoints_max == -1) || (npoints > npoints_max)) npoints_max = npoints;

        if ((npoints_min == -1) || (npoints < npoints_min)) npoints_min = npoints;

        npoints_sum += npoints;

        if (npoints > 0.8 * _numberOfPointsThreshold) ++nboxes_dense;

        if (npoints < 0.1 * _numberOfPointsThreshold) ++nboxes_sparse;

        ++nboxes_count[npoints / npoints_bin];
    }

    ss << "Maximum number of points per box: " << npoints_max << "\n";

    ss << "Minimum number of points per box: " << npoints_min << "\n";

    ss << "Total number of points in all boxes: " << npoints_sum << "\n";

    ss << "Average number of points per box: " << npoints_sum / getNumberOfBoxes() << "\n";

    ss << nboxes_dense << " boxes with npoints > 80\% of threshold " << "\n";

    ss << nboxes_sparse << " boxes with npoints < 10\% of threshold " << "\n";

    ss << "-----------------\n";

    ss << " NPoints  NBoxes\n";

    ss << "-----------------\n";

    for (auto i = 0; i < static_cast<int32_t>(nboxes_count.size()); i++)
    {
        ss << std::right << std::setfill(' ') << std::setw(4) << i * npoints_bin << "-";

        ss << std::left << std::setfill(' ') << std::setw(4) << (i + 1) * npoints_bin - 1 << ":";

        ss << std::right << std::setfill(' ') << std::setw(6) << nboxes_count[i] << "\n";
    }

    ss << "-----------------\n";

    return ss.str();
}

int32_t
CGridPartitioner::getTotalNumberOfGridPoints() const
{
    int32_t totalnpoints = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        totalnpoints += box->getNumberOfGridPoints();
    }

    return totalnpoints;
}

CMemBlock2D<double>
CGridPartitioner::getAllGridPoints() const
{
    CMemBlock2D<double> allgridpoints(getTotalNumberOfGridPoints(), 4);

    int32_t count = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        auto npoints = box->getNumberOfGridPoints();

        auto xcoords = box->getCoordinatesX();

        auto ycoords = box->getCoordinatesY();

        auto zcoords = box->getCoordinatesZ();

        auto weights = box->getWeights();

        std::memcpy(allgridpoints.data(0) + count, xcoords, npoints * sizeof(double));

        std::memcpy(allgridpoints.data(1) + count, ycoords, npoints * sizeof(double));

        std::memcpy(allgridpoints.data(2) + count, zcoords, npoints * sizeof(double));

        std::memcpy(allgridpoints.data(3) + count, weights, npoints * sizeof(double));

        count += npoints;
    }

    return allgridpoints;
}

CMemBlock<int32_t>
CGridPartitioner::getGridPointCounts() const
{
    CMemBlock<int32_t> counts(getNumberOfBoxes());

    int32_t index = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        counts.data(index)[0] = box->getNumberOfGridPoints();

        ++index;
    }

    return counts;
}

CMemBlock<int32_t>
CGridPartitioner::getGridPointDisplacements() const
{
    CMemBlock<int32_t> displacements(getNumberOfBoxes());

    int32_t index = 0, count = 0;

    for (auto box = _boxes.cbegin(); box != _boxes.cend(); ++box)
    {
        displacements.data(index)[0] = count;

        ++index;

        count += box->getNumberOfGridPoints();
    }

    return displacements;
}
