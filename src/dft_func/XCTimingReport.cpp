//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "XCTimingReport.hpp"

#include <algorithm>
#include <cstdio>
#include <map>

namespace xcprof {

auto
report(const std::string             &what,
       const CMultiTimer             &timer,
       const std::vector<CMultiTimer> &omptimers,
       const size_t                   boxes) -> void
{
    const auto nthreads = std::max<size_t>(omptimers.size(), 1);

    // NOTE: the labels in the order they were first started, which is the order the
    // phases run in, rather than sorted by name or by cost. A profile which reads in
    // the order of the work is easier to hold against the source.

    std::vector<std::string> order;

    std::map<std::string, double> totals;

    for (const auto &omptimer : omptimers)
    {
        for (const auto &[label, seconds] : omptimer.getTimings())
        {
            if (totals.find(label) == totals.end()) order.push_back(label);

            totals[label] += seconds;
        }
    }

    double total = 0.0;

    for (const auto &[label, seconds] : timer.getTimings())
    {
        if (label == "Total timing") total = seconds;
    }

    std::printf("XC %s on %zu threads, %zu boxes, %.3f s\n", what.c_str(), nthreads, boxes, total);

    double accounted = 0.0;

    for (const auto &label : order)
    {
        const auto shared = totals[label] / static_cast<double>(nthreads);

        accounted += shared;

        std::printf("XC   %-24s %8.3f s %6.1f %%\n", label.c_str(), shared,
                    (total > 0.0) ? 100.0 * shared / total : 0.0);
    }

    // NOTE: what the threads did not account for is the serial part -- everything
    // outside the parallel region -- plus whatever the threads spent waiting for
    // each other. The two are not separated here and the second is the interesting
    // one, so a large remainder is a question and not an answer.

    std::printf("XC   %-24s %8.3f s %6.1f %%\n", "serial and imbalance", total - accounted,
                (total > 0.0) ? 100.0 * (total - accounted) / total : 0.0);

    std::fflush(stdout);
}

auto
report_blocks(const size_t naos,
              const size_t boxes,
              const size_t points,
              const size_t ao_sum,
              const size_t ao_max,
              const size_t point_ao,
              const size_t point_ao_sq) -> void
{
    if ((boxes == 0) || (points == 0) || (naos == 0)) return;

    const auto plain = static_cast<double>(ao_sum) / static_cast<double>(boxes);

    const auto weighted = static_cast<double>(point_ao) / static_cast<double>(points);

    // NOTE: against what a calculation with no screening at all would do, which is
    // every point against every function twice over. One over this is how many times
    // more work a dense quadrature would be, and it is the number which has to stay
    // put as the molecule grows for the quadrature to be linear in it.
    const auto dense = static_cast<double>(points) * static_cast<double>(naos) * static_cast<double>(naos);

    const auto share = static_cast<double>(point_ao_sq) / dense;

    std::printf("XC   blocks: %zu of %zu functions kept on average, %.1f weighted by points, %zu at most\n",
                static_cast<size_t>(plain), naos, weighted, ao_max);

    std::printf("XC   blocks: %.1f %% of the functions, and %.2f %% of the work a dense quadrature would do\n",
                100.0 * weighted / static_cast<double>(naos), 100.0 * share);

    std::fflush(stdout);
}

}  // namespace xcprof
