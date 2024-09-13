#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import time as tm
import sys
import os


class Profiler:
    """
    Impements multifunctional profiler.

    :param settings:
        The dictionary of profiler settings.

    Instance variable
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_tracing: The flag for tracing memory allocation using
        - start_time: The starting time.
        - start_avail_mem: The starting available memory.
        - timing_dict: The dictionary containing information about timing.
        - timing_key: The key for the timing dictionary.
        - pr: The Profile object.
    """

    def __init__(self, settings=None):
        """
        Initializes multifunctional profiler.

        :param settings:
            The dictionary containing profiler settings.
        """

        self.timing = False
        self.profiling = False
        self.memory_tracing = False

        self.start_time = None
        self.start_avail_mem = None

        self.timing_dict = None
        self.timing_key = None

        self.pr = None

        if settings is not None:
            self.begin(settings)

    def begin(self, settings):
        """
        Starts the multifunctional profiler.

        :param settings:
            The dictionary containing profiler settings.
        """

        if 'timing' in settings:
            self.timing = settings['timing']
        if 'profiling' in settings:
            self.profiling = settings['profiling']
        if 'memory_tracing' in settings:
            self.memory_tracing = settings['memory_tracing']

        self.start_time = tm.time()

        self.timing_dict = {}

        if self.profiling:
            import cProfile
            self.pr = cProfile.Profile()
            self.pr.enable()

        if self.memory_tracing:
            import tracemalloc
            tracemalloc.start()

    def end(self, ostream, scf_flag=False):
        """
        Stops the profiler and print the output.

        :param ostream:
            The output stream.
        :param scf_flag:
            The flag for SCF.
        """

        self.print_timing(ostream, scf_flag)
        self.print_profiling_summary(ostream)
        self.print_memory_tracing(ostream)

    def print_profiling_summary(self, ostream):
        """
        Prints profiling summary.

        :param ostream:
            The output stream.
        """

        if self.profiling:
            import pstats
            import io
            self.pr.disable()
            s = io.StringIO()
            sortby = 'time'
            ps = pstats.Stats(self.pr, stream=s).sort_stats(sortby)
            ps.print_stats(50)

            valstr = 'Profiling summary'
            ostream.print_info('   ' + valstr)
            ostream.print_info('   ' + '-' * len(valstr))
            ostream.print_info('')

            for line in s.getvalue().splitlines():
                if os.sep in line[20:]:
                    text = line.split(os.sep, 1)[0]
                    fname, lineno = line.split(os.sep, 1)[1].split(':')
                    if 'numpy' in fname or 'h5py' in fname:
                        fname = os.sep.join(fname.split(os.sep)[-3:])
                    else:
                        fname = os.sep.join(fname.split(os.sep)[-2:])
                    line = '{:s} ...{:s}:{:s}'.format(text, fname, lineno)
                ostream.print_info(line)
            ostream.print_blank()

    def print_memory_tracing(self, ostream):
        """
        Prints memory tracing output.

        :param ostream:
            The output stream.
        """

        if self.memory_tracing:
            import tracemalloc
            mem_curr = tracemalloc.take_snapshot()
            cur_stats = mem_curr.statistics('lineno')
            ostream.print_info('[ Tracemalloc ]')
            for index, stat in enumerate(cur_stats[:15]):
                frame = stat.traceback[0]
                if 'numpy' in frame.filename or 'h5py' in frame.filename:
                    filename = os.sep.join(frame.filename.split(os.sep)[-3:])
                else:
                    filename = os.sep.join(frame.filename.split(os.sep)[-2:])
                text = '#{:<3d} ...{:s}:{:d}:'.format(index + 1, filename,
                                                      frame.lineno)
                ostream.print_info('{:<45s} {:s}'.format(
                    text, self.memory_to_string(stat.size)))
            ostream.print_blank()

    def set_timing_key(self, key):
        """
        Sets the timing key for the timing dictionary.
        """

        if self.timing:
            self.timing_key = key
            self.timing_dict[key] = {}

    def start_timer(self, label):
        """
        Starts the timer for a label.

        :param label:
            The label.
        """

        key = self.timing_key

        if self.timing and (key is not None) and (key in self.timing_dict):
            self.timing_dict[key][label] = tm.time()

    def stop_timer(self, label):
        """
        Stops the timer for a label.

        :param label:
            The label.
        """

        key = self.timing_key

        if self.timing and (key is not None) and (key in self.timing_dict):
            t0 = self.timing_dict[key][label]
            self.timing_dict[key][label] = tm.time() - t0

    def add_timing_info(self, label, dt):
        """
        Adds timing info for a label.

        :param label:
            The label.
        :param dt:
            The timing info to be added.
        """

        key = self.timing_key

        if self.timing and (key is not None) and (key in self.timing_dict):
            if label not in self.timing_dict[key]:
                self.timing_dict[key][label] = 0.0
            self.timing_dict[key][label] += dt

    def print_timing(self, ostream, scf_flag=False):
        """
        Prints timing.

        :param ostream:
            The output stream.
        :param scf:
            The flag for SCF.
        """

        if self.timing:
            if not self.timing_dict:
                return

            width = 92

            valstr = 'Timing (in sec)'
            ostream.print_header(valstr.ljust(width))
            ostream.print_header(('-' * len(valstr)).ljust(width))

            key_0 = list(self.timing_dict.keys())[0]
            labels = list(self.timing_dict[key_0].keys())

            valstr = '{:<18s}'.format('')
            for label in labels:
                valstr += ' {:>12s}'.format(label)
            ostream.print_header(valstr.ljust(width))

            for key, val in self.timing_dict.items():
                #if scf_flag and key.lower() == 'iteration 0':
                #    continue
                valstr = f'{key:<18s}'
                for label in labels:
                    if label in val:
                        valstr += f' {val[label]:12.2f}'
                    else:
                        valstr += ' {:12s}'.format('')
                ostream.print_header(valstr.ljust(width))

            ostream.print_blank()
            ostream.print_blank()

    def memory_to_string(self, memsize):
        """
        Gets memory size as text string.

        :param memsize:
            The memory size.
        :return:
            The text string for memory size.
        """

        units = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB']

        unit_index = 0
        while (memsize >= 1000 and unit_index < len(units) - 1):
            memsize /= 1000
            unit_index += 1

        return '{:.2f} {:s}'.format(float(memsize), units[unit_index])
