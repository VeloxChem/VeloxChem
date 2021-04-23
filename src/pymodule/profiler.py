#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
import psutil
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
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation using
        - start_time: The starting time.
        - start_avail_mem: The starting available memory.
        - timing_list: The list containing information about timing.
        - memory_usage: The list containing information about memory usage.
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
        self.memory_profiling = False
        self.memory_tracing = False

        self.start_time = None
        self.start_avail_mem = None

        self.timing_list = None
        self.memory_usage = None
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
        if 'memory_profiling' in settings:
            self.memory_profiling = settings['memory_profiling']
        if 'memory_tracing' in settings:
            self.memory_tracing = settings['memory_tracing']

        self.start_time = tm.time()
        self.start_avail_mem = psutil.virtual_memory().available

        self.timing_list = []
        self.memory_usage = []

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
        self.print_memory_usage(ostream, scf_flag)
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
            sortby = 'cumulative'
            ps = pstats.Stats(self.pr, stream=s).sort_stats(sortby)
            ps.print_stats(30)

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

    def start_timer(self, iteration, label):
        """
        Starts the timer for an iteration and a label.

        :param iteration:
            The iteration.
        :param label:
            The label.
        """

        if self.timing:
            if iteration >= len(self.timing_list):
                num = iteration - len(self.timing_list) + 1
                self.timing_list += [{}] * num
            self.timing_list[iteration][label] = tm.time()

    def stop_timer(self, iteration, label):
        """
        Stops the timer for an iteration and a label.

        :param iteration:
            The iteration.
        :param label:
            The label.
        """

        if self.timing:
            t0 = self.timing_list[iteration][label]
            self.timing_list[iteration][label] = tm.time() - t0

    def update_timer(self, iteration, timing_dict):
        """
        Updates the timer with a dictionary of timing information.

        :param iteration:
            The iteration.
        :param timing_dict:
            The dictionary containing timing information.
        """

        if self.timing:
            for key, val in timing_dict.items():
                self.timing_list[iteration][key] = val

    def print_timing(self, ostream, scf_flag=False):
        """
        Prints timing.

        :param ostream:
            The output stream.
        :param scf:
            The flag for SCF.
        """

        if self.timing:
            width = 92

            valstr = 'Timing (in sec)'
            ostream.print_header(valstr.ljust(width))
            ostream.print_header(('-' * len(valstr)).ljust(width))

            keys = list(self.timing_list[0].keys())

            valstr = '{:<18s}'.format('')
            for key in keys:
                valstr += ' {:>12s}'.format(key)
            ostream.print_header(valstr.ljust(width))

            for i, d in enumerate(self.timing_list):
                if scf_flag and i == 0:
                    continue
                index = i if scf_flag else i + 1
                valstr = '{:<18s}'.format('Iteration {:d}'.format(index))
                for key in keys:
                    if key in d:
                        valstr += ' {:12.2f}'.format(d[key])
                    else:
                        valstr += ' {:12s}'.format('')
                ostream.print_header(valstr.ljust(width))

            ostream.print_blank()
            ostream.print_blank()

    def check_memory_usage(self, remark=''):
        """
        Checks memory usage.

        :param remark:
            Descriptive text about the point of checking.
        """

        if self.memory_profiling:
            avail_mem = psutil.virtual_memory().available
            used_mem = max(0, self.start_avail_mem - avail_mem)
            self.memory_usage.append(
                (tm.time() - self.start_time, used_mem, remark))

    def print_memory_usage(self, ostream, scf_flag=False):
        """
        Prints memory usage.

        :param ostream:
            The output stream.
        :param scf_flag:
            The flag for SCF.
        """

        if self.memory_profiling:
            mem_str = 'Estimated memory usage'
            ostream.print_header(mem_str.ljust(92))
            ostream.print_header(('-' * len(mem_str)).ljust(92))

            mem_str = '{:20s}'.format('Elapsed Time')
            mem_str += ' {:20s}'.format('Memory Usage')
            mem_str += ' {:s}'.format('Remark')
            ostream.print_header(mem_str.ljust(92))

            for dt, mem, remark in self.memory_usage:
                if scf_flag and remark.lower() == 'iteration 0':
                    continue
                mem_str = '{:.2f} sec'.format(dt).ljust(20)
                mem_str += ' {:20s}'.format(self.memory_to_string(mem))
                mem_str += ' {:s}'.format(remark)
                ostream.print_header(mem_str.ljust(92))

            ostream.print_blank()
            ostream.print_blank()

    def comp_memory_object(self, obj, counted_ids=None):
        """
        Computes the memory usage of an object recursively.

        :param obj:
            The object.
        :param counted_ids:
            The list of id's of counted objects.
        :return:
            The memory usage in bytes.
        """

        memsize = 0
        if counted_ids is None:
            counted_ids = []
        if id(obj) not in counted_ids:
            memsize += sys.getsizeof(obj)
            counted_ids.append(id(obj))

        obj_is_dict = isinstance(obj, dict)
        if isinstance(obj, (dict, list, tuple, set, frozenset)):
            for x in obj:
                memsize += self.comp_memory_object(x, counted_ids)
                if obj_is_dict:
                    memsize += self.comp_memory_object(obj[x], counted_ids)

        return memsize

    def get_memory_object(self, obj):
        """
        Gets memory usage of an object as text string.

        :param obj:
            The object.
        :return:
            The amount of memory usage of the object as text string.
        """

        return self.memory_to_string(self.comp_memory_object(obj))

    def print_memory_subspace(self, d, ostream):
        """
        Prints memory usage in subspace solver.

        :param d:
            The dictionary containing the data used in subspace.
        :param ostream:
            The output stream.
        """

        mem_usage, mem_detail = self.get_memory_dictionary(d)
        mem_avail = self.get_available_memory()

        ostream.print_info(
            '{:s} of memory used for subspace procedure on the master node'.
            format(mem_usage))
        if self.memory_profiling:
            for m in mem_detail:
                ostream.print_info('  {:<15s} {:s}'.format(*m))
        ostream.print_info(
            '{:s} of memory available for the solver on the master node'.format(
                mem_avail))
        ostream.print_blank()

    def get_memory_dictionary(self, d):
        """
        Gets memory usage of a dictionary.

        :param d:
            The dictionary.
        :return:
            Memory usage of the dictionary and the list of memory usage of each
            item.
        """

        mem_usage = 0.0
        mem_detail = []

        for key, obj in d.items():
            flag = ''
            if isinstance(obj, np.ndarray) and not obj.flags.owndata:
                flag = '  (not accurate)'
            mem_obj = self.comp_memory_object(obj)
            mem_usage += mem_obj
            mem_detail.append((key, self.memory_to_string(mem_obj) + flag))

        return self.memory_to_string(mem_usage), mem_detail

    def get_available_memory(self):
        """
        Gets available memory as text string.

        :return:
            The amount of available memory as text string.
        """

        return self.memory_to_string(psutil.virtual_memory().available)

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
