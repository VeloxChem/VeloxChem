#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from os import environ
from pathlib import Path
import sys
import time as tm

from .errorhandler import assert_msg_critical


class OutputStream:
    """
    Implements the output stream.

    :param filename:
        Name of the output file (or sys.stdout).
    :param width:
        Width of the output.

    Instance variables
        - width: Width of the output.
        - buffer_lines: The buffered lines of output.
        - stream: The stream to which output is printed.
        - state: The flag for writing to the stream.
    """

    def __init__(self, stream=sys.stdout, width=122):
        """
        Initializes the output stream.
        """

        self.width = width
        self.buffer_lines = []

        if stream is None:
            self.stream = None
            self.state = False

        elif stream is sys.stdout:
            self.stream = sys.stdout
            self.state = True

        elif (isinstance(stream, str) and stream) or isinstance(stream, Path):
            fname = str(stream)
            try:
                self.stream = open(fname, 'w')
                self.state = True
            except OSError:
                self.state = False
            errio = f'OutputStream: cannot open output file {fname}'
            assert_msg_critical(self.state, errio)

        else:
            self.state = False
            errio = f'OutputStream: invalid argument {stream}'
            assert_msg_critical(self.state, errio)

        self._state_backup = None
        self._mute_level = 0

    def __del__(self):
        """
        Deletes the output stream.
        """

        self.close()

    def close(self):
        """
        Closes the output stream.
        """

        if self.state:
            if self.buffer_lines:
                self.flush()
            if self.stream is not sys.stdout:
                self.stream.close()
            self.state = False

    @property
    def is_muted(self):
        """
        Checks if the output stream is muted.
        """

        return (self._state_backup is not None and self._mute_level > 0)

    def mute(self):
        """
        Mutes the output stream.
        """

        self._mute_level += 1

        if self._state_backup is None:
            self._state_backup = self.state

        self.state = False

    def unmute(self):
        """
        "Unmutes" the output stream. Note that this only restores the original
        state of the output stream.
        """

        if self._mute_level > 0:
            self._mute_level -= 1

        if self._mute_level == 0 and self._state_backup is not None:
            self.state = self._state_backup
            self._state_backup = None

    def get_state(self):
        """
        :return:
            State of the output stream.
        """

        return self.state

    def flush(self):
        """
        Flushes the buffered output to stream.
        """

        if self.state:
            for line in self.buffer_lines:
                self.stream.write(line + '\n')
            self.stream.flush()
            del self.buffer_lines[:]

    @staticmethod
    def header(line, width):
        """
        Gets the header string.

        :param line:
            The line of text.
        :param width:
            Width of the output.

        :return:
            The header string.
        """

        length = len(line)
        left = (width - length) // 2
        right = width - length - left
        return ' ' * left + line + ' ' * right

    @staticmethod
    def title(line, width):
        """
        Gets the title string.

        :param line:
            The line of text.
        :param width:
            Width of the output.

        :return:
            The title string.
        """

        length = len(line)
        left = (width - 2 - length) // 2
        right = width - 2 - length - left
        return '!' + ' ' * left + line + ' ' * right + '!'

    @staticmethod
    def info(line, width):
        """
        Gets the information string.

        :param line:
            The line of text.
        :param width:
            Width of the output.

        :return:
            The information string.
        """

        length = len(line)
        left = 9
        right = width - length - left
        return '* Info * ' + line + ' ' * right

    @staticmethod
    def warning(line, width):
        """
        Gets the warning string.

        :param line:
            The line of text.
        :param width:
            Width of the output.

        :return:
            The warning string.
        """

        length = len(line)
        left = 12
        right = width - length - left
        return '* Warning * ' + line + ' ' * right

    @staticmethod
    def tsep(width):
        """
        Gets the separator string.

        :param width:
            Width of the output.

        :return:
            The separator string.
        """

        return '!' + '=' * (width - 2) + '!'

    def print_line(self, line):
        """
        Prints the line to stream.

        :param line:
            The line of text.
        """

        if not self.state:
            return
        self.buffer_lines.append(line)

    def print_blank(self):
        """
        Prints a blank line to stream.
        """

        if not self.state:
            return
        self.buffer_lines.append(' ' * self.width)

    def print_header(self, line):
        """
        Prints a header line to stream.

        :param line:
            The line of text.
        """

        if not self.state:
            return
        self.buffer_lines.append(self.header(line, self.width))

    def print_title(self, line):
        """
        Prints a title line to stream.

        :param line:
            The line of text.
        """

        if not self.state:
            return
        self.buffer_lines.append(self.title(line, self.width))

    def print_info(self, line):
        """
        Prints an information line to stream.

        :param line:
            The line of text.
        """

        if not self.state:
            return
        self.buffer_lines.append(self.info(line, self.width))

    def print_warning(self, line):
        """
        Prints a warning line to stream.

        :param line:
            The line of text.
        """

        if not self.state:
            return
        self.buffer_lines.append('*' * 11 + ' ' * (self.width - 11))
        self.buffer_lines.append(self.warning(line, self.width))
        self.buffer_lines.append('*' * 11 + ' ' * (self.width - 11))

    def print_separator(self):
        """
        Prints a separator line to stream.
        """

        if not self.state:
            return
        self.buffer_lines.append(self.tsep(self.width))

    def print_block(self, block_lines):
        """
        Prints a block of lines to stream.

        :param block_lines:
            The multiple lines of text.
        """

        if not self.state:
            return

        lines = block_lines.splitlines()
        for line in lines:
            self.print_header(line)

    def print_reference(self, line):
        """
        Prints reference to output stream.

        :param lines:
            The reference line.
        """

        if not self.state:
            return

        cur_line = line
        spaces = ' ' * 9

        while len(cur_line) > 100:
            index = cur_line[:100].strip().rfind(' ')
            self.buffer_lines.append(spaces + cur_line[:index].strip())
            cur_line = cur_line[index:]

        self.buffer_lines.append(spaces + cur_line.strip())

    def print_start_header(self, num_nodes):
        """
        Prints start header to output stream.

        :param num_nodes:
            The number of MPI processes.
        """

        if not self.state:
            return tm.time()

        start_time = tm.time()

        self.print_separator()
        self.print_title('')
        self.print_title('VELOXCHEM')
        self.print_title('AN ELECTRONIC STRUCTURE CODE')
        self.print_title('')
        self.print_title('Copyright (C) 2018-2024 VeloxChem developers.')
        self.print_title('All rights reserved.')
        self.print_separator()
        exec_str = 'VeloxChem execution started'
        if num_nodes > 1:
            exec_str += ' on ' + str(num_nodes) + ' compute nodes'
        exec_str += ' at ' + tm.asctime(tm.localtime(start_time)) + '.'
        self.print_title(exec_str)
        self.print_separator()
        self.print_blank()

        if 'OMP_NUM_THREADS' in environ:
            self.print_info('Using {} OpenMP threads per compute node.'.format(
                environ['OMP_NUM_THREADS']))
            self.print_blank()

        return start_time

    def print_finish_header(self, start_time):
        """
        Prints finish header to output stream.

        :param start_time:
            The start time of the computation.
        """

        if not self.state:
            return

        end_time = tm.time()

        self.print_separator()
        exec_str = 'VeloxChem execution completed at '
        exec_str += tm.asctime(tm.localtime(end_time)) + '.'
        self.print_title(exec_str)
        self.print_separator()
        exec_str = 'Total execution time is '
        exec_str += '{:.2f}'.format(end_time - start_time) + ' sec.'
        self.print_title(exec_str)
        self.print_separator()

        ref_str = 'Rinkevicius, Z.; Li, X.; Vahtras, O.; Ahmadzadeh, K.; '
        ref_str += 'Brand, M.; Ringholm, M.;'
        self.print_title(ref_str)
        ref_str = 'List, N. H.; Scheurer, M.; Scott, M.; Dreuw, A.; Norman, P.'
        self.print_title(ref_str)
        ref_str = 'VeloxChem: A Python-driven Density-functional Theory '
        ref_str += 'Program for Spectroscopy'
        self.print_title(ref_str)
        ref_str = 'Simulations in High-performance Computing Environments.'
        self.print_title(ref_str)
        ref_str = 'WIREs Comput Mol Sci 2020, 10 (5), e1457.'
        self.print_title(ref_str)
        self.print_separator()

        self.flush()
