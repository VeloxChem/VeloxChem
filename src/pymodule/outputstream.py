import sys
import time as tm
from .veloxchemlib import assert_msg_critical


class OutputStream:

    def __init__(self, filename=None, width=122):
        self.width = width
        self.buffer_lines = []

        # filename is...
        #   None or empty:  stream is None
        #   sys.stdout:     stream is sys.stdout
        #   string:         stream is file handle

        if not filename:
            self.stream = None
            self.state = False

        elif filename is sys.stdout:
            self.stream = sys.stdout
            self.state = True

        else:
            try:
                self.stream = open(filename, 'w')
            except OSError:
                errio = "OutputStream: cannot open output file %s" % filename
                assert_msg_critical(False, errio)
            self.state = True

    def __del__(self):
        if self.state:
            if self.buffer_lines:
                self.flush()
            if self.stream is not sys.stdout:
                self.stream.close()

    def get_state(self):
        return self.state

    def flush(self):
        if self.state:
            for line in self.buffer_lines:
                self.stream.write(line)
            self.stream.flush()
            del self.buffer_lines[:]

    @staticmethod
    def header(line, width):
        length = len(line)
        left = (width - length) // 2
        right = width - length - left
        return ' ' * left + line + ' ' * right

    @staticmethod
    def title(line, width):
        length = len(line)
        left = (width - 2 - length) // 2
        right = width - 2 - length - left
        return '!' + ' ' * left + line + ' ' * right + '!'

    @staticmethod
    def info(line, width):
        length = len(line)
        left = 9
        right = width - length - left
        return '* Info * ' + line + ' ' * right

    @staticmethod
    def tsep(width):
        return '!' + '=' * (width - 2) + '!'

    def print_line(self, line):
        if not self.state:
            return
        self.buffer_lines.append(line + '\n')

    def print_blank(self):
        if not self.state:
            return
        self.buffer_lines.append(' ' * self.width + '\n')

    def print_header(self, line):
        if not self.state:
            return
        self.buffer_lines.append(self.header(line, self.width) + '\n')

    def print_title(self, line):
        if not self.state:
            return
        self.buffer_lines.append(self.title(line, self.width) + '\n')

    def print_info(self, line):
        if not self.state:
            return
        self.buffer_lines.append(self.info(line, self.width) + '\n')

    def print_separator(self):
        if not self.state:
            return
        self.buffer_lines.append(self.tsep(self.width) + '\n')

    def print_block(self, block_lines):
        """
        Prints multiline text block to output stream.
        """
        if not self.state:
            return

        lines = block_lines.split('\n')
        for line in lines:
            self.print_header(line)

    def print_start_header(self, num_nodes):
        """
        Prints start header to output stream.
        """
        if not self.state:
            return tm.time()

        start_time = tm.time()

        self.print_separator()
        self.print_title("")
        self.print_title("VELOXCHEM (V.0.0 2019)")
        self.print_title("AN ELECTRONIC STRUCTURE CODE")
        self.print_title("")
        self.print_title("Copyright (c) 2019 VeloxChem developers.")
        self.print_title("All rights reserved.")
        self.print_separator()
        exec_str = "VeloxChem execution started"
        if num_nodes > 1:
            exec_str += " on " + str(num_nodes) + " compute nodes"
        exec_str += " at " + tm.asctime(tm.localtime(start_time)) + "."
        self.print_title(exec_str)
        self.print_separator()
        self.print_blank()

        return start_time

    def print_finish_header(self, start_time):
        """
        Prints start header to output stream.
        """
        if not self.state:
            return

        end_time = tm.time()

        self.print_separator()
        exec_str = "VeloxChem execution completed at "
        exec_str += tm.asctime(tm.localtime(end_time)) + "."
        self.print_title(exec_str)
        self.print_separator()
        exec_str = "Total execution time is "
        exec_str += "{:.2f}".format(end_time - start_time) + " sec."
        self.print_title(exec_str)
        self.print_separator()
        self.flush()
