from .veloxchemlib import OutputStream
import time as tm


def _print_start_header(self, num_nodes):
    """
    Prints start header to output stream.
    """

    start_time = tm.time()

    self.print_separator()
    self.print_title("")
    self.print_title("VELOX CHEM MP (V.0.0 2018)")
    self.print_title("AN ELECTRONIC STRUCTURE CODE FOR NANOSCALE")
    self.print_title("")
    self.print_title("Copyright (c) 2018 Velox Chem MP developers.")
    self.print_title("All rights reserved.")
    self.print_separator()
    exec_str = "VeloxChem MP execution started"
    if num_nodes > 1:
        exec_str += " on " + str(num_nodes) + " compute nodes"
    exec_str += " at " + tm.asctime(tm.localtime(start_time)) + "."
    self.print_title(exec_str)
    self.print_separator()
    self.print_blank()

    return start_time


def _print_finish_header(self, start_time):
    """
    Prints start header to output stream.
    """

    end_time = tm.time()

    self.print_separator()
    exec_str = "VeloxChem MP execution completed at "
    exec_str += tm.asctime(tm.localtime(end_time)) + "."
    self.print_title(exec_str)
    self.print_separator()
    exec_str = "Total execution time is "
    exec_str += "{:.2f}".format(end_time - start_time) + " sec."
    self.print_title(exec_str)
    self.print_separator()
    self.flush()


def _print_block(self, block_lines):
    """
    Prints multiline text block to output stream.
    """

    lines = block_lines.split('\n')
    self.print_header(lines.pop(0))
    for line in lines:
        self.print_line(line)


OutputStream.print_start_header = _print_start_header
OutputStream.print_finish_header = _print_finish_header
OutputStream.print_block = _print_block
