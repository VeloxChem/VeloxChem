from .VeloxChemLib import OutputStream
import time


def _print_start_header(self, num_nodes):
    """
    Prints start header to output stream.
    """
    # get start time

    start_time = time.time()

    # print start header

    self.put_separator()
    self.put_title("")
    self.put_title("VELOX CHEM MP (V.0.0 2018)")
    self.put_title("AN ELECTRONIC STRUCTURE CODE FOR NANOSCALE")
    self.put_title("")
    self.put_title("Copyright (c) 2018 Velox Chem MP developers.")
    self.put_title("All rights reserved.")
    self.put_separator()
    exec_str = "VeloxChem MP execution started"
    if num_nodes > 1:
        exec_str += " on " + str(num_nodes) + " compute nodes"
    exec_str += " at " + time.asctime(time.localtime(start_time)) + "."
    self.put_title(exec_str)
    self.put_separator()
    self.new_line()

    # return start time

    return start_time


def _print_finish_header(self, start_time):
    """
    Prints start header to output stream.
    """
    # get end time

    end_time = time.time()

    # print finish header

    self.put_separator()
    exec_str = "VeloxChem MP execution completed at "
    exec_str += time.asctime(time.localtime(end_time)) + "."
    self.put_title(exec_str)
    self.put_separator()
    exec_str = "Total execution time is "
    exec_str += str(int(end_time - start_time)) + " sec."
    self.put_title(exec_str)
    self.put_separator()
    self.flush()


OutputStream.print_start_header = _print_start_header
OutputStream.print_finish_header = _print_finish_header
