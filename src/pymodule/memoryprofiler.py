import numpy as np
import time as tm
import psutil
import sys


class MemoryProfiler:
    """
    Impements memory profiler.

    :param remark:
        Descriptive text about the task that is being monitored.

    Instance variable
        - remark: Descriptive text about the task that is being monitored.
        - time_start: The starting time.
        - avail_mem_start: The starting available memory.
        - memory_usage: The list containing information about memory
          usage.
    """

    def __init__(self, remark=''):
        """
        Initializes memory profiler.
        """

        self.remark = remark
        self.time_start = tm.time()
        self.avail_mem_start = psutil.virtual_memory().available
        self.memory_usage = [(0, 0, self.remark + ' start')]

    def check_memory_system(self, remark=''):
        """
        Checks total memory usage.

        :param remark:
            Descriptive text about the point of checking.
        """

        avail_mem = psutil.virtual_memory().available
        used_mem = max(0, self.avail_mem_start - avail_mem)
        self.memory_usage.append(
            (tm.time() - self.time_start, used_mem, remark))

    def print_memory_usage(self, ostream):
        """
        Prints total memory usage.

        :param ostream:
            The output stream.
        """

        mem_str = '*** {} memory usage ***'.format(self.remark)
        ostream.print_header(mem_str.ljust(92))
        ostream.print_blank()

        mem_str = 'Note: Memory usage is calculated based on the amount of'
        ostream.print_header(mem_str.ljust(92))
        mem_str = 'available memory on the master node. For informational'
        ostream.print_header(mem_str.ljust(92))
        mem_str = 'purposes only.'
        ostream.print_header(mem_str.ljust(92))
        ostream.print_blank()

        mem_str = '{:20s}'.format('Elapsed Time')
        mem_str += ' {:20s}'.format('Memory Usage')
        mem_str += ' {:s}'.format('Remark')
        ostream.print_header(mem_str.ljust(92))
        ostream.print_header(('-' * (len(mem_str) + 10)).ljust(92))

        for dt, mem, remark in self.memory_usage:
            mem_str = '{:.2f} sec'.format(dt).ljust(20)
            mem_str += ' {:20s}'.format(self.memory_string(mem))
            mem_str += ' {:s}'.format(remark)
            ostream.print_header(mem_str.ljust(92))
        ostream.print_blank()

        max_used_mem = max([mem for dt, mem, remark in self.memory_usage])
        mem_str = 'Maximum memory usage: {:s} out of {:s} ({:.1f}%)'.format(
            self.memory_string(max_used_mem),
            self.memory_string(self.avail_mem_start),
            max_used_mem / self.avail_mem_start * 100.0,
        )
        ostream.print_header(mem_str.ljust(92))
        ostream.print_blank()
        ostream.flush()

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

        return self.memory_string(self.comp_memory_object(obj))

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
            mem_detail.append((key, self.memory_string(mem_obj) + flag))

        return self.memory_string(mem_usage), mem_detail

    def get_available_memory(self):
        """
        Gets available memory as text string.

        :return:
            The amount of available memory as text string.
        """

        return self.memory_string(psutil.virtual_memory().available)

    def memory_string(self, memsize):
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
