import numpy as np
from sys import getsizeof
from psutil import virtual_memory

from .errorhandler import assert_msg_critical


def sizeof_numpy_arrays(objs, flag=None):
    """
    Gets the size of containers of numpy arrays.

    :param objs:
        The list of container objects.
    :param flag:
        The flag for printing memory size in bytes.
    """

    memsize = getsizeof(objs)

    for obj in objs:

        if isinstance(obj, (list, tuple, set)):
            for x in obj:
                assert_msg_critical(
                    isinstance(x, np.ndarray),
                    'sizeof_numpy_nparrays: incorrect type of object ' +
                    str(type(x)))
                memsize += getsizeof(x)

        elif isinstance(obj, dict):
            for key, x in obj.items():
                assert_msg_critical(
                    isinstance(x, np.ndarray),
                    'sizeof_numpy_arrays: incorrect type of object ' +
                    str(type(x)))
                memsize += getsizeof(x)

        else:
            assert_msg_critical(
                isinstance(obj, np.ndarray),
                'sizeof_numpy_arrays: incorrect type of object ' +
                str(type(obj)))
            memsize += getsizeof(obj)

    if flag == 'bytes':
        return memsize
    else:
        return mem_string(memsize)


def avail_mem(flag=None):
    """
    Returns the available system memory.

    :param flag:
        The flag for printing memory size in bytes.
    """

    if flag == 'bytes':
        return virtual_memory().available
    else:
        return mem_string(virtual_memory().available)


def mem_string(memsize):
    """
    Prints memory size in readable format.

    :param memsize:
        The memory size.
    """

    units = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB']

    unit_index = 0
    while memsize >= 1000 and unit_index < len(units) - 1:
        memsize /= 1000
        unit_index += 1

    return '{:.2f} {:s}'.format(float(memsize), units[unit_index])
