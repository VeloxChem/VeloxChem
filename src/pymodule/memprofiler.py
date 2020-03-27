from sys import getsizeof
from itertools import chain
from psutil import virtual_memory


def object_size(obj, flag=None):
    """
    Returns the approximate memory footprint of an object or list of objects
    and all of its contents along with the units of the value. Automatically
    finds the contents of the following builtin containers and their
    subclasses:  tuple, list, dict, set and frozenset.

    :param obj:
        The object.
    :param flag:
        The flag for printing memory size in bytes.
    """

    def dict_handler(d):
        return chain.from_iterable(d.items())

    all_handlers = {
        tuple: iter,
        list: iter,
        dict: dict_handler,
        set: iter,
        frozenset: iter,
    }
    seen = set()
    default_size = getsizeof(0)

    def sizeof(o):
        if id(o) in seen:
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    if flag == 'bytes':
        return sizeof(obj)
    else:
        return mem_string(sizeof(obj))


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
