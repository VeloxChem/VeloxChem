from sys import getsizeof
from itertools import chain
from psutil import virtual_memory

def object_size(obj):
    """
    Returns the approximate memory footprint of an object or list of objects
    and all of its contents along with the units of the value.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, dict, set and frozenset.

    :param obj:
        The object or list of objects.
    """

    units = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB']

    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
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

    size = sizeof(obj)

    unit_indx = 0
    unit = units[unit_indx]
    while size >= 1000:
        size = size / 1000
        unit_indx += 1
        unit = units[unit_indx]

    return size, unit


def avail_mem():
    """
    Returns the available system memory along with the units of the value.
    """

    units = ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB']

    virtmem = virtual_memory().available

    unit_indx = 0
    unit = units[unit_indx]
    while virtmem >= 1000:
        virtmem = virtmem / 1000
        unit_indx += 1
        unit = units[unit_indx]

    return virtmem, unit

