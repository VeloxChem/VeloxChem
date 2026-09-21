"""The threadpoolctl branch of node.py, which only runs where it is installed.

The laptop has no threadpoolctl and the node does, so the branch which matters on
the node is the one which cannot be exercised where the code is written. A stub
stands in for the library and checks the three things that branch is for:

  * the pool is **read back** rather than assumed, so a library which clamps the
    limit is recorded as what it gave and not as what it was asked for;
  * a library serving fewer threads than the rank was given is fatal, because it
    does not refuse them -- it warns once per thread and dies later somewhere else;
  * one serving enough of them is not.

    python test_node.py          # prints each check
    pytest benchmarks/scripts/test_node.py

It is not collected by the repository's pytest run, whose testpaths is `tests`.
"""
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

LIBS = []          # what the stubbed controller reports
LIMITS = []        # what limit() was asked for, so the call can be checked


class _Selected:
    def __init__(self, prefix):
        self.prefix = prefix

    def limit(self, limits):
        LIMITS.append((self.prefix, limits))
        # NOTE: the clamp this exists to catch -- a library built with a ceiling of
        # its own takes the limit and gives back less.
        for lib in LIBS:
            if Path(lib["filepath"]).name.startswith(self.prefix):
                lib["num_threads"] = min(limits, lib["ceiling"])


class _Controller:
    def info(self):
        return [dict(lib) for lib in LIBS]

    def select(self, prefix):
        return _Selected(prefix)


def _install_stub():
    """Puts the stub where node.py's lazy import will find it."""
    stub = types.ModuleType("threadpoolctl")
    stub.ThreadpoolController = _Controller
    sys.modules["threadpoolctl"] = stub


_install_stub()

import node  # noqa: E402  -- the stub has to be in place first


def _set(*libs):
    LIBS[:] = [dict(lib) for lib in libs]
    LIMITS.clear()


def _openblas(path, threads, ceiling):
    return {"filepath": path, "internal_api": "openblas",
            "num_threads": threads, "ceiling": ceiling}


def test_pool_is_sized_to_the_ranks_threads():
    """numpy's library is the only one limited, and to this rank's share."""
    _set(_openblas("/x/libscipy_openblas.so", 1, 64),
         {"filepath": "/x/libomp.dylib", "internal_api": "openmp",
          "num_threads": 14, "ceiling": 256})

    assert node.widen_for_numpy(32) == 32
    # the driver's own library is left alone: it binds per region from OMP_PLACES
    # and wants every thread it was given
    assert LIMITS == [("libscipy_openblas", 32)]


def test_pool_is_read_back_and_not_assumed():
    """A library with a ceiling of its own gives back less than it was asked."""
    _set(_openblas("/x/libscipy_openblas.so", 1, 16))

    assert node.widen_for_numpy(32) == 16


def test_a_library_below_the_ranks_threads_is_fatal():
    """The module's OpenBLAS is MAX_THREADS=48; a rank of 128 must not start."""
    _set(_openblas("/mod/libopenblas.so", 48, 48))

    fatal, warnings = node.check_blas_ceilings(128)

    assert len(fatal) == 1
    assert "libopenblas.so" in fatal[0]
    assert warnings == []


def test_a_library_which_serves_the_threads_is_not():
    """The rebuilt one, 256 threads, against a rank given 128."""
    _set(_openblas("/new/libopenblas.so", 256, 256))

    assert node.check_blas_ceilings(128) == ([], [])


def test_the_library_list_reaches_the_record():
    """Path and pool both, since which file answered is half of what a number is."""
    _set(_openblas("/new/libopenblas.so", 256, 256))

    libraries = node._blas_libraries()

    assert len(libraries) == 1
    assert libraries[0]["path"] == "/new/libopenblas.so"
    assert libraries[0]["threads"] == 256


def main():
    failures = []

    for name, case in sorted(globals().items()):
        if not name.startswith("test_") or not callable(case):
            continue
        try:
            case()
        except AssertionError as exc:
            failures.append(f"{name}: {exc}")
            print(f"  FAIL {name}")
        else:
            print(f"  ok   {name}")

    print("\nall checks passed" if not failures else "\nFAILURES:")
    for failure in failures:
        print("  " + failure)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
