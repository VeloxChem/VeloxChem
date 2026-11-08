import gc

import pytest

try:
    import matplotlib
except ImportError:
    matplotlib = None
else:
    # Use a non-interactive backend so plot tests run reliably in headless environments.
    matplotlib.use('Agg')


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "solvers: tests for response solvers",
    )
    config.addinivalue_line(
        "markers",
        "timeconsuming: time-consuming tests that are not included in CI",
    )


@pytest.fixture(autouse=True)
def _release_autocdam_output_streams(request):
    """Delete unreachable AutoCDAM output streams while pytest capture is active.

    A default ``OutputStream`` holds pytest's capture file and closes it on
    deletion unless it is the current ``sys.stdout``, so collect such streams
    before capture is suspended.
    """
    yield
    if request.module.__name__.startswith("test_autocdam_"):
        gc.collect()


@pytest.fixture
def ostream(tmp_path_factory):
    """A file-backed output stream for objects that outlive the test body.

    The log lives outside ``tmp_path`` so tests can assert no files are written.
    """
    from veloxchem.outputstream import OutputStream

    stream = OutputStream(str(tmp_path_factory.mktemp("ostream") / "autocdam.log"))
    yield stream
    stream.close()
