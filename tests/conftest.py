try:
    import matplotlib
except ImportError:
    matplotlib = None
else:
    # Use a non-interactive backend so plot tests run reliably in headless environments.
    matplotlib.use('Agg')


def pytest_make_parametrize_id(config, val, argname):
    if argname == 'ri' and isinstance(val, bool):
        return 'RI' if val else '4C'

    return None


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "solvers: tests for response solvers",
    )
    config.addinivalue_line(
        "markers",
        "timeconsuming: tests taking >= 3 s with 16 threads; excluded from CI",
    )
