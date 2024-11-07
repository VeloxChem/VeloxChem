def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "solvers: tests for response solvers",
    )
    config.addinivalue_line(
        "markers",
        "timeconsuming: time-consuming tests that are not included in CI",
    )
