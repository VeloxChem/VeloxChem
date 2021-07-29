def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "solvers: tests for response solvers (deselect with '-m \"not solvers\"')",
    )
