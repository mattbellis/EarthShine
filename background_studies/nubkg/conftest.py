"""Local pytest configuration so the package is self-contained.

Registers the `literature` marker without needing a pytest.ini next to the
tests, which matters because the notebooks invoke pytest via --pyargs from
whatever directory they happen to be opened in.

Run only the external cross-checks with:   pytest -m literature
Skip them with:                            pytest -m "not literature"
"""


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "literature: comparison against a published measurement")
