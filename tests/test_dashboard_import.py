"""Smoke tests for the analytics dashboard package."""

from gvmagdb.analytics.app import create_app


def test_create_app_smoke():
    app = create_app()
    assert app.title == "gvmagDB Analytics"
