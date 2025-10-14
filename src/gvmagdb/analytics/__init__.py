"""
Interactive analytics dashboards for the gvmagDB project.

This module exposes a Dash application (see ``gvmagdb.analytics.app``) and a small
set of helpers for querying the DuckDB catalog in a dashboard-friendly manner.
"""

from . import (
    config,  # re-export for convenience
    data_access,
)

__all__ = ["config", "data_access"]
