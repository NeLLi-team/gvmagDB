"""
Dashboard configuration helpers.

Default values can be adjusted via environment variables:

``GVMAGDB_PARQUET_GLOB``   override Parquet glob used for analytics
``GVMAGDB_DASH_HOST``      host interface for Dash server (default: 127.0.0.1)
``GVMAGDB_DASH_PORT``      server port (default: 8050)
``GVMAGDB_DASH_DEBUG``     enable Dash debug mode when set to ``true``
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DashboardSettings:
    """Resolved dashboard configuration."""

    parquet_glob: str
    host: str
    port: int
    debug: bool
    cache_dir: Path
    cache_enabled: bool


def _env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def load_settings() -> DashboardSettings:
    """Load configuration from environment variables."""

    parquet_glob = _env("GVMAGDB_PARQUET_GLOB", "artifacts/parquet/sequences/**/*.parquet")
    host = _env("GVMAGDB_DASH_HOST", "127.0.0.1")
    port = int(_env("GVMAGDB_DASH_PORT", "8050"))
    debug = _env("GVMAGDB_DASH_DEBUG", "false").lower() in {"1", "true", "yes", "on"}
    cache_dir = Path(_env("GVMAGDB_ANALYTICS_CACHE_DIR", "artifacts/analytics")).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_enabled = _env("GVMAGDB_ANALYTICS_CACHE", "true").lower() in {"1", "true", "yes", "on"}

    return DashboardSettings(
        parquet_glob=parquet_glob,
        host=host,
        port=port,
        debug=debug,
        cache_dir=cache_dir,
        cache_enabled=cache_enabled,
    )


SETTINGS = load_settings()
