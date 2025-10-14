"""Dash application entrypoint."""

from __future__ import annotations

from pathlib import Path

import dash
import dash_bootstrap_components as dbc

from . import layout
from .config import SETTINGS

ASSETS_PATH = Path(__file__).resolve().parents[3] / "dashboard" / "assets"


def create_app() -> dash.Dash:
    """Create and configure the Dash app."""

    app = dash.Dash(
        __name__,
        use_pages=True,
        suppress_callback_exceptions=True,
        external_stylesheets=[dbc.themes.LUX],
        title="gvmagDB Analytics",
        assets_folder=str(ASSETS_PATH) if ASSETS_PATH.exists() else None,
    )
    app.layout = layout.build_layout()
    return app


app = create_app()
server = app.server


def main() -> None:
    """CLI entrypoint used by pixi."""

    app.run_server(host=SETTINGS.host, port=SETTINGS.port, debug=SETTINGS.debug)


if __name__ == "__main__":
    main()
