"""Shared layout components for the analytics dashboard."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
from dash import html


def build_navigation() -> dbc.Nav:
    """Create navigation links from the registered Dash pages."""

    nav_links: list[dbc.NavLink] = []
    for page in sorted(dash.page_registry.values(), key=lambda p: p.get("order", 999)):
        nav_links.append(
            dbc.NavLink(
                page["name"],
                href=page["path"],
                active="exact",
                class_name="nav-link px-3",
            )
        )

    return dbc.Nav(nav_links, pills=True, class_name="flex-column flex-lg-row gap-2")


def build_layout() -> dbc.Container:
    """Return the root layout for the Dash app."""

    navigation = build_navigation()

    header = dbc.Navbar(
        dbc.Container(
            [
                html.Div(
                    [
                        html.H1("gvmagDB Analytics", className="mb-0"),
                        html.Small(
                            "Interactive exploration of GVClass, taxonomy, environments, and annotations",
                            className="text-muted",
                        ),
                    ],
                    className="d-flex flex-column",
                ),
            ]
        ),
        class_name="mb-4 shadow-sm",
        color="dark",
        dark=True,
    )

    return dbc.Container(
        [
            header,
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Navigation", className="card-title"),
                                    navigation,
                                ]
                            ),
                            class_name="mb-4",
                        ),
                        lg=3,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        html.Div(dash.page_container, id="page-container"),
                        lg=9,
                    ),
                ],
                class_name="g-4",
            ),
            html.Footer(
                [
                    html.Small(
                        [
                            "Built with Plotly Dash · ",
                            html.A(
                                "Documentation",
                                href="https://github.com/NeLLi-team/gvmagDB/tree/master/docs",
                                target="_blank",
                            ),
                        ],
                        className="text-muted",
                    )
                ],
                className="mt-5 pb-4",
            ),
        ],
        fluid=True,
        className="pt-4",
    )
