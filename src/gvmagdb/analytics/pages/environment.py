"""Environment and ecosystem visualisation page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import Input, Output, dcc, html

from .. import data_access

dash.register_page(__name__, path="/environment", name="Environment & Ecosystems", order=2)

DIMENSION_OPTIONS = [
    {"label": "Ecosystem", "value": "ecosystem"},
    {"label": "Ecosystem category", "value": "ecosystem_category"},
    {"label": "Ecosystem type", "value": "ecosystem_type"},
    {"label": "Ecosystem subtype", "value": "ecosystem_subtype"},
    {"label": "Habitat", "value": "habitat"},
    {"label": "Source", "value": "source"},
]


def _empty_figure(message: str) -> dict:
    fig = px.scatter()
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        annotations=[
            {
                "text": message,
                "xref": "paper",
                "yref": "paper",
                "showarrow": False,
                "font": {"size": 14},
            }
        ],
        xaxis={"visible": False},
        yaxis={"visible": False},
    )
    return fig


def layout(**kwargs) -> html.Div:
    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Select(
                            id="environment-dimension", options=DIMENSION_OPTIONS, value="ecosystem"
                        ),
                        md=4,
                        class_name="mb-3",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="environment-bar-chart"))),
                        lg=6,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="environment-pie-chart"))),
                        lg=6,
                        class_name="mb-4",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Environment breakdown", className="card-title"),
                                    dcc.Loading(html.Div(id="environment-table"), type="dot"),
                                ]
                            )
                        )
                    )
                ]
            ),
        ],
        className="environment-page",
    )


@dash.callback(
    Output("environment-bar-chart", "figure"),
    Output("environment-pie-chart", "figure"),
    Output("environment-table", "children"),
    Input("environment-dimension", "value"),
)
def update_environment_charts(dimension: str):
    df = data_access.fetch_environment_distribution(dimension)

    if df.empty:
        empty = _empty_figure("No environment data available")
        return empty, empty, html.P("No data to display.")

    bar = px.bar(
        df.sort_values("genomes"),
        x="genomes",
        y="label",
        orientation="h",
        title=f"Genomes per {dimension.replace('_', ' ')}",
        labels={"label": dimension.replace("_", " ").title(), "genomes": "Genomes"},
    )

    pie = px.pie(
        df,
        values="genomes",
        names="label",
        title=f"Share of genomes by {dimension.replace('_', ' ')}",
    )

    table = dbc.Table.from_dataframe(
        df,
        striped=True,
        bordered=False,
        hover=True,
        size="sm",
        columns=[
            {"name": dimension.replace("_", " ").title(), "id": "label"},
            {"name": "Genomes", "id": "genomes"},
            {"name": "Sequences", "id": "sequences"},
        ],
    )

    return bar, pie, table
