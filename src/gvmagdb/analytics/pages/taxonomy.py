"""Taxonomy explorer page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import Input, Output, dcc, html

from gvmagdb.analytics import data_access

dash.register_page(__name__, path="/taxonomy", name="Taxonomy Explorer", order=1)

LEVEL_OPTIONS = [
    {"label": label.title(), "value": value}
    for value, label in [
        ("domain", "domain"),
        ("phylum", "phylum"),
        ("class", "class"),
        ("order", "order"),
        ("family", "family"),
        ("genus", "genus"),
        ("species", "species"),
    ]
]

SOURCE_OPTIONS = [
    {"label": "GVClass taxonomy", "value": "gvclass"},
    {"label": "Phylogenomic taxonomy", "value": "phylo"},
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
                            id="taxonomy-level",
                            options=LEVEL_OPTIONS,
                            value="phylum",
                            class_name="mb-2",
                        ),
                        md=4,
                    ),
                    dbc.Col(
                        dbc.Select(
                            id="taxonomy-source",
                            options=SOURCE_OPTIONS,
                            value="gvclass",
                            class_name="mb-2",
                        ),
                        md=4,
                    ),
                ],
                class_name="align-items-end",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="taxonomy-distribution-graph"))),
                        class_name="mb-4",
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Top labels", className="card-title"),
                                    dcc.Loading(
                                        html.Div(id="taxonomy-table"),
                                        type="dot",
                                    ),
                                ]
                            )
                        )
                    )
                ]
            ),
        ],
        className="taxonomy-page",
    )


@dash.callback(
    Output("taxonomy-distribution-graph", "figure"),
    Output("taxonomy-table", "children"),
    Input("taxonomy-level", "value"),
    Input("taxonomy-source", "value"),
)
def update_taxonomy_graph(level: str, source: str):
    df = data_access.fetch_taxonomy_distribution(level, source)

    if df.empty:
        return _empty_figure("No taxonomy data available"), html.P("No data to display.")

    title_prefix = "GVClass" if source == "gvclass" else "Phylogenomic"
    fig = px.bar(
        df.sort_values("genomes"),
        x="genomes",
        y="label",
        orientation="h",
        title=f"{title_prefix} {level.title()} distribution",
        labels={"label": level.title(), "genomes": "Genomes"},
    )

    table = dbc.Table.from_dataframe(
        df.rename(columns={"label": level.title(), "genomes": "Genomes"}),
        striped=True,
        bordered=False,
        hover=True,
        size="sm",
    )

    return fig, table
