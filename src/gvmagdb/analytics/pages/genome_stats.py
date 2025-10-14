"""Genome statistics dashboard page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import Input, Output, dcc, html

from .. import data_access

dash.register_page(__name__, path="/genomes", name="Genome Statistics", order=3)


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
                        [
                            html.Label("Colour by taxonomy"),
                            dbc.Select(
                                id="genome-color-field",
                                options=[
                                    {"label": "GVClass taxonomy", "value": "taxonomy_majority"},
                                    {"label": "Ecosystem", "value": "ecosystem"},
                                ],
                                value="taxonomy_majority",
                            ),
                        ],
                        md=4,
                        class_name="mb-3",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="genome-length-hist"))),
                        lg=4,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="gc-vs-coding"))),
                        lg=4,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="completeness-violin"))),
                        lg=4,
                        class_name="mb-4",
                    ),
                ]
            ),
        ],
        className="genome-stats-page",
    )


@dash.callback(
    Output("genome-length-hist", "figure"),
    Output("gc-vs-coding", "figure"),
    Output("completeness-violin", "figure"),
    Input("genome-color-field", "value"),
)
def update_genome_charts(color_field: str):
    df = data_access.fetch_genome_statistics()

    if df.empty:
        empty = _empty_figure("No genome statistics available")
        return empty, empty, empty

    colour = color_field if color_field in df.columns else None
    group_field = colour or "taxonomy_majority"

    hist = px.histogram(
        df,
        x="genome_length",
        color=colour,
        nbins=40,
        title="Genome length distribution",
        labels={"genome_length": "Genome length (bp)"},
    )

    scatter = px.scatter(
        df,
        x="gc_percent",
        y="coding_percent",
        color=colour,
        hover_name="dataset_id",
        title="GC% vs coding density",
        labels={"gc_percent": "GC %", "coding_percent": "Coding %"},
    )

    violin = px.violin(
        df,
        x=group_field,
        y="order_completeness",
        color=group_field,
        title="Completeness by group",
        labels={"order_completeness": "Completeness"},
    )
    violin.update_xaxes(title_text="Group")

    return hist, scatter, violin
