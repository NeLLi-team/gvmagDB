"""Genome statistics dashboard page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import Input, Output, dcc, html

from gvmagdb.analytics import data_access

dash.register_page(__name__, path="/genomes", name="Genome Statistics", order=3)

TAXONOMY_OPTIONS = [
    {"label": "Domain", "value": "taxonomy_domain"},
    {"label": "Phylum", "value": "taxonomy_phylum"},
    {"label": "Class", "value": "taxonomy_class"},
    {"label": "Order", "value": "taxonomy_order"},
    {"label": "Family", "value": "taxonomy_family"},
    {"label": "Genus", "value": "taxonomy_genus"},
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
                        [
                            html.Label("Group genomes by"),
                            dbc.Select(
                                id="genome-group-field",
                                options=TAXONOMY_OPTIONS,
                                value="taxonomy_order",
                            ),
                        ],
                        md=4,
                        class_name="mb-3",
                    ),
                    dbc.Col(
                        [
                            html.Label("Colour points by"),
                            dbc.Select(
                                id="genome-color-field",
                                options=[
                                    {"label": "Ecosystem", "value": "ecosystem"},
                                    {"label": "GVClass Order", "value": "taxonomy_order"},
                                    {"label": "GVClass Family", "value": "taxonomy_family"},
                                ],
                                value="taxonomy_order",
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
                        dbc.Card(dbc.CardBody(dcc.Graph(id="completeness-box"))),
                        lg=4,
                        class_name="mb-4",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="contamination-box"))),
                        lg=6,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Summary by group", className="card-title"),
                                    dcc.Loading(html.Div(id="genome-summary-table"), type="dot"),
                                ]
                            )
                        ),
                        lg=6,
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
    Output("completeness-box", "figure"),
    Output("contamination-box", "figure"),
    Output("genome-summary-table", "children"),
    Input("genome-group-field", "value"),
    Input("genome-color-field", "value"),
)
def update_genome_charts(group_field: str, color_field: str):
    df = data_access.fetch_genome_statistics()

    if df.empty:
        empty = _empty_figure("No genome statistics available")
        return empty, empty, empty, empty, html.P("No data to display.")

    group_field = group_field if group_field in df.columns else "taxonomy_order"
    colour = color_field if color_field in df.columns else group_field

    hist = px.histogram(
        df,
        x="genome_length_mb",
        color=colour,
        nbins=40,
        title="Genome length distribution",
        labels={"genome_length_mb": "Genome length (Mb)"},
    )

    scatter = px.scatter(
        df,
        x="gc_percent",
        y="coding_percent",
        color=colour,
        hover_data={
            "dataset_id": True,
            "genome_length_mb": ":.2f",
            "gene_count": True,
            "completeness": ":.1f",
            "contamination": ":.2f",
        },
        title="GC% vs coding density",
        labels={"gc_percent": "GC %", "coding_percent": "Coding %"},
    )

    completeness_box = px.box(
        df,
        x=group_field,
        y="completeness",
        color=group_field,
        title="Completeness by group",
        labels={"completeness": "Completeness (%)", group_field: "Group"},
        points="outliers",
    )
    completeness_box.update_layout(xaxis_title="Group")

    contamination_box = px.box(
        df,
        x=group_field,
        y="contamination",
        color=group_field,
        title="Estimated contamination by group",
        labels={"contamination": "Contamination (%)", group_field: "Group"},
        points=False,
    )
    contamination_box.update_layout(xaxis_title="Group")

    summary = (
        df.groupby(group_field)
        .agg(
            genomes=("dataset_id", "nunique"),
            median_length_mb=("genome_length_mb", "median"),
            median_gc=("gc_percent", "median"),
            median_completeness=("completeness", "median"),
            median_contamination=("contamination", "median"),
        )
        .reset_index()
        .sort_values("genomes", ascending=False)
        .head(10)
    )
    summary = summary.round(
        {
            "median_length_mb": 3,
            "median_gc": 2,
            "median_completeness": 2,
            "median_contamination": 2,
        }
    )
    summary = summary.rename(
        columns={
            group_field: "Group",
            "genomes": "Genomes",
            "median_length_mb": "Median length (Mb)",
            "median_gc": "Median GC %",
            "median_completeness": "Median completeness",
            "median_contamination": "Median contamination",
        }
    )

    summary_table = dbc.Table.from_dataframe(
        summary,
        striped=True,
        bordered=False,
        hover=True,
        size="sm",
    )

    return hist, scatter, completeness_box, contamination_box, summary_table
