"""Overview dashboard page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import dcc, html

from .. import data_access

dash.register_page(__name__, path="/", name="Overview", order=0)


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


def _build_metric_cards(metrics: dict[str, float]) -> dbc.Row:
    card_info = [
        ("Total sequences", f"{metrics['total_sequences']:,.0f}"),
        ("Unique genomes", f"{metrics['unique_genomes']:,.0f}"),
        ("GVClass species", f"{metrics['gvclass_species']:,.0f}"),
        ("ANI clusters", f"{metrics['ani_clusters']:,.0f}"),
        ("Avg GC (NT)", f"{metrics['avg_gc_nt']:.2f}%"),
    ]

    cards = [
        dbc.Col(
            dbc.Card(
                dbc.CardBody(
                    [
                        html.H6(title, className="card-title text-uppercase text-muted mb-2"),
                        html.H3(value, className="card-text"),
                    ]
                ),
                class_name="shadow-sm",
            ),
            xs=12,
            md=6,
            lg=4,
            class_name="mb-3",
        )
        for title, value in card_info
    ]

    return dbc.Row(cards, class_name="gy-3")


def layout(**kwargs) -> html.Div:
    metrics = data_access.fetch_overview_metrics()
    taxonomy_df = data_access.fetch_taxonomy_distribution("phylum").head(10)
    env_df = data_access.fetch_environment_distribution("ecosystem").head(10)

    if not taxonomy_df.empty:
        taxonomy_fig = px.bar(
            taxonomy_df.sort_values("genomes"),
            x="genomes",
            y="label",
            orientation="h",
            title="Top GVClass phyla",
            labels={"label": "Phylum", "genomes": "Genomes"},
        )
    else:
        taxonomy_fig = _empty_figure("No taxonomy data available")

    if not env_df.empty:
        environment_fig = px.bar(
            env_df.sort_values("genomes"),
            x="genomes",
            y="label",
            orientation="h",
            title="Top ecosystems",
            labels={"label": "Ecosystem", "genomes": "Genomes"},
        )
    else:
        environment_fig = _empty_figure("No environment data available")

    genome_stats = data_access.fetch_genome_statistics()
    if not genome_stats.empty:
        genome_fig = px.scatter(
            genome_stats,
            x="genome_length",
            y="gc_percent",
            color="taxonomy_majority",
            title="Genome length vs GC%",
            labels={"genome_length": "Genome length (bp)", "gc_percent": "GC %"},
            hover_name="dataset_id",
        )
    else:
        genome_fig = _empty_figure("No genome statistics available")

    return html.Div(
        [
            _build_metric_cards(metrics),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(figure=taxonomy_fig))),
                        lg=6,
                        class_name="mb-4",
                    ),
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(figure=environment_fig))),
                        lg=6,
                        class_name="mb-4",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(figure=genome_fig))),
                        class_name="mb-4",
                    )
                ]
            ),
        ],
        className="overview-page",
    )
