"""ANI cluster summary page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import dcc, html

from .. import data_access

dash.register_page(__name__, path="/clusters", name="Clusters & Quality", order=5)


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
    df = data_access.fetch_cluster_summary()

    if df.empty:
        chart = _empty_figure("No cluster data available")
        table = html.P("No data to display.")
    else:
        chart = px.bar(
            df.sort_values("genomes"),
            x="genomes",
            y="cluster_id",
            orientation="h",
            title="Largest ANI clusters",
            labels={"cluster_id": "Cluster", "genomes": "Genomes"},
        )

        pretty_df = df.rename(
            columns={
                "cluster_id": "Cluster",
                "genomes": "Genomes",
                "representatives": "Representatives",
                "gvclass_order": "GVClass order",
            }
        )
        table = dbc.Table.from_dataframe(
            pretty_df,
            striped=True,
            bordered=False,
            hover=True,
            size="sm",
        )

    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(figure=chart))),
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
                                    html.H5("Cluster summary", className="card-title"),
                                    table,
                                ]
                            )
                        )
                    )
                ]
            ),
        ],
        className="clusters-page",
    )
