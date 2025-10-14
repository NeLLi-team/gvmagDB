"""Functional annotation explorer page."""

from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
import plotly.express as px
from dash import Input, Output, dcc, html

from .. import data_access

dash.register_page(__name__, path="/annotations", name="Annotations", order=4)

FIELD_OPTIONS = [
    {"label": "COG category", "value": "emapper_COG_category"},
    {"label": "KEGG pathway", "value": "emapper_KEGG_Pathway"},
    {"label": "PFAM domains", "value": "emapper_PFAMs"},
    {"label": "Descriptions", "value": "emapper_Description"},
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
                            id="annotation-field",
                            options=FIELD_OPTIONS,
                            value="emapper_COG_category",
                        ),
                        md=4,
                        class_name="mb-3",
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Card(dbc.CardBody(dcc.Graph(id="annotation-heatmap"))),
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
                                    html.H5("Top annotation hits", className="card-title"),
                                    dcc.Loading(html.Div(id="annotation-table"), type="dot"),
                                ]
                            )
                        )
                    )
                ]
            ),
        ],
        className="annotations-page",
    )


@dash.callback(
    Output("annotation-heatmap", "figure"),
    Output("annotation-table", "children"),
    Input("annotation-field", "value"),
)
def update_annotations(field: str):
    df = data_access.fetch_annotation_matrix(field)

    if df.empty:
        return _empty_figure("No annotation data available"), html.P("No data to display.")

    pivot = df.pivot_table(
        index="gvclass_order",
        columns="annotation_value",
        values="sequences",
        fill_value=0,
    )

    fig = px.imshow(
        pivot,
        color_continuous_scale="Viridis",
        aspect="auto",
        labels={"color": "Sequences"},
        title="Annotation enrichment across GVClass orders",
    )

    table_df = df.copy()
    table_df.rename(
        columns={
            "gvclass_order": "GVClass order",
            "annotation_value": "Annotation",
            "sequences": "Sequences",
        },
        inplace=True,
    )

    table = dbc.Table.from_dataframe(
        table_df,
        striped=True,
        bordered=False,
        hover=True,
        size="sm",
    )

    return fig, table
