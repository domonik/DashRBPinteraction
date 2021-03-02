import pandas as pd
import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import pickle
import dash  # (version 1.12.0) pip install dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_table
from dash.dependencies import Input, Output, State
import numpy as np
import plotly.io as pio

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

with open("./goa_intersection_table_no_prop_all.pickl", "rb") as f:
    intersect_table = pickle.load(f)

s = intersect_table.intersection.str.len().sort_values(ascending=False).index
godict = dict()


intersect_table = intersect_table.reindex(s)
intersect_table = intersect_table.iloc[:10]

for x in range(len(intersect_table)):
    data = intersect_table.iloc[x]["protein1_goterms"]
    godict = {**godict, **data}
intersect_table["intersection"] = intersect_table["intersection"].apply(lambda x: ";".join(list(godict[y] for y in x)))
icolumns = ["protein1", "protein2", "intersection"]


def get_data_table():
    p = dash_table.DataTable(id="table",
                             columns=[{"name": i, "id": i} for i in intersect_table[icolumns].columns],
                             data=intersect_table[icolumns].to_dict("records"),
                             style_cell={"color": "black"}
                             )
    return p


app.layout = html.Div([

    html.Div(html.H1("Protein Interaction Dasboard", style={'text-align': 'center'}), className="page-header"),
    html.Div(get_data_table()
             , className="databox")

], id="wrapper"
)

if __name__ == '__main__':
    app.run_server(debug=True, port=8080, host="0.0.0.0")
