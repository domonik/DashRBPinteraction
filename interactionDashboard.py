import pandas as pd
import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import pickle
import dash  # (version 1.12.0) pip install dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import numpy as np
import plotly.io as pio

app = dash.Dash(__name__, )

filter_columns = ["input_min", "input_max", "filter_value", "filter_size", "cutoff"]
PLOTLY_COLORS = px.colors.qualitative.Dark24
pio.templates["plotly_dark_custom"] = pio.templates["plotly_dark"]

templa = pio.templates["plotly_dark_custom"].update({"layout": {
    # e.g. you want to change the background to transparent
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': '#E5ECF6'
}})


def get_table(filter_columns):
    table = html.Table(
        [html.Tr([html.Td(col) for col in filter_columns])] +
        [html.Tr(
            [html.Td(dcc.Input(id="{}".format(inp), type="number", placeholder="{}".format(inp))) for inp in
             filter_columns] +
            [html.Button("Switch mode", id="switch_button", n_clicks=0)]
        )]
    )
    return table


with open("plotly_visualize_DATA_200mid.pickle", "rb") as f:
    df = pickle.load(f)
x_axis = np.array(list(range(len(df["normalized_interactions"].iloc[0])))) * 5
selection = sorted(df["protein1"].unique())
opts = []

colordict = dict()
ct = 0
for entry in selection:
    colordict[entry] = PLOTLY_COLORS[ct]
    opt = {"label": entry, "value": entry}
    opts.append(opt)
    ct += 1
    if ct >= len(PLOTLY_COLORS):
        ct = 0


app.layout = html.Div([

    html.H1("Web Application Dashboards with Dash", style={'text-align': 'center'}),

    dcc.Dropdown(id="slct_year",
                 options=opts,
                 multi=False,
                 value=selection[0],
                 style={'width': "70%"}
                 ),
    html.Div(id='output_container', children=[]),
    html.Br(),
    dcc.Graph(id='plotly_graph', figure={"layout": {"height": 700}}),
    html.Br(),
    get_table(filter_columns),

],
)


@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='plotly_graph', component_property='figure')],
    [Input(component_id='slct_year', component_property='value'),
     Input(component_id="input_min", component_property="value"),
     Input(component_id="input_max", component_property="value"),
     Input(component_id="filter_value", component_property="value"),
     Input(component_id="filter_size", component_property="value"),
     Input(component_id="cutoff", component_property="value"),
     Input(component_id="switch_button", component_property="n_clicks"),
     ]
)
def update_graph(option_slctd, inputmin, inputmax, filter_value, filter_size, cutoff, switch):
    container = "The protein chosen by user was: {}".format(option_slctd)
    viewmode = switch % 2
    if filter_value is None:
        filter_value = 0
    if inputmin is None:
        inputmin = 0
    if inputmax is None:
        inputmax = np.inf
    if filter_size is None or filter_size <= 0:
        filter_size = 1
    if cutoff is None or cutoff < 0:
        cutoff = 0
    dff = df.copy()
    dff = dff[dff["protein1"] == option_slctd]

    fig = go.Figure()
    to_plot = []
    for i in range(len(dff)):
        if viewmode == 0:
            y = dff["normalized_interactions"].iloc[i]
        else:
            y = dff["interactions"].iloc[i]
        num_i = dff["interactions"].iloc[i][-1]
        name = dff["protein2"].iloc[i]
        name_ext = dff["protein2"].iloc[i] + ":" + str(num_i)

        diff = dff["differences"].iloc[i]
        diff_window = get_sum(diff, filter_size)
        max_diff = max(diff_window[cutoff:])
        if inputmax >= num_i >= inputmin:
            if max_diff > filter_value:
                to_plot.append(go.Scatter(x=x_axis, y=y, mode="lines", name=name_ext, hovertext=name_ext,
                                          line={"width": 4, "color": colordict[name]}))
            else:
                fig.add_trace(go.Scatter(x=x_axis, y=y, mode="lines", name=name_ext, hovertext=name_ext,
                                         line={"color": "grey"}))
    for element in to_plot:
        fig.add_trace(element)
    # fig.layout.template = "plotly_dark_custom"
    fig.update_layout(hoverlabel=dict(namelength=0))
    return container, fig


def get_sum(arr, dif):
    result_array = np.zeros(len(arr))
    for x in range(len(arr)):
        window = (x - dif, x + dif + 1)
        if window[0] < 0:
            area = arr[0:window[1]]
        elif window[1] > len(arr):
            area = arr[window[0]:]
        else:
            area = arr[window[0]:window[1]]
        med = np.sum(area)
        result_array[x] = med
    return result_array


if __name__ == '__main__':
    app.run_server(debug=True)
