import pandas as pd
import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import pickle
import dash  # (version 1.12.0) pip install dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import numpy as np
import plotly.io as pio

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

STEPSIZE = 5
PLOTLY_COLORS = px.colors.qualitative.Light24
del PLOTLY_COLORS[12]

pio.templates["plotly_white"].update({"layout": {
    # e.g. you want to change the background to transparent
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': ' rgba(0,0,0,0)',
    "font": dict(color="white")
}})

# load interaction window data file
with open("pickled_full_interaction_table400.pckl", "rb") as f:
    df = pickle.load(f)
df = df.sort_values(by=['protein1', "protein2"])

# set x-axis
x_axis = np.array(list(range(len(df["interactions"].iloc[0])))) * STEPSIZE
selection = sorted(df["protein1"].unique())
opts = []
cellline_selection = sorted(df["protein1_cellline"].unique())
cellline_opts = [{"label": entry, "value": entry} for entry in cellline_selection]
region_opts = [{"label": entry, "value": entry} for entry in sorted(df["region"].unique())]

colordict = dict()
ct = 0
for entry in selection:
    opt = {"label": entry, "value": entry}
    opts.append(opt)
    if entry not in colordict:
        colordict[entry] = dict()
    for cellline in cellline_selection:
        colordict[entry][cellline] = PLOTLY_COLORS[ct]
        ct += 1
        if ct >= len(PLOTLY_COLORS):
            ct = 0

# load GO Analysis table
with open("./goa_mastertable_all.pickl", "rb") as f:
    goa_table = pickle.load(f)
    goa_table = goa_table[goa_table["p_fdr_bh"] < 0.05]
    goa_table = goa_table[goa_table["NS"] == "BP"]
    goa_table = goa_table[goa_table["enrichment"] == "e"]

# load GO intersection table
with open("./goa_intersection_table_all.pickl", "rb") as f:
    intersect_table = pickle.load(f)

s = intersect_table.intersection.str.len().sort_values().index

intersect_table = intersect_table.reindex(s)

prot_dict = dict()
for el in goa_table.protein.unique():
    prot_dict[el] = {}
    for cel in goa_table.cellline.unique():
        prot_dict[el][cel] = set()
for x in range(len(goa_table)):
    data = goa_table.iloc[x]
    protein = data["protein"]
    go_term = data["# GO"]
    id = data["id"]
    ns = data["NS"]
    cellline = data["cellline"]
    prot_dict[protein][cellline].add(go_term)


with open("id2numbindingsites.pckl", "rb") as f:
    id2numbs = pickle.load(f)
# values = []
# name = []
# for protein in prot_dict:
#     for cellline in prot_dict[protein]:
#
#         for slct_protein in prot_dict:
#             for slct_cellline in prot_dict[slct_protein]:
#                 if protein != slct_protein and cellline != slct_cellline:
#                     comp = prot_dict[protein][cellline]
#                     sel = prot_dict[slct_protein][slct_cellline]
#                     wlen = len(sel) + len(comp)
#                     val2 = len(sel.intersection(comp))
#                     if wlen != 0:
#                         val = len(sel.intersection(comp)) / wlen
#
#                     else:
#                         val = 0
#
#                     bla = "{0}-{1}:{2}-{3}".format(protein, cellline, slct_protein, slct_cellline)
#                     values.append((val2, bla))
#                     name.append(bla)
# arr = np.array(values, dtype=[("x", int), ("y", "S30")])
# arr.sort(order="x")
# x = 0


dropdown_style = {'width': "100%", "margin": "auto", }
tr_style = {'width': "50%", "margin": "auto", "padding": "5px"}


def get_table():
    filter_c = {
        "input_min": "minimum interactions at max nt",
        "input_max": "maximum interactions at max nt",
        "filter_value": "minimum slope within a window set at filter_size",
        "filter_size": "window used to calculate the minimum slope. 0 is default. Value multiplied by " + str(STEPSIZE),
        "cutoff": "filters ignore the first n * " + str(STEPSIZE) + " nucleotides"
    }
    table = html.Table(
        [html.Tr([html.Th(html.Div(col, className="table-header"), title=filter_c[col]) for col in filter_c])] +
        [html.Tr(
            [html.Td(html.Div(
                dcc.Input(id="{}".format(inp), type="number", placeholder="{}".format(inp), style={"width": "100%"})),
                className="filter_table_interactive") for inp in filter_c] +
            [html.Div(dbc.Button("Switch mode", id="switch_button", color="light", n_clicks=0,
                                 style={"white-space": "nowrap", "margin": "auto"}),
                      className="filter_table_interactive")]
        )], id="filter-table", style={"margin": "auto", "width": "80%"}
    )
    return table


def dropdown_menues():
    dropdowns = [
        dcc.Dropdown(id="slct_protein",
                     options=opts,
                     multi=False,
                     value="UPF1",
                     style=dropdown_style),
        dcc.Dropdown(id="slct_cellline",
                     options=cellline_opts,
                     value=cellline_opts[0]["value"],
                     style=dropdown_style
                     ),
        dcc.Dropdown(id="slct_region",
                     options=region_opts,
                     value=region_opts[0]["value"],
                     style=dropdown_style
                     )

    ]
    menu = html.Table(
        html.Tr(
            [html.Td(html.Div(col, style={"width": "100%"}), className="selection-table_interactive") for col in
             dropdowns]
        ),
        id="selection-table", style={"width": "80%", "margin": "auto"}

    )

    return menu


app.layout = html.Div([

    html.Div(html.H1("Protein Interaction Dasboard", style={'text-align': 'center'}), className="page-header"),
    html.Div(dropdown_menues(), className="databox"),

    html.Div([html.H2(id='output_container', children=[], style={"text-align": "center"}),
              dcc.Graph(id='plotly_graph', figure={"layout": {"height": 700}}),
              get_table()], className="plotly-graph"),
    html.Div(
        [
            html.H2(id='goa_output_container', children=[], style={"text-align": "center"}),
            html.H3(id='goa_output_container2', children=[], style={"text-align": "center"}),
            dcc.Graph(id='goa-graph', figure={"layout": {"width": "100%"}}),
            html.Div(
                dbc.Button("Switch mode", id="goa_switch_button", color="light", n_clicks=0,
                           style={"white-space": "nowrap", "margin": "auto"}),
                style={"display": "flex", "align-items": "center", "justify-content": "center"})
        ]
        , className="databox")

], id="wrapper"
)


@app.callback(
    [Output(component_id='output_container', component_property='children'),
     Output(component_id='plotly_graph', component_property='figure')],
    [Input(component_id='slct_protein', component_property='value'),
     Input(component_id="slct_cellline", component_property="value"),
     Input(component_id="slct_region", component_property="value"),
     Input(component_id="input_min", component_property="value"),
     Input(component_id="input_max", component_property="value"),
     Input(component_id="filter_value", component_property="value"),
     Input(component_id="filter_size", component_property="value"),
     Input(component_id="cutoff", component_property="value"),
     Input(component_id="switch_button", component_property="n_clicks"),
     ]
)
def update_graph(option_slctd, cellline_slctd, region_slctd, inputmin, inputmax, filter_value, filter_size, cutoff,
                 switch):
    container = "{} Interactions".format(option_slctd)
    viewmode = switch % 2
    if viewmode == 0:
        y_axis_label = "# of interactions / max # of interactions"
    else:
        y_axis_label = "# of interactions"
    if filter_value is None:
        filter_value = 0
    if inputmin is None:
        inputmin = 0
    if inputmax is None:
        inputmax = np.inf
    if filter_size is None or filter_size <= 0:
        filter_size = 1
    filter_size = int(filter_size)
    if cutoff is None or cutoff < 0:
        cutoff = 0
    dff = df.copy()
    dff = dff[dff["region"] == region_slctd]

    dff = dff[dff["protein1"] == option_slctd]
    dff = dff[dff["protein1_cellline"] == cellline_slctd]

    fig = go.Figure()
    to_plot = []
    if len(dff) == 0:
        container = "{} {} seems to be missing".format(option_slctd, cellline_slctd)
        fig.update_yaxes(range=[0, 1])
        fig.update_xaxes(range=[0, 400])
    for i in range(len(dff)):
        data = dff.iloc[i]
        interactions = data["interactions"]
        numbs = id2numbs[data["protein2_id"]]
        if viewmode == 0 and interactions[-1] != 0:
            y = interactions / interactions[-1]
        else:
            y = interactions / numbs
        max_i = data["interactions"][-1]
        name = data["protein2"]
        cellline = data["protein2_cellline"]
        name_ext = name + "-" + cellline + ":" + str(max_i)

        diff = np.diff(y)
        diff_window = get_sum(diff, filter_size)
        max_diff = max(diff_window[cutoff:])
        print(numbs)
        print(name)

        if inputmax >= max_i >= inputmin:
            if max_diff > filter_value:
                to_plot.append(go.Scatter(x=x_axis, y=y, mode="lines", name=name_ext, hovertext=name_ext,
                                          line={"width": 4, "color": colordict[name][cellline]}))
            else:
                fig.add_trace(go.Scatter(x=x_axis, y=y, mode="lines", name=name_ext, hovertext=name_ext,
                                         line={"color": "grey"}))
    for element in to_plot:
        fig.add_trace(element)
    fig.layout.template = "plotly_white"
    fig.update_layout(hoverlabel=dict(namelength=0),
                      xaxis_title="distance [nt]",
                      yaxis_title=y_axis_label,
                      )
    return container, fig


@app.callback(
    [
        Output(component_id='goa_output_container', component_property='children'),
        Output(component_id='goa_output_container2', component_property='children'),
        Output(component_id='goa-graph', component_property='figure'),
    ],
    [Input(component_id='slct_protein', component_property='value'),
     Input(component_id="slct_cellline", component_property="value"),
     Input(component_id="goa_switch_button", component_property="n_clicks"),

     ],
    [State("goa-graph", "relayoutData")]
)
def update_goanalysis(slct_protein, slct_cellline, nclicks, relayout_data):
    container = "{} GO-terms Intersection Enrichment".format(slct_protein)

    heatline = []
    names = []
    heatline2 = []
    if slct_protein in prot_dict:
        sel = prot_dict[slct_protein][slct_cellline]

        for protein in prot_dict:
            for cellline in prot_dict[protein]:
                if protein != slct_protein and cellline != slct_cellline:
                    comp = prot_dict[protein][cellline]
                    wlen = len(sel) + len(comp)
                    val2 = len(sel.intersection(comp))

                    if wlen != 0:
                        val = len(sel.intersection(comp)) / wlen
                    else:
                        val = 0
                    heatline.append(val)
                    heatline2.append(val2)
                    names.append(protein + "-" + cellline)
        if nclicks % 2 == 1:
            heatline = heatline2
            container2 = "Absolute Values"
        else:
            container2 = "Jaccard Index"

        heatline = [heatline]
        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(z=heatline, y=[slct_protein + "-" + slct_cellline], x=names, zmin=0, zauto=False,
                       colorscale="Hot", colorbar=go.heatmap.ColorBar(ticklen=100, tickwidth=100, tickangle=90),
                       zmax=np.max(heatline), ))



    else:
        heatline = [[0, 2, 1]]
        fig = px.imshow(heatline)
    # fig.update_traces(showscale=False)
    fig.layout.template = "plotly_white"

    if relayout_data:
        if 'xaxis.range[0]' in relayout_data:
            fig['layout']['xaxis']['range'] = [
                relayout_data['xaxis.range[0]'],
                relayout_data['xaxis.range[1]']
            ]
        if 'yaxis.range[0]' in relayout_data:
            fig['layout']['yaxis']['range'] = [
                relayout_data['yaxis.range[0]'],
                relayout_data['yaxis.range[1]']
            ]

    return container, container2, fig


def get_sum(arr, dif):
    result_array = np.zeros(len(arr))
    for x in range(len(arr)):
        window = (x, x + dif)
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
    app.run_server(debug=True, port=8080, host="0.0.0.0")
