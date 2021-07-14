import base64
import datetime
import io
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import numpy as np
import plotly.express as px
import pandas as pd

global dataframe
global pie_curves
pie_curves = {}
global pie_split
pie_split = {}

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                suppress_callback_exceptions=True)

complete_data = pd.ExcelFile('All_Data copy.xls')  # /home/ishangupta/mysite/data/All_Data copy.xls

sheet_names = complete_data.sheet_names  # global

main_df = complete_data.parse(sheet_names[0])

columns = main_df.columns
main_df.columns = ['Mass'] + ['2mm ' + x for x in columns[1:]]

for x in sheet_names[1:]:
    df = complete_data.parse(x)
    columns = df.columns
    df = df.drop('Mass', axis=1)
    df.columns = [x + ' ' + i for i in columns[1:]]
    main_df = main_df.join(df)

del df
columns = main_df.columns[1:]


# print(df.head(3))
# print(energy_values)
def reduce(sheet, threshold):
    transposed = sheet.T
    columns = transposed.columns
    delete = []
    for x in columns:
        if transposed[x][0] > 1000 or sum(transposed[x][1:]) < threshold:
            delete.append(x)
    transposed = transposed.drop(delete, axis=1)
    return transposed.T


# educed_df = reduce(main_df, 48)
# get_discrete works takes in 5 inputs and returns a spin off of the main df that has the clusters specific to the
# molecule ratios / constitution being looked for
def get_discrete(mass_1, number_1, error, mass_2=0, number_2=0):
    def process(x, mass_number, bounds, start=0):
        if (x - start) % mass_number > (mass_number - bounds):
            return ((x - start) % mass_number) - mass_number
        else:
            return (x - start) % mass_number

    mass_start_1 = mass_1 + 1
    mass_end_1 = number_1 * mass_1 + 1
    mass_start_2 = mass_end_1 + mass_2
    mass_end_2 = mass_start_2 + (number_2 - 1) * mass_2
    temp_df = main_df[(main_df['Mass'] < (mass_end_2 + error)) & (main_df['Mass'] > (mass_start_1 - error))]
    array_of_masses = [mass_start_1 + (mass_1 * i) for i in range(number_1)] + [mass_start_2 + (mass_2 * i) for i in
                                                                                range(number_2)]
    # print(array_of_masses)
    yo = [process(x, mass_1, error, start=mass_start_1) if x < (mass_end_1 + error) else process(x, mass_2, error,
                                                                                                 start=mass_start_2) for
          x in temp_df['Mass']]

    # print(temp_df['Mass'])
    # print(yo)
    # print(len(temp_df))
    # print(len(yo))
    yo = [x if error >= x >= -error else 0 for x in yo]
    a = []
    counter = 1
    for x in range(len(yo) - 1):
        if yo[x] != 0:
            a.append(counter)
            if yo[x] != 0 and yo[x + 1] == 0:
                counter = counter + 1
        else:
            a.append(0)
    a.append(counter)
    # print(a)
    # print(temp_df['Mass'])
    temp_df = temp_df.drop(['Mass'], axis=1)
    temp_df['Mass'] = a
    temp_df = temp_df.groupby('Mass').max()
    temp_df = temp_df.drop(labels=[0])
    # print(temp_df.columns)
    temp_df["Peaks_Mass"] = np.ones(len(temp_df)) * np.array(array_of_masses)
    # print(temp_df)
    return temp_df


app.layout = html.Div(style={'backgroundColor': '#FFFFEE'}, children=[
    html.Div(
        className='app-header',
        children=[
            html.H1('PHYSICAL CHEMISTRY BERKELEY LAB - MASS SPEC ANALYSIS',
                    style={'textAlign': 'center', 'color': '#026A7A'}),
            html.Br(),
            html.H3('Module 1: Mass Spec Graphs (COUNT VS M/Z)',
                    style={'textAlign': 'center', 'color': '#026A7A'}),
            html.Br(),
            html.P('This graph plots ToF Mass Spec graphs for all energies and distance values'
                   ' you can select multiple energy and distance pairs to plot superimposed Mass'
                   ' specs. Additionally, you can click on a set of points consecutively to create'
                   ' a trace of the peaks you want to focus on. Use the "Undo Trace" button to'
                   ' delete and entire trace and use the "Undo Last" to delete the last added'
                   ' point.',
                   style={'textAlign': 'center', 'color': '#026A7A'}),
            html.Div(children=[
                html.Br(),
                dcc.Dropdown(
                    id='energies',
                    options=[{'label': i, 'value': i} for i in columns],
                    value=[columns[0]],
                    multi=True
                ),
            ], style={'width': '100%', 'display': 'inline-block'}
            ),
            html.Br(),
            html.Br(),
            html.Div(children=[
                html.Button('Undo Trace', id='undo', n_clicks=0,
                            style={'marginLeft': '150px', 'marginRight': '150px'}),
                html.Button('Undo Last', id='undo-last', n_clicks=0,
                            style={'marginRight': '150px'}),
                dcc.Store(id='memory-storage')
            ]
            ),
            html.Div(children=[
                html.Br(),
                dcc.Graph(id='checklist-graph'),
                html.Br(),
            ], style={'padding': '10 10', 'width': '100%', 'display': 'inline-block'}
            ),
        ]),
    html.Div([
        html.Br(),
        html.H3('Module 2: Mass Specs With Select Cluster Composition',
                style={'textAlign': 'center', 'color': '#026A7A'}),
        html.Br(),
        html.P('The following Graph allows you to select a certain molecular composition for the '
               'clusters and plot the mass specs for that specific arrangement at all energy and distance pairs.'
               ' For example: if you want to focus on a cluster with 2 molecules of Ethanol and 40 '
               'water molecules with peaks having a spread of around 1 around the m/z ratios, you should input '
               '46, 2, 0.5, 18, 40. Additionally, you can click on any cluster peaks and get the PIE '
               'curve for that specific molecular composition at all distances',
               style={'textAlign': 'center', 'color': '#026A7A'}),
        html.Br(),
        dcc.Input(id="mol_mass1", type="number", placeholder="molecular mass",
                  style={'marginRight': '100px',
                         'marginLeft': '100px'}),
        dcc.Input(id="mol_amount1", type="number", placeholder="# of molecules", debounce=True,
                  style={'marginRight': '100px'}),
        dcc.Input(id="error", type="number", placeholder="error", debounce=True,
                  style={'marginRight': '100px'}),
        dcc.Input(id="mol_mass2", type="number", placeholder="molecular mass 2",
                  style={'marginRight': '100px'}),
        dcc.Input(id="mol_amount2", type="number", placeholder="# of molecules 2", debounce=True,
                  style={'marginRight': '100px'})
        # might need to add a dropdown for energy to shorten the stored datasets
    ], style={}),
    html.Br(),
    html.Div(children=[
        dcc.Dropdown(
            id='modular-cluster',
            options=[{'label': i, 'value': i} for i in columns],
            value=[columns[0]],
            multi=True
        ),
        html.Label(id='placeholder')
    ], style={'padding': '10 10', 'width': '100%', 'display': 'inline-block'}
    ),
    html.Div([
        html.Br(),
        html.Div([
            dcc.Graph(id='checklist-graph-3'),
        ], style={'padding': '100 100', 'width': '100%', 'display': 'inline-block'}
        ),
        html.Br(),
        html.H3('Module 3: A Closer Look at PIE Curves',
                style={'textAlign': 'center', 'color': '#026A7A'}),
        html.Br(),
        html.P('After the graph above has a plot displayed, click on any point to see the super-PIE'
               ' curve for the m/z ratio. A super-PIE curve all the PIE curves for that m/z stitched '
               'together based on a changing variable (in this case distance). Graph "3" will show '
               'one unique super-PIE curve at a time while Graph "4" will show a superimposition of '
               'super-PIE curves for all the points that are clicked. Both, "2" and "3" are meant for '
               'quick looks at the PIE data to see what should be focused on.',
               style={'textAlign': 'center', 'color': '#026A7A'}),
        html.Br(),
        html.P('Once you have specific m/z and distance values you want to compare in greater depth, '
               'you can select the m/z on Graph "2" and then from the drop down below choose the '
               'corresponding distance. If you wish to compare one m/z value over all distances, '
               'keep changing the distance chosen in the drop down!',
               style={'textAlign': 'center', 'color': '#026A7A'}),
        html.Br(),
        html.Div([dcc.Graph(id='pie-curve'),
                  ], style={'padding': '100 100', 'width': '40%', 'display': 'inline-block',
                            'marginLeft': '7.5%', 'marginRight': '5%'}
                 ),
        html.Div([dcc.Graph(id='pie-curve-superimposed'),
                  ], style={'padding': '100 100', 'width': '40%', 'display': 'inline-block',
                            'marginRight': '7.5%'}
                 ),
        html.Br(),
        html.Br(),
        dcc.Dropdown(
            id='select-pie',
            options=[{'label': i, 'value': i} for i in sheet_names],
            value=[sheet_names[0]],
            multi=False
        ),
        html.Br(),
        html.Br(),
        html.Div([dcc.Graph(id='final-graph'),
                  ], style={'padding': '100 100', 'width': '70%', 'display': 'inline-block',
                            'marginLeft': '15%'}
                 ),
    ]),

    html.Br(),
    html.Br(),
    html.Br(),
    html.P('Designed, Managed and Coded by Ishan Gupta :)',
           style={'textAlign': 'center', 'color': '#FF0000'}),
    html.P("ishan.gupta@berkeley.edu | ChemE and Data Science @ UC Berkeley '23",
           style={'textAlign': 'center', 'color': '#FF0000'})

])

undo_trace = 0
undo_last = 0
x_clicked, y_clicked = [], []


@app.callback(
    Output('checklist-graph', 'figure'),
    [Input('energies', 'value'), Input('memory-storage', 'data'), Input('undo', 'n_clicks'),
     Input('undo-last', 'n_clicks')])
def make_figure(selected_energies, data, click_trace, click_last):
    fig = go.Figure()
    x_axis = main_df['Mass']
    for x in selected_energies:
        y_axis = main_df[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    global undo_trace, undo_last
    global x_clicked
    global y_clicked
    if click_trace > undo_trace:
        undo_trace = click_trace
        x_clicked = []
        y_clicked = []
        # print(x_clicked)
        data = None
    if data != None:
        x_clicked.append(data[0])
        y_clicked.append(data[1])
        if click_last > undo_last:
            undo_last = click_last
            x_clicked = x_clicked[:-2]
            y_clicked = y_clicked[:-2]
            fig.add_trace(go.Scatter(x=x_clicked, y=y_clicked, name='Trend Line'))
        else:
            fig.add_trace(go.Scatter(x=x_clicked, y=y_clicked, name='Trend Line'))
    fig.update_layout(
        title='1: Raw Mass Specs',
        height=500,
        margin=dict(l=20, r=20, b=20, t=30, pad=10),
        paper_bgcolor="LightSteelBlue"
    )
    return fig


@app.callback(
    Output('memory-storage', 'data'),
    [Input('checklist-graph', 'clickData')])
def display_click_data(clickData):
    points_dict = clickData['points'][0]
    coordinates = [points_dict['x'], points_dict['y']]
    # print(coordinates)
    return coordinates


@app.callback(
    Output('placeholder', 'children'),
    [Input('mol_mass1', 'value'), Input('mol_amount1', 'value'), Input('error', 'value'),
     Input('mol_mass2', 'value'), Input('mol_amount2', 'value')])
def make_figure(mass_1, number_1, error, mass_2, number_2):
    global dataframe
    dataframe = get_discrete(mass_1, number_1, error, mass_2, number_2)
    # print(dataframe)
    return 'none'


@app.callback(
    Output('checklist-graph-3', 'figure'),
    [Input('modular-cluster', 'value')])
def make_figure(select_energy):
    global dataframe
    fig = go.Figure()
    x_axis = dataframe['Peaks_Mass']
    for x in select_energy:
        y_axis = dataframe[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    fig.update_layout(
        title='2: Cluster Mass Specs',
        height=500,
        margin=dict(l=20, r=20, b=20, t=30, pad=10),
        paper_bgcolor="LightSteelBlue"
    )
    return fig


@app.callback(
    Output('pie-curve', 'figure'),
    [Input('checklist-graph-3', 'clickData')])
def update_pie_curve(clickData):
    # print(df.head())
    global dataframe
    points_dict = clickData['points'][0]
    x_value = points_dict['x']
    transposed_cluster = dataframe.T
    transposed_cluster.columns = dataframe['Peaks_Mass']
    transposed_cluster = transposed_cluster.drop(['Peaks_Mass'])
    x_axis = transposed_cluster.index
    y_axis = transposed_cluster[x_value]

    fig = px.line(x=x_axis, y=y_axis)
    fig.update_layout(
        title="3: unique super-PIE curve",
        xaxis_title="energy distance pairs",
        yaxis_title="Intensity",
        height=500,
        margin=dict(l=20, r=20, b=20, t=30, pad=10),
        paper_bgcolor="LightSteelBlue"
    )
    return fig


@app.callback(
    Output('pie-curve-superimposed', 'figure'),
    [Input('checklist-graph-3', 'clickData')])
def update_pie_curve(clickData):
    # print(df.head())
    global dataframe
    global pie_curves
    points_dict = clickData['points'][0]
    x_value = points_dict['x']

    transposed_cluster = dataframe.T
    transposed_cluster.columns = dataframe['Peaks_Mass']
    transposed_cluster = transposed_cluster.drop(['Peaks_Mass'])

    x_axis = transposed_cluster.index
    pie_curves[x_value] = transposed_cluster[x_value]

    fig = go.Figure()

    for x in pie_curves:
        y_axis = pie_curves[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    fig.update_layout(
        title="4: multiple super-PIE curves",
        xaxis_title="energy distance pairs",
        yaxis_title="Intensity",
        height=500,
        margin=dict(l=20, r=20, b=20, t=30, pad=10),
        paper_bgcolor="LightSteelBlue"
    )
    return fig


@app.callback(
    Output('final-graph', 'figure'),
    [Input('checklist-graph-3', 'clickData'), Input('select-pie', 'value')])
def final_pie(clickData, distance):
    # print(df.head())
    global dataframe
    global pie_split
    points_dict = clickData['points'][0]
    x_value = points_dict['x']

    transposed_cluster = dataframe.T
    transposed_cluster.columns = dataframe['Peaks_Mass']
    transposed_cluster = transposed_cluster.drop(['Peaks_Mass'])
    energies_list = list(transposed_cluster.index)
    transposed_cluster = transposed_cluster[(np.char.find(energies_list, distance, start=0) >= 0)]
    x_axis = [x[len(distance):] for x in energies_list]
    pie_split[str(x_value) + 'm/z at ' + distance] = list(transposed_cluster[x_value])
    # print(pie_split)
    fig = go.Figure()
    for x in pie_split:
        y_axis = pie_split[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    fig.update_layout(
        title="5: PIE curve with select distance and m/z",
        xaxis_title="energy",
        yaxis_title="Intensity",
        height=500,
        margin=dict(l=20, r=20, b=20, t=30, pad=10),
        paper_bgcolor="LightSteelBlue"
    )
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)


