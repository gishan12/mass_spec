import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px

global dataframe

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

complete_data = pd.ExcelFile('All_Data copy.xls')

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


reduced_df = reduce(main_df, 48)


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

    print(temp_df['Mass'])
    print(yo)
    print(len(temp_df))
    print(len(yo))
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
    print(a)
    print(temp_df['Mass'])
    temp_df = temp_df.drop(['Mass'], axis=1)
    temp_df['Mass'] = a
    temp_df = temp_df.groupby('Mass').max()
    temp_df = temp_df.drop(labels=[0])
    print(temp_df.columns)
    temp_df["Peaks_Mass"] = np.ones(len(temp_df)) * np.array(array_of_masses)
    return temp_df


app.layout = html.Div([
    html.Div(children=[
        html.Label('Checklist mass'),
        html.Div(children=[
            dcc.Graph(id='checklist-graph'),
            dcc.Dropdown(
                id='energies',
                options=[{'label': i, 'value': i} for i in columns],
                value=[columns[0]],
                multi=True
            ),
        ], style={'padding': '10 10', 'width': '85%', 'display': 'inline-block'}
        ),
        html.Div(children=[
            html.Div(children=[
                html.Button('Undo Trace', id='undo', n_clicks=0)
            ]
            ),
            html.Div(children=[
                html.Button('Undo Last', id='undo-last', n_clicks=0),
                dcc.Store(id='memory-storage')
            ], style={'margin-top': '6vw'}
            )
        ], style={'padding': '10 10', 'display': 'inline-block', 'vertical-align': 'top', 'margin-top': '4vw'}
        ),
    ]),
    html.Div(children=[
        dcc.Graph(id='checklist-graph-2'),
        dcc.Dropdown(
            id='energies-new',
            options=[{'label': i, 'value': i} for i in columns],
            value=[columns[0]],
            multi=True
        ),
    ], style={'padding': '10 10', 'width': '85%', 'display': 'inline-block'}
    ),
    html.Div([
        html.I("Enter the Molecular Mass of species one followed by # of molecules"),
        html.Br(),
        dcc.Input(id="mol_mass1", type="number", placeholder="molecular mass", style={'marginRight': '10px'}),
        dcc.Input(id="mol_amount1", type="number", placeholder="# of molecules", debounce=True),
        dcc.Input(id="error", type="number", placeholder="error", debounce=True),
        dcc.Input(id="mol_mass2", type="number", placeholder="molecular mass 2", style={'marginRight': '10px'}),
        dcc.Input(id="mol_amount2", type="number", placeholder="# of molecules 2", debounce=True),
        html.P(id='placeholder')
        # might need to add a dropdown for energy to shorten the stored datasets
    ]),
    html.Div(children=[
        dcc.Graph(id='checklist-graph-3'),
        dcc.Dropdown(
            id='modular-cluster',
            options=[{'label': i, 'value': i} for i in columns],
            value=[columns[0]],
            multi=True
        ),
    ], style={'padding': '10 10', 'width': '85%', 'display': 'inline-block'}
    )
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

    return fig


@app.callback(
    Output('checklist-graph-2', 'figure'),
    [Input('energies-new', 'value')])
def make_figure(selected_energies):
    fig = go.Figure()
    x_axis = reduced_df['Mass']
    for x in selected_energies:
        y_axis = reduced_df[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    return fig


@app.callback(
    Output('placeholder', 'children'),
    [Input('mol_mass1', 'value'), Input('mol_amount1', 'value'), Input('error', 'value'),
     Input('mol_mass2', 'value'), Input('mol_amount2', 'value')])
def make_figure(mass_1, number_1, error, mass_2, number_2):
    global dataframe
    dataframe = get_discrete(mass_1, number_1, error, mass_2, number_2)
    return 'none'


@app.callback(
    Output('checklist-graph-3', 'figure'),
    [Input('modular-cluster', 'value')])
def make_figure(selected_energies):
    global dataframe
    fig = go.Figure()
    x_axis = dataframe['Peaks_Mass']
    for x in selected_energies:
        y_axis = dataframe[x]
        fig.add_trace(go.Scatter(x=x_axis, y=y_axis, name=x))
    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
