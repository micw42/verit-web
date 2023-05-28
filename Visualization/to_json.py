import pandas as pd
import numpy as np
import seaborn as sns
import math
import pickle

from .layeredConcentric import get_xy


def clean_nodes():
    nodes_df = pd.read_csv("query_nodes.csv", header=0)

    # Duplicate the depth column in nodes table
    nodes_df["rank"] = nodes_df["depth"]

    # Convert rank values so small values are large and vice versa
    def get_rank_value(rank):
        if rank == 0:
            return 1000000000
        elif rank == 1:
            return 1000000
        elif rank == 2:
            return 10000
        elif rank == 3:
            return 1000
        elif rank == 4:
            return 1

    nodes_df["rank"] = nodes_df["rank"].apply(get_rank_value)

    # Convert depth to color in nodes table
    def get_node_color(depth):
        if depth == 0:
            return "#fc0800"
        elif depth == 1:
            return "#f1c9f2"
        elif depth == 2:
            return "#c9ddf2"
        elif depth == 3:
            return "#d7f2c9"
        else:
            return "#f7d4ab"

    nodes_df["color"] = nodes_df["depth"].apply(get_node_color)

    return nodes_df


def clean_edges():
    edges_df = pd.read_csv("query_edges.csv", header=0)

    def get_width(x):
        if x < 50:
            return x+10
        else:
            return math.log(x, 10)+60

    # So edge thicknesses aren't too big
    edges_df["edge_width"] = edges_df["thickness"].apply(get_width)

    #Convert the color col into hex color strings
    def convert_col(color_val, palette):
        c_segments = np.linspace(-1, 1, len(palette))
        c_i = np.argmin((c_segments - color_val) ** 2)

        return palette[c_i]

    pal = list(sns.color_palette("coolwarm_r", as_cmap=False, n_colors=25).as_hex())
    edges_df["color"]=edges_df["color"].apply(convert_col, args=(pal,))

    return edges_df


# Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df, sq=False):
    elements = []

    # Construct nodes datatable
    if not sq:
        for i in range(len(nodes_df)):
            ndrow = nodes_df.iloc[i]
            node_dict = {"data": {"id": ndrow["Id"],
                                  "label": ndrow.Label,
                                  "color": ndrow.color,
                                  "KB": ndrow.KB,
                                  "display_id": ndrow.display_id,
                                  "rank": int(ndrow["rank"])
                                  }}

    if sq:
        # Sort the nodes by thickness to be arranged polarly
        query_id = nodes_df[nodes_df["depth"] == 0].iloc[0]["Id"]

        nq1 = edges_df[edges_df.source == query_id][["target", "thickness"]].rename(columns={"target": "Id"})
        nq2 = edges_df[edges_df.target == query_id][["source", "thickness"]].rename(columns={"source": "Id"})
        
        nq_df = pd.concat([nq1, nq2])
        nq_df = nq_df.groupby("Id").max().reset_index()
        nodes_df = nodes_df.merge(nq_df, on="Id", how="left")
        nodes_df = nodes_df.sort_values(["depth", "thickness"], ascending=[True, False])


        # Calculate x, y coordinates for each node
        Xs, Ys, _, __ = get_xy(len(nodes_df)-1, n_fl_co=20, r=1050)
        Xs = [0] + list(Xs)    # First index is 0 because it's the query node
        Ys = [0] + list(Ys)

        for i in range(len(nodes_df)):
            ndrow = nodes_df.iloc[i]

            node_dict = {"data": {"id": ndrow["Id"],
                                  "label": ndrow.Label,
                                  "color": ndrow.color,
                                  "KB": ndrow.KB,
                                  "display_id": ndrow.display_id,
                                  "rank": int(ndrow["rank"])
                                  }}
            node_dict["position"] = {"x": Xs[i], "y": Ys[i]}

            elements.append(node_dict)

    # Construct edges datatable
    # Does the order of edges and nodes have to be the same?
    edges_df["edge_id"] = edges_df.source.str.cat(edges_df.target)
    for i in range(len(edges_df)):
        erow = edges_df.iloc[i]
        edge_dict = {"data": {"id": erow.edge_id,
                              "source": erow.source, "target": erow.target,
                              "weight": float(erow.edge_width),
                              "color": erow.color,
                              "files": erow.files,
                              "thickness": int(erow.thickness)}}
        elements.append(edge_dict)

    return elements


def clean(sq=False):
    '''
    sq: ad hoc change for singlequery to manually add x and y values
    '''
    nodes_df = clean_nodes()
    edges_df = clean_edges()
    elements = convert(nodes_df, edges_df, sq=sq)
    
    return elements
