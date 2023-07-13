import pandas as pd
import numpy as np
import seaborn as sns
import math
import pickle

from .layeredConcentric import get_xy


def clean_nodes(nodes_df, layer):

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
    nodes_df["layer"] = layer
    if layer=="reach":
        nodes_df["display"] = "element"
    else:
        nodes_df["display"] = "none"
    nodes_df["border_width"] = 2
    nodes_df["border_color"] = "#0000FFFF"
    return nodes_df


def get_width(x):
        if x < 50:
            return x+10
        else:
            return math.log(x, 10)+60


def clean_edges(edges_df, layer):
    # So edge thicknesses aren't too big
    edges_df["edge_width"] = edges_df["thickness"].apply(get_width)

    #Convert the color col into hex color strings
    def convert_col(color_val, palette):
        c_segments = np.linspace(-1, 1, len(palette))
        c_i = np.argmin((c_segments - color_val) ** 2)

        return palette[c_i]

    pal = list(sns.color_palette("coolwarm_r", as_cmap=False, n_colors=25).as_hex())
    edges_df["color"] = edges_df["color"].apply(convert_col, args=(pal,))
    edges_df["layer"] = layer
    if layer=="reach":
        edges_df["display"] = "element"
    else:
        edges_df["display"] = "none"    
    return edges_df


def clean_union(nodes_df, edges_df):
    ## Nodes
    gb_size = nodes_df.groupby("Id").size()
    union_ids = gb_size[gb_size == 2].index

    union_nodes_df = nodes_df.copy()
    union_nodes_df = union_nodes_df.drop_duplicates(subset="Id")

    outline_vec = union_nodes_df["layer"].copy().replace({
        "reach": "#9d49f2",
        "biogrid": "#77ed40"
    })
    outline_vec[union_nodes_df["display_id"].isin(union_ids)] = "#42a7f5"

    union_nodes_df["layer"] = "union"

    union_nodes_df["border_color"] = outline_vec
    union_nodes_df["border_width"] = 50
    union_nodes_df["display"] = "none"
    
    nodes_df = pd.concat([nodes_df, union_nodes_df])
    
    ## Edges
    union_edges_df = edges_df.copy()
    # Sum the thicknesses together between reach and BIOGRID
    union_thickness = union_edges_df.groupby(["source", "target"]).apply(lambda x: x.thickness.sum())
    union_thickness = union_thickness.reset_index()

    # Format layer, thickness, and edge_width columns appropriately
    union_edges_df = union_edges_df.drop_duplicates(subset=["source", "target"])
    union_edges_df = union_edges_df.merge(
        union_thickness).drop(
        columns="thickness").rename(
        columns={0: "thickness"})

    union_edges_df["edge_width"] = union_edges_df["thickness"].apply(get_width)
    
    union_edges_df["layer"] = "union"
    union_edges_df["display"] = "none"
    edges_df = pd.concat([edges_df, union_edges_df])

    return nodes_df, edges_df


# Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):
    elements = []
    layers = list(nodes_df.layer.unique())

    for layer in layers:
        nodes_layer = nodes_df[nodes_df.layer == layer].copy()
        edges_layer = edges_df[edges_df.layer == layer].copy()
        
        # Sort the nodes by thickness to be arranged polarly
        query_id = nodes_layer[nodes_layer["depth"] == 0].iloc[0]["Id"]

        nq1 = edges_layer[edges_layer.source == query_id][["target", "thickness"]].rename(columns={"target": "Id"})
        nq2 = edges_layer[edges_layer.target == query_id][["source", "thickness"]].rename(columns={"source": "Id"})

        nq_df = pd.concat([nq1, nq2])
        nq_df = nq_df.groupby("Id").max().reset_index()
        nodes_layer = nodes_layer.merge(nq_df, on="Id", how="left")
        nodes_layer = nodes_layer.sort_values(["depth", "thickness"], ascending=[True, False])


        # Calculate x, y coordinates for each node
        Xs, Ys, _, __ = get_xy(len(nodes_layer)-1, n_fl_co=20, r=1050)
        Xs = [0] + list(Xs)    # First index is 0 because it's the query node
        Ys = [0] + list(Ys)

        for i in range(len(nodes_layer)):
            ndrow = nodes_layer.iloc[i]

            node_dict = {"data": {"id": ndrow["Id"]+ndrow["layer"],
                                  "label": ndrow.Label,
                                  "color": ndrow.color,
                                  "KB": ndrow.KB,
                                  "display_id": ndrow.display_id,
                                  "syn": ndrow["name"],
                                  "rank": int(ndrow["rank"]),
                                  "layer": layer,
                                  "display":ndrow.display,
                                  "border_color":ndrow["border_color"], "border_width":int(ndrow["border_width"]),
                                 "depth":int(ndrow.depth)
                                  }}
            node_dict["position"] = {"x": Xs[i], "y": Ys[i]}

            elements.append(node_dict)

        # Construct edges datatable
        # Does the order of edges and nodes have to be the same?
        edges_layer["edge_id"] = edges_layer.source.str.cat(edges_layer.target)
        for i in range(len(edges_layer)):
            erow = edges_layer.iloc[i]
            edge_dict = {"data": {"id": erow.edge_id+erow.layer,
                                  "source": erow.source+erow.layer, "target": erow.target+erow.layer,
                                  "weight": float(erow.edge_width),
                                  "color": erow.color,
                                  "files": erow.files,
                                  "thickness": int(erow.thickness),
                                 "layer": layer,
                                 "display":erow.display}}
            elements.append(edge_dict)

    return elements


def clean(nodes_df_reach, edges_df_reach, nodes_df_bg=None, edges_df_bg=None, biogrid=False):
    if biogrid:
        nodes_df_reach = clean_nodes(nodes_df_reach, layer="reach")
        edges_df_reach = clean_edges(edges_df_reach, layer="reach")

        nodes_df_bg = clean_nodes(nodes_df_bg, layer="biogrid")
        edges_df_bg = clean_edges(edges_df_bg, layer="biogrid")
        edges_df_bg.to_csv("edges_df_bg.csv", index=False)
        edges_df_bg_switched = edges_df_bg.rename(columns={"source":"target", "target":"source", "source_id":"target_id", "target_id":"source_id"})   #Bidirectional BG edges
        edges_df_bg = pd.concat([edges_df_bg, edges_df_bg_switched]).drop_duplicates()

        nodes_df = pd.concat([nodes_df_reach, nodes_df_bg])
        edges_df = pd.concat([edges_df_reach, edges_df_bg])
        
        nodes_df, edges_df = clean_union(nodes_df, edges_df)
        
        nodes_df.to_csv("qnodes_union.csv", index=False)
        edges_df.to_csv("qedges_union.csv", index=False)
    else:
        nodes_df = clean_nodes(nodes_df_reach, layer="reach")
        edges_df = clean_edges(edges_df_reach, layer="reach")
    
    elements = convert(nodes_df, edges_df)
    
    with open("elements.pkl", "wb") as p:
        pickle.dump(elements, p)
    
    return elements
