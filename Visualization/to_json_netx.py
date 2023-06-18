import pandas as pd
import numpy as np
import math
import seaborn as sns
from operator import itemgetter
import networkx as nx

from Visualization import layeredConcentric


def filter_graph():
    nodes_df=pd.read_csv("query_nodes.csv",header = 0)
    edges_df = pd.read_csv("query_edges.csv",header = 0)
    
    edges_df = edges_df[edges_df["thickness"] > 20]
    all_ids = list(np.union1d(edges_df["source"], edges_df["target"]))
    nodes_df = nodes_df[nodes_df["Id"].isin(all_ids)]
    
    edges_df.to_csv("query_edges.csv", index=False)
    nodes_df.to_csv("query_nodes.csv", index=False)


def clean_nodes(nodes_df, layer):
    def get_color(query_type):
        if query_type == "Query":
            return "#fc0800"
        else:
            return "#99beff"

    def get_class(query_type):
        if query_type == "Query":
            return "level3"
        else:
            return "level0"

    def get_display(query_type):
        if query_type=="Direct":
            return "none"
        else:
            return "element"

    nodes_df["Color"] = nodes_df["Type"].apply(get_color)
    nodes_df["class"] = nodes_df["Type"].apply(get_class)
    nodes_df["display"] = nodes_df["Type"].apply(get_display)
    nodes_df["layer"] = layer

    return nodes_df


def clean_edges(nodes_df, edges_df, layer):
    direct_nodes = nodes_df[nodes_df["Type"] == "Direct"]["Id"].tolist()

    def get_width(x):
        if x < 50:
            return x + 10
        else:
            return math.log(x, 10) + 60

    #Take square root of thickness column so values aren't too big
    edges_df["edge_width"] = edges_df["thickness"].apply(get_width)

    #Convert the color col into hex color strings
    def convert_col(color_val, palette):
        c_segments = np.linspace(-1, 1, len(palette))
        c_i = np.argmin((c_segments - color_val) ** 2)
        
        return palette[c_i]

    pal = list(sns.color_palette("coolwarm_r", as_cmap=False, n_colors=25).as_hex())
    edges_df["color"] = edges_df["color"].apply(convert_col, args=(pal,))

    def get_display(id1, id2, direct_nodes):
        if id1 in direct_nodes or id2 in direct_nodes:
            return "none"
        else:
            return "element"

    edges_df['display'] = edges_df.apply(lambda x: get_display(x.source, x.target, direct_nodes=direct_nodes), axis=1)
    edges_df["layer"] = layer
    
    return edges_df


#Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):
    elements=[]
    coord_dict = dict()
    layers = list(nodes_df.layer.unique())

    for layer in layers:
        nodes_layer = nodes_df[nodes_df.layer == layer].copy()
        edges_layer = edges_df[edges_df.layer == layer].copy()

        vc_nodes = nodes_layer.Type.value_counts()

        n_query = vc_nodes["Query"]
        n_links = vc_nodes.sum() - vc_nodes["Query"]

        # Compute X and Y for concentric layout
        Xs = []; Ys = []

        r1 = 500
        Xs_q, Ys_q, R_arr_q, n_arr_q = layeredConcentric.get_xy(n_query, r=r1)
        Xs.extend(Xs_q); Ys.extend(Ys_q)

        r2 = 100
        n_fl_co_d = 2 * np.pi * (R_arr_q[-1] + 3*r1) / (2 * r2)
        Xs_d, Ys_d, R_arr_d, n_arr_d = layeredConcentric.get_xy(n_links, n_fl_co_d, r=r2)
        Xs.extend(Xs_d); Ys.extend(Ys_d)


        # Construct nodes datatable
        for i in range(len(nodes_layer)):
            ndrow = nodes_layer.iloc[i]
            node_dict={"data":{"id":ndrow["Id"], "label":ndrow.Label, "KB":ndrow.KB, "type":ndrow["Type"],
                               "syn":ndrow["name"], "color":ndrow.Color, "classes":ndrow["class"],
                               "display":ndrow.display, "orig_display":ndrow.display, "display_id":ndrow.display_id, 
                               "layer":ndrow.layer, "lcX": Xs[i], "lcY": Ys[i]
                              }}

            elements.append(node_dict)

        # Construct edges datatable
        edges_layer["edge_id"] = edges_layer.source.str.cat(edges_layer.target)
        for i in range(len(edges_layer)):
            erow = edges_layer.iloc[i]
            edge_dict = {"data": {"id": erow.edge_id,
                                  "source": erow.source, "target": erow.target,
                                  "weight": float(erow.edge_width),
                                  "color": erow.color,
                                  "files": erow.files,
                                  "thickness": int(erow.thickness),
                                  "layer": erow.layer
                                 }}
            elements.append(edge_dict)

    return elements


def clean(biogrid=False):
    if biogrid:
        nodes_df_reach = pd.read_csv("query_nodes.csv", header=0)
        edges_df_reach = pd.read_csv("query_edges.csv", header=0)

        nodes_df_reach = clean_nodes(nodes_df_reach, layer="reach")
        edges_df_reach = clean_edges(nodes_df_reach, edges_df_reach, layer="reach")

        nodes_df_bg = pd.read_csv("query_nodes_BIOGRID.csv", header=0)
        edges_df_bg = pd.read_csv("query_edges_BIOGRID.csv", header=0)

        nodes_df_bg = clean_nodes(nodes_df_bg, layer="biogrid")
        edges_df_bg = clean_edges(nodes_df_bg, edges_df_bg, layer="biogrid")

        nodes_df = pd.concat([nodes_df_reach, nodes_df_bg])
        edges_df = pd.concat([edges_df_reach, edges_df_bg])

    else:
        nodes_df = pd.read_csv("query_nodes.csv", header=0)
        edges_df = pd.read_csv("query_edges.csv", header=0)

        nodes_df = clean_nodes(nodes_df, layer="reach")
        edges_df = clean_edges(nodes_df, edges_df, layer="reach")
    
    elements = convert(nodes_df, edges_df)

    return elements

