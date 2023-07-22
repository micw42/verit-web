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

    def get_display(query_type, layer):
        if query_type=="Direct" or layer=="biogrid":
            return "none"
        else:
            return "element"

    nodes_df["Color"] = nodes_df["Type"].apply(get_color)
    nodes_df["class"] = nodes_df["Type"].apply(get_class)
    nodes_df["display"] = nodes_df["Type"].apply(get_display, layer=layer)
    nodes_df["layer"] = layer
    nodes_df["border_width"] = 2
    nodes_df["border_color"] = "#0000FFFF"

    return nodes_df


def get_width(x):
    if x < 50:
        return x
    else:
        return math.log(x, 10) + 50


def clean_edges(nodes_df, edges_df, layer):
    direct_nodes = nodes_df[nodes_df["Type"] == "Direct"]["Id"].tolist()
    
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
    
    def get_type(id1, id2, direct_nodes):
        if id1 in direct_nodes or id2 in direct_nodes:
            return "Direct"
        else:
            return "Query"
        
    edges_df['type'] = edges_df.apply(lambda x: get_type(x.source, x.target, direct_nodes=direct_nodes), axis=1)   
    edges_df["layer"] = layer
    return edges_df


def clean_union(nodes_df, edges_df_reach, edges_df_bg):
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
    union_nodes_df["border_width"] = 16
    union_nodes_df["display"] = "none"
    
    nodes_df = pd.concat([nodes_df, union_nodes_df])
        
    ## Edges
    edges_df_bg_switched = edges_df_bg.rename(columns={"source":"target", "target":"source", "source_id":"target_id", "target_id":"source_id"})   #Bidirectional BG edges
    edges_df_bg_switched["files"] = edges_df_bg_switched["source"] + "_" + edges_df_bg_switched["target"] + ".txt"
    edges_df_bg_all=pd.concat([edges_df_bg, edges_df_bg_switched]).drop_duplicates()
    union_edges_df = pd.concat([edges_df_reach, edges_df_bg_all])
    
    # Retain thickness for each layer separately
    bg_thickness = union_edges_df[union_edges_df["layer"] == "biogrid"][["source", "target", "thickness"]].rename(columns={"thickness": "thickness_bg"})
    r_thickness = union_edges_df[union_edges_df["layer"] == "reach"][["source", "target", "thickness"]].rename(columns={"thickness": "thickness_r"})
    
    # Sum the thicknesses together between reach and BIOGRID
    union_thickness = union_edges_df.groupby(["source_id", "target_id"]).apply(lambda x: x.thickness.sum())
    union_thickness = union_thickness.reset_index()

    # Format layer, thickness, and edge_width columns appropriately
    union_edges_df = union_edges_df.drop_duplicates(subset=["source_id", "target_id"])
    union_edges_df = union_edges_df.merge(
        union_thickness).drop(
        columns="thickness").rename(
        columns={0: "thickness"})

    # Keep separated thicknesses for union layer to separately filter
    union_edges_df = union_edges_df.merge(r_thickness, on=["source", "target"], how="outer")
    union_edges_df = union_edges_df.merge(bg_thickness, on=["source", "target"], how="outer")

    # Create edge width from sum of thickness
    union_edges_df["edge_width"] = union_edges_df["thickness"].apply(get_width)
    
    union_edges_df["layer"] = "union"

    edges_df = pd.concat([edges_df_reach, edges_df_bg, union_edges_df])
    edges_df = edges_df.replace({np.nan: 0})

    return nodes_df, edges_df


#Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):
    elements=[]
    layers = list(nodes_df.layer.unique())

    for layer in layers:
        nodes_layer = nodes_df[nodes_df.layer == layer].copy()
        nodes_layer = nodes_layer.sort_values('Type', ascending=False)
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
            node_dict={"data":{"id":ndrow["Id"]+ndrow["layer"], "label":ndrow.Label, "KB":ndrow.KB, "type":ndrow["Type"],
                               "syn":ndrow["name"], "color":ndrow.Color, "classes":ndrow["class"],
                               "display":ndrow.display, "orig_display":ndrow.display, "display_id":ndrow.display_id, 
                               "layer":ndrow.layer, "lcX": Xs[i], "lcY": Ys[i], "border_color":ndrow["border_color"], "border_width":int(ndrow["border_width"])
                              }}
            
            elements.append(node_dict)
            
        # Construct edges datatable
        edges_layer["edge_id"] = edges_layer.source.str.cat(edges_layer.target)
        if layer == "union":
            for i in range(len(edges_layer)):
                erow = edges_layer.iloc[i]
                edge_dict = {"data": {"id": erow.edge_id+erow.layer,
                                      "source": erow.source+erow.layer, "target": erow.target+erow.layer,
                                      "weight": float(erow.edge_width) * (erow.layer == "reach" or erow.layer=="union") + 10,
                                      "color": erow.color,
                                      "files": erow.files,
                                      "thickness": int(erow.thickness),
                                      "thickness_bg": int(erow.thickness_bg),
                                      "thickness_r": int(erow.thickness_r),
                                      "layer": layer,
                                      "type":erow.type}}
                elements.append(edge_dict)

        else:
            for i in range(len(edges_layer)):
                erow = edges_layer.iloc[i]
                edge_dict = {"data": {"id": erow.edge_id+erow.layer,
                                      "source": erow.source+erow.layer, "target": erow.target+erow.layer,
                                      "weight": float(erow.edge_width) * (erow.layer == "reach" or erow.layer=="union") + 10,
                                      "color": erow.color,
                                      "files": erow.files,
                                      "thickness": int(erow.thickness),
                                      "layer": layer,
                                      "type":erow.type}}
                elements.append(edge_dict)

    return elements


def clean(nodes_df_reach, edges_df_reach, nodes_df_bg=None, edges_df_bg=None, biogrid=False):
    if biogrid:
        nodes_df_reach = clean_nodes(nodes_df_reach, layer="reach")
        edges_df_reach = clean_edges(nodes_df_reach, edges_df_reach, layer="reach")

        nodes_df_bg = clean_nodes(nodes_df_bg, layer="biogrid")
        edges_df_bg = clean_edges(nodes_df_bg, edges_df_bg, layer="biogrid")
        
        nodes_df = pd.concat([nodes_df_reach, nodes_df_bg])
        nodes_df, edges_df = clean_union(nodes_df, edges_df_reach, edges_df_bg)

    else:
        nodes_df = clean_nodes(nodes_df_reach, layer="reach")
        edges_df = clean_edges(nodes_df_reach, edges_df_reach, layer="reach")
    
    elements = convert(nodes_df, edges_df)

    return elements
