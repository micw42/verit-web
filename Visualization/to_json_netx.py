import pandas as pd
import numpy as np
import math
import seaborn as sns
from operator import itemgetter
import networkx as nx

from .layeredConcentric import get_xy


def filter_graph():
    nodes_df=pd.read_csv("query_nodes.csv",header = 0)
    edges_df = pd.read_csv("query_edges.csv",header = 0)
    
    edges_df = edges_df[edges_df["thickness"] > 20]
    all_ids = list(np.union1d(edges_df["source"], edges_df["target"]))
    nodes_df = nodes_df[nodes_df["Id"].isin(all_ids)]
    
    edges_df.to_csv("query_edges.csv", index=False)
    nodes_df.to_csv("query_nodes.csv", index=False)


def clean_nodes():

    nodes_df=pd.read_csv("query_nodes.csv",header = 0)

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

    return nodes_df


def clean_edges():
    nodes_df = pd.read_csv("query_nodes.csv",header = 0)
    direct_nodes = nodes_df[nodes_df["Type"]=="Direct"]["Id"].tolist()
    edges_df=pd.read_csv("query_edges.csv",header = 0)

    def get_width(x):
        if x<50:
            return x+10
        else:
            return math.log(x,10)+60

    #Take square root of thickness column so values aren't too big
    edges_df["edge_width"]=edges_df["thickness"].apply(get_width)

    #Convert the color col into hex color strings
    def convert_col(color_val, palette):
        c_segments = np.linspace(-1, 1, len(palette))
        c_i = np.argmin((c_segments - color_val) ** 2)
        
        return palette[c_i]

    pal = list(sns.color_palette("coolwarm_r", as_cmap=False, n_colors=25).as_hex())
    edges_df["color"]=edges_df["color"].apply(convert_col, args=(pal,))

    def get_display(id1, id2, direct_nodes):
        if id1 in direct_nodes or id2 in direct_nodes:
            return "none"
        else:
            return "element"

    edges_df['display'] = edges_df.apply(lambda x: get_display(x.source, x.target, direct_nodes=direct_nodes), axis=1)

    return edges_df


def get_square_clusters():
    nodes_df=pd.read_csv("query_nodes.csv",header = 0)
    edges_df = pd.read_csv("query_edges.csv", header=0)

    is_query = nodes_df[nodes_df["Type"] == "Query"]["Id"].tolist()
    is_connected = list(set(edges_df["source"].tolist()) | set(edges_df["target"].tolist()))  #All non-orphans
    query_orphans = list(set(is_query) - set(is_connected))  # Non-connected query nodes
    query_conn = list(set(is_query) - set(query_orphans))    # Connected query nodes

    linkers = nodes_df[nodes_df["Type"] == "Linker"]["Id"].tolist()
    direct = nodes_df[nodes_df["Type"] == "Direct"]["Id"].tolist()
    n_linkers = int(math.sqrt(len(linkers)))

    # Align connected query nodes
    if len(linkers) == 0:
        align = [{"nodeId": query_conn[i], "position": {"x": 1000*(i%2), "y": 1000*(i//2)}} for i in range(len(query_conn))]
    else:
        align = [{"nodeId": query_conn[i], "position": {"x": (n_linkers*500+5000)*(i%2), "y": (n_linkers)*100*(i//2+1)}} for i in range(len(query_conn))]

    # Align linker nodes
    y_coord = 0
    x_coord = 2500
    for i in range(len(linkers)):
        align.append({"nodeId":linkers[i], "position":{"x":x_coord, "y":y_coord}})
        x_coord += 500
        if i%n_linkers == 0:
            x_coord = 2500
            y_coord += 500

    # Align orphan query nodes
    y_coord = -500
    x_coord = -500
    n_orphans = int(math.sqrt(len(query_orphans)))
    for i in range(len(query_orphans)):
        align.append({"nodeId":query_orphans[i], "position":{"x":x_coord, "y":y_coord}})
        x_coord -= 500
        if (i+1)%n_orphans == 0:
            x_coord = -500
            y_coord -= 500

    return align


#original function
def get_orig_clusters():
    nodes_df=pd.read_csv("query_nodes.csv",header = 0)
    edges_df = pd.read_csv("query_edges.csv", header=0)

    is_query = nodes_df[nodes_df["Type"] == "Query"]["Id"].tolist()
    is_connected = list(set(edges_df["source"].tolist()) | set(edges_df["target"].tolist()))  #All non-orphans
    query_orphans = list(set(is_query) - set(is_connected))  # Non-connected query nodes
    query_conn = list(set(is_query) - set(query_orphans))    # Connected query nodes

    linkers = nodes_df[nodes_df["Type"] == "Linker"]["Id"].tolist()
    direct = nodes_df[nodes_df["Type"] == "Direct"]["Id"].tolist()
    if len(linkers) == 0:
        align = [{"nodeId": query_conn[i], "position": {"x": 1000*(i%2), "y": 1000*(i//2)}} for i in range(len(query_conn))]
    else:
        align = [{"nodeId": query_conn[i], "position": {"x": 10000*(i%2), "y": 10000*(i//2)}} for i in range(len(query_conn))]

    orph_align = [{"nodeId": query_orphans[i],
                   "position": {"x": (-500-500*(i%2)), "y": -500*(i//2)}} for i in range(len(query_orphans))]
    align.extend(orph_align)
    
    align_linkers_n1 = [{"left": query_conn[0], "right": x, "gap": 2500} for x in linkers]
    align_linkers_n2 = [{"left": x, "right": query_conn[1], "gap": 2500} for x in linkers]
    align_linkers = align_linkers_n1 + align_linkers_n2
    
    return (align, align_linkers)


#Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):
    nodes_df = nodes_df.sort_values("Type", ascending=False)

    n_query = len(nodes_df[nodes_df.Type == "Query"])
    n_direct = len(nodes_df[nodes_df.Type == "Direct"])

    elements=[]

    # Construct nodes datatable
    for i in range(len(nodes_df)):
        ndrow = nodes_df.iloc[i]
        node_dict={"data":{"id":ndrow["Id"], "label":ndrow.Label, "KB":ndrow.KB, "type":ndrow["Type"],
                           "syn":ndrow["name"], "color":ndrow.Color, "classes":ndrow["class"],
                           "display":ndrow.display, "orig_display":ndrow.display, "display_id":ndrow.display_id}}

        elements.append(node_dict)

    # Construct edges datatable
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

    return elements, n_query, n_direct


def clean():
    nodes_df=clean_nodes()
    edges_df=clean_edges()
    elements, n_query, n_direct = convert(nodes_df, edges_df)

    return elements, n_query, n_direct


def get_mc_clust():
    nodes_df = pd.read_csv("query_nodes.csv")
    q_nodes = nodes_df[nodes_df["Type"]=="Query"]["Id"].tolist()
    
    edges_df = pd.read_csv("query_edges.csv")
    edges_df = edges_df[edges_df["source"].isin(q_nodes) & edges_df["target"].isin(q_nodes)]
    
    G = nx.from_pandas_edgelist(edges_df, edge_attr=True, source="source", target="target", create_using=nx.DiGraph())
    
    mat = nx.to_scipy_sparse_matrix(G)

    result = mc.run_mcl(mat, inflation=1.2)           # run MCL with default parameters
    clusters = mc.get_clusters(result)

    node_list = list(G.nodes())
    
    clust_list = [itemgetter(*list(x))(node_list) for x in clusters]
    clust_list = [list(x) if type(x) is tuple else [x] for x in clust_list]
    
    small = list(filter(lambda x: len(x) < 5, clust_list))
    small = [item for sublist in small for item in sublist]    #Combine all small clusters into one
    
    clust_list = list(filter(lambda x: len(x) >= 5, clust_list))
    clust_list.append(small)
    
    return clust_list
