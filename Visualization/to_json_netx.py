import pandas as pd
import numpy as np
import math
import seaborn as sns

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
        if float(color_val) < -0.8:
            return palette[0]
        elif float(color_val)< -0.6:
            return palette[1]
        elif float(color_val)< -0.4:
            return palette[2]
        elif float(color_val)< -0.2:
            return palette[3]
        elif float(color_val)< 0:
            return palette[4]
        elif float(color_val)< 0.2:
            return palette[5]
        elif float(color_val)< 0.4:
            return palette[6]
        elif float(color_val)< 0.6:
            return palette[7]
        elif float(color_val)< 0.8:
            return palette[8]
        return palette[9]
    
    pal = list(sns.color_palette("RdBu", 10).as_hex())
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
        
    orph_align = [{"nodeId": query_orphans[i], "position": {"x": (-500-500*(i%2)), "y": -500*(i//2)}} for i in range(len(query_orphans))]
    align.extend(orph_align)
    align_linkers_n1 = [{"left": query_conn[0], "right": x, "gap": 2500} for x in linkers]
    align_linkers_n2 = [{"left": x, "right": query_conn[1], "gap": 2500} for x in linkers]
    align_linkers = align_linkers_n1 + align_linkers_n2
    return (align, align_linkers)   

#Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):

    nodes=[]
    for _, row in nodes_df.iterrows():
        parts=row.values.tolist()
        nodes.append(parts)

    elements=[]
    for node in nodes:
        node_dict={"data":{"id":node[0], "label":node[1], "type":node[3], "syn":node[2],
                           "color":node[4], "classes":node[5], "display":node[6], "orig_display":node[6]}}
        elements.append(node_dict)

    edges=[]
    for _, row in edges_df.iterrows():
        parts=row.values.tolist()
        edges.append(parts)
    
    for edge in edges:
        edge_id=edge[3]+edge[4]
        edge_dict={"data":{"id":edge_id, "source":edge[3], "target":edge[4], 
                           "weight":edge[5], "color":edge[0],  
                           "display":edge[6], "orig_display":edge[6],
                          "thickness":edge[1], "files":edge[2]}}
        elements.append(edge_dict)

    return elements

def clean():
    nodes_df=clean_nodes()
    edges_df=clean_edges()
    return convert(nodes_df, edges_df)





