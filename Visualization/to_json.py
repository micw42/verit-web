import pandas as pd
import numpy as np
import seaborn as sns
import math

def clean_nodes():

    nodes_df=pd.read_csv("query_nodes.csv",header = 0)

    #Duplicate the depth column in nodes table
    nodes_df["rank"]=nodes_df["depth"]

    #Convert rank values so small values are large and vice versa
    def get_rank_value(rank):
        if rank==0:
            return 1000000000
        elif rank==1:
            return 1000000
        elif rank==2:
            return 10000
        elif rank==3:
            return 1000
        elif rank==4:
            return 1

    nodes_df["rank"]=nodes_df["rank"].apply(get_rank_value)


    #Convert depth to color in nodes table
    def get_node_color(depth):
        if depth==0:
            return "#fc0800"
        elif depth==1:
            return "#f1c9f2"
        elif depth==2:
            return "#c9ddf2"
        elif depth==3:
            return "#d7f2c9"
        else:
            return "#f7d4ab"

    nodes_df["depth"]=nodes_df["depth"].apply(get_node_color)

    return nodes_df

def clean_edges():

    edges_df=pd.read_csv("query_edges.csv",header = 0)
 
    def get_width(x):
        if x<50:
            return x+10
        else:
            return math.log(x,10)+60

    # So edge thicknesses aren't too big
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

    return edges_df

#Convert nodes and edges tables into one json-style list
def convert(nodes_df, edges_df):

    nodes=[]
    for _, row in nodes_df.iterrows():
        parts=row.values.tolist()
        nodes.append(parts)

    elements=[]
    for node in nodes:
        node_dict={"data":{"id":node[0], "label":node[1], "color":node[2], "rank":int(node[3])}}
        elements.append(node_dict)

    edges=[]
    for _, row in edges_df.iterrows():
        parts=row.values.tolist()
        edges.append(parts)

    for edge in edges:
        edge_id=edge[3]+edge[4]
        edge_dict={"data":{"id":edge_id, "source":edge[3], "target":edge[4], 
                           "weight":edge[5], "color":edge[0], 
                           "files":edge[2], "thickness":edge[1]}}
        elements.append(edge_dict)


    return elements

def clean():
    nodes_df=clean_nodes()
    edges_df=clean_edges()
    elements=convert(nodes_df, edges_df)
    return elements











