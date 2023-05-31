import pandas as pd
import copy
import math
import random
import numpy as np

def makeQuery(edges_df, query_list, thresh=20):
    subset_df = edges_df[edges_df["source"].isin(query_list) | edges_df["target"].isin(query_list)]
    subset_df = subset_df[subset_df["thickness"]>thresh]
    return subset_df

def filter_nodes(full_df, prev_query, max_nodes=100):
    all_nodes = full_df["target"].tolist() + full_df["source"].tolist()
    all_thickness = (full_df["thickness"].tolist()) * 2
    ranked_df = pd.DataFrame.from_dict({"node":all_nodes, "thickness":all_thickness})
    ranked_df = ranked_df[~ranked_df["node"].isin(prev_query)]   #Remove already-found nodes
    ranked_df = ranked_df.groupby(['node']).agg(max).sort_values(by = 'thickness', ascending=False) #Get maximum thickness edge that node has
    if len(ranked_df.index) > max_nodes:     
        passing_nodes = (ranked_df.index.tolist())[:max_nodes] 
        #Keep only top X nodes
        full_df = full_df[(full_df["source"].isin(prev_query) & (full_df["target"].isin(passing_nodes) | full_df["target"].isin(prev_query))) | (full_df["source"].isin(passing_nodes) & (full_df["target"].isin(prev_query)))]
    return full_df

def query(G, edges_df, nodes_df, db_df, query, depth):
    user_query = list(query.keys())[0]
    # For name queries (accounts for the user selecting multiple IDs)
    if user_query != "QUERY_ID":
        query_list = list(query.values())
        query_list = [item for sublist in query_list for item in sublist]
    # For ID queries
    else:
        query_list = list(query.values())
        

    if depth==1:
        full_df = makeQuery(edges_df, query_list, thresh=0)
        full_df["depth"]=1
    else:
        # Recursively find target nodes
        queried = []   #Already-queried nodes
        targets = copy.deepcopy(query_list)
        full_df_cols = list(edges_df.columns)
        full_df_cols.append("depth")
        full_df = pd.DataFrame(columns = full_df_cols)
        for i in range(depth):
            query_df = makeQuery(edges_df, targets, thresh=20)  #pre filter to avoid large networks
            queried.extend(targets)
            query_df["depth"] = i+1
            nodes = np.union1d(query_df["target"].tolist(), query_df["source"].tolist())
            targets = np.setdiff1d(nodes, queried)    #Get new nodes to use in the next query
            full_df = full_df.append(query_df)

    #Bidirectional edges
    opp_df = full_df.merge(edges_df, left_on=["source", "target"], right_on=["target","source"])
    opp_df = opp_df.drop(labels=["source_x","target_x","color_x", "thickness_x"], axis=1).rename(columns={"source_y":"source", "target_y":"target", "color_y":"color", "thickness_y":"thickness"})
    full_df = pd.concat([full_df, opp_df]).drop_duplicates(subset=["source", "target"])
    
    # Make query nodes depth = 0, so they're in the center of the visualization
    df_dict = {"target":query_list, "depth":[0]*len(query_list)}
    zero_rows = pd.DataFrame.from_dict(df_dict)
    # Only need the targets, since every node (except for some query nodes) are a target at least once
    targets = full_df[["target", "depth"]]
    nodes = zero_rows.append(targets)
    sources = full_df[["source", "depth"]]
    sources = sources.rename(columns={'source': 'target'})
    nodes = nodes.append(sources)
    nodes = nodes.sort_values('depth').drop_duplicates('target').sort_index()
    nodes_df = nodes_df.drop_duplicates(subset='Id', keep="first")
    nodes = pd.merge(nodes, nodes_df, left_on = "target", right_on = "Id", how="inner")
    
    # Add synonyms to the nodes
    nodes = nodes.merge(db_df, left_on="Id", right_on="id", how="left")
    nodes["name"] = nodes["name"].fillna(nodes["Label"])
    nodes["id"] = nodes["id"].fillna(nodes["Id"])
    nodes["name"] = nodes["name"].fillna("NAN")
    syn_concat = lambda x: "%%".join(x)  # Separate each synonym with %%
    aggregation_functions = {'Id': 'first', 'Label':"first", "depth":"first", "KB":"first", "PR":"first", "name":syn_concat}
    nodes = nodes.groupby('id').aggregate(aggregation_functions)
    nodes.to_csv("intermediate_nodes.csv", index=False)
    
    # If the user selected multiple IDs, merge them all into one node
    if user_query != "QUERY_ID":
        nodes.loc[(nodes.Id.isin(query_list)),'display_id'] = ", ".join(query_list)
        nodes["display_id"] = nodes["display_id"].fillna(nodes["Id"])
        nodes.loc[(nodes.Id.isin(query_list)),'Label']=user_query
        nodes.loc[(nodes.Id.isin(query_list)),'Id']=user_query
        nodes = nodes.drop_duplicates(subset="Id")
        full_df["orig_source"] = full_df["source"]
        full_df["orig_target"] = full_df["target"]
        full_df.loc[(full_df.source.isin(query_list)),'source']=user_query
        full_df.loc[(full_df.target.isin(query_list)),'target']=user_query
    else:
        nodes["display_id"] = nodes["Id"]
        full_df["orig_source"] = full_df["source"]
        full_df["orig_target"] = full_df["target"]
        
    syn_concat = lambda x: "%%".join(x)
    aggregation_functions = { "depth":"first", "PR":"first", "Label":"first", "KB":"first", "display_id":"first", 'name': syn_concat}
    nodes = nodes.groupby("Id").aggregate(aggregation_functions).reset_index()

    full_df = full_df[["color", "thickness", "source", "target", "orig_source", "orig_target"]]
    full_df["color2"] = full_df["color"] * full_df["thickness"]
    id_concat = lambda x: "%%".join(x) # Concat all source and target IDs of merged nodes
    aggregation_functions = {'color2': 'sum','thickness': 'sum', "orig_source":id_concat, "orig_target":id_concat}
    full_df = full_df.groupby(["source", "target"]).aggregate(aggregation_functions).reset_index()
    full_df["color"] = full_df["color2"]/full_df["thickness"]
    def formatter(sources, targets):
        source_list = sources.split("%%")
        target_list = targets.split("%%")
        file_list = [f"{source_list[i]}_{target_list[i]}.txt" for i in range(len(source_list))]
        files = "%%".join(file_list)
        return files
    full_df["files"] = full_df.apply(lambda x: formatter(x.orig_source, x.orig_target), axis=1)
    
    nodes = nodes[["Id", "Label", "depth", "KB", "display_id", "name"]]
    full_df = full_df[["color", "thickness", 
                     "files", "source", "target"]]
    
    
    full_df.to_csv("query_edges.csv", index=False)
    nodes.to_csv("query_nodes.csv", index=False)
