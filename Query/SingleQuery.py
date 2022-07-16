import pandas as pd
import copy

def query(G, edges_df, nodes_df, query, depth):
    user_query = list(query.keys())[0]
    # For name queries (accounts for the user selecting multiple IDs)
    if user_query != "QUERY_ID":
        query_list = list(query.values())
        query_list = [item for sublist in query_list for item in sublist]
    # For ID queries
    else:
        query_list = list(query.values())

    # Recursively find target nodes
    # Go one depth further than selected depth to find any bidirectional connections of last layer
    targets = copy.deepcopy(query_list)
    full_df_cols = list(edges_df.columns)
    full_df_cols.append("depth")
    full_df = pd.DataFrame(columns = full_df_cols)
    for i in range(depth+1):
        query_df = edges_df[edges_df["source"].isin(targets)]
        if len(query_df.index) > 10000 and i>0:
            query_df = query_df[query_df["thickness"]>500]
        query_df["depth"] = i+1
        targets = query_df["target"].drop_duplicates().tolist()
        full_df = full_df.append(query_df)
    
    # Drop the nodes in the outermost layer that are not connected to any inner nodes
    # Leaves only connections from last layer nodes to inner layer nodes
    inner_edges = full_df[full_df["depth"] != depth+1]
    outer_edges = full_df[full_df["depth"] == depth+1]
    inner_source = full_df[full_df["depth"] == depth]["source"].drop_duplicates().tolist()
    outer_edges = outer_edges[outer_edges["target"].isin(inner_source)]
    full_df = inner_edges.append(outer_edges)

    # Make query nodes depth = 0, so they're in the center of the visualization
    df_dict = {"target":query_list, "depth":[0]*len(query_list)}
    zero_rows = pd.DataFrame.from_dict(df_dict)
    # Only need the targets, since every node (except for some query nodes) are a target at least once
    unique_targets = full_df[["target", "depth"]].drop_duplicates(subset = "target", keep="first")
    nodes = zero_rows.append(unique_targets)
    nodes = nodes.drop_duplicates(subset = "target", keep="first") 
    nodes_df = nodes_df.drop_duplicates(subset='Id', keep="first")
    nodes = pd.merge(nodes, nodes_df, left_on = "target", right_on = "Id", how="inner")
    
    # If the user selected multiple IDs, merge them all into one node
    if user_query != "QUERY_ID":
        nodes.loc[(nodes.Id.isin(query_list)),'Label']=user_query
        nodes.loc[(nodes.Id.isin(query_list)),'Id']=user_query
        nodes = nodes.drop_duplicates(subset="Id")
        full_df["orig_source"] = full_df["source"]
        full_df["orig_target"] = full_df["target"]
        full_df.loc[(full_df.source.isin(query_list)),'source']=user_query
        full_df.loc[(full_df.target.isin(query_list)),'target']=user_query
    else:
        full_df["orig_source"] = full_df["source"]
        full_df["orig_target"] = full_df["target"]
        
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
    
    nodes = nodes[["Id", "Label", "depth"]]
    full_df = full_df[["color", "thickness", 
                     "files", "source", "target"]]

    full_df.to_csv("query_edges.csv", index=False)
    nodes.to_csv("query_nodes.csv", index=False)