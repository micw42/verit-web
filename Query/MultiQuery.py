import networkx as nx
import pandas as pd
import math
import itertools as it
import numpy as np
import copy
import boto3
import time
from multiprocessing import Pool
from Query import GetEv
from functools import partial
from joblib import Parallel, delayed


def run_nx(query_pairs, G, qtype, max_linkers):
    sources = []; targets = []
    for query_pair in query_pairs:
        source, target = query_pair
        print("Running", query_pair)
        if qtype == "all_simple_paths":
            try:
                path = list(nx.all_simple_paths(G, source, target, cutoff=max_linkers))
                # Loop through interaction pairs in a path
                for ind in path:
                    for n1, n2 in zip(ind, ind[1:]):
                        sources.append(n1)
                        targets.append(n2)

            except (nx.NetworkXNoPath, nx.NodeNotFound) as e:
                pass

        elif qtype == "all_shortest_paths":
            try:
                path = list(nx.all_shortest_paths(G, source, target))
                path = [x for x in path if (len(x)-1)<=max_linkers]
                for ind in path:
                    for n1, n2 in zip(ind, ind[1:]):
                        sources.append(n1)
                        targets.append(n2)

            except (nx.NetworkXNoPath, nx.NodeNotFound) as e:
                pass

    return sources, targets

# Get number of valid (non-direct) paths (for shortQuery allShortestPaths)
def get_valid_paths(path_list, direct_paths):
    valid = list(set(path_list) - set(direct_paths))
    valid = [x for x in valid if len(set(x))>1]
    n_valid = len(valid)
    return n_valid

def shortQuery(edges_df, query_list, q_combinations, q_type, cutoff):
    # All edges to and from query nodes
    subset_df = edges_df[(edges_df["source"].isin(query_list)) | (edges_df["target"].isin(query_list))]
    # Direct connections between query nodes
    direct_edges = subset_df[(subset_df["source"].isin(query_list)) & (subset_df["target"].isin(query_list))]
    
    # Return just direct connections if max linkers is 1
    if cutoff==1:
        return direct_edges
    
    direct_path_list = list(direct_edges[["source", "target"]].itertuples(index=False, name=None))
    # List of direct paths that don't exist
    no_direct = list(set(q_combinations) - set(direct_path_list))
    # List of query nodes with no direct path
    no_direct = list(sum(no_direct, ()))

    # Edges between query node and non-query node
    paths_df = subset_df[(~subset_df["source"].isin(query_list)) | (~subset_df["target"].isin(query_list))]
    # List of nodes that target each non-query node
    target_df = paths_df[["source", "target"]][~paths_df["target"].isin(query_list)].groupby('target')['source'].apply(list).reset_index(name='sources')
    # List of nodes that each non-query node targets
    source_df = paths_df[["source", "target"]][~paths_df["source"].isin(query_list)].groupby('source')['target'].apply(list).reset_index(name='targets')
    merged_df = target_df.merge(source_df, left_on="target", right_on="source")
    # Number of sources and targets each non-query node has
    merged_df["source_len"] = merged_df["sources"].apply(lambda x: len(x))
    merged_df["target_len"] = merged_df["targets"].apply(lambda x: len(x))
    # Paths that each non-query node acts as a linker in 
    merged_df["paths"] = merged_df.apply(lambda x: list(it.product(x.sources, x.targets)), axis=1)
    # Number of paths that aren't already a direct connection that each non-query node acts as a linker in
    merged_df["paths"] = merged_df["paths"].apply(lambda x: get_valid_paths(x, direct_path_list))

    if q_type == "all_simple_paths":
        linkers = merged_df[((merged_df["source_len"]>1) & (merged_df["target_len"]>0)) | ((merged_df["target_len"]>1) & (merged_df["source_len"]>0)) | (merged_df["sources"]!=merged_df["targets"])]["target"].tolist()
        link_edges = subset_df[((subset_df["source"].isin(linkers)) | (subset_df["target"].isin(linkers)))]

    elif q_type == "all_shortest_paths":
        # Only keep linker paths between query nodes w/o direct connection
        merged_df = merged_df[merged_df["paths"]!=0]  
        linkers = merged_df[((merged_df["source_len"]>1) & (merged_df["target_len"]>0)) | ((merged_df["target_len"]>1) & (merged_df["source_len"]>0)) | (merged_df["sources"]!=merged_df["targets"])]["target"].tolist()
        link_edges = subset_df[((subset_df["source"].isin(linkers)) | (subset_df["target"].isin(linkers))) & (subset_df["source"].isin(no_direct) | (subset_df["target"].isin(no_direct)))]

    result_df = pd.concat([link_edges, direct_edges])
    return result_df



# queries_id: dict of format {query_name:[list of selected IDs]} (if name query)
# or {"QUERY_ID":"comma separated string of query IDs"} (if ID query)

# query_type: name or ID

#qtype: all_simple_paths or all_shortest_paths

def query(G, edges_df, nodes_df, queries_id, max_linkers, qtype, query_type, get_direct_linkers, db_df, access_key, secret_key, bucket="all-abstract-ev", parallel_threshold=40):   
    nodes_df = nodes_df.drop_duplicates(subset='Id', keep="first")
    edges_df = edges_df.drop_duplicates(subset=['source', 'target'], keep="first")
    
    if query_type == "name":
        query_list = list(queries_id.values())
        query_lengths = set([len(x) for x in query_list])
        # Get all possible pairs of queries (note: each query is a list of IDs)
        li_perm = list(it.permutations(query_list, 2)) 
        # Get all possible pairs of element from 1st list with element from 2nd list
        q_combinations = [list(it.product(*sub_li)) for sub_li in li_perm] 
        # Unroll to list of tuples
        q_combinations = [item for sublist in q_combinations for item in sublist]
        # Unroll list of query IDs
        query_list = [item for sublist in query_list for item in sublist]
        
    elif query_type == "id":
        query_list = queries_id["QUERY_ID"].split(",")
        q_combinations = list(it.permutations(query_list, 2))
        query_lengths = {1}
        
    
    if max_linkers < 3 and query_lengths=={1}:
        # Pandas path finding
        print("Running short query")
        rel_df = shortQuery(edges_df, query_list, q_combinations, qtype, cutoff=max_linkers)
        
    else:
        # Networkx path finding
        # Emperically tested threshold for parallelization
        print("Running nx query")
        if len(q_combinations) >= parallel_threshold:
            start = time.time()
            source_targets = Parallel(n_jobs=4, require="sharedmem", verbose=10)(delayed(run_nx)(pair_chunk, G, qtype, max_linkers)
                                                for pair_chunk in np.array_split(np.array(q_combinations), len(q_combinations)))
            print("Finished querying in", time.time()-start, "sec.")

            sources = [x[0] for x in source_targets]
            targets = [x[1] for x in source_targets]

            sources = list(it.chain.from_iterable(sources))
            targets = list(it.chain.from_iterable(targets))

        elif len(q_combinations) < parallel_threshold:
            sources, targets = run_nx(q_combinations, G, qtype, max_linkers)

        st_dict = {"source": sources, "target": targets}
        st_df = pd.DataFrame(st_dict).drop_duplicates()

        # Making edges
        rel_df = st_df.merge(edges_df, on=["source", "target"], how="left")

        # Bidirectional edges
        opp_df = rel_df.merge(edges_df, left_on=["source", "target"], right_on=["target","source"])
        opp_df = opp_df.drop(labels=["source_x","target_x","color_x", "thickness_x"], axis=1).rename(columns={"source_y":"source", "target_y":"target", "color_y":"color", "thickness_y":"thickness"})
        rel_df = pd.concat([rel_df, opp_df]).drop_duplicates(subset=["source", "target"])
    
    # Create nodes df
    nodes = list(it.chain(*q_combinations)) # List of all query IDs
    # Add found nodes
    nodes.extend(rel_df["source"].tolist())
    nodes.extend(rel_df["target"].tolist())
    nodes = list(set(nodes))
    nodes = pd.DataFrame({"Id": nodes})
    nodes = nodes.merge(nodes_df, on="Id", how="inner")[["Id", "Label"]]
    
    # IDs that were found (not in original query list)
    found_ids = set(rel_df["source"].tolist()) | set(rel_df["target"].tolist()) - set(query_list)
    
    # Get all direct connections to query nodes that were not already found by query 
    # (user has option to show them in the visualization)
    links = edges_df[((edges_df["source"].isin(query_list)) & ~(edges_df["target"].isin(found_ids))) | 
                     ((edges_df["target"].isin(query_list)) & ~(edges_df["source"].isin(found_ids)))]
    
    # Filter thickness of direct connections so visualization isn't crowded (arbitrarily selected threshold)
    links = links[links["thickness"] > 20]  

    # Nodes directly targeted by query nodes
    targets = links[(links["source"].isin(query_list)) & ~(links["target"].isin(query_list))]
    targets = targets[["target"]]
    # Nodes that directly target query nodes
    sources = links[(links["target"].isin(query_list)) & ~(links["source"].isin(query_list))]
    sources = sources[["source"]]
    targets = targets.rename(columns = {"target":"Id"})
    sources = sources.rename(columns = {"source":"Id"})
    full_nodes = pd.concat([targets, sources]).drop_duplicates(subset = ["Id"])
    full_nodes = full_nodes.merge(nodes_df, on = "Id", how = "inner")
    full_nodes = full_nodes[["Id", "Label"]]
    
    # Add direct connection nodes to the rest of the nodes
    nodes = pd.concat([nodes, full_nodes])

    # Add direct connection edges to the rest of the edges
    rel_df = pd.concat([rel_df, links]).drop_duplicates(subset = ["source", "target"])
    # Add evidence column
    sources = rel_df["source"].tolist()
    targets = rel_df["target"].tolist()
    rel_df = rel_df[["source", "target", "color", "thickness"]]
        
    # Add synonyms to the nodes
    nodes = nodes.merge(db_df, left_on="Id", right_on="id", how="left")
    nodes["name"] = nodes["name"].fillna(nodes["Label"])
    nodes["id"] = nodes["id"].fillna(nodes["Id"])
    syn_concat = lambda x: "%%".join(x)  # Separate each synonym with %%
    aggregation_functions = {'Id': 'first', 'Label':"first", "name":syn_concat}
    nodes = nodes.groupby('id').aggregate(aggregation_functions)
    
    # Fix the node labels to account for combined IDs (can ignore)
    if query_type == "name":
        # Make df with user queries and corresponding IDs
        name_df = pd.DataFrame([(key, var) for (key, L) in queries_id.items() for var in L], 
                 columns=['key', 'variable'])
        
        # Fix edges table
        src = rel_df[["source"]]
        merged_src = pd.merge(src, name_df, how="left", left_on="source", right_on = "variable")["key"].tolist()   
        tar = rel_df[["target"]]
        merged_tar = pd.merge(tar, name_df, how="left", left_on="target", right_on = "variable")["key"].tolist()
        rel_df["merged_source"] = merged_src
        rel_df["merged_target"] = merged_tar      
        rel_df.merged_source.fillna(rel_df.source, inplace=True)
        rel_df.merged_target.fillna(rel_df.target, inplace=True)
        rename_dict1 = {"source":"orig_source", "target":"orig_target"}
        rename_dict2 = {"merged_source":"source", "merged_target":"target"}
        rel_df = rel_df.rename(columns=rename_dict1).rename(columns=rename_dict2)
        
        # Fix nodes table
        id_converted = pd.merge(nodes, name_df, how="left", left_on="Id", right_on="variable")["key"]
        nodes["Id2"] = id_converted.tolist()
        nodes["Label2"] = id_converted.tolist()
        nodes.Id2.fillna(nodes.Id, inplace=True)
        nodes.Label2.fillna(nodes.Label, inplace=True)
        nodes = nodes.drop(labels = ["Label", "Id"], axis=1)
        rename_dict = {"Id2":"Id", 
                      "Label2":"Label"}
        nodes = nodes.rename(columns= rename_dict)
    else:
        rel_df["orig_source"] = rel_df["source"]
        rel_df["orig_target"] = rel_df["target"]
        
    # Combine synonyms again for nodes that were just fixed
    syn_concat = lambda x: "%%".join(x)
    aggregation_functions = {"Label":"first", 'name': syn_concat}
    nodes = nodes.groupby("Id").aggregate(aggregation_functions).reset_index()
  
    # For edges connecting nodes that correspond to multiple IDs, 
    # take average of all color values (weighted by thickness)
    
    rel_df = rel_df[["color", "thickness", "source", "target", "orig_source", "orig_target"]]
    rel_df["color2"] = rel_df["color"] * rel_df["thickness"]
    id_concat = lambda x: "%%".join(x) # Concat all source and target IDs of merged nodes
    aggregation_functions = {'color2': 'sum','thickness': 'sum', "orig_source":id_concat, "orig_target":id_concat}
    rel_df = rel_df.groupby(["source", "target"]).aggregate(aggregation_functions).reset_index()
    rel_df["color"] = rel_df["color2"]/rel_df["thickness"]
    def formatter(sources, targets):
        source_list = sources.split("%%")
        target_list = targets.split("%%")
        file_list = [f"{source_list[i]}_{target_list[i]}.txt" for i in range(len(source_list))]
        files = "%%".join(file_list)
        return files
    rel_df["files"] = rel_df.apply(lambda x: formatter(x.orig_source, x.orig_target), axis=1)
    
    #For consistency: source and target columns disappear if dataframe happens to be empty
    if len(rel_df.index)==0:
        rel_df = pd.DataFrame(columns=["color",
                                       "thickness", 
                                       "files"
                                       "source", "target"])
        
    rel_df = rel_df[["color", "thickness", 
                     "files", "source", "target"]]
    nodes = nodes[["Id", "Label", "name"]]
    
    # Indicate whether each node is query (part of user query list), 
    # linker (found by Netx algorithm), or direct (direct connection to query node not found by Netx)
    # For visualization 
    if query_type == "name":
        nodes["Type"] = ["Query" if x in queries_id.keys() else ("Linker" if x in found_ids else "Direct") for x in nodes['Id']]   
    elif query_type == "id":
        nodes["Type"] = ["Query" if x in query_list else ("Linker" if x in found_ids else "Direct") for x in nodes['Id']]
        
    nodes["Label"] = nodes["Label"].str.replace("SPACE", " ")
    
    # Remove recursive edges
    rel_df = rel_df[rel_df["source"] != rel_df["target"]]

    # Write edges and nodes
    nodes.to_csv("query_nodes.csv", index=False)
    rel_df.to_csv("query_edges.csv", index=False)