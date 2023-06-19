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
from memory_profiler import profile, LogFile


def run_nx(query_pairs, G, qtype, max_linkers):
    sources = []; targets = []
    for query_pair in query_pairs:
        source, target = query_pair
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

        elif qtype == "all_shortest_paths" and max_linkers < 3:
            try:
                path = list(nx.all_simple_paths(G, source, target, cutoff=max_linkers))
                if len(path) > 0:
                    path_lengths = [len(x) for x in path]
                    shortest_length = min(path_lengths)
                    for ind in path:
                        if len(ind)==shortest_length:
                            for n1, n2 in zip(ind, ind[1:]):
                                sources.append(n1)
                                targets.append(n2)

            except (nx.NetworkXNoPath, nx.NodeNotFound) as e:
                pass

        elif qtype == "all_shortest_paths" and max_linkers >= 3:
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


# queries_id: dict of format {query_name:[list of selected IDs]} (if name query)
# or {"QUERY_ID":"comma separated string of query IDs"} (if ID query)

# query_type: name or ID

#qtype: all_simple_paths or all_shortest_paths
fp=open('memory_profiler_RB_E2F.log','w')
@profile(stream=fp)
def query(G, edges_df, nodes_df, queries_id, max_linkers, qtype, query_type, get_direct_linkers, db_df, access_key, secret_key, bucket="all-abstract-ev", parallel_threshold=40):
    nodes_df = nodes_df.drop_duplicates(subset='Id', keep="first")
    edges_df = edges_df.drop_duplicates(subset=['source', 'target'], keep="first")
    if query_type == "name":
        query_list = list(queries_id.values())    # Expecting queries_id to be a dictionary. Each element is a list.
        print("Len query_list:", len(query_list))
        
        # Get all possible pairs of queries (note: each query is a list of IDs)
        li_perm = list(it.permutations(query_list, 2))
        # Get all possible pairs of element from 1st list with element from 2nd list
        q_combinations = [list(it.product(*sub_li)) for sub_li in li_perm]
        # Unroll to list of tuples
        q_combinations = [item for sublist in q_combinations for item in sublist]
        # Unroll list of query IDs
        query_list = [item for sublist in query_list for item in sublist]
        print("Got query list")

    elif query_type == "id":
        query_list = queries_id["QUERY_ID"].split(",")
        q_combinations = list(it.permutations(query_list, 2))

    if len(q_combinations) >= parallel_threshold:
        start = time.time()
        source_targets = Parallel(n_jobs=4, require="sharedmem", verbose=10)(delayed(run_nx)(pair_chunk, G, qtype, max_linkers)
                                            for pair_chunk in np.array_split(np.array(q_combinations), len(q_combinations)))

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
    nodes = nodes.merge(nodes_df, on="Id", how="inner")[["Id", "Label", "KB"]]

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
    full_nodes = full_nodes[["Id", "Label", "KB"]]

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
    nodes["name"] = nodes["name"].fillna("NAN")
    syn_concat = lambda x: "%%".join(x)  # Separate each synonym with %%
    aggregation_functions = {'Id': 'first', 'Label':"first", "KB":"first", "name":syn_concat}
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
    aggregation_functions = {"Label":"first", "KB":"first", 'name': syn_concat}
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
        
    name_df["variable"] = name_df.groupby(['key'])['variable'].transform(lambda x: ', '.join(x))
    name_df = name_df.drop_duplicates()
    source_ids = pd.merge(rel_df, name_df, how="left", left_on="source", right_on = "key")["variable"].tolist()
    target_ids = pd.merge(rel_df, name_df, how="left", left_on="target", right_on = "key")["variable"].tolist()
    rel_df["source_id"] = source_ids
    rel_df["target_id"] = target_ids
    rel_df.source_id.fillna(rel_df.source, inplace=True)
    rel_df.target_id.fillna(rel_df.target, inplace=True)
    rel_df = rel_df[["color", "thickness",
                     "files", "source", "target", "source_id", "target_id"]]
    
    display_ids = pd.merge(nodes, name_df, how="left", left_on="Id", right_on="key")["variable"].tolist()
    nodes["display_id"] = display_ids
    nodes.display_id.fillna(nodes.Id, inplace=True)
    nodes = nodes[["Id", "Label", "KB", "name", "display_id"]]
    
    # Indicate whether each node is query (part of user query list),
    # linker (found by Netx algorithm), or direct (direct connection to query node not found by Netx)
    # For visualization
    if query_type == "name":
        nodes["Type"] = ["Query" if x in queries_id.keys() else ("Linker" if x in found_ids else "Direct") for x in nodes['Id']]
    elif query_type == "id":
        nodes["Type"] = ["Query" if x in query_list else ("Linker" if x in found_ids else "Direct") for x in nodes['Id']]

    nodes["Label"] = nodes["Label"].str.replace("SPACE", " ")
    nodes = nodes[["Id", "Label", "KB", "name", "Type", "display_id"]]

    # Remove recursive edges
    rel_df = rel_df[rel_df["source"] != rel_df["target"]]

    # Write edges and nodes
    nodes.to_csv("query_nodes.csv", index=False)
    rel_df.to_csv("query_edges.csv", index=False)
    
    nodes['name'] = nodes['name'].str.replace('%%',', ')
    nodes = nodes.rename(columns={"name": "Synonyms"})
    rel_df = rel_df.merge(nodes, left_on="source", right_on="Id").drop(labels="Id", axis=1).rename(columns={"Label":"source_name"})
    rel_df = rel_df.merge(nodes, left_on="target", right_on="Id").drop(labels="Id", axis=1).rename(columns={"Label":"target_name"})
    rel_df = rel_df.rename(columns={"color":"score", "thickness":"evidence_count"})
    rel_df = rel_df[["source_id", "source_name", "target_id", "target_name", "evidence_count", "score"]]
    nodes = nodes.drop(labels="Id", axis=1).rename(columns={"display_id":"Id"})
    nodes = nodes[["Id", "Label", "KB", "Synonyms", "Type"]]
    
    nodes.to_csv("query_nodes_cleaned.csv", index=False)
    rel_df.to_csv("query_edges_cleaned.csv", index=False)


def BIOGRID_query(G, edges_df, nodes_df, queries_id,
          max_linkers, qtype, query_type, get_direct_linkers,
          db_df, parallel_threshold=40):
    
    ## Permutations of queries
    query_list = list(queries_id.values())    # Expecting queries_id to be a dictionary. Each element is a list.
    # Get all possible pairs of queries (note: each query is a list of IDs)
    li_perm = list(it.permutations(query_list, 2))
    # Get all possible pairs of element from 1st list with element from 2nd list
    q_combinations = [list(it.product(*sub_li)) for sub_li in li_perm]
    # Unroll to list of tuples
    q_combinations = [item for sublist in q_combinations for item in sublist]
    # Unroll list of query IDs
    query_list = [item for sublist in query_list for item in sublist]


    ## Run network querying
    if len(q_combinations) >= parallel_threshold:
        start = time.time()

        source_targets = Parallel(n_jobs=4, require="sharedmem", verbose=10)(
            delayed(run_nx)(
                pair_chunk, G, qtype, max_linkers) for pair_chunk in np.array_split(
                np.array(q_combinations), len(q_combinations)))

        print("Finished querying in", time.time()-start, "sec.")

        sources = [x[0] for x in source_targets]
        targets = [x[1] for x in source_targets]

        sources = list(it.chain.from_iterable(sources))
        targets = list(it.chain.from_iterable(targets))

    elif len(q_combinations) < parallel_threshold:
        sources, targets = run_nx(q_combinations, G, qtype, max_linkers)

    st_dict = {"source_id": sources, "target_id": targets}
    st_df = pd.DataFrame(st_dict).drop_duplicates()


    ## Construct nodes dataframe
    # Initialize
    qnode_ids = list(set(pd.concat([st_df.source_id, st_df.target_id])))

    # Direct linkers
    links = edges_df[((edges_df["source"].isin(query_list)) &
                      ~(edges_df["target"].isin(set(qnode_ids) - set(query_list)))) |
                     ((edges_df["target"].isin(query_list)) &
                      ~(edges_df["source"].isin(set(qnode_ids) - set(query_list))))]
    links.columns = ["source_id", "target_id"]

    # Keep direct-linkers to keep track
    direct = list(set(pd.concat([links.source_id, links.target_id])))
    qnode_ids.extend(direct)
    st_df = pd.concat([st_df, links]).reset_index(drop=True)

    qnodes_df = pd.DataFrame.from_dict({"Id": qnode_ids})
    qnodes_df = qnodes_df.merge(
        nodes_df, on="Id").drop_duplicates(
        subset="Id").rename(
        columns={"name": "Label"}).reset_index(
        drop=True)

    # Dummy knowledgebase column
    qnodes_df["KB"] = "BIOGRID"

    # Distinguish which nodes were which
    ntype = pd.Series(["Linker"] * len(qnodes_df))
    ntype[qnodes_df["Id"].isin(direct)] = "Direct"
    ntype[qnodes_df["Id"].isin(query_list)] = "Query"
    qnodes_df["Type"] = ntype

    # Get the synonyms
    qnodes_syn = nodes_df.merge(qnodes_df, on=["Id"])
    qnodes_syn = qnodes_syn.groupby("Id").agg({'name': lambda x: "%%".join(x)}).reset_index()
    qnodes_syn.columns = ["Id", "name"]
    qnodes_df = qnodes_df.merge(qnodes_syn, on="Id")
    qnodes_df["display_id"] = qnodes_df["Id"].copy()    # Equal because BG only used with gene id query

    qnodes_df.to_csv("./query_nodes_BIOGRID.csv", index=False)


    ## Construct edges dataframe
    qedges_df = st_df.merge(qnodes_df[["Id", "Label"]], left_on="source_id", right_on="Id").drop(
        columns="Id"
    ).rename(columns={"Label": "source"})

    qedges_df = qedges_df.merge(qnodes_df[["Id", "Label"]], left_on="target_id", right_on="Id").drop(
        columns="Id"
    ).rename(columns={"Label": "target"})
    
    qedges_df["source"] = qedges_df["source_id"]
    qedges_df["target"] = qedges_df["target_id"]

    # Dummy columns for consistency. Files are None
    qedges_df["color"] = 0
    qedges_df["thickness"] = 1e6
    qedges_df["files"] = None

    qedges_df.to_csv("./query_edges_BIOGRID.csv", index=False)
    