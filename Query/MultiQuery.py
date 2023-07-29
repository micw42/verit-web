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
import pickle

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
def query(G, edges_df, nodes_df, queries_id, max_linkers, qtype, query_type, get_direct_linkers, db_df,
          access_key, secret_key, bucket="all-abstract-ev",
          parallel_threshold=40, bg_edges=None, min_thickness=20):
    
    # QA steps
    nodes_df = nodes_df.drop_duplicates(subset='Id', keep="first")
    edges_df = edges_df.drop_duplicates(subset=['source', 'target'], keep="first")

    if (query_type == "name") | (query_type == "gene"):
        query_list = list(queries_id.values())    # Expecting queries_id to be a dictionary. Each element is a list.    
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

    # Parallel processing if large number of queries
    if len(q_combinations) >= parallel_threshold:
        start = time.time()
        source_targets = Parallel(n_jobs=4, require="sharedmem", verbose=10)(
            delayed(run_nx)(
                pair_chunk, G, qtype, max_linkers) for pair_chunk in np.array_split(
                np.array(q_combinations), len(q_combinations)))

        sources = [x[0] for x in source_targets]
        targets = [x[1] for x in source_targets]

        sources = list(it.chain.from_iterable(sources))
        targets = list(it.chain.from_iterable(targets))

    elif len(q_combinations) < parallel_threshold:
        sources, targets = run_nx(q_combinations, G, qtype, max_linkers)

    #### ####  #### ####

    # Construct qedges_df predecessor
    st_dict = {"source": sources, "target": targets}
    qedges_df = pd.DataFrame(st_dict).drop_duplicates()

    # Get reverse direction also (not all are expected to exist in edges_df)
    qedges_df = pd.concat([qedges_df, qedges_df.rename(columns={"source": "target", "target": "source"})])

    # Get autoregulation also (not all expected to exist in edges_df)
    autoreg_df = pd.DataFrame.from_dict({
        "source": query_list,
        "target": query_list
    })
    qedges_df = pd.concat([qedges_df, autoreg_df])

    # only now retrieve the edges that exist and their attributes
    qedges_df = qedges_df.merge(edges_df, how="inner")

    # Find direct neighbors to the query list
    # This will include duplicates to the previous steps
    found_ids = set(pd.concat([qedges_df.source, qedges_df.target])) - set(query_list)
    links = edges_df[(edges_df["source"].isin(query_list)) |
                     (edges_df["target"].isin(query_list))]

    # Add direct linkers to output edges
    qedges_df = pd.concat([qedges_df, links])

    # We expect duplicates to exist
    qedges_df = qedges_df.drop_duplicates()

    # Construct files column
    qedges_df["files"] = qedges_df.apply(lambda x: f"{x.source}_{x.target}.txt", axis=1)

    # Retrieve nodes
    qnodes_df = pd.DataFrame(pd.concat([qedges_df.source, qedges_df.target]).drop_duplicates(), columns=["Id"])
    qnodes_df = qnodes_df.merge(nodes_df, on = "Id", how = "inner")

    # Retrieve node synonyms
    qnodes_df = qnodes_df.merge(db_df, left_on="Id", right_on="id", how="left")
    qnodes_df["name"] = qnodes_df["name"].fillna(qnodes_df["Label"])
    qnodes_df["id"] = qnodes_df["id"].fillna(qnodes_df["Id"])

    # Aggregate synonyms into one column
    syn_concat = lambda x: "%%".join(x)  # Separate each synonym with %%
    aggregation_functions = {'Id': 'first', 'Label':"first", "KB":"first", "name":syn_concat}
    qnodes_df = qnodes_df.groupby('id').aggregate(aggregation_functions)


    #### Unchanged ####
    # Fix the node labels to account for combined IDs (can ignore)
    if query_type == "name":
        # Make df with user queries and corresponding IDs
        name_df = pd.DataFrame([(key, var) for (key, L) in queries_id.items() for var in L],
                 columns=['key', 'variable'])
        # Fix edges table
        src = qedges_df[["source"]]
        merged_src = pd.merge(src, name_df, how="left", left_on="source", right_on = "variable")["key"].tolist()
        tar = qedges_df[["target"]]
        merged_tar = pd.merge(tar, name_df, how="left", left_on="target", right_on = "variable")["key"].tolist()

        qedges_df["merged_source"] = merged_src
        qedges_df["merged_target"] = merged_tar
        qedges_df.merged_source.fillna(qedges_df.source, inplace=True)
        qedges_df.merged_target.fillna(qedges_df.target, inplace=True)
        rename_dict1 = {"source":"orig_source", "target":"orig_target"}
        rename_dict2 = {"merged_source":"source", "merged_target":"target"}
        qedges_df = qedges_df.rename(columns=rename_dict1).rename(columns=rename_dict2)

        # Fix nodes table
        id_converted = pd.merge(qnodes_df, name_df, how="left", left_on="Id", right_on="variable")["key"]
        qnodes_df["Id2"] = id_converted.tolist()
        qnodes_df["Label2"] = id_converted.tolist()
        qnodes_df.Id2.fillna(qnodes_df.Id, inplace=True)
        qnodes_df.Label2.fillna(qnodes_df.Label, inplace=True)
        qnodes_df = nodes.drop(labels = ["Label", "Id"], axis=1)
        rename_dict = {"Id2":"Id",
                      "Label2":"Label"}
        qnodes_df = qnodes_df.rename(columns= rename_dict)

        # Combine synonyms again for nodes that were just fixed
        syn_concat = lambda x: "%%".join(x)
        aggregation_functions = {"Label":"first", "KB":"first", 'name': syn_concat}
        qnodes_df = qnodes_df.groupby("Id").aggregate(aggregation_functions).reset_index()

        # For edges connecting nodes that correspond to multiple IDs,
        # take average of all color values (weighted by thickness)
        qedges_df = qedges_df[["color", "thickness", "source", "target", "orig_source", "orig_target", "files"]]
        rel_df["color2"] = qedges_df["color"] * qedges_df["thickness"]
        id_concat = lambda x: "%%".join(x) # Concat all source and target IDs of merged nodes
        aggregation_functions = {'color2': 'sum','thickness': 'sum', "orig_source":id_concat, "orig_target":id_concat, "files":id_concat}
        qedges_df = rel_df.groupby(["source", "target"]).aggregate(aggregation_functions).reset_index()
        qedges_df["color"] = qedges_df["color2"]/qedges_df["thickness"]

        name_df["variable"] = name_df.groupby(['key'])['variable'].transform(lambda x: ', '.join(x))
        name_df = name_df.drop_duplicates()
        source_ids = pd.merge(qedges_df, name_df, how="left", left_on="source", right_on = "key")["variable"].tolist()
        target_ids = pd.merge(qedges_df, name_df, how="left", left_on="target", right_on = "key")["variable"].tolist()
        qedges_df["source_id"] = source_ids
        qedges_df["target_id"] = target_ids
        qedges_df.source_id.fillna(rel_df.source, inplace=True)
        qedges_df.target_id.fillna(rel_df.target, inplace=True)
        qedges_df = qedges_df[["color", "thickness",
                         "files", "source", "target", "source_id", "target_id"]]
        display_ids = pd.merge(qnodes_df, name_df, how="left", left_on="Id", right_on="key")["variable"].tolist()
        qnodes_df["display_id"] = display_ids
        qnodes_df.display_id.fillna(qnodes_df.Id, inplace=True)
        qnodes_df = qnodes_df[["Id", "Label", "KB", "name", "display_id"]]

    else:
        # New columns for source and target (Not sure what it does)
        qedges_df["source_id"] = qedges_df["source"]
        qedges_df["target_id"] = qedges_df["target"]

        qnodes_df["display_id"] = qnodes_df["Id"]


    # Indicate whether each node is query (part of user query list),
    # linker (found by Netx algorithm), or direct (direct connection to query node not found by Netx)
    # For visualization
    if (query_type == "name") :
        qnodes_df["Type"] = ["Query" if x in queries_id.keys() 
                         else ("Linker" if x in found_ids 
                               else "Direct") 
                         for x in qnodes_df['Id']]

    elif (query_type == "id") | (query_type == "gene"):
        qnodes_df["Type"] = ["Query" if x in query_list 
                         else ("Linker" if x in found_ids 
                               else "Direct") 
                         for x in qnodes_df['Id']]


    # nodes["Label"] = nodes["Label"].str.replace("SPACE", " ")
    # nodes = nodes[["Id", "Label", "KB", "name", "Type", "display_id"]]

    #### Unchanged ####
    nodes_cleaned = qnodes_df.copy()
    edges_cleaned = qedges_df.copy()

    nodes_cleaned['name'] = nodes_cleaned['name'].str.replace('%%',', ')
    nodes_cleaned = nodes_cleaned.rename(columns={"name": "Synonyms"})

    edges_cleaned = edges_cleaned.merge(nodes_cleaned, left_on="source", right_on="Id").drop(labels="Id", axis=1).rename(columns={"Label":"source_name"})
    edges_cleaned = edges_cleaned.merge(nodes_cleaned, left_on="target", right_on="Id").drop(labels="Id", axis=1).rename(columns={"Label":"target_name"})
    edges_cleaned = edges_cleaned.rename(columns={"color":"score", "thickness":"evidence_count"})
    edges_cleaned = edges_cleaned[["source_id", "source_name", "target_id", "target_name", "evidence_count", "score"]]

    nodes_cleaned = nodes_cleaned.drop(labels="Id", axis=1).rename(columns={"display_id":"Id"})
    nodes_cleaned = nodes_cleaned[["Id", "Label", "KB", "Synonyms", "Type"]]
    
    return qnodes_df, qedges_df, nodes_cleaned, edges_cleaned    


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

    # Unidirectional
    st_df = st_df.merge(
        edges_df, left_on=["source_id", "target_id"], right_on=["source", "target"]
    )[["source_id", "target_id", "thickness"]]

    # auto-regulation
    auto_edges = edges_df[(edges_df["source"].isin(query_list)) & (edges_df["target"].isin(query_list)) & (edges_df["source"]==edges_df["target"])]
    auto_edges=auto_edges.rename(columns={"source":"source_id", "target":"target_id"})
    st_df = pd.concat([st_df, auto_edges]).drop_duplicates(subset=["source_id", "target_id"])
    
    ## Construct nodes dataframe
    # Initialize
    qnode_ids = list(set(pd.concat([st_df.source_id, st_df.target_id])))

    # Direct linkers
    links = edges_df[((edges_df["source"].isin(query_list)) &
                      ~(edges_df["target"].isin(set(qnode_ids) - set(query_list)))) |
                     ((edges_df["target"].isin(query_list)) &
                      ~(edges_df["source"].isin(set(qnode_ids) - set(query_list))))]
    print("links columns:", links.columns)
    links.columns = ["source_id", "target_id", "thickness"]

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

    qedges_df["files"] = qedges_df["source_id"] + "_" + qedges_df["target_id"] + ".txt"

    # Dummy columns for consistency.
    qedges_df["color"] = 0
    
    qedges_df = qedges_df.drop_duplicates()

    qedges_df.to_csv("./query_edges_BIOGRID.csv", index=False)
    
    return qnodes_df, qedges_df
    