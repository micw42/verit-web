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
        #Concatenate ids of compound nodes
        name_df["combined_id"] = name_df.groupby(['key'])['variable'].transform(lambda x: ', '.join(x))
        name_df = name_df.drop_duplicates()
        #For concatenating rows of compound nodes/edges
        concat = lambda x: "%%".join(x) 

        #Fix edges table
        qedges_df = qedges_df.rename(columns={"source":"orig_source", "target":"orig_target"})
        #source_id, target_id used for display
        #For multi-ID nodes, source_id, target_id correspond to concatenated IDs (i.e. "ID_1, ID_2")
        #For everything else, source_id, target_id are equal to ID
        qedges_df = pd.merge(qedges_df, name_df, how="left", left_on="orig_source", right_on = "variable").rename(columns={"key":"source", "combined_id":"source_id"}).drop(columns="variable")
        qedges_df = pd.merge(qedges_df, name_df, how="left", left_on="orig_target", right_on = "variable").rename(columns={"key":"target", "combined_id":"target_id"}).drop(columns="variable")
        qedges_df.source.fillna(qedges_df.orig_source, inplace=True)   #Replace IDs of non-query nodes
        qedges_df.target.fillna(qedges_df.orig_target, inplace=True)
        qedges_df.source_id.fillna(qedges_df.orig_source, inplace=True)
        qedges_df.target_id.fillna(qedges_df.orig_target, inplace=True)
        qedges_df["color2"] = qedges_df["color"] * qedges_df["thickness"]   #for weighted avg
        aggregation_functions = {'color2': 'sum','thickness': 'sum', "files":concat, "source_id":"first", "target_id":"first"}
        qedges_df = qedges_df.groupby(["source", "target"]).aggregate(aggregation_functions).reset_index() #Group compound edges
        qedges_df["color"] = qedges_df["color2"]/qedges_df["thickness"]    #Weighted avg of color for compound edges
        qedges_df.drop(columns="color2", inplace=True)

        # Fix nodes table
        qnodes_df = pd.merge(qnodes_df, name_df, how="left", left_on="Id", right_on="variable").drop(columns="variable").rename(columns={"key":"Id2", "combined_id":"display_id"})
        qnodes_df["Label2"] = qnodes_df["Id2"]     #Label query nodes as what the user queried
        qnodes_df.Id2.fillna(qnodes_df.Id, inplace=True)  #Replace non-query ids, labels
        qnodes_df.Label2.fillna(qnodes_df.Label, inplace=True)
        qnodes_df.display_id.fillna(qnodes_df.Id, inplace=True)
        qnodes_df = qnodes_df.drop(labels = ["Label", "Id"], axis=1).rename(columns= {"Id2":"Id", "Label2":"Label"})
        aggregation_functions = {"Label":"first", "KB":"first", 'name':concat, "display_id":"first"}
        qnodes_df = qnodes_df.groupby("Id").aggregate(aggregation_functions).reset_index()  #Combine multi-ID nodes into one node
    else:
        # Dummy columns for compatibility with visualization
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

    #Add labels to edges table For export as tsv
    qnodes_df["Label"] = qnodes_df["Label"].str.replace("SPACE", " ")
    qedges_df = pd.merge(qedges_df, qnodes_df, left_on="source", right_on="Id", how="left")\
    .rename(columns={"Label":"source_lab"}).drop(columns=["Id", "KB", "name", "display_id"])
    qedges_df = pd.merge(qedges_df, qnodes_df, left_on="target", right_on="Id", how="left")\
    .rename(columns={"Label":"target_lab"}).drop(columns=["Id", "KB", "name", "display_id"])
    return qnodes_df, qedges_df


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
    qedges_df["source_lab"] = qedges_df.merge(qnodes_df, left_on="source", right_on="Id", how="left")["Label"].tolist()
    qedges_df["target_lab"] = qedges_df.merge(qnodes_df, left_on="target", right_on="Id", how="left")["Label"].tolist()
    
    qedges_df = qedges_df.drop_duplicates()

    qedges_df.to_csv("./query_edges_BIOGRID.csv", index=False)
    
    return qnodes_df, qedges_df
    
