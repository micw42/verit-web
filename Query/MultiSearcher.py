import pandas as pd
import numpy as np
import pickle


def name_query(name, nodes, full_df, string_type):
    if string_type == "starts_with":
        name_query = full_df.query('name.str.startswith(@name, na=False)', engine='python').reset_index(drop=True)

    elif string_type == "ends_with":
        name_query = full_df.query('name.str.endswith(@name, na=False)', engine='python').reset_index(drop=True)

    elif string_type == "contains":
        name_query = full_df.query('name.str.contains(@name, na=False)', engine='python').reset_index(drop=True)

    # The following nodes are in the network, but REACH likely grounds to more
    in_net = pd.merge(nodes, name_query, left_on="Id", right_on="id", how="inner").drop_duplicates().reset_index(drop=True)
    in_net.to_csv("in_net.csv", index=False)

    # Reverse search to collect all aliases
    u_ids = list(in_net["id"].drop_duplicates())
    in_net = full_df.query("id == @u_ids").iloc[::-1].reset_index(drop=True)

    in_net["user_query"] = name
    in_net = pd.merge(in_net, nodes, left_on="id", right_on = "Id", how="inner")
    in_net = in_net.drop_duplicates(subset = ["name", "id"])

    return in_net


def multi_query(query_list, nodes, full_df, bg_nodes, string_type):
    # Gene querying can be streamlined
    if string_type == "gene":
        # Get found and unfound queries
        # full_df in this case is uniprot_df. uniprot_df is has 1-1 correspondence between Id and Label at this stage.
        full_df["Label"] = full_df["Label"].str.upper()
        mapped_ids = full_df[full_df["Label"].isin(query_list)].drop_duplicates(subset=["Id", "Label"])
        
        print(bg_nodes.columns)
        
        # Only line pertaining to BIOGRID network for searching.
        bg_in_net = bg_nodes[bg_nodes["name"].isin(query_list)].drop_duplicates(subset=["Id"])
        bg_in_net["user_query"] = bg_in_net["name"]
        bg_in_net = bg_in_net.drop(columns="name").set_index("user_query")
        bg_in_net = bg_in_net.T.to_dict(orient="list")

        unmapped = np.setdiff1d(query_list, list(full_df["Label"]))    # Unused

        # Get only IDs and PR values
        nodes_dd = nodes[["Id", "PR"]].drop_duplicates()

        # Unique IDs within VERIT network
        in_net = nodes_dd.merge(mapped_ids, on="Id", how="inner")

        # Reconstruct expected format
        in_net.columns = ["id", "PR", "name"]
        in_net["user_query"] = in_net["name"]

        return in_net, bg_in_net

    else:
        # initialize empty dataframe
        full_query = pd.DataFrame()
        for query in query_list:
            full_query = full_query.append(name_query(query, nodes, full_df, string_type))
        return full_query, None


def query(query_list, nodes, full_df, uniprot_df, bg_nodes, string_type):
    if string_type == "gene":
        query_list = [x.upper() for x in query_list]
        result = multi_query(query_list, nodes, uniprot_df, bg_nodes, string_type)
    else:
        # Turns all entries into lower case
        case_proof = [x.lower() for x in query_list]
        result = multi_query(case_proof, nodes, full_df, bg_nodes, string_type)
    return result
