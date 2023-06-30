import pandas as pd

def name_query(name, nodes, full_df, string_type):
    if string_type == "starts_with":
        name_query = full_df.query('name.str.startswith(@name, na=False)', engine='python').reset_index(drop=True)
    
    elif string_type == "ends_with":
        name_query = full_df.query('name.str.endswith(@name, na=False)', engine='python').reset_index(drop=True)
    
    elif string_type == "contains":
        name_query = full_df.query('name.str.contains(@name, na=False)', engine='python').reset_index(drop=True)
    # The following nodes are in the network, but REACH likely grounds to more
    #in_net = pd.merge(nodes, name_query, left_on="Id", right_on="id", how="inner").iloc[:, 3:6].drop_duplicates().reset_index(drop=True)
    in_net = pd.merge(nodes, name_query, left_on="Id", right_on="id", how="inner").drop_duplicates().reset_index(drop=True)
    u_ids = list(in_net["id"].drop_duplicates())
    in_net = full_df.query("id == @u_ids").iloc[::-1].reset_index(drop=True)
    #in_net = pd.merge(nodes, in_net, left_on="Id", right_on="id", how="inner").iloc[:, 3:6].drop_duplicates().reset_index(drop=True)
    in_net = pd.merge(nodes, in_net, left_on="Id", right_on="id", how="inner").drop_duplicates().reset_index(drop=True)
    in_net = in_net.drop_duplicates(subset = ["name", "id"])
    in_net = in_net.sort_values(by = "PR", ascending = False)
    
    return in_net

def single_query(query, nodes, full_df):
    # Get found and unfound queries
    mapped_ids = full_df[full_df["Label"]==query].drop_duplicates(subset="Label")

    # Get only IDs and PR values
    nodes_dd = nodes[["Id", "PR"]].drop_duplicates()

    # Unique IDs within VERIT network
    in_net = nodes_dd.merge(mapped_ids, on="Id", how="inner")

    # Reconstruct expected format
    in_net.columns = ["id", "PR", "name"]
    in_net["user_query"] = in_net["name"]

    return in_net

def query(name, nodes, full_df, uniprot_df, string_type):
    if string_type=="gene":
        name = name.upper()
        result = single_query(name, nodes, uniprot_df)
    else:
        name = name.lower()
        result = name_query(name, nodes, full_df, string_type)
    return result

