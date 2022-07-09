import pandas as pd

def name_query(name, nodes, full_df, string_type):
    if string_type == "starts_with":
        name_query = full_df.query('name.str.startswith(@name, na=False)', engine='python').reset_index(drop=True)
    
    elif string_type == "ends_with":
        name_query = full_df.query('name.str.endswith(@name, na=False)', engine='python').reset_index(drop=True)
    
    elif string_type == "contains":
        name_query = full_df.query('name.str.contains(@name, na=False)', engine='python').reset_index(drop=True)

    # The following nodes are in the network, but REACH likely grounds to more
    in_net = pd.merge(nodes, name_query, left_on="Id", right_on="id", how="inner").iloc[:, 3:6].drop_duplicates().reset_index(drop=True)
    print("COLS:", in_net.columns)
    u_ids = list(in_net["id"].drop_duplicates())
    in_net = full_df.query("id == @u_ids").iloc[::-1].reset_index(drop=True)
    #in_net = pd.merge(nodes, in_net, left_on="Id", right_on="id", how="inner").iloc[:, 3:6].drop_duplicates().reset_index(drop=True)
    in_net = pd.merge(nodes, in_net, left_on="Id", right_on="id", how="inner").drop_duplicates().reset_index(drop=True)
    in_net.to_csv("in_net.csv", index=False)
    in_net = in_net.drop_duplicates(subset = ["name", "id"])
    in_net = in_net.sort_values(by = "PR", ascending = False)
    
    return in_net


def query(name, nodes, full_df, string_type):
    name = name.lower()
    print("NAME:", name)
    result = name_query(name, nodes, full_df, string_type)
    result.to_csv("singleSearchOut.csv", index=False)

