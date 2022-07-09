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
    
    # Reverse search to collect all aliases
    u_ids = list(in_net["id"].drop_duplicates())
    in_net = full_df.query("id == @u_ids").iloc[::-1].reset_index(drop=True)
    
    in_net["user_query"] = name
    in_net = pd.merge(in_net, nodes, left_on="id", right_on = "Id", how="inner")
    in_net = in_net.drop_duplicates(subset = ["name", "id"])
    
    return in_net


def multi_query(query_list, nodes, full_df, string_type):
    # initialize empty dataframe
    full_query = pd.DataFrame()

    for query in query_list:
        full_query = full_query.append(name_query(query, nodes, full_df, string_type))

    return full_query



def query(query_list, nodes, full_df, string_type):
    # Turns all entries into lower case
    case_proof = [x.lower() for x in query_list]
    multi_query(case_proof, nodes, full_df, string_type).to_csv("multiSearchOut.csv", index=False)


