import pandas as pd
import numpy as np
import json

#Returns a dictionary with the ids that were found and the ids that weren't found
def check(edges_table, queries_id, query_type):
    '''
    Args:
        - edges_table (pd.DataFrame): table of connections between species
        - queries_id (array): array of user-input queried IDs
    '''
    if query_type=="single":
        queries_id = [queries_id]
        
    all_specs = pd.concat([edges_table["source"], edges_table["target"]]).unique()
    
    present = np.intersect1d(all_specs, queries_id)
    not_in = np.setdiff1d(queries_id, all_specs)
    
    result_dict = {"not_in": list(not_in), "present":list(present)}
    
    with open("result_dict.json", "w") as outfile:
        json.dump(result_dict, outfile)

    query = ",".join(queries_id)
    query=query.replace(" ","SPACE")
    query_dict = {"QUERY_ID":query}
    with open("query_dict.json", "w") as outfile:
        json.dump(query_dict, outfile)

