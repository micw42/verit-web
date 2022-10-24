import pandas as pd
import numpy as np

#Returns a dictionary with the ids that were found and the ids that weren't found
def check(edges_table, queries_id):
    '''
    Args:
        - edges_table (pd.DataFrame): table of connections between species
        - queries_id (array): array of user-input queried IDs
    '''
    
    all_specs = pd.concat([edges_table["source"], edges_table["target"]]).unique()
    
    present = np.intersect1d(all_specs, queries_id)
    not_in = np.setdiff1d(queries_id, all_specs)
    
    result_dict = {"not_in": list(not_in), "present":list(present)}
    
    return result_dict
