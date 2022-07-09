import pandas as pd

#Returns a dictionary with the ids that were found and the ids that weren't found
def check(edges_table, queries_id):
    sources = edges_table["source"].tolist()
    targets = edges_table["target"].tolist()

    all_species = list(set(sources + targets))
    not_in=[]

    for query in queries_id:
        print("DictChecker Query:", query)
        if query not in all_species:
            not_in.append(query)
        else:
            continue
    
    present = list(set(queries_id) - set(not_in))    #ids that were found
    result_dict = {"not_in":not_in, "present":present}
    return result_dict


