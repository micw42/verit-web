import pandas as pd
import time
import numpy as np
import json

def single_convert():
    results=pd.read_csv("singleSearchOut.csv")
    ids=results["id"].tolist()
    names=results["name"].tolist()
    pr_vals = results["PR"].tolist()
    print("PR vals:", pr_vals)
    result_dict={}
    for i in range(len(ids)):
        if ids[i] not in result_dict:
            result_dict[ids[i]]={"names":[]}
        result_dict[ids[i]]["names"].append(names[i])
        result_dict[ids[i]]["PR"] = pr_vals[i]

    if result_dict:
        result_dict["max_PR"] = ids[0]
    return result_dict

def multi_convert():
    results=pd.read_csv("multiSearchOut.csv")

    ids=results["id"].tolist()
    names=results["name"].tolist()
    pr_vals = results["PR"].tolist()
    queries=results["user_query"].tolist()
    query_display=results["user_query"].tolist()  #The name to display, but not use as a key, b/c flask doesn't like spaces in request.form[]

    s = time.time()
    print("Starting convert loop 1.")
    result_dict={}
    for i in range(len(ids)):
        queries[i]=queries[i].replace(" ", "SPACE")
        if queries[i] not in result_dict:
            result_dict[queries[i]]={}
        if ids[i] not in result_dict[queries[i]]:
            result_dict[queries[i]][ids[i]]=[]
        name_tuple = (names[i], pr_vals[i])
        result_dict[queries[i]][ids[i]].append(name_tuple)
    print("Finished CV1 in", time.time()-s)

    s = time.time()
    print("Starting convert loop 2.")
    for query in result_dict:
        ids = list(result_dict[query].items())
        sorted_ids = sorted(ids, key=lambda x:x[1][0][1], reverse = True)
        max_id = sorted_ids[0][0]
        id_dict = dict(sorted_ids)
        result_dict[query] = id_dict
        result_dict[query]["Display"] = query.replace("SPACE", " ")
        result_dict[query]["max_PR"] = max_id
    print("Finished CV2 in", time.time()-s)

    return result_dict

def get_missing(query_list, result_dict, string_type):
    if string_type=="gene":
        query_list = [x.upper() for x in query_list]
    in_net = list(result_dict.keys())
    in_net = [x.replace("SPACE", " ") for x in in_net]
    not_in = np.setdiff1d(query_list, in_net)
    result_dict = {"not_in": list(not_in), "present":in_net, "query":query_list}
    with open("result_dict.json", "w") as outfile:
        json.dump(result_dict, outfile)
    print("Finished get missing")



