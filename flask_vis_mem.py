from flask import Flask, render_template, url_for, request, redirect, jsonify, Response, make_response
from Query import DictChecker, SingleQuery, SingleSearcher, MultiSearcher, ConvertSearch, MultiQuery, GeneConvert, GetEv
from Visualization import to_json, to_json_netx, layeredConcentric
import pandas as pd
import colorama as clr
import networkx as nx
import time
import copy
import json
import os
import numpy as np
from werkzeug.utils import secure_filename
from os.path import expanduser
import random
from memory_profiler import profile, LogFile
import sys


def validate(query, query_type="name", string_type="contains", default_PR=True, max_linkers=1, qtype="all_simple_paths"):
    if string_type == "id":
        DictChecker.check(edges_df, query, query_type=query_type)
        return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
    else:
        MultiSearcher.query(query, nodes_df, full_df, uniprot_df, string_type)
        result = pd.read_csv("multiSearchOut.csv")
        result_dict=ConvertSearch.multi_convert()
        if string_type=="gene" or default_PR:
            query_dict={}
            for query_key in result_dict:
                query_dict[query_key] = [result_dict[query_key]["max_PR"]]
            with open("query_dict.json", "w") as outfile:
                json.dump(query_dict, outfile)
            ConvertSearch.get_missing(query, query_dict, string_type)
            bfs_query_result(max_linkers=max_linkers, qtype=qtype, query_type=query_type, get_direct_linkers=True)
        else:
            with open("query_dict.json", "w") as outfile:
                json.dump(result_dict, outfile)
            query = ",".join(query)
            query=query.replace(" ","SPACE")
            return redirect(url_for("pick_query", query=query, query_type = "name"))

def bfs_query_result(max_linkers=1, qtype="all_simple_paths", query_type="name", get_direct_linkers=True):
    global edges_df
    global nodes_df
    global ev_df
    global G
    global access_key
    global secret_key

    with open('query_dict.json') as json_file:
        query = json.load(json_file)
    MultiQuery.query(G, edges_df, nodes_df, query, max_linkers, qtype, query_type, get_direct_linkers = get_direct_linkers, db_df = full_df, access_key=access_key, secret_key=secret_key, bucket=bucket)

    query_edges = pd.read_csv("query_edges.csv",header = 0)
    elements, n_query, n_direct = to_json_netx.clean()

    # Compute X and Y for concentric layout
    r1 = 500

    Xs = []; Ys = []
    Xs_q, Ys_q, R_arr_q, n_arr_q = layeredConcentric.get_xy(n_query, r=r1)
    Xs.extend(Xs_q); Ys.extend(Ys_q)

    r2 = 100
    n_fl_co_d = 2 * np.pi * (R_arr_q[-1] + 3*r1) / (2 * r2)
    Xs_d, Ys_d, R_arr_d, n_arr_d = layeredConcentric.get_xy(n_direct, n_fl_co_d, r=r2)
    Xs.extend(Xs_d); Ys.extend(Ys_d)

if __name__=="__main__":
    with open("./settings.txt") as file:
        settings = [x.strip("\n") for x in file.readlines()]

    pickle_path = settings[0]
    aws_path = settings[1]

    UPLOAD_FOLDER = "uploads/"

    bucket = settings[2]

    app = Flask(__name__)
    app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

    with open(f"{aws_path}aws.txt", "r") as f:
        creds = f.read().split(",")
    access_key = creds[0]
    secret_key = creds[1]

    start = time.time()

    print("Reading edges...", end="\x1b[1K\r")
    edges_df=pd.read_pickle(f"{pickle_path}edges.pkl")
    print("Loaded edges.")

    print("Reading nodes...", end="\x1b[1K\r")
    nodes_df=pd.read_pickle(f"{pickle_path}nodes_withKB.pkl")
    print("Loaded nodes.")

    print("Reading databases...", end="\x1b[1K\r")
    full_df=pd.read_pickle(f"{pickle_path}combinedDBs.pkl")
    uniprot_df=pd.read_pickle(f"{pickle_path}FiltUniProtMappings_human.pkl")

    print("Loaded databases.")

    print(f"{clr.Fore.GREEN}Loaded pickles in {round(time.time() - start, 3)}s.{clr.Style.RESET_ALL}")

    G = nx.from_pandas_edgelist(edges_df, edge_attr=True, source="source", target="target", create_using=nx.DiGraph())
    query=["e2f", "rb", "cyclin e", "cyclind", "myc"]
    
    validate(query)
    

