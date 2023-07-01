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
from memory_profiler import profile
from flask_caching import Cache
import uuid

with open("./settings.txt") as file:
    settings = [x.strip("\n") for x in file.readlines()]
print(settings)

pickle_path = settings[0]
aws_path = settings[1]

UPLOAD_FOLDER = "uploads/"

bucket = settings[2]

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config["CACHE_TYPE"] = "SimpleCache" # better not use this type w. gunicorn
cache = Cache(app)

user_id = str(uuid.uuid4())
print("Starting new session with user id:", user_id)

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
bg_nodes = pd.read_pickle(f"{pickle_path}BIOGRID_nodes.pkl")
bg_edges = pd.read_pickle(f"{pickle_path}BIOGRID_edges.pkl")
bg_G = pd.read_pickle(f"{pickle_path}BIOGRID_graph.pkl")

print("Loaded databases.")

print(f"{clr.Fore.GREEN}Loaded pickles in {round(time.time() - start, 3)}s.{clr.Style.RESET_ALL}")

G = nx.from_pandas_edgelist(edges_df, edge_attr=True, source="source", target="target", create_using=nx.DiGraph())

@profile
@app.route('/_get_evidence')
def get_evidence():
    global access_key
    global secret_key
    files = request.args.get('files', 0, type=str)
    files = files.split("%%")
    ev = [GetEv.get_ev(file, access_key=access_key, secret_key=secret_key, bucket=bucket) for file in files]
    ev = "%%".join(ev)
    return jsonify(result=ev)


@app.route('/')
def home():
    return render_template("home.html")


@app.route('/readme')
def readme():
    return render_template("README.html")


@app.route("/selectQuery", methods=["POST","GET"])
def select_query():
    if request.method=="POST":
        query_type=request.form["query_type"]
        if query_type=="Make a Multi Query":
            query_type="dijkstra"
            return redirect(url_for("validate", query_type=query_type))
        else:
            query_type="single"
            return redirect(url_for("validate", query_type=query_type))
    else:
        return render_template("select_query.html")

    
@app.route('/validate/<query_type>', methods=["POST","GET"])
def validate(query_type):
    global edges_df
    if request.method=="POST":
        string_type=request.form["queryType"]
        default_PR = request.form.get("default_PR")
        
        #Only multi query accepts file input
        if query_type == "dijkstra":
            file = request.files['file']
            #Handles file input
            if file:
                file_contents = file.read().decode("utf-8")
                query = file_contents.splitlines()
            else:
                query=request.form["query"].split(",")
                query = [x.strip() for x in query]

            if string_type == "id":
                query_dict, result_dict = DictChecker.check(edges_df, query, query_type=query_type)
                cache.set(f"query_dict_{user_id}", query_dict)
                cache.set(f"result_dict_{user_id}", result_dict)
                return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
            #If the query is a text box format with name contains/begins/ends
            else:
                results = MultiSearcher.query(query, nodes_df, full_df, uniprot_df, bg_nodes, string_type)
                result_dict=ConvertSearch.multi_convert(results)
                if string_type == "gene" or default_PR:
                    query_dict = {}
                    for query_key in result_dict:
                        query_dict[query_key] = [result_dict[query_key]["max_PR"]]
                    result_dict = ConvertSearch.get_missing(query, query_dict, string_type)
                    cache.set("query_dict", query_dict)
                    cache.set("result_dict", result_dict)
                    return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
                else:
                    cache.set("result_dict", result_dict)
                    query = ",".join(query)
                    query=query.replace(" ","SPACE")
                    return redirect(url_for("pick_query", query=query, query_type = "name"))

        else:
            query=request.form["query"].strip()
            if string_type=="id":
                query_dict, result_dict = DictChecker.check(edges_df, query, query_type=query_type)
                cache.set("query_dict", query_dict)
                cache.set("result_dict", result_dict)
                return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
            else:
                results = SingleSearcher.query(query, nodes_df, full_df, uniprot_df, string_type)
                result_dict=ConvertSearch.single_convert(results)
                if default_PR or string_type=="gene":
                    if result_dict:
                        query_dict = {"QUERY_ID":result_dict["max_PR"]}
                        result_dict = {"not_in":[], "present":[query]}
                    else:
                        query_dict = {"QUERY_ID":[]}
                        result_dict = {"not_in":[query], "present":[]}
                    cache.set("query_dict", query_dict)
                    cache.set("result_dict", result_dict)
                    return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
                else:
                    cache.set("result_dict", result_dict)
                    query=query.replace(" ","SPACE")
                    return redirect(url_for("pick_query", query=query, query_type=query_type))

    else:
        if query_type=="single":
            return render_template("validate_single.html")
        else:
            return render_template("validate_bfs.html")


@app.route('/pick_query/<query>/<query_type>', methods=["POST","GET"])
def pick_query(query, query_type):
    global nodes_df
    global full_df

    user_query=query.replace("SPACE"," ")
    user_query = user_query.split(",")

    result_dict = cache.get("result_dict")

    if request.method=="POST":
        if query_type=="single":
            if request.form["query"]=="Try another query":
                    return redirect(url_for("select_query"))
        else:
            if request.form["submit"]=="Try another query":
                    return redirect(url_for("select_query"))

        if query_type=="single":
            query=request.form.getlist("query")
            query_dict = {user_query[0]:query}   #Query[0] since there's only 1 query, and no commas
            cache.set("query_dict", query_dict)
            return redirect(url_for("make_single_query", query_type="name"))
        else:
            query_dict={}
            for query in result_dict:
                query_dict[query] = request.form.getlist(query)
            cache.set("query_dict", query_dict)
            return redirect(url_for("make_bfs_query", string_type="name", query_type="name"))

    else:
        if query_type=="single":
            return render_template("pick_results_single-test.html", result_dict=result_dict)
        else:
            not_in=list(set([q.lower() for q in user_query])-set([k["Display"].lower() for k in result_dict.values()]))
            return render_template("pick_results_multi.html", query=query, result_dict=result_dict, not_in=not_in)


@app.route('/options/<query_type>/<string_type>', methods=["POST","GET"])
def display_options(query_type, string_type):
    query = cache.get(f"query_dict_{user_id}")
    result_dict = cache.get(f"result_dict_{user_id}")
    not_in = result_dict["not_in"]
    present = result_dict["present"]
    if request.method=="POST" or len(not_in)==0:
        if query_type=="dijkstra":
            if string_type == "id":
                return redirect(url_for("make_bfs_query", string_type=string_type, query_type="id"))
            else:
                return redirect(url_for("make_bfs_query", string_type=string_type, query_type="name"))
        elif query_type=="single":
            return redirect(url_for("make_single_query", query_type = "id"))
    else:
        return render_template("validate_result.html", not_in=not_in, query_type=query_type, present=present)


@app.route('/bfs/<string_type>/<query_type>', methods=["POST","GET"])
def make_bfs_query(string_type, query_type):
    query = cache.get(f"query_dict_{user_id}")
    q_len = len(query)
    if q_len > 100:
        return redirect(url_for("bfs_query_result",
                                max_linkers=1,
                                qtype="all_shortest_paths",
                                string_type=string_type,
                                query_type=query_type,
                                get_direct_linkers = True)
                       )

    if request.method=="POST":
        max_linkers=int(request.form["max_linkers"])
        qtype = request.form["qtype"]
        return redirect(
            url_for("bfs_query_result",
                    max_linkers=max_linkers,
                    string_type=string_type,
                    qtype=qtype,
                    query_type=query_type,
                    get_direct_linkers=True)
        )

    else:
        return render_template("bfs_search.html")


@app.route('/single_query', methods=["POST","GET"])
@app.route('/single_query/<query_type>', methods=["POST","GET"])
def make_single_query(query_type):
    if request.method=="POST":
        depth=request.form["depth"]
        return redirect(url_for("single_query_result", depth=depth, query_type = query_type))
    else:
        return render_template("single_search.html")


@app.route('/bfsresult/<max_linkers>/<qtype>/<string_type>/<query_type>/<get_direct_linkers>')
def bfs_query_result(max_linkers, qtype, string_type, query_type, get_direct_linkers):
    global edges_df
    global nodes_df
    global ev_df
    global G
    global access_key
    global secret_key

    query = cache.get(f"query_dict_{user_id}")
    print(query)
        
    print(len(query))
    max_linkers=int(max_linkers)
    if get_direct_linkers == "True":
        get_direct_linkers = True
    else:
        get_direct_linkers = False

    
    # -- BIOGRID portion --
    if string_type == "gene":
        # Run biogrid query
        queries_id = pd.read_pickle("bg_multiSearchOut.pkl")
        query_bg_nodes, query_bg_edges = MultiQuery.BIOGRID_query(
            bg_G,
            bg_edges,
            bg_nodes,
            queries_id,
            max_linkers,
            qtype,
            query_type,
            get_direct_linkers=get_direct_linkers,
            db_df=full_df
        )

        #Run REACH query
        query_nodes, query_edges, nodes_cleaned, edges_cleaned = MultiQuery.query(
            G,
            edges_df,
            nodes_df,
            query,
            max_linkers,
            qtype,
            query_type,
            get_direct_linkers=get_direct_linkers,
            db_df=full_df,
            access_key=access_key,
            secret_key=secret_key,
            bucket=bucket,
            bg_edges=query_bg_edges
        )
        
        #Change IDs back to uniprot for compatibility with biogrid       
        query_nodes["Id"] = query_nodes["display_id"]
        query_edges["source"] = query_edges["source_id"]
        query_edges["target"] = query_edges["target_id"]
    else:
        query_bg_nodes = None
        query_bg_edges = None
        query_nodes, query_edges, nodes_cleaned, edges_cleaned = MultiQuery.query(
            G,
            edges_df,
            nodes_df,
            query,
            max_linkers,
            qtype,
            query_type,
            get_direct_linkers=get_direct_linkers,
            db_df=full_df,
            access_key=access_key,
            secret_key=secret_key,
            bucket=bucket
        )
    n_edges = len(query_edges.index)
    filtered = False
    if n_edges > 5000:
        to_json_netx.filter_graph()
        filtered = True
    
    elements = to_json_netx.clean(query_nodes, query_edges, query_bg_nodes, query_bg_edges, biogrid=string_type=="gene")        
    
    cache.set(f"nodes_cleaned_{user_id}", nodes_cleaned)
    cache.set(f"edges_cleaned_{user_id}", edges_cleaned)
    return render_template(
        "bfs_result.html",
        elements = elements,
        filtered = filtered,
        string_type=string_type
    )


@app.route('/singleresult/<depth>/<query_type>')
def single_query_result(depth, query_type, methods=["GET"]):
    global edges_df
    global nodes_df
    global G

    query = cache.get("query_dict")
    depth=int(depth)

    query_nodes, query_edges = SingleQuery.query(G, edges_df, nodes_df, full_df, query, depth)
    elements=to_json.clean(query_nodes, query_edges, sq=True)
    return render_template("single_query_result.html", elements=elements)


@app.route("/process", methods=["POST"])
def process_data():
    query=request.form["next_query"]
    query_dict = {"QUERY_ID":query}
    cache.set("query_dict", query_dict)
    return redirect(url_for("make_single_query", query=query_dict, query_type = "id"))


@app.route("/go_home", methods=["POST"])
def go_home():
    return redirect(url_for("home"))


@app.route("/getNodes")
def getNodes():
    # with open("outputs/Adjacency.csv") as fp:
    #     csv = fp.read()
    df = cache.get(f"nodes_cleaned_{user_id}")
    resp = make_response(df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=query_nodes.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp


@app.route("/getEdges")
def getEdges():
    # with open("outputs/Adjacency.csv") as fp:
    #     csv = fp.read()
    df = cache.get(f"edges_cleaned_{user_id}")
    resp = make_response(df.to_csv())
    resp.headers["Content-Disposition"] = "attachment; filename=query_edges.csv"
    resp.headers["Content-Type"] = "text/csv"
    return resp


if __name__=="__main__":
    print(f'\n{clr.Fore.BLUE}*****  IGNORE ADDRESS GIVEN BELOW. RUN USING localhost:5001  *****\n')
    app.run(host='0.0.0.0', port=5001, debug=True)
