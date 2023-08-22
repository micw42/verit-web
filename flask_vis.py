from flask import Flask, render_template, url_for, request, redirect, jsonify, Response, make_response, session
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
import pickle


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

app.secret_key = 'verit_sk'


with open(f"{aws_path}aws.txt", "r") as f:
    creds = f.read().split(",")
access_key = creds[0]
secret_key = creds[1]

start = time.time()

print("Reading edges...", end="\x1b[1K\r")
edges_df=pd.read_pickle(f"{pickle_path}edges.pkl")
print("Loaded edges.")

print("Reading nodes...", end="\x1b[1K\r")
nodes_df=pd.read_pickle(f"{pickle_path}nodes.pkl").dropna()
print("Loaded nodes.")

print("Reading databases...", end="\x1b[1K\r")
full_df=pd.read_pickle(f"{pickle_path}combinedDBs.pkl")
uniprot_df=pd.read_pickle(f"{pickle_path}uniprot_nodes.pkl")
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

@app.route('/rerunCLC', methods=["POST"])
def rerunCLC():
    subset_nodes = request.form['subset_nodes']
    subset_edges = request.form["subset_edges"]
    layout = request.form["layout"]
    subset_nodes = json.loads(subset_nodes)
    subset_edges = json.loads(subset_edges)
    node_df = pd.DataFrame.from_dict({"Id":[node["data"]["id"] for node in subset_nodes], 
                                 "Type": [node["data"]["type"] for node in subset_nodes]})
    edge_df = pd.DataFrame.from_dict({"source":[edge["data"]["source"] for edge in subset_edges],
                                     "target":[edge["data"]["target"] for edge in subset_edges],
                                     "thickness":[edge["data"]["thickness"] for edge in subset_edges]})
    if layout=="clc":
        icp = int(request.form["icp"])
        coord_df = layeredConcentric.cluster_layered_concentric(node_df, edge_df, icp=icp, arr_family=True).set_index('Id')
    elif layout=="lc":
        coord_df = layeredConcentric.layered_concentric(node_df).set_index('Id')
    coord_dict = coord_df.to_dict(orient="index")
    return jsonify(coords=coord_dict)

@app.route('/')
def home():
    resp = make_response(render_template('home.html'))
    user_id = str(uuid.uuid4())
    resp.set_cookie('user_id', user_id)
    return resp
    


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
        user_id = request.cookies.get('user_id')
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
                if string_type == "gene" or default_PR:
                    results, result_bg = MultiSearcher.query(query, nodes_df, full_df, uniprot_df, bg_nodes, string_type)
                    if result_bg is not None:
                        cache.set(f"bg_multi_{user_id}", result_bg)
                    result_dict=ConvertSearch.multi_convert(results)
                    query_dict = {}
                    for query_key in result_dict:
                        query_dict[query_key] = [result_dict[query_key]["max_PR"]]
                    result_dict = ConvertSearch.get_missing(query, query_dict, string_type)
                    
                    cache.set(f"query_dict_{user_id}", query_dict)
                    cache.set(f"result_dict_{user_id}", result_dict)
                    
                    return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
                
                else:
                    results, _ = MultiSearcher.query(query, nodes_df, full_df, uniprot_df, bg_nodes, string_type)
                    result_dict=ConvertSearch.multi_convert(results)
                    cache.set(f"result_dict_{user_id}", result_dict)
                    query = ",".join(query)
                    query=query.replace(" ","SPACE")
                    return redirect(url_for("pick_query", query=query, query_type = "name"))

        else:
            query=request.form["query"].strip()
            if string_type=="id":
                query_dict, result_dict = DictChecker.check(edges_df, query, query_type=query_type)
                cache.set(f"query_dict_{user_id}", query_dict)
                cache.set(f"result_dict_{user_id}", result_dict)
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
                    cache.set(f"query_dict_{user_id}", query_dict)
                    cache.set(f"result_dict_{user_id}", result_dict)
                    return redirect(url_for("display_options", query_type=query_type, string_type=string_type))
                else:
                    cache.set(f"result_dict_{user_id}", result_dict)
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
    
    user_id = request.cookies.get('user_id')
    user_query=query.replace("SPACE"," ")
    user_query = user_query.split(",")

    result_dict = cache.get(f"result_dict_{user_id}")

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
            cache.set(f"query_dict_{user_id}", query_dict)
            return redirect(url_for("make_single_query", query_type="name", string_type="name"))
        else:
            query_dict={}
            for query in result_dict:
                query_dict[query] = request.form.getlist(query)
            cache.set(f"query_dict_{user_id}", query_dict)
            return redirect(url_for("make_bfs_query", string_type="name", query_type="name"))

    else:
        if query_type=="single":
            return render_template("pick_results_single.html", result_dict=result_dict)
        else:
            not_in=list(set([q.lower() for q in user_query])-set([k["Display"].lower() for k in result_dict.values()]))
            return render_template("pick_results_multi.html", query=query, result_dict=result_dict, not_in=not_in)


@app.route('/options/<query_type>/<string_type>', methods=["POST","GET"])
def display_options(query_type, string_type):
    user_id = request.cookies.get('user_id')
    query = cache.get(f"query_dict_{user_id}")
    result_dict = cache.get(f"result_dict_{user_id}")
    not_in = result_dict["not_in"]
    present = result_dict["present"]
    if request.method=="POST" or len(not_in)==0:
        if query_type=="dijkstra":
            if string_type == "id":
                return redirect(url_for("make_bfs_query", string_type=string_type, query_type="id"))
            elif string_type == "gene":
                return redirect(url_for("make_bfs_query", string_type=string_type, query_type="gene"))
            else:
                return redirect(url_for("make_bfs_query", string_type=string_type, query_type="name"))
        elif query_type=="single":
            return redirect(url_for("make_single_query", query_type="id", string_type=string_type))
    else:
        return render_template("validate_result.html", not_in=not_in, query_type=query_type, present=present)


@app.route('/bfs/<string_type>/<query_type>', methods=["POST","GET"])
def make_bfs_query(string_type, query_type):
    user_id = request.cookies.get('user_id')
    query = cache.get(f"query_dict_{user_id}")

    if request.method=="POST":
        max_linkers=int(request.form["max_linkers"])
        min_thickness=int(request.form["min_thickness"])
        qtype = request.form["qtype"]
        return redirect(
            url_for("bfs_query_result",
                    max_linkers=max_linkers,
                    string_type=string_type,
                    qtype=qtype,
                    query_type=query_type,
                    min_thickness=min_thickness,
                    get_direct_linkers=True)
        )

    else:
        return render_template("bfs_search.html")


@app.route('/single_query', methods=["POST","GET"])
@app.route('/single_query/<query_type>/<string_type>', methods=["POST","GET"])
def make_single_query(query_type, string_type):
    if request.method=="POST":
        depth=request.form["depth"]
        return redirect(url_for("single_query_result", depth=depth, query_type=query_type, string_type=string_type))
    else:
        return render_template("single_search.html")


@app.route('/bfsresult/<max_linkers>/<qtype>/<string_type>/<query_type>/<min_thickness>/<get_direct_linkers>')
def bfs_query_result(max_linkers, qtype, string_type, query_type, min_thickness, get_direct_linkers):
    global edges_df
    global nodes_df
    global ev_df
    global G
    global access_key
    global secret_key

    user_id = request.cookies.get('user_id')
    query = cache.get(f"query_dict_{user_id}")

    max_linkers=int(max_linkers)
    min_thickness=int(min_thickness)
    if get_direct_linkers == "True":
        get_direct_linkers = True
    else:
        get_direct_linkers = False

    
    # -- BIOGRID portion --
    if string_type == "gene":
        # Run biogrid query
        queries_id = cache.get(f"bg_multi_{user_id}")
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
        query_nodes, query_edges = MultiQuery.query(
            G,
            edges_df,
            nodes_df,
            query,
            max_linkers,
            qtype,
            query_type="gene",
            get_direct_linkers=get_direct_linkers,
            db_df=full_df,
            access_key=access_key,
            secret_key=secret_key,
            bucket=bucket,
            bg_edges=query_bg_edges,
            min_thickness=min_thickness
        )
        
        #Change IDs back to uniprot for compatibility with biogrid       
        query_nodes["Id"] = query_nodes["display_id"]
        query_edges["source"] = query_edges["source_id"]
        query_edges["target"] = query_edges["target_id"]

    else:
        query_bg_nodes = None
        query_bg_edges = None
        query_nodes, query_edges = MultiQuery.query(
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
            min_thickness=min_thickness
        )
    
    elements = to_json_netx.clean(query_nodes, query_edges, query_bg_nodes, query_bg_edges, biogrid=string_type=="gene")        
    
    return render_template(
        "bfs_result.html",
        elements = elements,
        string_type=string_type
    )


@app.route('/singleresult/<depth>/<query_type>/<string_type>')
def single_query_result(depth, query_type, string_type, methods=["GET"]):
    global edges_df
    global nodes_df
    global G
    
    
    user_id = request.cookies.get('user_id')
    q = cache.get(f"query_dict_{user_id}")
    print(f"Query: {q}")
    depth=int(depth)

    query_nodes, query_edges = SingleQuery.query(
        G, 
        edges_df, 
        nodes_df,
        full_df,
        q,
        depth
    )
    
    if string_type == "gene":
        query_bg_nodes, query_bg_edges = SingleQuery.BIOGRID_query(
            bg_G,
            bg_edges,
            bg_nodes,
            q,
            depth
        )
    else:
        query_bg_nodes = None
        query_bg_edges = None
    
    elements = to_json.clean(
        query_nodes,
        query_edges,
        query_bg_nodes,
        query_bg_edges
    )
    
    return render_template("single_query_result.html", elements=elements, string_type=string_type)


@app.route("/process", methods=["POST"])
def process_data():
    user_id = request.cookies.get('user_id')
    query=request.form["next_query"]
    query_dict = {"QUERY_ID":query}
    cache.set(f"query_dict_{user_id}", query_dict)
    return redirect(url_for("make_single_query", query_type = "id", string_type="id"))


@app.route("/go_home", methods=["POST"])
def go_home():
    return redirect(url_for("home"))


@app.route("/getNodes", methods = ['POST'])
def getNodes():
    # with open("outputs/Adjacency.csv") as fp:
    #     csv = fp.read()
    subset_nodes = request.form['subset_nodes']
    subset_nodes = json.loads(subset_nodes)
    ids = [node["data"]["display_id"] for node in subset_nodes]
    labels = [node["data"]["label"] for node in subset_nodes]
    node_dict = {"Id":ids, "Label":labels}
    df = pd.DataFrame.from_dict(node_dict).to_csv(sep='\t')
    return jsonify(nodes_df=df)



@app.route("/getEdges", methods = ['POST'])
def getEdges():
    # with open("outputs/Adjacency.csv") as fp:
    #     csv = fp.read()
    subset_edges = request.form['subset_edges']
    subset_edges = json.loads(subset_edges)
    source_ids = [edge["data"]["source_DI"] for edge in subset_edges]
    target_ids = [edge["data"]["target_DI"] for edge in subset_edges]
    source_labs = [edge["data"]["source_lab"] for edge in subset_edges]
    target_labs = [edge["data"]["target_lab"] for edge in subset_edges]
    thicknesses = [edge["data"]["thickness"] for edge in subset_edges]
    colors = [edge["data"]["color"] for edge in subset_edges]
    edge_dict = {"source_id":source_ids, "source_name":source_labs,
                 "target_id":target_ids, "target_name":target_labs,
                "thickness":thicknesses, "color":colors}
    df = pd.DataFrame.from_dict(edge_dict).to_csv(sep='\t')
    return jsonify(edges_df=df)


if __name__=="__main__":
    print(f'\n{clr.Fore.BLUE}*****  IGNORE ADDRESS GIVEN BELOW. RUN USING localhost:5001  *****\n')
    app.run(host='0.0.0.0', port=5001, debug=True)