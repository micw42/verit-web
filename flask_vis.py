from flask import Flask, render_template, url_for, request, redirect, jsonify
from Query import DictChecker, SingleQuery, SingleSearcher, MultiSearcher, ConvertSearch, MultiQuery, GeneConvert, GetEv
from Visualization import to_json, to_json_netx
import pandas as pd
import colorama as clr
import networkx as nx
import time
import copy
import json
import os
from werkzeug.utils import secure_filename
from os.path import expanduser

pickle_path = "./fixed_ev_pickles/"

UPLOAD_FOLDER = "uploads/"

bucket = "updated-abstract-evidence"

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

with open("aws.txt", "r") as f:
    creds = f.read().split(",")
access_key = creds[0]
secret_key = creds[1]

start = time.time()

print("Reading edges...", end="\x1b[1K\r")
edges_df=pd.read_pickle(f"{pickle_path}edges.pkl")
print("Loaded edges.")

print("Reading nodes...", end="\x1b[1K\r")
nodes_df=pd.read_pickle(f"{pickle_path}nodes.pkl")
print("Loaded nodes.")

print("Reading databases...", end="\x1b[1K\r")
full_df=pd.read_pickle(f"{pickle_path}combinedDBs.pkl")

print("Loaded databases.")

print(f"{clr.Fore.GREEN}Loaded pickles in {round(time.time() - start, 3)}s.{clr.Style.RESET_ALL}")

G = nx.from_pandas_edgelist(edges_df, edge_attr=True, source="source", target="target", create_using=nx.DiGraph())

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
        
        #Only multi query accepts file input
        if query_type == "dijkstra":
            file = request.files['file']
            
            #Handles file input
            if file:
                #filename = secure_filename(file.filename)
                #file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                #file.save(file_path)
                #with open(file_path) as file:
                #    file_contents = file.read()
                file_contents = file.read().decode("utf-8") 
                query = file_contents.splitlines()
                string_type=request.form["queryType"]
                if string_type == "id":
                    result_dict=DictChecker.check(edges_df, query)
                    result_dict = json.dumps(result_dict)
                    query = ",".join(query)
                    query=query.replace(" ","SPACE")
                    query_dict = {"QUERY_ID":query}
                    query_dict = json.dumps(query_dict)
                    return redirect(url_for("display_options", query_type=query_type, query=query_dict, result_dict=result_dict))
                else:
                    query = ",".join(query)
                    query=query.replace(" ","SPACE")
                    return redirect(url_for("pick_query", query_type=query_type, query=query, string_type=string_type))

        if request.form["queryType"]=="id":
            query=request.form["query"].split(",")
            result_dict=DictChecker.check(edges_df, query)
            print("DictChecker result:", result_dict)
            result_dict = json.dumps(result_dict)
            query=request.form["query"]
            query_dict = {"QUERY_ID":query}
            query_dict = json.dumps(query_dict)
            #Id query with text box goes directly to page that shows which ones were found
            return redirect(url_for("display_options", result_dict=result_dict, query_type=query_type, query=query_dict))
        
        #If the query is a text box format with name contains/begins/ends
        else:
            query=request.form["query"]
            string_type=request.form["queryType"]
            default_PR = request.form.get("default_PR")
            print("Default PR:", default_PR)
            if query_type=="dijkstra":
                default_PR = request.form.get("default_PR")
                if default_PR:
                    query=query.split(",")
                    MultiSearcher.query(query, nodes_df, full_df, string_type)
                    result = pd.read_csv("multiSearchOut.csv")
                    result_dict=ConvertSearch.multi_convert()
                    query_dict={}
                    for query in result_dict:
                        query_dict[query] = [result_dict[query]["max_PR"]]
                    query_dict = json.dumps(query_dict)
                    return redirect(url_for("make_bfs_query", query=query_dict, query_type = "name"))
                    
            if query_type=="single":
                if default_PR:
                    print("Default PR")
                    SingleSearcher.query(query, nodes_df, full_df, string_type)
                    result_dict=ConvertSearch.single_convert()
                    query = {"QUERY_ID":result_dict["max_PR"]}
                    query = json.dumps(query)
                    return redirect(url_for("make_single_query", query=query, query_type=query_type))
                query=request.form["query"]
                if string_type=="gene":
                    conv_dict = GeneConvert.convert_genes([query])
                    query = conv_dict[query]
                    query = {"QUERY_ID":query}
                    query = json.dumps(query)
                    return redirect(url_for("make_single_query", query=query, query_type=query_type))
                #Remove spaces from query entries
                
            query=query.replace(" ","SPACE")
            return redirect(url_for("pick_query", query_type=query_type, query=query, string_type=string_type))
            
    else:
        if query_type=="single":
            return render_template("validate_single.html")
        else:
            return render_template("validate_bfs.html")

@app.route('/pick_query/<query_type>/<query>/<string_type>', methods=["POST","GET"])
def pick_query(query_type, query, string_type):
    global nodes_df
    global full_df
    
    user_query=query.replace("SPACE"," ")
    user_query = user_query.split(",")
    
    if query_type=="single":
        SingleSearcher.query(user_query[0], nodes_df, full_df, string_type)
        result_dict=ConvertSearch.single_convert()
        print("Pick query result dict:", result_dict)
       
        
    else:
            #Handles gene ids    
            if string_type == "gene":
                id_dict = GeneConvert.convert_genes(user_query)
                conv_genes = list(id_dict.keys())    #the genes that were able to be converted to uniprot id
                conv_ids = list(id_dict.values())   #converted uniprot ids
                not_in = list(set(user_query) - set(conv_genes))
                no_id=DictChecker.check(edges_df, conv_ids)   #The uniprot ids not found in the network
                no_id_gene = [conv_genes[i] for i in range(len(conv_ids)) if conv_ids[i] in no_id["not_in"]]   #genes w/ uniprot id not in network
                not_in.extend(no_id_gene)
                present = list(set(conv_genes) - set(no_id_gene)) #genes w/ uniprot id present in network
                result_dict = {"not_in":not_in, "present":present}
                result_dict=json.dumps(result_dict)
                query=",".join(conv_ids)    #new query consists just of uniprot ids found in network
                query_dict = {"QUERY_ID":query}
                query_dict = json.dumps(query_dict)
                return redirect(url_for("display_options", result_dict = result_dict, query_type=query_type, query=query_dict))
            
            #Handles name queries
            else:  
                MultiSearcher.query(user_query, nodes_df, full_df, string_type)
                result = pd.read_csv("multiSearchOut.csv")
                result_dict=ConvertSearch.multi_convert()

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
            query_dict = json.dumps(query_dict)
            return redirect(url_for("make_single_query", query=query_dict, query_type="name"))
        
        else:
            query_dict={}
            default_PR = request.form.get("default_PR")
            if default_PR:
                for query in result_dict:
                    query_dict[query] = [result_dict[query]["max_PR"]]
            else:
                for query in result_dict:
                    query_dict[query] = request.form.getlist(query)
            query_dict = json.dumps(query_dict)
            return redirect(url_for("make_bfs_query", query=query_dict, query_type = "name"))
          
    else:
        if query_type=="single":
            print("Result dict type:", type(result_dict))
            return render_template("pick_results_single-test.html", result_dict=result_dict)
        else:
            not_in=list(set([q.lower() for q in user_query])-set([k["Display"].lower() for k in result_dict.values()]))
            return render_template("pick_results_multi.html", query=query, result_dict=result_dict, not_in=not_in)
        
    
@app.route('/options/<query_type>/<result_dict>/<query>', methods=["POST","GET"])
def display_options(result_dict, query_type, query):
    result_dict = json.loads(result_dict)
    not_in = result_dict["not_in"]
    present = result_dict["present"]
    query_parts=(json.loads(query))["QUERY_ID"].split(",")
   
    if request.method=="POST":
        choice=request.form["choice"]
        if choice=="Try another query":
            return redirect(url_for("select_query"))
        else:
            if query_type=="dijkstra":
                return redirect(url_for("make_bfs_query", query=query, query_type = "id"))
            elif query_type=="single":
                return redirect(url_for("make_single_query", query=query, query_type = "id"))
    else:
        return render_template("validate_result.html", not_in=not_in, query_type=query_type, present=present, query=query)

@app.route('/bfs/<query>/<query_type>', methods=["POST","GET"])
def make_bfs_query(query, query_type):
    if request.method=="POST":
        query_string=query
        max_linkers=int(request.form["max_linkers"])
        qtype = request.form["qtype"]
        return redirect(url_for("bfs_query_result",query_string=query_string, max_linkers=max_linkers, qtype = qtype, query_type = query_type, get_direct_linkers = True))
    else:
        return render_template("bfs_search.html", query=query)

@app.route('/single_query', methods=["POST","GET"])
@app.route('/single_query/<query>/<query_type>', methods=["POST","GET"])
def make_single_query(query, query_type):
    if request.method=="POST":
        depth=request.form["depth"]
        return redirect(url_for("single_query_result",query=query, depth=depth, query_type = query_type))
    else:
        return render_template("single_search.html")


@app.route('/bfsresult/<query_string>/<max_linkers>/<qtype>/<query_type>/<get_direct_linkers>')
def bfs_query_result(query_string, max_linkers, qtype, query_type, get_direct_linkers):
    global edges_df
    global nodes_df
    global ev_df
    global G
    global access_key
    global secret_key
    
    query=json.loads(query_string)
    max_linkers=int(max_linkers)
    if get_direct_linkers == "True":
        get_direct_linkers = True
    else:
        get_direct_linkers = False
    MultiQuery.query(G, edges_df, nodes_df, query, max_linkers, qtype, query_type, get_direct_linkers = get_direct_linkers, db_df = full_df, access_key=access_key, secret_key=secret_key, bucket=bucket)

    elements = to_json_netx.clean()
    sq_align = to_json_netx.get_square_clusters()
    orig_align_q, orig_align_l = to_json_netx.get_orig_clusters()
    return render_template("bfs_result.html", elements = elements, sq_align = sq_align, orig_align_q = orig_align_q, orig_align_l = orig_align_l)



@app.route('/singleresult/<query>/<depth>/<query_type>')
def single_query_result(query, depth, query_type, methods=["GET"]):
    global edges_df
    global nodes_df
    global G
    
    query = json.loads(query)
    depth=int(depth)

    SingleQuery.query(G, edges_df, nodes_df, query, depth)
    elements=to_json.clean()
    return render_template("single_query_result.html", elements=elements)

@app.route("/process", methods=["POST"])
def process_data():
    query=request.form["next_query"]
    query_dict = {"QUERY_ID":query}
    query_dict = json.dumps(query_dict)
    return redirect(url_for("make_single_query", query=query_dict, query_type = "id"))

@app.route("/go_home", methods=["POST"])
def go_home():
    return redirect(url_for("home"))

if __name__=="__main__":
    print(f'\n{clr.Fore.BLUE}*****  IGNORE ADDRESS GIVEN BELOW. RUN USING localhost:5001  *****\n')
    app.run(host='0.0.0.0', port=5001, debug=True)
