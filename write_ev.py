import pandas as pd
import os
from glob import glob

#Read in evidence df
r_df = pd.read_csv("/xdisk/guangyao/michellewei/PMC_OA_Processed/evidence.csv")
b_df = pd.read_pickle("/xdisk/guangyao/michellewei/BIOGRID_pickles/BIOGRID_evidence.pkl")
full_df = r_df.merge(b_df, how="outer", on=["source", "target"], suffixes=["_r", "_b"])
full_df = full_df.fillna("")
full_df["evidence"] = full_df["evidence_r"] + "&&&" + full_df["evidence_b"]
#full_df.to_pickle("/xdisk/guangyao/michellewei/BIOGRID_pickles/concat_ev.pkl")

#Get lists of ev, target ID, source ID
ev_list = full_df["evidence"].tolist()
targets = full_df["target"].tolist()
sources = full_df["source"].tolist()

n_dir=25

for i in range(len(ev_list)):
    delim = f"{sources[i]}_{targets[i]}.txt"
    dir_num=i%n_dir
    dir_name = f"/xdisk/guangyao/michellewei/PMC_OA_BIOGRID/ev_{dir_num}"
    if not os.path.exists(dir_name):   
        os.mkdir(dir_name)
    filename = f"{dir_name}/{delim}"
    
    #Append evidence to file
    ev = ev_list[i]
    with open(filename, "w") as file:
        file.write(ev)       
