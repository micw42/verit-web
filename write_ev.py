import pandas as pd
import os
from glob import glob

#Read in evidence df
df = pd.read_csv("/xdisk/guangyao/michellewei/PMC_OA_Processed/evidence_BIOGRID.csv")

#Get lists of ev, target ID, source ID
ev_list = df["evidence"].tolist()
targets = df["target"].tolist()
sources = df["source"].tolist()

for i in range(len(ev_list)):
    delim = f"{sources[i]}_{targets[i]}.txt"
    #See if file already exists
    file_list = glob(f'/xdisk/guangyao/michellewei/PMC_OA_Processed/ev_*/{delim}', recursive=True)
    #If file doesn't exist put file in the "no_reach" dir
    if len(file_list)==0:
        dir_name = "/xdisk/guangyao/michellewei/PMC_OA_Processed/no_reach"
        if not os.path.exists(dir_name):   
            os.mkdir(dir_name)
        filename = f"{dir_name}/{delim}"
    #If file exists
    else:
        filename = file_list[0]
    
    #Append to file
    ev = "&&&"+ev_list[i]
    with open(filename, "a") as file:
        file.write(ev)       
