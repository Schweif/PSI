# -*- coding: utf-8 -*-
"""
Returns maximum value of a Spectra data file
Created on Wed Feb 14 11:18:08 2024

@author: just_d
"""

import glob
import os
import json



##### Configuration

specific_file = None #leave None if latest file should be analyzed
base_path = r"C:\Users\just_d\Desktop\tmp\Spectra" #path where to look for latest file
ending = "\*.json"

##### Configuration


def returnNewestFile(base_path, ending):
    list_of_files = glob.glob(base_path + ending) # * means all if need specific format then *.csv
    latest_file = max(list_of_files, key=os.path.getctime)
    return(latest_file)

if __name__ == '__main__':
    #assign latest or specific file to task
    if not specific_file:
        file = returnNewestFile(base_path, ending)
    else:
        file = specific_file
    
    with open(file, "r") as f: #opens json file
        data = json.load(f)
    
    distance = data["Input"]["Configurations"]["Distance from the Source (m)"]
    title = data["Output"]["titles"][2]
    unit = data["Output"]["units"][2]
    max_value = max(data["Output"]["data"][2])
    
    print( "The maximum " + title + " is " + str(max_value) + " " + unit + " at " + str(distance) + " m distance.")
        
