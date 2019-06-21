# -*- coding:utf-8 -*-
import pandas as pd
from pandas import Series,DataFrame
import numpy as np
import os


folder="sdf"
files=os.listdir(folder)

atom_list={"H":1,"C":6,"N":7,"O":8,}

for fname in files:
    print(fname)
    if fname[-3:-1]+fname[-1]!="sdf":
        break
    file_name=fname[0:-4]

    molecules_file=open("sdf/"+file_name+".sdf","r")

    pos=True
    after_name=0
    atom_num=0
    data={"1x":[],"2y":[],"3z":[],"4label":[]}
    for line in molecules_file.readlines():
        if pos:
            name=line[0:-1]
            pos=not pos
            after_name=1
            continue
        elif(after_name<3):
            after_name+=1
            continue
        elif(after_name==3):
            points_num=int(line.split()[0])
            A=np.zeros((points_num,3))
            label=[0]*points_num
            after_name+=1
            continue
        elif(atom_num<points_num):
            hoge=line.split()
            data["1x"].append(float(hoge[0]))
            data["2y"].append(float(hoge[1]))
            data["3z"].append(float(hoge[2]))
            data["4label"].append(atom_list[hoge[3]])
            atom_num+=1
            continue
        elif(line=="$$$$\n"):
            pos=True
            after_name=0
            atom_num=0
    data["1x"].insert(0,points_num)
    data["2y"].insert(0,points_num)
    data["3z"].insert(0,points_num)
    data["4label"].insert(0,points_num)
    frame=DataFrame(data)
    frame.to_csv("csv/"+file_name+".csv",header=False,index=False)