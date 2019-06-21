import pandas as pd
from pandas import Series,DataFrame
import numpy as np
import os

f=open("hoge.txt","r")

ene=[]
for line in f.readlines():
    hoge=line.split()
    ene.append(hoge[1])

s=Series(ene)
s.to_csv("energy.csv",header=False,index=False)