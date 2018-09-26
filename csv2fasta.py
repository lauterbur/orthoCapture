#!/bin/python

import sys
import pandas

with open(sys.argv[1], "r") as f:
        name = sys.argv[1].split(".")[0]+".fa"
        with open(name, "a") as o:
                s = pandas.read_csv(f)
                for i in range(len(s)):
                        if not pandas.isnull(s.iloc[i][0]):
                                if pandas.isnull(s.iloc[i][2]):
#                                       print(type(s.iloc[i][1]))
#                                       print(type(s.iloc[i][3]))
                                        o.write(">"+str(s.iloc[i][1])+"_"+str(s.iloc[i][3])+"\n")
                                        o.write(s.iloc[i][4]+"\n")
                                else:
                                        o.write(">"+str(s.iloc[i][1])+"_"+str(s.iloc[i][2])+"_"+str(s.iloc[i][3])+"\n")
                                        o.write(str(s.iloc[i][4])+"\n")