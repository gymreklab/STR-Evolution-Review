#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
import csv
import statistics
from statistics import mean

file = open(sys.argv[1],"r")
strcount= [] # counts of # STRs for each repeated units
while True:
    line = file.readline()
    if line == "":
        break
    else:
        strcount.append(int(line.strip().split()[1]))
        
file = open(sys.argv[2],"r")
f = file.readline().replace("\n","")
if f == "":
    genomelength = 0
else:
    genomelength= int(f)

        
file = open(sys.argv[3],"r")
f = file.readline().replace("\n","")
if f == "":
    strbase = 0
else:
    strbase= int(f)

species_name= sys.argv[4]
filepath = sys.argv[5]

# In[ ]:
strnumsum = sum(strcount)

#percentage for each repeated unites
strpercent = []
for i in range(0,len(strcount)):
    if strnumsum == 0:
        strpercent.append(0)
    else:
        strpercent.append(strcount[i]/float(strnumsum))

# In[ ]:
final = []
row = []
for i in strcount:
    row.append(i)
for k in strpercent:
    row.append(k)
row.append(genomelength)
row.append(strbase)
if genomelength == 0:
    row.append(0)
else:
    row.append(strbase/float(genomelength))
row.append(strnumsum)
if genomelength == 0:
    row.append(0)
else:
    row.append(strnumsum/float(genomelength))

#for a in repeatlength:
 #   row.append(a)

final.append(species_name)

for n in row:
    final.append(str(n))
# In[ ]:

outf = open(filepath, "w")
outf.write(",".join([str(item) for item in final])+"\n")
outf.close()
