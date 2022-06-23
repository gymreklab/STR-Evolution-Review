#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import sys
import csv


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plt

infile = sys.argv[1]
out_period_perc = sys.argv[2]
out_genome_perc = sys.argv[3]
out_density = sys.argv[4]
out_1 = sys.argv[5]
out_2 = sys.argv[6]
out_3 = sys.argv[7]
out_4 = sys.argv[8]
out_5 = sys.argv[9]
out_6 = sys.argv[10]
out_str_num = sys.argv[11]



# Load ordered species list
species_list = [item.strip() for item in open("STR_SpeciesList.txt", "r").readlines()]
print(species_list)
#species_list = [item.replace("\n","") for item in open("SpeciesList.txt", "r").readlines()]
species_order = dict(zip(species_list, range(1, len(species_list)+1)))


# Load data and order by species order

#data = pd.read_txt(sys.argv[1])

col_name = ["name","homopolymers","dinucleotide","trinucleotide","tetranucleotide","pentanucleotide","hexanucleotide","per1.perc","per2.perc","per3.perc","per4.perc","per5.perc","per6.perc","genome_length","str_base","str.perc.bp","str.num","str.density","ave_homopolymers","ave_dinucleotide","ave_trinucleotide","ave_tetranucleotide","ave_pentanucleotide","ave_hexanucleotide"]
data = pd.read_csv(infile,header = None, names = col_name)
data["order"] = data["name"].apply(lambda x: species_order[x])
data = data.sort_values("order")

#full full names into dictionary
reader = csv.reader(open('fullname.csv'))

result = {}
for row in reader:
    key = row[0]
    if key in result:
        # implement your duplicate row handling here
        pass
    result[key] = row[1]

data["fullname"] = data["name"].apply(lambda x: result[x])

# Remove species we don't want. e.g. Ebola
#rmspecies = ["Ebola_virus"]
#data = data[~data["name"].isin(rmspecies)]


# In[8]:



data


# In[35]:


# Plot repeat percentages

fig = plt.figure()
fig.set_size_inches((20, 12))
ax = fig.add_subplot(111)

bottoms = np.array([0]*data.shape[0])

percolors = ["gray","red","gold","blue","green","purple"]

for period in range(1, 7):
    values = np.array(data["per%s.perc"%period])
    ax.bar(range(len(values)), values, color=percolors[period-1], bottom=bottoms)
    bottoms = bottoms+values
    
ax.set_xticks(range(len(bottoms)))
ax.set_xticklabels(data["fullname"], rotation=90);

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
fig.savefig(out_period_perc,bbox_inches = 'tight')


# In[36]:


# Plot repeat density

fig = plt.figure()
fig.set_size_inches((20, 12))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["str.density"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_density)

# Plot total str number
fig = plt.figure()
fig.set_size_inches((20, 12))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["str.num"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_str_num)


fig = plt.figure()
fig.set_size_inches((20, 12))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["str.perc.bp"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_genome_perc)

fig = plt.figure()
fig.set_size_inches((20, 12))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_homopolymers"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_1)

# In[44]:


fig = plt.figure()
fig.set_size_inches((15, 8))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_dinucleotide"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_2)


# In[45]:


fig = plt.figure()
fig.set_size_inches((15, 8))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_trinucleotide"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_3)


# In[46]:


fig = plt.figure()
fig.set_size_inches((15, 8))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_tetranucleotide"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_4)


# In[47]:


fig = plt.figure()
fig.set_size_inches((15, 8))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_pentanucleotide"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_5)


# In[48]:


fig = plt.figure()
fig.set_size_inches((15, 8))
ax = fig.add_subplot(111)
ax.bar(range(data.shape[0]), data["ave_hexanucleotide"], color="black", edgecolor="white");
ax.set_xticks(range(data.shape[0]))
ax.set_xticklabels(data["name"], rotation=90);
fig.savefig(out_6)
