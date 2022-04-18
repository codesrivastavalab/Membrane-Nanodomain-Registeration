"""
Created on Sat Dec 25 09:29:21 2021

@author: saroj
"""
import re
import string
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec

#Main#
#------------------------------------------------------------------------------#

filename = "PSstate.txt"
file_data = np.loadtxt(filename)
    
vs = file_data[:, 2]
aa = file_data[:, 3]
bb = file_data[:, 4]
ab = file_data[:, 5]
aap = file_data[:, 6]
bbp = file_data[:, 7]
abp = file_data[:, 8]
state = file_data[:, 9]

fig = plt.figure(figsize=(8, 12), dpi=300)

# Bins for enthalpic V'
b = np.array([2.0, 14.0, 26.0, 38.0])

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(aa[i])
    elif state[i]==5:
        ppstemp.append(aa[i])
    elif state[i]==4:
        npstemp.append(aa[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 1)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{AA}$")

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(bb[i])
    elif state[i]==5:
        ppstemp.append(bb[i])
    elif state[i]==4:
        npstemp.append(bb[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 2)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{BB}$")
ax.set_ylabel("No. of Parameter Space Points")

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(ab[i])
    elif state[i]==5:
        ppstemp.append(ab[i])
    elif state[i]==4:
        npstemp.append(ab[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 3)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{AB}$")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 3.525),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D3411P'.png",bbox_inches='tight')

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fig = plt.figure(figsize=(8, 12), dpi=300)

# Bins for enthalpic V
b = np.array([2.5, 5.5, 8.5, 11.5, 14.5, 17.5])

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(aap[i])
    elif state[i]==5:
        ppstemp.append(aap[i])
    elif state[i]==4:
        npstemp.append(aap[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 1)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{AA}$")

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(bbp[i])
    elif state[i]==5:
        ppstemp.append(bbp[i])
    elif state[i]==4:
        npstemp.append(bbp[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 2)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{BB}$")
ax.set_ylabel("No. of Parameter Space Points")

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(abp[i])
    elif state[i]==5:
        ppstemp.append(abp[i])
    elif state[i]==4:
        npstemp.append(abp[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

ax = fig.add_subplot(3, 1, 3)
ax.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{AB}$")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 3.525),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D3411P.png",bbox_inches='tight')

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fig = plt.figure(figsize=(8, 4), dpi=300)

# Bins for enthalpic V
b = np.array([0, 0.0125, 0.0275, 0.0425, 0.0575])

pstemp = []
ppstemp = []
npstemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        pstemp.append(vs[i])
    elif state[i]==5:
        ppstemp.append(vs[i])
    elif state[i]==4:
        npstemp.append(vs[i])
ps = np.array(pstemp)
pps = np.array(ppstemp)
nps = np.array(npstemp)

#ax = fig.add_subplot(3, 1, 1)
plt.hist([ps,pps,nps], b, label = ['PS','PPS','NPS'])
plt.xticks([0.005, 0.020, 0.035, 0.050])
plt.xlabel("$V'_{S}$")
plt.ylabel("No. of Parameter Space Points")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 1.11),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D3411Ps.png",bbox_inches='tight')