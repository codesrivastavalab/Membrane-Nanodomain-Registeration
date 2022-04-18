"""
Created on Sat Dec 25 11:23:50 2021

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


filename = "Regstate.txt"
file_data = np.loadtxt(filename)
    
vs = file_data[:, 2]
aa = file_data[:, 3]
bb = file_data[:, 4]
ab = file_data[:, 5]
aap = file_data[:, 6]
bbp = file_data[:, 7]
abp = file_data[:, 8]
state = file_data[:, 10]

fig = plt.figure(figsize=(8, 12), dpi=300)

# Bins for enthalpic V'
b = np.array([2.0, 14.0, 26.0, 38.0])

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(aa[i])
    elif state[i]==5:
        prtemp.append(aa[i])
    elif state[i]==4:
        urtemp.append(aa[i])
    elif state[i]==3:
        partemp.append(aa[i])
    elif state[i]==2:
        artemp.append(aa[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 1)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{AA}$")

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(bb[i])
    elif state[i]==5:
        prtemp.append(bb[i])
    elif state[i]==4:
        urtemp.append(bb[i])
    elif state[i]==3:
        partemp.append(bb[i])
    elif state[i]==2:
        artemp.append(bb[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 2)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{BB}$")
ax.set_ylabel("No. of Parameter Space Points")

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(ab[i])
    elif state[i]==5:
        prtemp.append(ab[i])
    elif state[i]==4:
        urtemp.append(ab[i])
    elif state[i]==3:
        partemp.append(ab[i])
    elif state[i]==2:
        artemp.append(ab[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 3)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([8.0, 20.0, 32.0])
ax.set_xlabel("$V'_{AB}$")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 3.525),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D2311R'.png",bbox_inches='tight')

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fig = plt.figure(figsize=(8, 12), dpi=300)

# Bins for enthalpic V
b = np.array([2.5, 5.5, 8.5, 11.5, 14.5, 17.5])

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(aap[i])
    elif state[i]==5:
        prtemp.append(aap[i])
    elif state[i]==4:
        urtemp.append(aap[i])
    elif state[i]==3:
        partemp.append(aap[i])
    elif state[i]==2:
        artemp.append(aap[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 1)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{AA}$")

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(bbp[i])
    elif state[i]==5:
        prtemp.append(bbp[i])
    elif state[i]==4:
        urtemp.append(bbp[i])
    elif state[i]==3:
        partemp.append(bbp[i])
    elif state[i]==2:
        artemp.append(bbp[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 2)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{BB}$")
ax.set_ylabel("No. of Parameter Space Points")

rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(abp[i])
    elif state[i]==5:
        prtemp.append(abp[i])
    elif state[i]==4:
        urtemp.append(abp[i])
    elif state[i]==3:
        partemp.append(abp[i])
    elif state[i]==2:
        artemp.append(abp[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

ax = fig.add_subplot(3, 1, 3)
ax.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
ax.set_xticks([4.0, 7.0, 10.0, 13.0, 16.0])
ax.set_xlabel("$V_{AB}$")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 3.525),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D2311R.png",bbox_inches='tight')

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fig = plt.figure(figsize=(8, 4), dpi=300)

# Bins for enthalpic V
b = np.array([0, 0.0125, 0.0275, 0.0425, 0.0575])


rtemp = []
prtemp = []
urtemp = []
partemp = []
artemp = []
for i in range(0, len(vs)):
    if state[i]==6:
        rtemp.append(vs[i])
    elif state[i]==5:
        prtemp.append(vs[i])
    elif state[i]==4:
        urtemp.append(vs[i])
    elif state[i]==3:
        partemp.append(vs[i])
    elif state[i]==2:
        artemp.append(vs[i])
r = np.array(rtemp)
pr = np.array(prtemp)
ur = np.array(urtemp)
par = np.array(partemp)
ar = np.array(artemp)

#ax = fig.add_subplot(3, 1, 2)
plt.hist([r,pr,ur,par,ar], b, label = ['R','PR','UR','PAR','AR'])
plt.xticks([0.005, 0.020, 0.035, 0.050])
plt.xlabel("$V'_{S}$")
plt.ylabel("No. of Parameter Space Points")
plt.legend(loc='upper center',bbox_to_anchor=(0.48, 1.11),fancybox=True,fontsize=10,ncol = 5)
plt.savefig("D2311Rs.png",bbox_inches='tight')