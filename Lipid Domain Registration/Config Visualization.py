import re
import string
import math
import numpy as np
import matplotlib.pyplot as plt
import cv2
import glob
#Parameter Declarations#
#----------------------------------------------------------------------------- #

rows = 100
columns = 100

#Main Code#
#----------------------------------------------------------------------------- #    

input = ["1.txt","2.txt","3.txt","4.txt"]
t=0
for file in input:
    t+=1
    f=open(file,"r").readlines()
    u = np.zeros(100)
    l = np.zeros(100)
    # loop to move from one configuration to the next in the input file
    for i in range(2,len(f),101):
        upper, lower = [],[]
        # reading one configuration in the file into 'upper' and 'lower' 2D arrays
        for j in range(0,100):
            if j==0:
                for k in range(0,100):
                   upper.append(int(f[i+j][k]))
                   lower.append(int(f[i+j][k+101]))
            else:
                for k in range(0,100):
                    u[k] = f[i+j][k]
                    l[k] = f[i+j][k+101]
                upper = np.vstack((upper,u))                    
                lower = np.vstack((lower,l))
        diff = np.zeros((100,100))
        for j in range(0,100):
            for k in range(0,100):
                diff[j][k] = abs(upper[j][k] - lower[j][k])
    
        plt.figure(figsize=(10,3), dpi=300)
        plt.subplot(1,3,1)
        plt.imshow(upper, cmap = 'viridis')
        plt.title('Upper Leaflet')
        plt.subplot(1,3,2)
        plt.imshow(lower, cmap = 'viridis')
        plt.title('Lower Leaflet')
        plt.subplot(1,3,3)
        plt.imshow(diff, cmap = 'summer')
        plt.title('Difference')
        image = str(t)+"move_"+str(i)+"_diff.png"
        plt.savefig(image)

            