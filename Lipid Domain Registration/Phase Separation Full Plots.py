import re
import string
import math
import numpy as np
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec

class PiecewiseNormalize(Normalize):
    def __init__(self, xvalues, cvalues):
        self.xvalues = xvalues
        self.cvalues = cvalues

        Normalize.__init__(self)

    def __call__(self, value, clip=None):
        if self.xvalues is not None:
            x, y = self.xvalues, self.cvalues
            return np.ma.masked_array(np.interp(value, x, y))
        else:
            return Normalize.__call__(self, value, clip)


#Main#
#------------------------------------------------------------------------------#


filename = "PSstate.txt"
file_data = np.loadtxt(filename)
    
vs = file_data[:, 2]
x = file_data[:, 3]
y = file_data[:, 4]
z = file_data[:, 5]
xp = file_data[:, 6]
yp = file_data[:, 7]
zp = file_data[:, 8]
state = file_data[:, 9]

fig = plt.figure(figsize=(21, 28), dpi=300)

outer = gridspec.GridSpec(4, 3, left = 0.05, bottom = 0.05, right = 0.95, top = 0.95, wspace=0.3, hspace=0.3)
for i in range(0, len(vs), 315):
    vsi = vs[i:i + 314]
    xi = x[i:i + 314]
    yi = y[i:i + 314]
    zi = z[i:i + 314]
    xpi = xp[i:i + 314]
    ypi = yp[i:i + 314]
    zpi = zp[i:i + 314]
    statei = state[i:i + 314]
    
    ind = np.int64(i/315)
    inner = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer[ind], wspace=0.5, hspace=0.5)

    for j in range(0, len(xi),105):
        vsij = vsi[j:j + 104]
        xij = xi[j:j + 104]
        yij = yi[j:j + 104]
        zij = zi[j:j + 104]
        xpij = xpi[j:j + 104]
        ypij = ypi[j:j + 104]
        zpij = zpi[j:j + 104]
        stateij = statei[j:j + 104]
        for k in range(0, len(zij),35):
            vsijk = vsij[k:k + 34]
            xijk = xij[k:k + 34]
            yijk = yij[k:k + 34]
            zijk = zij[k:k + 34]
            xpjk = xpij[k:k + 34]
            ypjk = ypij[k:k + 34]
            zpjk = zpij[k:k + 34]
            statejk = stateij[k:k + 34]
                
            index = k/35 + 3*j/105 + 1
            index2 = np.int64(k/35)
            index1 = np.int64(j/105)
            ax = fig.add_subplot(inner[index1,index2], projection='3d')
            norm = PiecewiseNormalize([4,5,6],[0.25, 0.575, 0.9])
            
            p = ax.scatter3D(xpjk,ypjk,zpjk,alpha = 1.0, norm = norm, s = 20, c=statejk, cmap = 'viridis');

            ax.set_xlim(3, 17)
            ax.set_ylim(3, 17)
            ax.set_zlim(3, 17)
            ax.set_xticks((4,7,10,13,16))
            ax.set_yticks((4,7,10,13,16))
            ax.set_zticks((4,7,10,13,16))
            ax.set_xlabel('$Vp_{AA}$')
            ax.set_ylabel('$Vp_{BB}$')
            ax.set_zlabel('$Vp_{AB}$')

            # Customize the view angle so it's easier to see the scatter points
            ax.view_init(elev=40.0, azim=-45.0)


fig.tight_layout()

plt.savefig("D23P.png")
#plt.savefig(output)
#plt.show()
