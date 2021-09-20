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


filename = "Regstate.txt"
file_data = np.loadtxt(filename)
    
vsi = file_data[:, 2]
xi = file_data[:, 3]
yi = file_data[:, 4]
zi = file_data[:, 5]
xp = file_data[:, 6]
yp = file_data[:, 7]
zp = file_data[:, 8]
state = file_data[:, 10]
fig = plt.figure(figsize=(21, 28), dpi=300)

outer = gridspec.GridSpec(4, 3, left = 0.05, bottom = 0.05, right = 0.95, top = 0.95, wspace=0.3, hspace=0.3)
for j in range(0, len(vsi), 315):
    vsij = vsi[j:j + 314]
    xij = xi[j:j + 314]
    yij = yi[j:j + 314]
    zij = zi[j:j + 314]
    xpj = xp[j:j + 314]
    ypj = yp[j:j + 314]
    zpj = zp[j:j + 314]
    statej = state[j:j + 314]
    
    ind = np.int64(j/315)
    inner = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer[ind], wspace=0.5, hspace=0.5)

    for k in range(0, len(xij),105):
        vsijk = vsij[k:k + 104]
        xijk = xij[k:k + 104]
        yijk = yij[k:k + 104]
        zijk = zij[k:k + 104]
        xpjk = xpj[k:k + 104]
        ypjk = ypj[k:k + 104]
        zpjk = zpj[k:k + 104]
        statejk = statej[k:k + 104]
        for l in range(0, len(zijk),35):
            vsijkl = vsijk[l:l + 34]
            xijkl = xijk[l:l + 34]
            yijkl = yijk[l:l + 34]
            zijkl = zijk[l:l + 34]
            xpjkl = xpjk[l:l + 34]
            ypjkl = ypjk[l:l + 34]
            zpjkl = zpjk[l:l + 34]
            statejkl = statejk[l:l + 34]
                
            index = l/35 + 3*k/105 + 1
            index2 = np.int64(l/35)
            index1 = np.int64(k/105)
            ax = fig.add_subplot(inner[index1,index2], projection='3d')
            norm = PiecewiseNormalize([2,3,4,5,6],[0.15, 0.325, 0.50, 0.675, 0.85])
            
            p = ax.scatter3D(xpjkl,ypjkl,zpjkl,alpha = 1.0, norm = norm, s = 20, c=statejkl, cmap = 'viridis');

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

plt.savefig("D23R.png")
#plt.savefig(output)
#plt.show()
