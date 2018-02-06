#!/Users/yuzhang/anaconda/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
import pdb

data = []
for dirname in os.listdir("./"):
    if 'h' in dirname and 'prod' in dirname:
        os.chdir(dirname)
        for filename in os.listdir("./"):
            if filename == 'NumberDensity.txt':
                print dirname
                if 'pzc' in dirname:
                    data.insert(0, np.loadtxt(filename))
                else:
                    data.append(np.loadtxt(filename))
        os.chdir("../")
data = np.array(data)
x = data[0][:,0]
CAT = []
AN = []
for i in data:
    if len(i) == len(x):
        CAT.append(i[:,1])
        AN.append(i[:,2])
CAT = np.array(CAT)
AN = np.array(AN)
y = np.arange(len(CAT))
X, Y = np.meshgrid(x, y)
fig = plt.figure(figsize=(16, 12),dpi=1000)
ax1 = fig.add_subplot(211, projection='3d')
ax2 = fig.add_subplot(212, projection='3d')

surf1 = ax1.plot_surface(X, Y, CAT ,rstride =2 , cstride =8, alpha=0.5, color='blue', label= 'CAT')
surf2 = ax2.plot_surface(X, Y, AN ,rstride =2 , cstride =8, alpha=0.5, color='red', label='Anion')
#cset = ax.contour(X, Y, Z, zdir='z', offset=-140, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
#ax.plot_wireframe(X, Y, Z ,rstride=1, cstride=1)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
#fig.colorbar(surf, shrink=0.5, aspect=5)
ax1.set_xlabel('X')
#ax1.set_xlim(0,5)
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax2.set_xlabel('X')
#ax2.set_xlim(0,5)
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')

def select_view():
    ''' Select from different angles of view'''
    for ii in xrange(0,360,30):
        ax1.view_init(elev=10., azim=ii)
        ax2.view_init(elev=10., azim=ii)
        plt.tight_layout()
        plt.savefig('blar%d.pdf' %ii)

proxy = [plt.Rectangle((0,0),1,1, fc = 'blue', alpha =0.5)]
ax1.legend(proxy, ['CAT'],frameon= False)
proxy = [plt.Rectangle((0,0),1,1, fc = 'red', alpha =0.5)]
ax2.legend(proxy, ['ANION'], frameon= False)

ax1.view_init(elev=15., azim=-85)
ax2.view_init(elev=15., azim=-85)
plt.tight_layout()
plt.savefig('blar.pdf')

ax1.view_init(elev=5., azim=-90)
ax2.view_init(elev=5., azim=-90)
plt.tight_layout()
plt.savefig('blar90.pdf')
