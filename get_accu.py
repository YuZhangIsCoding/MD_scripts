#!/Users/yuzhang/anaconda/bin/python
import numpy as np

pos = np.loadtxt('NDinSlit_POS.txt')
neg = np.loadtxt('NDinSlit_NEG.txt')

n_p = (pos[1, 0]-pos[0, 0])*sum(pos[:,1:])
n_n = (neg[1, 0]-neg[0, 0])*sum(neg[:, 1:])
print n_p[0], n_p[1], n_n[0], n_n[1]
