import numpy as np

pos = np.loadtxt('NDinSlit_POS.txt')
neg = np.loadtxt('NDinSlit_NEG.txt')

print sum(pos[:,1:3]*(pos[1,0]-pos[0,0])),sum(neg[:,1:3]*(neg[1,0]-neg[0,0]))
