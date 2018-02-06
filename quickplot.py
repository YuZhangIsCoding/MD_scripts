#!/Users/yuzhang/anaconda/envs/py3/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import calc_common as comm
import os, sys

args = comm.args

if args.filename.endswith('.txt'):
    try:
        data = np.loadtxt(args.filename)
    except:
        df = pd.read_cvs(args.filename, header = None)
        data = df.values
elif args.filename.endswith('.xvg'):
    data = comm.load_xvg(args.filename)
else:
    print('Currently only support .xvg and .txt files')
    sys.exit()

if args.bound != None:
    bound = args.bound
else:
    bound = [min(data[:, 0]), max(data[:, 0])]

if args.sub == True:
    fig, axes = plt.subplots(len(data[0])-1)
    for i in range(1, len(data[0])):
        axes[i-1].plot(data[:,0],data[:,i], label = 'column'+str(i))
        axes[i-1].legend(loc = 'best')
        axes[i-1].grid()
        axes[i-1].set_xlim((bound[0], bound[1]))
else:
    for i in range(1, len(data[0])):
        plt.plot(data[:,0],data[:,i], label = 'column'+str(i))
        plt.legend(loc = 'best')
        plt.xlim([bound[0], bound[1]])
        plt.grid()

plt.savefig('quicklook.pdf')
os.system('open quicklook.pdf')
