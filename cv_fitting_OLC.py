#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: cv_fitting_raw.py
# Description:  This is a python script to fit the sigma-potential curve with polynomial
#               This script is built on the base of cv_fitting.py but uses different data structrue
#               The input data should start from pzc to charged system:
#               -0.10 ...
#               -0.09 ...
#                 ...
#                   0 ...
#                 ...
#                0.09 ...
#                0.10 ...
# Date: 11-28-2017 Created

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
fsize = 28
lwidth = 4.0
msize = 12


parser = argparse.ArgumentParser(description = 'User specified filenames, etc.')
parser.add_argument('-i', '--input', dest = 'filename', default = 'cv.txt', help = 'Filename for sigma and potential data')
parser.add_argument('--suffix', dest = 'suffix', help = 'suffix for output name')
parser.add_argument('-n', '--number', dest = 'nrun', type = int, help = 'Number of runs')
parser.add_argument('-p', '--potential', dest = 'pot', nargs = '*', type = float, default = [-0.5, 0.5], help = 'Potential range for integral capacitance')
parser.add_argument('-t', '--type', dest = 'type', type = int, default = 2, help = 'The type of cv file')
args = parser.parse_args()

if args.suffix == None:
    outname = 'DC'
else:
    outname = 'DC_'+args.suffix

df = pd.read_csv(args.filename, names = ['sigma', 'pot'], sep = '\s+')
sigma = df.sigma
pot = df.pot-df.pot[df.sigma == 0].values

fig, (ax1, ax2) = plt.subplots(2, sharex = True, figsize = (16, 12), dpi = 1000) 
ax1.plot(pot, sigma, 'o', markersize = msize, label = 'raw data')

p_fit = np.linspace(min(pot)-0.1, max(pot)+0.1, 50)

dc = []
for i in range(3, 7):
    z = np.polyfit(pot, sigma, i)
    f = np.poly1d(np.polyfit(pot, sigma, i))
    f = np.poly1d(z)
    df = np.poly1d(z[:-1]*np.arange(len(z)-1, 0, -1))
    sigma_fit = f(p_fit)
    dc.append(100*df(np.sort(pot)))
    ax1.plot(p_fit, sigma_fit, '-',linewidth = lwidth, label = str(i)+' fitting')
    ax2.plot(np.sort(pot), dc[-1], '-o', linewidth = lwidth, markersize = msize,label = str(i)+' fitting')
    for j in range(int(len(args.pot)/2)):
        print('Potential from %.2f to %.2f V:' %(args.pot[2*j], args.pot[2*j+1]), 100*(f(args.pot[2*j+1])-f(args.pot[2*j]))/(args.pot[2*j+1]-args.pot[2*j]))
ax1.legend(fontsize = 28, frameon = False, numpoints = 1, bbox_to_anchor = [1, 0.7])
ax1.set_ylabel('Surface charge density\n($\mathregular{C/m^2\!}$)', fontsize = fsize)
ax2.legend(fontsize = 28, frameon = False, numpoints = 1, bbox_to_anchor = [1, 0.5])
ax2.set_xlabel('Potential (V)', fontsize = fsize)
ax2.set_ylabel('Differential capacitance\n($\mathregular{\mu F/cm^2\!}$)', fontsize = fsize)
for axis in fig.axes:
    axis.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axis.spines[direction].set_linewidth(2)
        axis.grid(b = True, which = 'major', color = 'grey', linestyle = '--')
plt.tight_layout()
fig.savefig(outname+'.pdf')

fmt = ['%12.4f' for row in range(1+len(dc))]
dc = np.transpose(np.array(dc))
np.savetxt(outname+'.txt', np.column_stack((np.sort(pot), dc)), fmt = fmt)

if args.nrun != None:
    fig = plt.figure(figsize = (16, 12), dpi = 1000) 
    pot_avg = []
    dc_avg = []
    ll = len(pot)/args.nrun/2
    for i in range(ll):
        pot_avg.append(np.mean([pot[i+j*ll] for j in range(args.nrun)]))
    pot_avg.append(0)
    for i in range(ll):
        pot_avg.append(np.mean([pot[ll*args.nrun+i+j*ll+1] for j in range(args.nrun)]))
    for i in range(3, 7):   
        z = np.polyfit(pot, sigma, i)
        f = np.poly1d(np.polyfit(pot, sigma, i))
        f = np.poly1d(z)
        df = np.poly1d(z[:-1]*np.arange(len(z)-1, 0, -1))
        sigma_fit = f(p_fit)
        dc_avg.append(100*df(pot_avg))
        plt.plot(pot_avg, dc_avg[-1], '-o', linewidth = lwidth, markersize = msize,label = str(i)+' fitting')
    for axis in fig.axes:
        axis.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
        for direction in ['top', 'bottom', 'left', 'right']:
            axis.spines[direction].set_linewidth(2)
            axis.grid(b = True, which = 'major', color = 'grey', linestyle = '--')
    plt.legend(fontsize = 28, frameon = False, numpoints = 1, bbox_to_anchor = [1, 0.7])
    plt.savefig(outname+'_avg.pdf')
    fmt = ['%12.4f' for row in range(1+len(dc_avg))]
    dc_avg = np.transpose(np.array(dc_avg))
    np.savetxt(outname+'_avg.txt', np.column_stack((pot_avg, dc_avg)), fmt = fmt)
