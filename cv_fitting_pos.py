#!/Users/yuzhang/anaconda/bin/python
# Filename: cv_fitting.py
# Description: This is a python script to fit the sigma-potential curve with polynomial
# Date: 02-17-2016 Created

import argparse
import numpy as np
import matplotlib.pyplot as plt

fsize = 28
lwidth = 4.0
msize = 12


parser = argparse.ArgumentParser(description = 'User specified filenames, etc.')
parser.add_argument('-i', '--input', dest = 'filename', default = 'cv.txt', help = 'Filename for sigma and potential data')
parser.add_argument('--suffix', dest = 'suffix', help = 'suffix for output name')
parser.add_argument('-n', '--number', dest = 'nrun', type = int, help = 'Number of runs')
args = parser.parse_args()

if args.suffix == None:
    outname = 'DC'
else:
    outname = 'DC_'+args.suffix

data = np.loadtxt(args.filename)
sigma = np.concatenate(( [0], data[:, 0]))
pot = np.concatenate(([0], data[:, 1]))

fig, (ax1, ax2) = plt.subplots(2, sharex = True, figsize = (16, 12), dpi = 1000) 
ax1.plot(pot, sigma, 'o', markersize = msize, label = 'raw data')

p_fit = np.linspace(pot[0]-0.1, pot[-1]+0.1, 50)

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
    print 100*(f(pot[-1])-f(pot[2]))/(pot[-1]-pot[2])
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
