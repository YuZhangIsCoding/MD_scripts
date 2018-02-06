#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_potential.py
# Description: This is a python script to calculate the electrical potential distribution along the channel for armchair edge
# Date: 09-23-2016 Created

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

parser = argparse.ArgumentParser(description = 'User specified density profile, surface charge density, etc.')
parser.add_argument('-i', '--input', dest = 'filename', type = str, default = 'Density.txt', help = 'File that contains charge density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('-ms', '--multsigma', dest = 'ms', nargs = '*', type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('--exb', dest = 'exb', nargs = '*', help = 'The location of the two charge plane relative to the substrate')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', help = 'Boundary of the distribution')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
parser.add_argument('--sp', dest = 'sp', type = int, default = 20, help = 'Data points to be split')
parser.add_argument('-c', '--config', dest = 'config', choices = ['Q', 'H2Q', 'Q-0.25', 'Q-0.5', 'Q-0.75'], help = 'Configuration of quinone')
args = parser.parse_args()

if args.suffix != None:
    outname = 'pd_'+args.suffix
else:
    outname = 'pd'

if args.exb == None:
    if args.config == 'Q':
        exb_1 = 0
        exb_2 = 0.343583333
    elif args.config == 'H2Q':
        exb_1 = 0
        exb_2 = 0.346857142857
    elif args.config == 'Q-0.25':
        exb_1 = 0
        exb_2 = 0.348592592593
    elif args.config == 'Q-0.5':
        exb_1 = 0
        exb_2 = 0.348
    else:
        exb_1 = 0
        exb_2 = 0.35758
else:
    exb_1 = args.exb[0]
    exb_2 = args.exb[1]

filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)

if args.ms == None:
    if args.sigma != None:
        sigma = args.sigma
    else:
        sigma = float(raw_input('Please input the surface charge density (C/m^2):'))
    sigma1 = 0
    sigma2 = sigma
else:
    sigma1 = args.ms[0]
    sigma2 = args.ms[1]

data = []
if filename.endswith('.txt'):
    data = np.loadtxt(filename)
    if args.bound != None:
        bound = [float(i) for i in args.bound]
    else:
        bound = [0, max(data[:, 0])]
    rho = data[:,2]
elif filename.endswith('.xvg'):
    if args.bound != None:
        bound = args.bound
    else:
        bound = raw_input('Please input the boundary (nm):').split()
    bound = [float(i) for i in bound]
    file_in = open(filename, 'r')
    for line in file_in:
        if line[0] != '#' and line[0] != '@':
            temp = [float(i) for i in line.split()]
            if temp[0] < bound[0]+exb_1:
                continue
            elif temp[0] > bound[1]-exb_1:
                break
            else:
                data.append([float(i) for i in line.split()])
    file_in.close()
    data = np.array(data)
    rho = data[:,1]
z = data[:,0]
print sigma,
for seq, myrho in enumerate([rho, rho[::-1]]):
    phi_1 = np.zeros(len(z))
    phi_2 = np.zeros(len(z))
    const = 1000/8.85418784762/6.242
    for i, u in enumerate(z):
        if z[i] < exb_1:
            phi_1[i] = 0
        elif z[i] > bound[1]-exb_1:
            phi_1[i] = phi_1[i-1]
        else:
            phi_1[i] = -np.trapz((u-z[:i+1])*myrho[:i+1], z[:i+1])*const-sigma1*(z[i]-exb_1)*const/0.160217657
    for i, u in enumerate(z):
        if z[i] < exb_2:
            phi_2[i] = 0
        elif z[i] > bound[1]-exb_2:
            phi_2[i] = phi_2[i-1]
        else:
            phi_2[i] = -sigma2*(-1)**seq*(z[i]-exb_2)*const/0.160217657
    phi = phi_1+phi_2
    bulk = np.mean(phi[len(z)/2-50:len(z)/2+50])

    for i, item in enumerate(z):
        if item >= exb_2:
            wall = phi[i-1]+(phi[i]-phi[i-1])*(exb_2-z[i-1])/(item-z[i-1])
            break
    print wall-bulk, 

    hl = len(phi)/4*3
    if seq == 0:
        plt.plot(z[:hl], phi[:hl], color = 'pink', linewidth = lwidth)
        plt.plot([0, z.max()],[bulk, bulk],'--', color = 'red')#, linewidth = lwidth)
        plt.plot([exb_2, exb_2], [min(phi[:hl]), max(phi[:hl])], '--', color = 'grey')
    else:
        plt.plot(z[-hl:][::-1], phi[:hl], color = 'cornflowerblue', linewidth = lwidth)
        plt.plot([0, z.max()],[bulk, bulk],'--', color = 'blue')
        plt.plot([bound[1]-exb_2, bound[1]-exb_2], [min(phi[-hl:]), max(phi[-hl:])], '--', color = 'grey')

for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.xlim([0, z.max()])
plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.tight_layout()
plt.savefig(outname+'.pdf')

np.savetxt(outname+'.txt', np.column_stack((z, phi)), fmt = ['%12.4f', '%12.4f'])
