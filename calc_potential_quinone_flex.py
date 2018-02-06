#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_potential.py
# Description: This is a python script to calculate the electrical potential distribution along the channel for armchair edge
# Date: 09-23-2016 Created

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

if comm.args.suffix != None:
    outname = 'pd_'+comm.args.suffix
else:
    outname = 'pd'

filename = comm.args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)

if comm.args.sc != None:
    sigma = comm.args.sc
else:
    sigma = float(raw_input('Please input the surface charge density (C/m^2):'))

sigma1 = 0
sigma2 = sigma

## Instead of a predetermined number, the exb will be determined by the charge density distribution of the electrode.
ele = comm.load_xvg('density_electrode.xvg')
poe1 = sum(ele[:, 0]*ele[:, 1])/sum(ele[:, 1])
poe2 = sum(ele[:, 0]*ele[:, 2])/sum(ele[:, 2])
#if sigma == 0:
#    poe1, poe2 = 0, 0
#else:
#    ele = comm.load_xvg('density_electrode.xvg')
#    poe1 = sum(ele[:, 0]*ele[:, 1])/sum(ele[:, 1])
#    poe2 = sum(ele[:, 0]*ele[:, 2])/sum(ele[:, 2])
print 'The location of electrodes: ', poe1, poe2

if filename.endswith('.txt'):
    data = np.loadtxt(filename)
    if comm.args.bound:
        bound = [float(i) for i in comm.args.bound]
    else:
        bound = [0, max(data[:, 0])]
    rho = data[:,2]
elif filename.endswith('.xvg'):
    if comm.args.bound:
        bound = comm.args.bound
    else:
        bound = raw_input('Please input the boundary (nm):').split()
    bound = [float(i) for i in bound]
    data = comm.load_xvg(filename)
    for i, item in enumerate(data):
        if item[0] >= bound[1]:
            break
    data = data[:i]
    rho = data[:,1]
z = data[:,0]
print sigma,
exb = [poe1, bound[1]-poe2]
for seq, myrho in enumerate([rho, rho[::-1]]):
    exb_1 = 0
    exb_2 = exb[seq]
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
