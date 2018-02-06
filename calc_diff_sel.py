#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_diff_sel.py
# Description:  This is a python script to calculate the diffusion coefficients using Einstein Equation on selected particles.
#               In QENS experiments, they can only detech particles that "mobile" enough, which is judged roughly
#               by their jump length in a certain time. For example, the jump length for OMIm is about 1.57 angstrom to 15.7 angstrom
#               This code tries to use similar way of selecting particles, using the final MSD of each molecule to determine whether
#               a particle is mobile or immobile.
#               If the MSD of a particle is higher than (pi/Q)^2, then it is considered mobile, and the MSD will be averaged.
# Date  :   03-28-2017  Created

import matplotlib.pyplot as plt
import numpy as np
import sys
import calc_common as comm

fig = plt.figure(figsize = (16, 12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

if comm.args.dim == None:
    dim = 6
    #print 'Diffusion coefficients in all 3 directions!'
else:
    dim = int(comm.args.dim[0])*2
    #print 'Diffusion coefficients in %d directions!' %(dim/2)

if comm.args.filename.endswith('.txt'):
    msd_all = np.loadtxt(comm.args.filename)
elif comm.args.filename.endswith('.xvg'):
    msd_all = comm.load_xvg(comm.args.filename)
else:
    #print 'File', comm.args.traj_name, 'not suitable for this script!'
    sys.exit()

if comm.args.suffix == None:
    outname = 'Diff'
else:
    outname = 'Diff_'+comm.args.suffix

#for i, item in enumerate(msd_all[:, 0]):
#    if item >= 400:
#        mark = i
#        break


mobile = []
for i in range(1, len(msd_all[0])):
    if msd_all[-1][i] >= (np.pi/comm.args.q/10)**2:
        mobile.append(i)
#print len(mobile), 'out of', len(msd_all[0])-1, 'selected'

msd = np.column_stack((msd_all[:, 0], np.mean(msd_all[:, mobile], axis = 1)))

if comm.args.scale == None:
    scale = [0.2, 0.9]
elif len(comm.args.scale) == 1:
    scale = [comm.args.scale[0], 0.9]
else:
    scale = comm.args.scale
#print 'Fitting from %f to %f' %(scale[0], scale[1])

pt_start = int(scale[0]*len(msd))
pt_end = int(scale[1]*len(msd))

label_name = ['OMIm', 'TFSI', 'Solvent']
#label_name = ['TEA', 'BF4', 'ACN']
colors = ['red', 'blue', 'green']
for i in range(1, len(msd[0])):
    coef = np.polyfit(msd[pt_start:pt_end, 0], msd[pt_start:pt_end, i], 1)
    #print 'Diffusion coefficient for group %d is %.4E m^2/s' %(i, coef[0]/dim/1e6)
    fit = np.poly1d(coef)
    plt.plot(msd[:, 0]/1000, msd[:, i], linewidth = lwidth, label = label_name[i-1], color = colors[i-1])
#    plt.plot(msd[:, 0]/1000, fit(msd[:, 0]), '--', linewidth = lwidth, label = 'fitting '+label_name[i-1])

print comm.args.q, len(mobile), len(msd_all[0])-1, '%.4E ' %(coef[0]/dim/1e6)

plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ns)', fontsize = fsize)
plt.ylabel('MSD ($\mathregular{nm^2\!}$)', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig(outname+'.pdf')
np.savetxt(outname+'.txt', msd, fmt = ['%12.4f' for _ in msd[0]])
