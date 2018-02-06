#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_diff_vcf.py
# Description:  This is a python script to calculate the diffusion coefficients by integrating the velocity autocorrelation functions
#               D = 1/3 * integral(<v(0)v(t)>dt) from 0 --> infinity
# Date:     03-23-2017 Created

import matplotlib.pyplot as plt
import numpy as np
import pdb, sys
import calc_common as comm

fig = plt.figure(figsize = (16, 12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

if comm.args.suffix == None:
    outname = 'Diff_vcf'
else:
    outname = 'Diff_vcf_'+comm.args.suffix

if comm.args.filename.endswith('.txt'):
    vcf = np.loadtxt(comm.args.filename)
else:
    print 'File', comm.args.traj_name, 'not suitable for this script!'
    sys.exit()

if comm.args.dim == None:
    dim = 3
    print 'Diffusion coefficients in all 3 directions!'
else:
    dim = int(comm.args.dim[0])
    print 'Diffusion coefficients in %d directions!' %dim

diff = np.zeros(vcf.shape)
diff[:, 0] = vcf[:, 0]
pt = int(0.9*len(diff))

for i, item in enumerate(vcf):
    for j in range(1, len(item)):
        diff[i, j] = np.trapz(vcf[:i+1, j], vcf[:i+1, 0])/dim*1e-12
for i in range(1, len(diff[0])):
    plt.plot(diff[:, 0], diff[:, i], lw = lwidth, label = str(i))
    print 'Diffusion coefficient for group %d (m^2/s):' %i, np.mean(diff[pt:, i]), '+/-', np.std(diff[pt:, i])



plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ps)', fontsize = fsize)
plt.ylabel('Diffusion Coefficients ($\mathregular{m^2\!/s}$)', fontsize = fsize)
plt.yscale('log')
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.grid()
plt.tight_layout()
plt.savefig(outname+'.pdf')
np.savetxt(outname+'.txt', diff, fmt = ['%14.6E' for _ in diff[0]])
