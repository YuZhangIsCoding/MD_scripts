#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_diff_nvcf.py
# Description:  This is a python script to calculate the diffusion coefficients by integrating the normalized velocity autocorrelation functions
#               C(t) = <v(0)*v(t)>/<v(0)*v(0)>
#               Using the fact that <v(0)*v(0)> = 3kT/m, the diffusion coefficients could be calculated by
#               D =  kT/m * integral(C(t)dt) from 0 to infinity
#                   k: Boltzmann constant 1.38064852 J/K, J = kg*m^2/s^2
#                   T: Temperature
#                   m: molecular mass
# Date:     03-23-2017 Created

import matplotlib.pyplot as plt
import numpy as np
import pdb, sys
import calc_common as comm

fig = plt.figure(figsize = (16, 12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

if comm.args.filename.endswith('.txt'):
    nvcf = np.loadtxt(comm.args.filename)
else:
    print 'File', comm.args.traj_name, 'not suitable for this script!'
    sys.exit()

if comm.args.mmass == None:
    mmass = [float(i) for i in raw_input('Please input the molar mass for %d group(s):' %(len(nvcf[0])-1)).split()]
else:
    mmass = comm.args.mmass
if len(mmass) != len(nvcf[0])-1:
    sys.exit('Exit: please input the exact number of molmasses for the groups selected!')

if comm.args.suffix == None:
    outname = 'Diff_nvcf'
else:
    outname = 'Diff_nvcf_'+comm.args.suffix

const = 1.38064852*1000*6.022*comm.args.temperature
print 'Calculate the diffusion coefficients at %f K' %comm.args.temperature

diff = np.zeros(nvcf.shape)
diff[:, 0] = nvcf[:, 0]
pt = int(0.9*len(diff))

for i, item in enumerate(nvcf):
    for j in range(1, len(item)):
        diff[i, j] = np.trapz(nvcf[:i+1, j], nvcf[:i+1, 0])*1e-12*const/mmass[j-1]
for i in range(1, len(diff[0])):
    plt.plot(diff[:, 0], diff[:, i], lw = lwidth, label = str(i))
    print 'Diffusion coefficient for group %d (m^2/s):' %i, np.mean(diff[pt:, i]), '+/-', np.std(diff[pt:, i])

plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ps)', fontsize = fsize)
plt.ylabel('Diffusion Coefficients ($\mathregular{m^2/s}$)', fontsize = fsize)
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
