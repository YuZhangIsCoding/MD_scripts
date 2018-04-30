#!/Users/yuzhang/anaconda3/bin/python
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
else:
    dim = int(comm.args.dim[0])*2

if comm.args.filename.endswith('.txt'): msd = np.loadtxt(comm.args.filename)
elif comm.args.filename.endswith('.xvg'):
    msd = comm.load_xvg(comm.args.filename)
else:
    print('File', comm.args.filename, 'not suitable for this script!')
    sys.exit()

if comm.args.scale == None:
    scale = [0.2, 0.9]
elif len(comm.args.scale) == 1:
    scale = [comm.args.scale[0], 0.9]
else:
    scale = comm.args.scale

if comm.args.v:
    print('Diffusion coefficients in %d directions!' %(dim//2))
    print('Fitting from %f to %f' %(scale[0], scale[1]))

pt_start = int(scale[0]*len(msd))
pt_end = int(scale[1]*len(msd))

#label_name = ['OMIm', 'TFSI', 'Solvent']
label_name = ['TEA', 'BF4', 'ACN']
colors = ['red', 'blue', 'green']

for i in range(1, len(msd[0])):
    coef = np.polyfit(msd[pt_start:pt_end, 0], msd[pt_start:pt_end, i], 1)
    if comm.args.v:
        print('Diffusion coefficient for group %d is %.4E m^2/s' %(i, coef[0]/dim/1e6))
    else:
        print(coef[0]/dim/1e6, end = ' ')
    fit = np.poly1d(coef)
    plt.plot(msd[:, 0]/1000, msd[:, i], linewidth = lwidth, label = label_name[i-1], color = colors[i-1])
#    plt.plot(msd[:, 0]/1000, fit(msd[:, 0]), '--', linewidth = lwidth, label = 'fitting '+label_name[i-1])
print('')
plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ns)', fontsize = fsize)
plt.ylabel('MSD ($\mathregular{nm^2\!}$)', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('Diff.pdf')
