import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('input.txt')

fig = plt.figure(figsize = (4,4), dpi = 300)
## fontsize, linewidth, markersize
fs = 12 
lw = 1.0
ms = 6.0

## setting color cycle 
#mycmp = plt.cm.gist_ncar
#plt.gca().set_color_cycle([mycmp(i) for i in np.linspace(0, 0.9, 10)])

plt.plot(data[:, 0], data[:, 1], lw = lw, ms = ms, mfc = 'w', mec = 'k', label = ' ')

plt.xlabel(' ', fontsize = fs)
plt.ylabel(' ', fontsize = fs)
plt.legend(loc = 'best', fontsize = fs, numpoints = 1)

## setting up ticks format
#for ax in fig.axes:
#    ax.tick_params(labelsize = fs, width = 0.5, length = 4, pad = 6)
#    ax.tick_params(which = 'minor', length = 1, width = 0.25)
#    for direction in ['top', 'bottom', 'left', 'right']:
#        ax.spines[direction].set_linewidth(1)

plt.tight_layout(pad = 0.4, w_pad = 0.5, h_pad = 1.0)
#fig.set_tight_layout({'pad' : 0.4, 'w_pad' : 0.5, 'h_pad' : 1.0})
plt.savefig('outname.pdf')
#plt.savefig('outname.pdf', bbox_inches='tight')
