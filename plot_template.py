import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0
# setting color cycle for the plots automatically
#mycmp = plt.cm.gist_ncar
#plt.gca().set_color_cycle([mycmp(i) for i in np.linspace(0, 0.9, 10)])

plt.plot(x, y, '-o', color = 'red', markeredgecolor = 'none', linewidth = lwidth, label = 'blar')
plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel(' ', fontsize = fsize)
plt.ylabel(' ', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig()
