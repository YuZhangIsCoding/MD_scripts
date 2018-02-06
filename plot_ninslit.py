import matplotlib.pyplot as plt
import numpy as np
import pdb

fig, (ax2, ax1) = plt.subplots(2, figsize = (12,12), dpi = 1000)
fsize = 28
lwidth = 4.0
msize = 16.0
data = np.loadtxt('NumINslit.txt')

sec_1 = 50
sec_2 = 500

ax1.plot(data[:sec_2, 0], data[:sec_2, 1], 'r-', linewidth = lwidth, label = 'OMI - anode')
ax1.plot(data[:sec_2, 0], data[:sec_2, 2], 'b-', linewidth = lwidth, label = 'TFSI - anode')
ax1.plot(data[:sec_2, 0], data[:sec_2, 3], '--', color = 'pink', linewidth = lwidth, label = 'OMI - cathode')
ax1.plot(data[:sec_2, 0], data[:sec_2, 4], '--', color = 'cornflowerblue', linewidth = lwidth, label = 'TFSI - cathode')

ax1.plot(data[0, 0], data[0, 1], 'o', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'red', markeredgewidth = 2)
ax1.plot(data[0, 0], data[0, 2], 'o', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'blue', markeredgewidth = 2)
ax1.plot(data[0, 0], data[0, 3], 's', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'pink', markeredgewidth = 2)
ax1.plot(data[0, 0], data[0, 4], 's', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'cornflowerblue', markeredgewidth = 2)


ax2.plot(data[:sec_1, 0], data[:sec_1, 1], 'r-', linewidth = lwidth, label = 'OMI - anode')
ax2.plot(data[:sec_1, 0], data[:sec_1, 2], 'b-', linewidth = lwidth, label = 'TFSI - anode')
ax2.plot(data[:sec_1, 0], data[:sec_1, 3], '--', color = 'pink', linewidth = lwidth, label = 'OMI - cathode')
ax2.plot(data[:sec_1, 0], data[:sec_1, 4], '--', color = 'cornflowerblue', linewidth = lwidth, label = 'TFSI - cathode')
ax2.plot(data[0, 0], data[0, 1], 'o', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'red', markeredgewidth = 2)
ax2.plot(data[0, 0], data[0, 2], 'o', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'blue', markeredgewidth = 2)
ax2.plot(data[0, 0], data[0, 3], 's', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'pink', markeredgewidth = 2)
ax2.plot(data[0, 0], data[0, 4], 's', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'cornflowerblue', markeredgewidth = 2)

ax2.legend(loc = 'best', fontsize = 24, frameon = True, numpoints = 1, ncol = 2)
ax1.set_xlim(-data[sec_2, 0]/100)
ax2.set_xlim(-data[sec_1, 0]/100)

for axes in fig.axes:
    axes.set_xlabel('Time (ps)', fontsize = fsize)
    axes.set_ylabel('Number in the slit', fontsize = fsize)
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('NIS_all.pdf')

fig, (ax2, ax1) = plt.subplots(2, figsize = (12, 12), dpi = 1000)
ax1.plot(data[:sec_2, 0], 0.0114442*55*26*2+data[:sec_2, 1]-data[:sec_2, 2], 'r-', linewidth = lwidth, label = 'charge - anode')
ax1.plot(data[:sec_2, 0], -0.0114442*55*26*2+data[:sec_2, 3]-data[:sec_2, 4], 'b--', linewidth = lwidth, label = 'charge - cathode')
ax1.plot(data[0, 0], -0.0114442*55*26*2+data[0, 1]-data[0, 2], 'ro', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'red', markeredgewidth = 2)
ax1.plot(data[0, 0], 0.0114442*55*26*2+data[0, 3]-data[0, 4], 'bs', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'blue', markeredgewidth = 2)
ax2.plot(data[:sec_1, 0], 0.0114442*55*26*2+data[:sec_1, 1]-data[:sec_1, 2], 'r-', linewidth = lwidth, label = 'charge - anode')
ax2.plot(data[:sec_1, 0], -0.0114442*55*26*2+data[:sec_1, 3]-data[:sec_1, 4], 'b--', linewidth = lwidth, label = 'charge - cathode')
ax2.plot(data[0, 0], -0.0114442*55*26*2+data[0, 1]-data[0, 2], 'ro', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'red', markeredgewidth = 2)
ax2.plot(data[0, 0], 0.0114442*55*26*2+data[0, 3]-data[0, 4], 'bs', markersize = msize, markerfacecolor = 'white', markeredgecolor = 'blue', markeredgewidth = 2)
ax1.set_xlim(-data[sec_2, 0]/100)
ax2.set_xlim(-data[sec_1, 0]/100)
for axes in fig.axes:
    axes.set_xlabel('Time (ps)', fontsize = fsize)
    axes.set_ylabel('Charge in the slit', fontsize = fsize)
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('CIS_all.pdf')
