#!/Users/yuzhang/anaconda/bin/python
import matplotlib.pyplot as plt
import numpy as np

hdist = np.loadtxt('hbond_ha_dist.txt')

vdw_hf = [0.1475, 0.0000]
es_hf = [0.416, -0.16]
fe = 138.9354859 #kJ*nm/mol/e**2
kt = 2.479 # kJ/mol #1.38064852*10**-23*298*6.02214*10**23
bins = np.arange(0.15, 0.31, 0.005)
counts, _ = np.histogram(hdist, bins = bins, normed = True)

fig = plt.figure(figsize = (4,4), dpi = 300)
fsize = 12
lwidth = 1.0
msize = 6.0

vdw = vdw_hf
es = es_hf

f = lambda r: 4*vdw[1]*((vdw[0]/r)**12-(vdw[0]/r)**6)
g = lambda r: fe*es[0]*es[1]/r
#plt.plot(bins, f(bins)/kt, 'r--', lw = lwidth, label = 'VDW')
plt.plot(bins, g(bins)/kt, 'b-.', lw = lwidth, label = 'Electrostatic')
plt.plot([0.2, 0.3], [0, 0], ':', color = 'gray')
#plt.plot(bins, (f(bins)+g(bins))/kt, 'k--', lw = lwidth, label = 'Total')
plt.xlim([0.2,0.3])
plt.xticks(np.arange(0.15, 0.31, 0.05))
plt.scatter(bins, (f(bins)+g(bins))/kt, s = counts*msize, label = 'Total', alpha = 0.5, color = 'k', edgecolors = 'k')
plt.scatter(np.mean(hdist), np.mean(f(hdist)+g(hdist))/kt, s = msize*20, color = 'violet', edgecolors = 'violet', label = 'Average', alpha = 0.7)
plt.ylim([-25, 0])
pos = (np.mean(hdist), np.mean(f(hdist)+g(hdist))/kt)
plt.annotate('%.3f nm, %.2f kT' %pos, xy = pos, xytext = (0.25, -22), arrowprops = dict(arrowstyle="->", edgecolor = 'black', connectionstyle="arc3", shrinkB = 8), ha = 'center', va = 'bottom')

plt.xlabel('H-acceptor distance (nm)', fontsize = fsize)
plt.ylabel('Energy (kT)', fontsize = fsize)
plt.legend(fontsize = fsize, numpoints = 1, scatterpoints = 1)

## setting the tick formats
for axes in fig.axes:
    axes.tick_params(pad = 6)
#    axes.tick_params(which = 'minor', length = 1, width = 0.25)
#    for direction in ['top', 'bottom', 'left', 'right']:
#        axes.spines[direction].set_linewidth(1)
#plt.grid()
plt.tight_layout(pad = 0.4, w_pad = 0.5, h_pad = 1.0)
plt.savefig('Hbonds_energy_hydroxyl_hf.pdf')
print 'average interaction energy is (kT):', np.mean(f(hdist)+g(hdist))/kt
