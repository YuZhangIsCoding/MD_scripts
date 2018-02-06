#!/Users/yuzhang/anaconda/envs/py3/bin/python
import matplotlib.pyplot as plt
import numpy as np
import calc_common as comm

if 'msd' not in comm.args.filename:
    filename = 'msd.txt'
else:
    filename = comm.args.filename

if comm.args.scale == None:
    scale = 2
else:
    scale = comm.args.scale

data = np.loadtxt(filename)[1:]

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

labelname = ['EMIm', 'TFSI', 'Sol']
color = ['red', 'blue', 'green']
subcolor = ['pink', 'cornflowerblue', 'lightgreen']
#for i in range(2):
#    plt.plot(data[:, 0]/1000, data[:, i+1], color = color[i], markeredgecolor = 'none', linewidth = lwidth, label = labelname[i])

t = data[:, 0]
k = [[] for row in range(len(data[0])-1)]

lstyle = ['--', '-.']
for i in range(1, len(data[0])):
    plt.plot(data[:, 0]/1000, data[:, i], color = color[i-1], markeredgecolor = 'none', linewidth = lwidth, label = labelname[i-1])
    k0 = np.log(data[1, i]/data[0, i])/np.log(t[1]/t[0])
    b0 = data[0, i]/t[0]**k0
#    plt.plot(t/1000, b0*t**k0, '--', color = subcolor[i-1], linewidth = lwidth/2, label = 'order of %.3f' %k0)
    k1 = np.log(data[-1, i]/data[-2, i])/np.log(t[-1]/t[-2])
    b1 = data[-1, i]/t[-1]**k1
#    plt.plot(t/1000, b1*t**k1, '-.', color = subcolor[i-1], linewidth = lwidth/2, label = 'order of %.3f' %k1)
    for j in range(len(data)-1):
        k[i-1].append(np.log(data[j+1, i]/data[j, i])/np.log(t[j+1]/t[j]))
    b = data[-1, i]/data[-2, 0]
    plt.plot(t[-len(t)/2:]/1000, b*t[-len(t)/2:], 'k--', lw = lwidth/2)

pt = [data[-1, 0]/1000, max(data[-1, 1:])/scale]
b = pt[1]/pt[0]
plt.plot(t[-len(t)/2:]/1000, b*t[-len(t)/2:]/1000, 'k--', lw = lwidth)


plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ns)', fontsize = fsize)
plt.ylabel('MSD (nm^2)', fontsize = fsize)
for axes in fig.axes:
    axes.set_xscale('log', nonposx = 'clip')
    axes.set_yscale('log', nonposx = 'clip')
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
#plt.xlim(t[0]/1000, int(t[-1]/10000*1.5)*10)
plt.grid(which = 'both')
plt.tight_layout()
plt.savefig('diff_log.pdf')

fig = plt.figure(figsize = (16,12), dpi = 300)

for i in range(2):
    plt.plot(t[:-1]/1000, k[i], color = color[i], linewidth = lwidth/2, label = labelname[i])

plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('time (ns)', fontsize = fsize)
plt.ylabel('k', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
#plt.xlim(10)
plt.grid(which = 'both')
plt.tight_layout()
plt.savefig('diff_alpha.pdf')
