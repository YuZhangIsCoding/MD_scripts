#!/Users/yuzhang/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description = 'User specified filenames, boundaries, etc.')
#parser.add_argument('-c', '--choice', dest = 'choice', choices = ['surface', 'bulk'], help = 'select the region, surface or bulk')
parser.add_argument('-b', '--begin', dest = 'begin', default = 0, type = float, help = 'input the first frame')
parser.add_argument('-d', '--dimension', dest = 'dim', nargs = '*', type = float, help = 'input the dimension of the simulaition box')
parser.add_argument('-t', '--timestep', dest = 'timestep', default = 0.1, type = float, help = 'The timestep for each frame')
parser.add_argument('--suffix', type = str, help = 'The suffix for output')
args = parser.parse_args()

if args.dim == None:
    print 'No boxsize specified, gonna use the default value: (38, 3.19, 6.1)'
    dim = [38.0, 3.19, 6.1]
else:
    dim = args.dim
if args.suffix == None:
    suffix = ''
else:
    suffix = '_'+args.suffix

nid = [int(dim[0]*5-110), int(dim[1]*10)+1]

fsize = 28
lwidth = 4.0
msize = 16.0

surf = [[[] for row in range(2*nid[0])] for j in range(2)]
surf_sel = [[[] for row in range(2*20)] for j in range(2)]

bulk = [[] for row in range(int(dim[0]*10))]
st = [[] for row in range(4)]
st_err = [[] for row in range(4)]
bt = [[] for row in range(len(bulk))]
bt_err = []
myfile = open('potential_surf_kr.dat', 'r')
chunk = 0
mark = 0
for line in myfile:
    if line == '\n' or 'gfeng' in line:
        if mark == 1:
            for i in range(len(bulk)):
                bt[i].append(np.mean(bulk[i][-nid[1]:]))
        count = 0
        mark = 0
        chunk += 1
    elif chunk >= args.begin:
        temp = [float(i) for i in line.split()]
        if count < nid[0]*nid[1]*4:
            label = [count%2, count/2/nid[1]]
            surf[label[0]][label[1]].append(sum(temp))
            if count >= 30*nid[1]*2  and count < 50*nid[1]*2:
                label_sel = [count%2, count/2/nid[1]-30]
                surf_sel[label_sel[0]][label_sel[1]].append(sum(temp))
            elif count >= (30+80)*nid[1]*2  and count < (50+80)*nid[1]*2:
                label_sel = [count%2, count/2/nid[1]-30-80+20]
                surf_sel[label_sel[0]][label_sel[1]].append(sum(temp))
        else:
            if mark == 0:
                for i in range(4):
                    temp = []
                    for j in range(nid[0]):
                        temp.append(np.mean(surf[i/2][j+i%2*nid[0]][-nid[1]:])*0.010364272)
                    st[i].append(np.mean(temp))
                    st_err[i].append(np.std(temp))
                mark = 1
            bulk[count/nid[1]-nid[0]*4].append(sum(temp))
        count += 1
surf_avg = [[], []]
surf_err = [[], []]
bulk_avg = []
bulk_err = []
bt_avg = []

surf_sel_avg = [[], []]
surf_sel_err = [[], []]
#### selected region ####
for i, item in enumerate(surf_sel):
    for j in item:
        surf_sel_avg[i].append(np.mean(j)*0.010364272)
        surf_sel_err[i].append(np.std(j)*0.010364272)

fig = plt.figure(figsize = (16, 12), dpi = 1000)
plt.plot(np.arange(20)/10.0, surf_sel_avg[0][:20], linewidth = lwidth, color = 'r', label = 'bottom left')
plt.plot(np.arange(20)/10.0, surf_sel_avg[0][20:], linewidth = lwidth, color = 'b', label = 'bottom right')
plt.plot(np.arange(20)/10.0, surf_sel_avg[1][:20], linewidth = lwidth, color = 'pink', label = 'top left')
plt.plot(np.arange(20)/10.0, surf_sel_avg[1][20:], linewidth = lwidth, color = 'cornflowerblue', label = 'top right')
#(_, cap1, _) = plt.errorbar(range(2), [surf_sel_avg[0][:20][10*i] for i in range(2)], [surf_err[0][:20][10*i] for i in range(2)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'r')
#(_, cap2, _) = plt.errorbar(range(2), [surf_sel_avg[0][20:][10*i] for i in range(2)], [surf_err[0][20:][10*i] for i in range(2)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'b')
#(_, cap3, _) = plt.errorbar(range(2), [surf_sel_avg[1][:20][10*i] for i in range(2)], [surf_err[1][:20][10*i] for i in range(2)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'pink')
#(_, cap4, _) = plt.errorbar(range(2), [surf_sel_avg[1][20:][10*i] for i in range(2)], [surf_err[1][20:][10*i] for i in range(2)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'cornflowerblue')
#
#for caps in [cap1, cap2, cap3, cap4]:
#    for cap in caps:
#        cap.set_markeredgewidth(lwidth/2)

plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.legend(loc = 'best', fontsize = 24, ncol = 4)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('pt_sel_surf'+suffix+'.pdf')


for i, item in enumerate(surf):
    for j in item:
        surf_avg[i].append(np.mean(j)*0.010364272)
        surf_err[i].append(np.std(j)*0.010364272)
for item in bulk:
    bulk_avg.append(np.mean(item)*0.010364272)
    bulk_err.append(np.std(item)*0.010364272)
for item in bt:
    bt_avg.append(np.mean(item)*0.010364272)
    bt_err.append(np.std(item)*0.010364272)
### const: 0.010364272
print 'Total chunks:', chunk

fig = plt.figure(figsize = (16, 12), dpi = 1000)
plt.plot(np.arange(len(bt_avg))/10.0, bt_avg, linewidth = lwidth, color = 'k', label = 'potential in IL')
(_, caps, _) = plt.errorbar(range(len(bt_err)/10), [bt_avg[10*i] for i in range(len(bt_err)/10)], [bt_err[10*i] for i in range(len(bulk_err)/10)], capsize = 10, elinewidth= lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'k')

for cap in caps:
    cap.set_markeredgewidth(lwidth/2)
plt.xlim([0, dim[0]])
plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.legend(loc = 'best', fontsize = 24, ncol = 4)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('pt_bulk_t'+suffix+'.pdf')

fig = plt.figure(figsize = (16, 12), dpi = 1000)

plt.plot(np.arange(len(bulk_avg))/10.0, bulk_avg, linewidth = lwidth, color = 'k', label = 'potential in IL')
(_, caps, _) = plt.errorbar(range(len(bulk_err)/10), [bulk_avg[10*i] for i in range(len(bulk_err)/10)], [bulk_err[10*i] for i in range(len(bulk_err)/10)], capsize = 10, elinewidth= lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'k')

for cap in caps:
    cap.set_markeredgewidth(lwidth/2)
plt.xlim([0, dim[0]])
plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.legend(loc = 'best', fontsize = 24, ncol = 4)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('pt_bulk'+suffix+'.pdf')

fig = plt.figure(figsize = (16, 12), dpi = 1000)
plt.plot(np.arange(nid[0])/10.0, surf_avg[0][:nid[0]], linewidth = lwidth, color = 'r', label = 'bottom left')
plt.plot(np.arange(nid[0])/10.0, surf_avg[0][nid[0]:], linewidth = lwidth, color = 'b', label = 'bottom right')
plt.plot(np.arange(nid[0])/10.0, surf_avg[1][:nid[0]], linewidth = lwidth, color = 'pink', label = 'top left')
plt.plot(np.arange(nid[0])/10.0, surf_avg[1][nid[0]:], linewidth = lwidth, color = 'cornflowerblue', label = 'top right')
(_, cap1, _) = plt.errorbar(range(nid[0]/10), [surf_avg[0][:nid[0]][10*i] for i in range(nid[0]/10)], [surf_err[0][:nid[0]][10*i] for i in range(nid[0]/10)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'r')
(_, cap2, _) = plt.errorbar(range(nid[0]/10), [surf_avg[0][nid[0]:][10*i] for i in range(nid[0]/10)], [surf_err[0][nid[0]:][10*i] for i in range(nid[0]/10)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'b')
(_, cap3, _) = plt.errorbar(range(nid[0]/10), [surf_avg[1][:nid[0]][10*i] for i in range(nid[0]/10)], [surf_err[1][:nid[0]][10*i] for i in range(nid[0]/10)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'pink')
(_, cap4, _) = plt.errorbar(range(nid[0]/10), [surf_avg[1][nid[0]:][10*i] for i in range(nid[0]/10)], [surf_err[1][nid[0]:][10*i] for i in range(nid[0]/10)], capsize = 10, elinewidth = lwidth/2, fmt = 'o', markeredgecolor = 'none', linewidth = lwidth, color = 'cornflowerblue')

for caps in [cap1, cap2, cap3, cap4]:
    for cap in caps:
        cap.set_markeredgewidth(lwidth/2)

plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.legend(loc = 'upper center', fontsize = 24, ncol = 4)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('pt_surf'+suffix+'.pdf')

fig = plt.figure(figsize = (16, 12), dpi = 1000)

plt.plot(np.arange(len(st[0]))*args.timestep, st[0], linewidth = lwidth, color = 'r', label = 'bottom left')
plt.plot(np.arange(len(st[0]))*args.timestep, st[1], linewidth = lwidth, color = 'b', label = 'bottom right')
plt.plot(np.arange(len(st[0]))*args.timestep, st[2], linewidth = lwidth, color = 'pink', label = 'top left')
plt.plot(np.arange(len(st[0]))*args.timestep, st[3], linewidth = lwidth, color = 'cornflowerblue', label = 'top right')

#plt.errorbar(np.arange(len(st[0]))*args.timestep, st[0], st_err[0], linewidth = lwidth, color = 'r', label = 'bottom left')
#plt.errorbar(np.arange(len(st[0]))*args.timestep, st[1], st_err[1], linewidth = lwidth, color = 'b', label = 'bottom right')
#plt.errorbar(np.arange(len(st[0]))*args.timestep, st[2], st_err[2], linewidth = lwidth, color = 'pink', label = 'top left')
#plt.errorbar(np.arange(len(st[0]))*args.timestep, st[3], st_err[3], linewidth = lwidth, color = 'cornflowerblue', label = 'top right')

plt.xlabel('Time (ns)', fontsize = fsize)
plt.ylabel('Potential (V)', fontsize = fsize)
plt.legend(loc = 'center', fontsize = 24, ncol = 4)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('pt_surf_t'+suffix+'.pdf')

print 'Surface potential:', np.mean(surf_avg[0][:nid[0]]), np.mean(surf_avg[1][:nid[0]]), np.mean(surf_avg[0][nid[0]:]), np.mean(surf_avg[1][nid[0]:])
print 'Surface potential in selected region:', np.mean(surf_sel_avg[0][:20]), np.mean(surf_sel_avg[1][:20]), np.mean(surf_sel_avg[0][20:]), np.mean(surf_sel_avg[1][20:])
