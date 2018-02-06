#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_eia.py
# This is a script to calculate the EIA (Effective Ion Accumulation) profile
# Date created : 08-04-2015

import numpy as np
import argparse
import matplotlib.pyplot as plt

fig, (ax1, ax3) = plt.subplots(2, sharex = True, figsize = (12, 16), dpi = 1000)
parser = argparse.ArgumentParser(description = 'Specify the charge density.')
parser.add_argument('-i', '--input', dest = 'filename', default = 'NumberDensity.txt', help = 'filename for number density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', type = float, help = 'Surface charge density')
parser.add_argument('--pzc', dest = 'pzc', action = 'store_true', help = 'charged surface or not')
parser.add_argument('--funct', dest = 'funct', action = 'store_true', help = 'functionalized or not')
parser.add_argument('--sec', dest = 'sec', action = 'store_true', help = 'whether choose second layer')
args = parser.parse_args()

if args.sigma == None:
    sigma = float(raw_input('Please input surface charge density (C/m^2):\n'))
else:
    sigma = args.sigma

numd = np.loadtxt(args.filename)
if args.pzc == False:
    if args.funct == True:
        ref_1 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_OH/prodRun/charged_1/TEABF4_OH_test_9/NumberDensity.txt')
        ref_2 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_OH/prodRun/charged_2/TEABF4_OH_test_9_2/NumberDensity.txt')
        ref_3 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_OH/prodRun/charged_3/TEABF4_OH_test_9_3/NumberDensity.txt')
    else:
        ref_1 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_pristine/prodrun/charged_pristine_1/TEABF4_pristine_pzc_3/NumberDensity.txt')
        ref_2 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_pristine/prodrun/charged_pristine_2/TEABF4_pristine_pzc_3_2/NumberDensity.txt')
        ref_3 = np.loadtxt('/Volumes/Macintosh_HD_2/mystore/TEABF4/TEABF4_pristine/prodrun/charged_pristine_3/TEABF4_pristine_pzc_3_2/NumberDensity.txt')
    ref_1[:, 1:] = (ref_1[:, 1:]+ref_2[::-1][:, 1:])/2
    ref_2[:, 1:] = (ref_2[:, 1:]+ref_2[::-1][:, 1:])/2
    ref_3[:, 1:] = (ref_3[:, 1:]+ref_3[::-1][:, 1:])/2
    ref = (ref_1+ref_2+ref_3)/3
#    numd[:, 1:] -= ref[:, 1:]
else:
    numd[:, 1:] = (numd[::-1][:, 1:]+numd[:, 1:])/2

dz = numd[1][0]-numd[0][0]
mysum = np.zeros(len(numd))
zlist = []
nlist = [0]
nacc = [[], [],[]]
ratio = [0]
if args.pzc == True:
    for i in range(len(numd)-1):
        mysum[i] = np.trapz(numd[:i+1, 2]-numd[:i+1, 1], numd[:i+1, 0])
        if i > 0:
            temp1 = mysum[i]
            temp2 = -mysum[i-1]
            if temp1*temp2 > 0:
                ratio.append(temp2/(temp1+temp2))
                nlist.append(i)
                zlist.append(numd[i-1][0]+temp2/(temp1+temp2)*dz)
else:
    for i in range(len(numd)-1):
        mysum[i] = np.trapz(numd[:i+1, 2]-numd[:i+1, 1], numd[:i+1, 0])
        mysum[i] /= sigma*10/1.602176565
        if i > 0:
            temp1 = mysum[i]-1
            temp2 = 1-mysum[i-1]
            if temp1*temp2 > 0:
                ratio.append(temp2/(temp1+temp2))
                nlist.append(i)
                zlist.append(numd[i-1][0]+temp2/(temp1+temp2)*dz)
nlist.append(-1)

np.savetxt('EIA.txt',np.column_stack((numd[:, 0],mysum)),fmt=['%12.4f','%12.4f'])
ax1.plot(numd[:, 0], mysum, 'r-', label = 'EIA', linewidth = 2)
ax3.plot(numd[:, 0], numd[:, 1], color = 'purple', label = 'cation')
ax3.plot(numd[:, 0], numd[:, 2], color = 'cyan', label = 'anion')
for i in zlist:
    ax3.plot([i,i],[0,10], '--', color = 'red')
ax1.plot([0,numd[-1][0]],[1,1], 'k--')
ax1.set_ylim([np.min(mysum)*1.1, np.max(mysum)*1.1])
ax3.set_xlim([0, np.max(numd[:, 0])])
for ax in fig.axes:
    ax.legend()
plt.savefig('EIA.pdf')

zeff = [[], [], []]
for j in range(1, 4):
    area_1 = 0
    t_1 = 0
    for i in range(1, len(nlist)):
        if i == len(nlist)-1:
            area_2 = 0
            t_2 = 0
            nacc[j-1].append(np.trapz(numd[nlist[i-1]:, j], numd[nlist[i-1]:, 0])+area_1)
            if nacc[j-1][-1] == 0:
                zeff[j-1].append(0)
            else:
                zeff[j-1].append((np.trapz(numd[nlist[i-1]:, j]*numd[nlist[i-1]:, 0], numd[nlist[i-1]:, 0])+t_1)/nacc[j-1][-1])
        else:
            area_2 = (numd[nlist[i]-1, j]+0.5*ratio[i]*(numd[nlist[i], j]-numd[nlist[i]-1, j]))*ratio[i]*dz
            #nacc[j-1].append(np.trapz(numd[nlist[i-1]:nlist[i], j], numd[nlist[i-1]:nlist[i], 0])+area_1+area_2)
            nacc[j-1].append(np.trapz(numd[0:nlist[i], j], numd[0: nlist[i], 0])+area_2)
            area_1 = (numd[nlist[i]-1, j]+numd[nlist[i], j])/2*dz-area_2
            if nacc[j-1][-1] == 0:
                zeff[j-1].append(0)
            else:
                mdpoint = (numd[nlist[i]-1, 0]+ratio[i]*dz)*(numd[nlist[i]-1, j]+ratio[i]*(numd[nlist[i], j]-numd[nlist[i]-1, j]))
                t_2 = (numd[nlist[i]-1, j]*numd[nlist[i]-1, 0]+mdpoint)/2*ratio[i]*dz
                #zeff[j-1].append((np.trapz(numd[nlist[i-1]:nlist[i], j]*numd[nlist[i-1]:nlist[i], 0], numd[nlist[i-1]:nlist[i], 0])+t_1+t_2)/nacc[j-1][-1])
                zeff[j-1].append((np.trapz(numd[0: nlist[i], j]*numd[0: nlist[i], 0], numd[0: nlist[i], 0])+t_2)/nacc[j-1][-1])
                t_1 = (numd[nlist[i], j]*numd[nlist[i], 0]+mdpoint)/2*(1-ratio[i])*dz

zeff1 = [[], [], []]
nacc1 = [[], [], []]
numd[:, 1:] = numd[::-1][:, 1:]
#print zlist
nlist = [len(numd)-i for i in nlist]
nlist.reverse()
for j in range(1, 4):
    area_1 = 0
    t_1 = 0
    for i in range(1, len(nlist)):
        if i == len(nlist)-1:
            area_2 = 0
            t_2 = 0
            nacc1[j-1].append(np.trapz(numd[nlist[i-1]:, j], numd[nlist[i-1]:, 0])+area_1)
            if nacc1[j-1][-1] == 0:
                zeff1[j-1].append(0)
            else:
                zeff1[j-1].append((np.trapz(numd[nlist[i-1]:, j]*numd[nlist[i-1]:, 0], numd[nlist[i-1]:, 0])+t_1)/nacc1[j-1][-1])
        else:
            area_2 = (numd[nlist[i]-1, j]+0.5*ratio[i]*(numd[nlist[i], j]-numd[nlist[i]-1, j]))*ratio[i]*dz
            #nacc[j-1].append(np.trapz(numd[nlist[i-1]:nlist[i], j], numd[nlist[i-1]:nlist[i], 0])+area_1+area_2)
            nacc1[j-1].append(np.trapz(numd[0:nlist[i], j], numd[0: nlist[i], 0])+area_2)
            area_1 = (numd[nlist[i]-1, j]+numd[nlist[i], j])/2*dz-area_2
            if nacc1[j-1][-1] == 0:
                zeff1[j-1].append(0)
            else:
                mdpoint = (numd[nlist[i]-1, 0]+ratio[i]*dz)*(numd[nlist[i]-1, j]+ratio[i]*(numd[nlist[i], j]-numd[nlist[i]-1, j]))
                t_2 = (numd[nlist[i]-1, j]*numd[nlist[i]-1, 0]+mdpoint)/2*ratio[i]*dz
                #zeff[j-1].append((np.trapz(numd[nlist[i-1]:nlist[i], j]*numd[nlist[i-1]:nlist[i], 0], numd[nlist[i-1]:nlist[i], 0])+t_1+t_2)/nacc[j-1][-1])
                zeff1[j-1].append((np.trapz(numd[0: nlist[i], j]*numd[0: nlist[i], 0], numd[0: nlist[i], 0])+t_2)/nacc1[j-1][-1])
                t_1 = (numd[nlist[i], j]*numd[nlist[i], 0]+mdpoint)/2*(1-ratio[i])*dz
def print_nd(nacc):
    if args.sec == True:
        #print nacc[0][0], nacc[1][0], nacc[2][0], sum(nacc[0][-2:]), sum(nacc[1][-2:]), sum(nacc[2][-2:])
        print sum(nacc[0][:2]), sum(nacc[1][:2]), sum(nacc[2][:2]), sum(nacc[0][-3:]), sum(nacc[1][-3:]), sum(nacc[2][-3:])
    else:
        #print nacc[0][0], nacc[1][0], nacc[2][0], nacc[0][-1], nacc[1][-1], nacc[2][-1]
        #print sum(nacc[1][:2])-sum(nacc[0][:2]), sum(nacc[0][-2:])-sum(nacc[1][-2:])
        print sum(nacc[0][:2]), sum(nacc[1][:2]), sum(nacc[2][:2]),  sum(nacc[0][-2:]), sum(nacc[1][-2:]), sum(nacc[2][-2:])
def print_z(zeff):
    if args.sec == True:
        #results = [zeff[0][0]-0.682, zeff[1][0]-0.682, zeff[2][0]-0.682, zeff1[0][1]-0.682, zeff1[1][1]-0.682, zeff1[2][1]-0.682]
        results = [zeff[0][1]-0.682, zeff[1][1]-0.682, zeff[2][1]-0.682, zeff1[0][2]-0.682, zeff1[1][2]-0.682, zeff1[2][2]-0.682]
    else:
        #results = [zeff[0][0]-0.682, zeff[1][0]-0.682, zeff[2][0]-0.682, zeff1[0][0]-0.682, zeff1[1][0]-0.682, zeff1[2][0]-0.682]
        results = [zeff[0][1]-0.682, zeff[1][1]-0.682, zeff[2][1]-0.682, zeff1[0][1]-0.682, zeff1[1][1]-0.682, zeff1[2][1]-0.682]
    for i in results:
        if i <= 0:
            print "%12.5f" %0,
        elif i > 2:
            print "%12.5f" %(6.364-0.682-i),
        else:
            print "%12.5f" %i,

print_z(zeff)
#print zeff
#print '\n', zlist
