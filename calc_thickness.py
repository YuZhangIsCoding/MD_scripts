import numpy as np
import pdb
import matplotlib.pyplot as plt

numd = np.loadtxt('NumberDensity.txt')
numd[:,0] -= 0.682+0.07
bound = max(numd[:,0])-0.682-0.07
dz = numd[1][0]-numd[0][0]
mark = 0
sum_z = 0
sum_n = 0
e0 = 1/6.242

choice = raw_input('PZC?(Y/N)\n')
if choice == 'Y' or choice == 'y':
    cri = 4
    for ion in [1, 2]:
        pz1 = 0
        pz2 = 0
        pn1 = 0
        pn2 = 0
        for i in range(len(numd)):
            if numd[i, ion] > cri:
                pz1 += numd[i,0]*numd[i,ion]
                pn1 += numd[i,ion]
                mark = 1
            elif numd[i, ion] < cri and mark == 1:
                mark = 0
                break
        for i in range(len(numd)):
            if numd[len(numd)-i-1, ion] > cri:
                pz2 += numd[len(numd)-i-1,0]*numd[len(numd)-i-1, ion]
                pn2 += numd[len(numd)-i-1,ion]
                mark = 1
            elif numd[len(numd)-i-1, ion] < cri and mark == 1:
                mark = 0
                break
        thickness1 = pz1/pn1
        thickness2 = bound-pz2/pn2
        print thickness1, thickness2,(thickness1+thickness2)/2
else:
    surfc = input('Please input the surface charge density (C/m^2):\n')
    eia = np.zeros(len(numd))
    temp = 0
    loc = [[],[]]
    for i, item in enumerate(numd):
        temp += (item[2]-item[1])*dz*e0/surfc
        if temp < -1 and mark == 0:
            loc[0].append(i)
            mark = -1
        elif temp > -1 and mark == -1:
            loc[0].append(i)
            mark = 0
        elif temp > 1 and mark == 0:
            loc[1].append(i)
            mark = 1
        elif temp < 1 and mark == 1:
            loc[1].append(i)
            mark = 0
        eia[i] = temp
    sum_z = 0
    sum_n = 0
    for i in range(loc[0][0], loc[0][1]):
        sum_z += numd[i,0]*numd[i,2]
        sum_n += numd[i,2]
    thickness1 = sum_z/sum_n
    sum_z = 0
    sum_n = 0
    for i in range(loc[1][-2], loc[1][-1]):
        sum_z += numd[i,0]*numd[i,1]
        sum_n += numd[i,1]
    thickness2 = bound-sum_z/sum_n
    print 'Thickness (nm):'
    print thickness1, thickness2
    ## figure output ##
    fsize = 36
    lwidth = 4.0
    fig = plt.figure(figsize=(16,12), dpi=1000)
    plt.plot(numd[:,0], eia, color = 'black', linewidth = lwidth)
    plt.plot([0, bound], [0,0], '--',color = 'black')
    plt.plot([0, bound], [-1,-1], '--', color = 'blue')
    plt.plot([0, bound], [1,1], '--', color = 'red')
    for i in loc[0]:
        plt.plot([numd[i,0], numd[i,0]], [0, -1], '--', color = 'blue')
    for i in loc[1]:
        plt.plot([numd[i,0], numd[i,0]], [0, 1], '--',color = 'red')
    
    plt.plot([thickness1, thickness1], [0, -10],'--', color = 'black', linewidth = lwidth)
    plt.plot([bound-thickness2, bound-thickness2], [0, 10],'--', color = 'black', linewidth = lwidth)
    plt.xlim([0, bound])
    plt.xlabel('Distance', fontsize = fsize) #, labelpad = 20)
    plt.ylabel('EIA', fontsize = fsize) #, labelpad = 20)
    plt.gca().tick_params(labelsize = fsize, width = 2, length = 6, pad = 20)
    for direction in ['top', 'bottom', 'left', 'right']:
        plt.gca().spines[direction].set_linewidth(2)
    plt.tight_layout()
    plt.savefig('eia.pdf')
