#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_potential.py
# Description: This is a python script to calculate the electrical potential distribution along the channel for armchair edge
# Date: 09-23-2016 Created

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = 'User specified density profile, surface charge density, etc.')
parser.add_argument('-i', '--input', dest = 'filename', type = str, default = 'Density.txt', help = 'File that contains charge density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('--exb', dest = 'exb', nargs = '*', help = 'The location of the two charge plane')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', help = 'Boundary of the distribution')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
parser.add_argument('--sp', dest = 'sp', type = int, default = 20, help = 'Data points to be split')
args = parser.parse_args()

#if parser.filename != None:
#    filename = args.filename
#else:
#    filename = 'density.xvg'

if args.suffix != None:
    outname = 'pd_'+args.suffix
else:
    outname = 'pd'

if args.exb == None:
#    exb_1 = 1.352
#    exb_2 = 1.494
    exb_1 = 1.078
    exb_2 = 1.171
else:
    exb_1 = float(args.exb[0])
    exb_2 = float(args.exb[1])
filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)


sigma1 = -0.1268*24*15/5.112/5.115/6.2415
sigma2 = -sigma1
if args.sigma != None:
    sigma = args.sigma
else:
    sigma = float(raw_input('Please input the surface charge density (C/m^2):'))

data = []
if filename.endswith('.txt'):
    data = np.loadtxt(filename)
    if args.bound != None:
        bound = [float(i) for i in args.bound]
    else:
        bound = [0, max(data[:, 0])]
#    for i in temp:
#        if i[0] >= bound[0]+exb and i[0] <= bound[1]-exb:
#            data.append(i)
#    data = np.array(data)
    rho = data[:,2]
elif filename.endswith('.xvg'):
    if args.bound != None:
        bound = args.bound
    else:
        bound = raw_input('Please input the boundary (nm):').split()
    bound = [float(i) for i in bound]
    file_in = open(filename, 'r')
    for line in file_in:
        if line[0] != '#' and line[0] != '@':
            data.append([float(i) for i in line.split()])
    file_in.close()
    data = np.array(data)
    rho = data[:,1]
z = data[:,0]
for seq, myrho in enumerate([rho, rho[::-1]]):
    phi_1 = np.zeros(len(z))
    phi_2 = np.zeros(len(z))
    const = 1000/8.85418784762/6.242
    for i, u in enumerate(z):
        if z[i] < exb_1:
            phi_1[i] = 0
        elif z[i] > bound[1]-exb_1:
            phi_1[i] = phi_1[i-1]
        else:
            phi_1[i] = -np.trapz((u-z[:i+1])*myrho[:i+1], z[:i+1])*const-sigma1*(z[i]-exb_1)*const/0.160217657
    for i, u in enumerate(z):
        if z[i] < exb_2:
            phi_2[i] = 0
        elif z[i] > bound[1]-exb_2:
            phi_2[i] = phi_2[i-1]
        else:
            phi_2[i] = -(sigma2+(-1)**seq*sigma)*(z[i]-exb_2)*const/0.160217657
    
    phi = phi_1+phi_2
    bulk = np.mean(phi[len(z)/2-50:len(z)/2+50])
    #print 'potential (V), left, bulk, right'
    #print "surface charge density (C/m^2):", sigma
    if sigma == 0:
        print bulk,
    else:
        print bulk,
    hl = len(phi)/2
    if seq == 0:
        plt.plot(z[:hl], phi[:hl], color = 'pink')
        plt.plot([0, z.max()],[bulk, bulk],'--', color = 'red')
    else:
        plt.plot(z[-hl:][::-1], phi[:hl], color = 'cornflowerblue')
        plt.plot([0, z.max()],[bulk, bulk],'--', color = 'blue')
#plt.plot(z, phi_1, color = 'red')
#plt.plot(z, phi_2, color = 'blue')
plt.xlim([0, z.max()])
plt.xlabel('Distance (nm)')
plt.ylabel('Potential (V)')
plt.savefig(outname+'.pdf')

np.savetxt(outname+'.txt', np.column_stack((z, phi)), fmt = ['%12.4f', '%12.4f'])
