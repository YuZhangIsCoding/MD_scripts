#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_potential.py
# Description: This is a python script to calculate the electrical potential distribution along the channel
# Date: 02-15-2016 Make this file an executable
#       02-17-2016 Add argparser module to get parameters from command line input

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
parser.add_argument('--fitting', dest = 'fitting', type = float, nargs = 2, default = [0.4, 0.6], help = 'Define the fitting range')
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
    exb_1 = 1.352
    exb_2 = 1.494
#    exb_1 = 1.244
#    exb_2 = 1.386
else:
    exb_1 = float(args.exb[0])
    exb_2 = float(args.exb[1])
filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)

if args.sigma != None:
    sigma = args.sigma/2
else:
    sigma = float(raw_input('Please input the surface charge density (C/m^2):'))/2

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
rho_1 = rho.copy()
rho_2 = rho.copy()
for i in range(len(rho)):
    if i >= args.sp and i <= len(rho)-args.sp:
        rho_2[i] = 0
    else:
        rho_1[i] = 0
phi_1 = np.zeros(len(z))
phi_2 = np.zeros(len(z))
const = 1000/8.85418784762/6.242
for i, u in enumerate(z):
    if z[i] < exb_1:
        phi_1[i] = 0
    elif z[i] > bound[1]-exb_1:
        phi_1[i] = phi_1[i-1]
    else:
        phi_1[i] = -np.trapz((u-z[:i+1])*rho_1[:i+1], z[:i+1])*const-sigma*(z[i]-exb_1)*const/0.160217657
for i, u in enumerate(z):
    if z[i] < exb_2:
        phi_2[i] = 0
    elif z[i] > bound[1]-exb_2:
        phi_2[i] = phi_2[i-1]
    else:
        phi_2[i] = -np.trapz((u-z[:i+1])*rho_2[:i+1], z[:i+1])*const-sigma*(z[i]-exb_2)*const/0.160217657

phi = phi_1+phi_2
bulk = np.mean(phi[int(len(z)*args.fitting[0]):int(len(z)*args.fitting[1])+1])
#print 'potential (V), left, bulk, right'
#print "surface charge density (C/m^2):", sigma
if sigma == 0:
    print sigma, -bulk, np.mean(phi[-5:])-bulk
else:
    print sigma*2, phi[0]-bulk, phi[-1]-bulk

plt.plot(z,phi, color = 'black')
plt.plot([0, z.max()],[bulk, bulk],'--', color = 'red')
plt.xlim([0, z.max()])
plt.axvspan(z[int(len(z)*args.fitting[0])], z[int(len(z)*args.fitting[1])], facecolor = 'pink', alpha = 0.5, edgecolor = 'pink')
plt.xlabel('Distance (nm)')
plt.ylabel('Potential (V)')
plt.savefig(outname+'.pdf')

np.savetxt(outname+'.txt', np.column_stack((z, phi)), fmt = ['%12.4f', '%12.4f'])
