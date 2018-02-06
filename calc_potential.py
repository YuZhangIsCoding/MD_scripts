#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_potential.py
# Description: This is a python script to calculate the electrical potential distribution along the channel
# Date: 02-15-2016 Make this file an executable #       02-17-2016 Add argparser module to get parameters from command line input

import sys, os, argparse, warnings
import numpy as np
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description = 'User specified density profile, surface charge density, etc.')
parser.add_argument('-i', '--input', dest = 'filename', type = str, default = 'Density.txt', help = 'File that contains charge density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('--exb', dest = 'exb', type = float, default = 0.682,help = 'The length of the electrode')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', help = 'Boundary of the distribution')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
args = parser.parse_args()


if args.suffix != None:
    outname = 'pd_'+args.suffix
else:
    outname = 'pd'

exb = args.exb
filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)

if args.sigma != None:
    sigma = args.sigma
else:
    sigma = float(raw_input('Please input the surface charge density (C/m^2):'))

data = []
if filename.endswith('.txt'):
    temp = np.loadtxt(filename)
#    ref = np.loadtxt('Density_ACN.txt')
#    temp[:, 2] -= ref[:, 2]
    bound = [0, max(temp[:, 0])]
    for i in temp:
        if i[0] >= bound[0]+exb and i[0] <= bound[1]-exb:
            data.append(i)
    data = np.array(data)
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
            temp = [float(i) for i in line.split()]
            if temp[0] < bound[0]+exb:
                continue
            elif temp[0] > bound[1]-exb:
                break
            else:
                data.append([float(i) for i in line.split()])
    file_in.close()
    data = np.array(data)
    rho = data[:,1]
rho[0] = 0
z = data[:,0]-exb
phi = np.zeros(len(z))
const = 1000/8.85418784762/6.242
for i, u in enumerate(z):
    phi[i] = -np.trapz((u-z[:i+1])*rho[:i+1], z[:i+1])*const-sigma*z[i]*const/0.160217657

bulk = np.mean(phi[len(z)/2-50:len(z)/2+50])
loc = 0
ref = phi[loc]
if sigma == 0:
    print sigma, ref-bulk, phi[-1-loc]-bulk
else:
    print sigma, ref-bulk, phi[-1-loc]-bulk

### Test quinone
#plt.plot([0.384, 0.384], [min(phi), max(phi)], '--', color = 'grey')
#plt.plot([5-0.384, 5-0.384], [min(phi), max(phi)], '--', color = 'grey')

plt.plot(z,phi, color = 'black')
plt.plot([0, z.max()],[bulk, bulk],'--', color = 'red')
plt.xlim([0, z.max()])
plt.xlabel('Distance (nm)')
plt.ylabel('Potential (V)')
plt.savefig(outname+'.pdf')

np.savetxt(outname+'.txt', np.column_stack((z, phi)), fmt = ['%12.4f', '%12.4f'])
