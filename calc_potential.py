#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_potential.py
# Description:  This is a python script to calculate the electrical potential 
#               distribution along the channel. The potential distribution is
#               obtained by integrating Poisson equation in one direction.
#               This script works for surfaces that partial charge is evenly
#               distributed on the surface. Initially only for planar surface.
#               Potential distribution along the radial direction of OLC 
# 		supercapacitor can also be calculated.The integral of poisson 
# 		equation for spherical electrode is:
#               phi(r) = -1/epsilon*(sum{R->r}((1-u/r)u*pho(u)+sigma*R(1-R/r)))
# Date: 02-15-2016  Make this file an executable.
#       02-17-2016  Add argparser module to get parameters from command line.
#       02-23-2018  Add spherical surface.

import sys, os, argparse, warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description = 'User specified density profile, surface charge density, etc.')
parser.add_argument('-i', '--input', dest = 'filename', type = str, default = 'Density.txt', \
        help = 'File that contains charge density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', default = 0, type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('-e', dest = 'e', type = float, default = 0, help = 'surface net charge (e)')
parser.add_argument('--exb', dest = 'exb', type = float, default = 0.682,help = 'The length of the electrode')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', type = float, help = 'Boundary of the electrodes')
parser.add_argument('-r', '--region', nargs = '*', default = (0.45, 0.55), type = float, help = 'Relative positive of bulk region')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
parser.add_argument('-g', '--geometry', choices = ('plane', 'sphere'), default = 'plane', help = 'surface geometry')
parser.add_argument('-v', action = 'store_true', help = 'show details')
args = parser.parse_args()

######### functions ############
def load_txt(filename, args):
    '''Read txt file and choose the data according to electrode positions'''
    data = np.loadtxt(filename)[:, [0, 2]]
    data = data[(data[:, 0] >= args.exb) & (data[:, 0] <= data[:, 0].max()-args.exb)]
    return data[:, 0], data[:, 1]

def load_xvg(filename, args):
    '''Read xvg file, append the data within boundaries'''
    if args.bound:
        bound = args.bound
    else:
        bound = [float(item) for item in input('Please input the boundaries (nm):').split()]
    file_in = open(filename, 'r')
    signs = set(('#', '@'))
    for line in file_in:
        if line[0] not in signs:
            temp = [float(i) for i in line.split()]
            if temp[0] > bound[1]-args.exb:
                break
            elif temp[0] >= bound[0]+args.exb:
                data.append(temp)
    file_in.close()
    data = np.array(data)
    return data[:, 0], data[:, 1]


########### initialize ################
filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: not density profile named %s found in current directory!' %filename)

if args.suffix:
    outname = 'pd_'+args.suffix
else:
    outname = 'pd'

if args.geometry == 'plane':
    sigma = args.sigma
else:
    sigma = args.e/6.242/4/np.pi/args.exb**2 

if args.v:
    print('Reading from %s ...' %filename)
    print('Surface charge density is:', sigma, 'C/m^2')
    print('Thickness of electrode:', args.exb, 'nm')
    print('Calculating the bulk potential from %.2f to %.2f' %(args.region[0], args.region[1]))
    print('Writing to file %s.txt and %s.pdf ...' %(outname, outname))
    print('Surface charge density, anode potential, cathode potential:')

if filename.endswith('.txt'):
    z, rho = load_txt(filename, args)
elif filename.endswith('.txt'):
    z, rho = load_xvg(filename, args)
else:
    raise NotImplementedError()

const = 1000/8.85418784762/6.242
phi = np.zeros(len(z))
if args.geometry == 'plane':
    z -= args.exb
    for i, u in enumerate(z):
        phi[i] = -np.trapz((u-z[:i+1])*rho[:i+1], z[:i+1])*const-sigma*u*const/0.160217657
else:
    for i, r in enumerate(z):
        phi[i] = -np.trapz((1-z[:i+1]/r)*z[:i+1]*rho[:i+1], z[:i+1])*const-sigma*args.exb*(1-args.exb/r)*const/0.160217657

bulk = phi[int(args.region[0]*len(phi)): int(args.region[1]*len(phi))].mean()
loc = 0
ref = phi[loc]
if args.geometry == 'plane':
    print(sigma, phi[loc]-bulk, phi[-1-loc]-bulk)
else:
    print(sigma, ref-bulk)
    plt.axvline(x = args.exb, ls = ':', color = 'gray')

np.savetxt(outname+'.txt', np.column_stack((z, phi)), fmt = ['%12.4f']*2)
plt.plot(z, phi, color = 'black')
plt.plot((0, z.max()), (bulk, bulk), '--', color = 'red')
plt.xlim((0, z.max()))
plt.xlabel('Distance (nm)')
plt.ylabel('Potential (V)')
plt.savefig(outname+'.pdf')
