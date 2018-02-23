#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_potential.py
# Description:  This is a python script to calculate the electrical potential 
#               distribution along the radial direction of OLC supercapacitor
#               The integral of poisson equation for spherical electrode is:
#               phi(r) = -1/epsilon*(sum{R->r}((1-u/r)u*pho(u)+sigma*R(1-R/r)))
#               The density profile is loaded as dataframe using pandas 
# Date: 02-15-2016  Created based on calc_potential.py

import sys, os, argparse, warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description = 'User specified density profile, surface charge density, etc.')
parser.add_argument('-i', '--input', dest = 'filename', type = str, default = 'Density.txt', help = 'File that contains charge density profile')
parser.add_argument('-s', '--sigma', dest = 'sigma', type = float, help = 'Surface charge density (C/m^2)')
parser.add_argument('-e', dest = 'e', type = float, help = 'surface net charge (e)')
parser.add_argument('--exb', dest = 'exb', type = float, default = 0.682,help = 'The radius of the OLC')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', help = 'Boundary of the distribution')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
parser.add_argument('-v', action = 'store_true', help = 'Print out details')
args = parser.parse_args()


if args.suffix != None:
    outname = 'pd_'+args.suffix
else:
    outname = 'pd'

exb = args.exb
if args.v:
    print('The radius of the OLC is %f nm' %exb)
filename = args.filename
if not os.path.isfile(filename):
    sys.exit('Exit: no density profile named %s found in current directory!' %filename)

if args.sigma != None:
    sigma = args.sigma
elif args.e != None:
    sigma = args.e/6.242/4/np.pi/exb**2 
else:
    sigma = float(raw_input('Please input the surface charge density (C/m^2):'))


data = []
if filename.endswith('.txt'):
    df = pd.read_csv(filename, sep = '\s+', names = ['z', 'charge'], usecols = [0, 2])
    df = df[df.z >= exb]
#    df.z = df.z-exb

phi = np.zeros(len(df))
const = 1000/8.85418784762/6.242
for i, r in enumerate(df.z):
    phi[i] = -np.trapz((1-df.z.iloc[:i+1]/r)*df.z.iloc[:i+1]*df.charge.iloc[:i+1], df.z.iloc[:i+1])*const-sigma*exb*(1-exb/r)*const/0.160217657

bulk = np.mean(phi[int(0.8*len(df.z)):int(len(df.z)*0.9)])
loc = 0
ref = phi[loc]
print(sigma, ref-bulk)

plt.plot(df.z, phi, color = 'black')
plt.plot([0, df.z.max()], [bulk, bulk],'--', color = 'red')
plt.axvline(x=exb, ls = ':', color = 'gray')
plt.xlim([0, df.z.max()])
plt.xlabel('Distance (nm)')
plt.ylabel('Potential (V)')
plt.savefig(outname+'.pdf')

np.savetxt(outname+'.txt', np.column_stack((df.z, phi)), fmt = ['%12.4f', '%12.4f'])
