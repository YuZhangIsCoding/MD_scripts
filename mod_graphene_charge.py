#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: mod_graphene_charge.py
# Description:  This is a python script that modifies the surface charge density
#               of graphene sheets. It contains 3 types: GPH, GPO and GNE, each
#               representing surface with no charge, positive charge and negative
#               charge.
# Date: 07-06-2017 Created

import argparse

parser = argparse.ArgumentParser(description = 'specify inputs')
parser.add_argument('-sc', dest = 'sc', default = 0, type = float, help = 'surface charge density (C/m^2)')
args = parser.parse_args()

bondlen = 0.142
delta = args.sc*bondlen**2*3*3**(0.5)/4*10/1.602

myfile = open('graphene_charge.itp', 'w')
myfile.write('[ moleculetype ]\n\
GPH   1\n\
[ atoms ]\n\
1   CG    1  GPH   CG    1   0.0000000  12.011\n\n')

myfile.write('[ moleculetype ]\n\
GPO   1\n\
[ atoms ]\n\
1   CG    1  GPO   CG    1   %12.8f  12.011\n\n' %delta)

myfile.write('[ moleculetype ]\n\
GNE   1\n\
[ atoms ]\n\
1   CG    1  GNE   CG    1   %12.8f  12.011\n\n' %(-delta))
myfile.close()
