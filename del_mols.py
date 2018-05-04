#!/Users/yuzhang/anaconda3/bin/python
# Filename: del_mols.py
# Description:  This is a python script that deletes some molecules from .gro files.
#               A molecule is recognized by its residue name and residue index.
#               Both forward deletion and backward deletion are suported.
#               This is useful when tuning the number of molecules in the the system.
# Dates:        05-03-2018  Created

import argparse, sys, queue, os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'Input file')
parser.add_argument('-o', '--output', default = 'out.gro', help = 'Output filename')
parser.add_argument('-m', '--mols', nargs = '*', help = 'molecules to be deleted')
parser.add_argument('-n', '--nums', nargs = '*', type = int, help = 'number of molecules to be deleted')
parser.add_argument('-b', '--backward', action = 'store_true', help = 'delete from back if specified')
parser.add_argument('-v', action = 'store_true', help = 'Show details')

try:
    __IPYTHON__
    args = parser.parse_args('-i begin.gro -m TPA BF4 -n 2 2 -b'.split())
except NameError:
    args = parser.parse_args()
    
if len(args.mols) != len(args.nums):
    sys.exit('molecule name and number does not match!')

outname = args.output
while os.path.isfile(outname):
    outname = outname+'.dup...'

if args.v:
    print('Reading input file %s...\nWriting to file %s' %(args.input, outname))
    
    
mol_dict = {}
q_dict = {}
for i, j in zip(args.mols, args.nums):
    mol_dict[i] = j
    if args.backward:
        q_dict[i] =  queue.Queue()

def forward_delete(filename, outname):
    '''Read from line to line and delele the first n molecules encountered,
    where n is the number of current molecules you want to delete.
    '''
    myfile = open(filename, 'r')
    outfile = open(outname, 'w')
    rec_dict = {}
    pre = None
    for i, line in enumerate(myfile):
        if i > 1:
            temp = line[5: 10].replace(' ', '')
            if temp in mol_dict and rec_dict.get(temp, 0) <= mol_dict[temp]:
                resid = int(line[:5])
                if resid != pre:
                    pre = resid
                    rec_dict[temp] =  rec_dict.get(temp, 0)+1
                if rec_dict[temp] <= mol_dict[temp]:
                    continue
        outfile.write(line)
    outfile.close()
    myfile.close()

# Use python queue module for the FIFO data structure
def backward_delete(filename, outname):
    '''Delete the molecules from backward.
    Read from line to line, if a molecule is found first put it in
    a queue. When the number of molecules in the queue reaches n, which
    is the number of current molecules you want to delete, we begin to 
    push molecules from the queue, so that the total molecules remained
    in the queue equals n.
    '''
    myfile = open(filename, 'r')
    outfile = open(outname, 'w')
    rec_dict = {}
    pre = None
    for i, line in enumerate(myfile):
        if i > 1:
            temp = line[5: 10].replace(' ', '')
            if temp in mol_dict:
                q_dict[temp].put(line)
                resid = int(line[:5])
                if resid != pre:
                    pre = resid
                    rec_dict[temp] =  rec_dict.get(temp, 0)+1
                if rec_dict[temp] > mol_dict[temp]:
                    outfile.write(q_dict[temp].get())
                    continue
        outfile.write(line)
    outfile.close()
    myfile.close()

if args.backward:
    backward_delete(args.input, outname)
else:
    forward_delete(args.input, outname)
