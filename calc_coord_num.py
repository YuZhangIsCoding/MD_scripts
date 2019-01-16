#!/Users/yuzhang/anaconda3/bin/python
# Filename calc_coord_num.py
# Description:  This is a python script that calculates the coordination 
#               number distribution given the cutoff distance
# Date: 10-23-2018  Created

import mdtraj as md
import matplotlib.pyplot as plt
import pandas as pd
import argparse

import seaborn as sns
sns.set_style('darkgrid')

parser = argparse.ArgumentParser()
parser.add_argument('-chunk', '--chunksize', type = int, default=100,
                    help='chunksize for loading trajectory')
parser.add_argument('-cut', '--cutoff', type=float,
                    help='cutoff distance for coordination shell')
parser.add_argument('-debug', '--debug', action='store_true', help='flag when debug')
parser.add_argument('-i', '--input', default='traj.trr',
                    help='Input trajectory file: .trr, .xtc, etc.')
parser.add_argument('-o', '--output', default='CoordNum.txt',
                    help='Output file for coordination number distribution')
parser.add_argument('-ref', help='reference atom name')
parser.add_argument('-sel', help='selected atom name')
parser.add_argument('-top', '--topology', default='begin.gro', help='topology file')
parser.add_argument('-fig', help='figure name for plot')

try:
    __IPYTHON__
    inputs = '-i test.xtc -ref Mg -sel OW -cut 0.25 -fig CoordNum.pdf -chunk 10 -debug'
    args = parser.parse_args(inputs.split())
except NameError:
    args = parser.parse_args()

assert args.cutoff, 'Pleae specify the cutoff distance!'

traj0 = md.load_frame(args.input, 0, top=args.topology)
ref_index = traj0.topology.select('name %s' %(args.ref))
assert len(ref_index), 'The system constains no atom named %s' %(args.ref)
sel_index = traj0.topology.select('name %s' %(args.sel))
assert len(sel_index), 'The system constains no atom named %s' %(args.sel)
pairs = traj0.topology.select_pairs(ref_index, sel_index)
df = pd.DataFrame(pairs, columns=['ref', 'sel'])

for chunk_index, traj in enumerate(md.iterload(args.input, chunk = args.chunksize,
                                               top = args.topology)):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*args.chunksize+sub_ind
        df['dist'] = md.compute_distances(frame, pairs)[0]
        temp = df[df['dist'] < args.cutoff].groupby('ref')[['dist']].count()
        current_cn = temp.reset_index().groupby('dist').count()
        try:
            cns = cns.add(current_cn, fill_value = 0)
        except NameError:
            cns = current_cn
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*args.chunksize+sub_ind+1)
    if args.debug:
        break

cns /= (frame_ind+1)
if (len(ref_index)-cns['ref'].sum()) >= 0.1:
    splm = pd.DataFrame([[0, len(ref_index)-cns['ref'].sum()]], columns = ['dist', 'ref'])
    splm = splm.set_index('dist')
    cns = splm.append(cns)
cns.to_csv(args.output, header = None, sep = ' ')

if args.fig:
    fig, ax = plt.subplots(1, figsize = (4, 3), dpi = 300)
    cns.plot(kind = 'bar', color = 'lightskyblue', ax = ax, legend = False)
    plt.xlabel('Coordination Number')
    plt.ylabel('Count per frame')
    ax.tick_params(axis = 'x', rotation = 0)
    plt.tight_layout()
    plt.savefig(args.fig)
