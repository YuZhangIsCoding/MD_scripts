#!/Users/yuzhang/anaconda3/bin/python
# Description:  This is a python script that generates itp files for free 
#               standing graphene. In current code, the graphene is cut to
#               have zigzag edge.

from itp_common import Atom, Compound
import itertools, argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--size', nargs = 2, type = int, 
                    help = 'size of the graphene sheet (carbon atoms in each lateral direction)')
parser.add_argument('-o', '--output', default = 'graphene.itp', help = 'name of the output file')
try:
    __IPYTHON__
    args = parser.parse_args(['-s', '4','4'])
except NameError:
    args = parser.parse_args()
class graphene(Compound):
    def __init__(self, size, name = 'GPH'):
        Atom.clear_index()
        self.size = size
        self.frame = {}
        super(graphene, self).__init__(name)
    def add_atoms(self):
        for i, j in itertools.product(*[range(_) for _ in self.size]):
            atom = Atom('CG')
            atom.assign_cgroup()
            self.atoms[atom.index] = atom
            self.frame[(i, j)] = atom
            # vertical bonds
            if j == self.size[1]-1:
                self.add_bond(atom, self.frame[(i, 0)])
            if j != 0:
                self.add_bond(atom, self.frame[(i, j-1)])
            # horizontal bonds
            if i != 0 and ((i%2 == j%2 == 0) or (i%2 == j%2 == 1)):
                    self.add_bond(atom, self.frame[(i-1, j)])
    def _add_improper_dihedrals(self):
        for atom in self.atoms.values():
            if len(atom.neighbors) == 3:
                self.improp.append(tuple([atom.index]+atom.neighbors))
    def gen_itp(self):
        self.add_atoms()
        self.add_angles()
        self.add_dihedrals(improper = True)
    
temp = graphene(size = args.size)
temp.write_itp(para = {'CG': (0.0153333, 12.011)}, filename = args.output)
