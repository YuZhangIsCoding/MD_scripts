#!/Users/yuzhang/anaconda3/bin/python
# Filename:     itp_common.itp
# Description:  This is a python module that includes commonly used classes 
#               and functions for generating itp files for GROMACS

import itertools

class Atom(object):
    index = 0
    cgroup = 0
    def __init__(self, atomname):
        '''Initialize an atom with a unique index and empty neighbor list and no charge group'''
        self.name = atomname
        Atom.index += 1
        self.index = Atom.index
        self.neighbors = []
        self.cgroup = None
    def assign_cgroup(self, n = None):
        if n is None:
            Atom.cgroup += 1
            self.cgroup = Atom.cgroup
        else:
            self.cgroup = n
    def add_neighbor(self, other):
        '''Add neighbor to an atom'''
        self.neighbors.append(other.index)
        other.neighbors.append(self.index)
    def add_neighbors(self, iterator):
        '''Add neighbors from an iterator object'''
        for item in iterator:
            self.add_neighbor(item)
    def get_neighbors(self):
        for item in self.neighbors:
            yield item
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name
    def __lt__(self, other):
        return self.index < other.index
    @staticmethod
    def clear_index():
        '''Clean the atom index'''
        Atom.index = 0
        Atom.cgroup = 0

class Compound(object):
    def __init__(self, name = 'Compound'):
        self.name = name
        self.atoms = {}
        self.bonds = set()
        self.angles = []
        self.dihedrals = []
        self.improp = []
    def get_number(self):
        return len(self.atoms)
    def add_bond(self, atom1, atom2):
        '''Manually build your molecule with neighbors and bonds'''
        atom1.add_neighbor(atom2)
        self.bonds.add((min(atom1.index, atom2.index), max(atom1.index, atom2.index)))
    def add_bonds(self):
        '''Build bonds given neighbor lists of atoms'''
        for atom in self.atoms.values():
            for item in atom.get_neighbors:
                self.bonds.add((min(atom.index, item), max(atom.index, item)))
    def add_angles(self):
        '''Add angles'''
        for atom in self.atoms.values():
            if len(atom.neighbors) > 1:
                for index1, index2 in itertools.combinations(atom.neighbors, 2):
                    self.angles.append((min(index1, index2), atom.index, max(index1, index2)))
    def add_dihedrals(self, improper = False):
        for bond in self.bonds:
            for i, j in itertools.product(self.atoms[bond[0]].neighbors, self.atoms[bond[1]].neighbors):                                                             
                if i != bond[1] and j != bond[0] and i != j: # avoid triangle loop like cyclopropane
                    self.dihedrals.append((i, bond[0], bond[1], j))
        if improper:
            self._add_improper_dihedrals()
    def _add_improper_dihedrals(self):
        raise NotImplementedError('Subclass should implement this!')
    def gen_itp(self):
        raise NotImplementedError('Subclass should implement this!')
    def write_itp(self, para, filename = 'out.itp'):
        '''
        para: dictionary maps atom name with its corresponding partial charge and mass
        '''
        self.gen_itp()
        compound_name = self.name[:3].upper()
        myfile = open(filename, 'w')
        myfile.write('[ moleculetype ]\n')
        myfile.write('%s\t3\n' %compound_name)
        
        myfile.write('[ atoms ]\n')
        for i in range(1, self.get_number()+1):
            atomname = self.atoms[i].name
            myfile.write('%5d%5s%5d%5s%5s%5d%12.8f%12.8f\n' %(i, atomname, 1, compound_name, atomname,
                                        self.atoms[i].cgroup, para[atomname][0], para[atomname][1]))
        
        myfile.write('[ bonds ]\n')
        for bond in self.bonds:
            myfile.write('%5d%5d%5d\n' %(bond+(1,)))
        
        myfile.write('[ angles ]\n')
        for angle in self.angles:
            myfile.write('%5d%5d%5d%5d\n' %(angle+(1,)))
            
        myfile.write('[ dihedrals ]\n')
        for dih in self.dihedrals:
            ## RD dihedrals
            myfile.write('%5d%5d%5d%5d%5d\n' %(dih+(1,))) 
        
        if self.improp:
            myfile.write('; improper dihedrals\n')
            for improp in self.improp:
                myfile.write('%5d%5d%5d%5d%5d\n' %(improp+(4,))) 
        myfile.close()

