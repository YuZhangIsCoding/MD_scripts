#!/Users/yuzhang/anaconda3/bin/python
# Filename:     gro_common.py
# Description:  This is a python module that includes some common functions
#               to deal with .gro files, including reading coordinates, etc.
# Date: 06-14-2018  Created
#       09-26-2018  Add get_atomtypes

class Gro(object):
    def __init__(self):
        self.box = None
        self.natoms = None
        self.nmols = None
        self.info = []
        self.atomtypes = {}
        self.itp = []
        self.mass = {}
        self.lj = {}
    def load_gro(self, filename):
        '''read *.gro file, save atom info and box size'''
        with open(filename) as myfile:
            for i, line in enumerate(myfile):
                if i == 1:
                    self.natoms = int(line)
                elif i > 1 and i < self.natoms+2:
                    self.info.append(Gro.read_line(line))
            self.box = [float(_) for _ in line.split()]
            # box vector for any coordinate system
            while len(self.box) < 9:
                self.box.append(0)
        self.check_itp_gro()
        self.get_atomtypes()
    def load_itp(self, filename):
        '''Read corresponding *.itp file and save charge and molar mass
        into self.itp and self.mass.
        Note that one atom type could have different charges but could
        only have one molar mass'''
        with open(filename) as myfile:
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'atoms' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    temp = line.split()
                    if temp[1] not in self.mass:
                        self.mass[temp[1]] = float(temp[-1])
                    self.itp.append([float(_) for _ in temp[-2:]])
        self.check_itp_gro()
    def load_nonbond(self, filename):
        '''Read gromacs topology file and find the Lennard-Jones parameter
        and write it with atomtypes'''
        with open(filename) as myfile:
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'atomtypes' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    temp = line.split()
                    if temp[0] in self.atomtypes:
                        self.lj[temp[0]] = [float(_) for _ in temp[-2:]]
    def get_atomtypes(self):
        # Note if we want to deal with multiple input gro
        # file, probably need to define a unique atomtype
        # number under Gro, something like self.typenum
        cnt = 1
        for item in self.info:
            temp = item[2].strip()
            if temp not in self.atomtypes:
                self.atomtypes[temp] = cnt
                cnt += 1
    def check_itp_gro(self):
        '''Check if the number of atoms in gro file is a multiple of
        the number of atoms in itp file'''
        if self.natoms and self.itp and self.natoms%len(self.itp):
            raise Exception('gro file does not match itp file!')
    @staticmethod
    def read_line(string):
        results = []
        #gro file format '%5d%5s%5s%5d%8.3f%8.3f%8.3f'
        for i in range(4):
            results.append(string[i*5:(i+1)*5])
        results.extend([float(_) for _ in string[20:44].split()])
        return results
