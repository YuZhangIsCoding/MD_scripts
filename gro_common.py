#!/Users/yuzhang/anaconda3/bin/python
# Filename:     gro_common.py
# Description:  This is a python module that includes some common functions
#               to deal with .gro files, including reading coordinates, etc.
# Date: 06-14-2018  Created
#       09-26-2018  Add get_atomtypes

class Gro(object):
    def __init__(self):
        '''
        info:   list of parsed line in .gro file, based on gro file format, 
                including: residue number, residue name, atom name, 
                atom number position(x, y, z)
        atomtypes:  dictionary of atomtype and its corresponding index
        mq:     list of atom mass and charge
        mass:   dictionary of atomtypes and their corresponding mass
        lj:     dictionary of atomtypes and their corresponding lj parameters
                the order is [sigma, epsilon]
        '''
        self.box = None
        self.natoms = None
        self.nmols = None
        self.info = []
        self.atomtypes = {}
        self.mq = []
        self.mass = {}
        self.lj = {}
        self.bondtypes = {}
        self.angletypes = {}
        self.dihedraltypes = {}
        self.bonds = []
        self.angles = []
        self.dihedrals = []
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
    def load_ff(self, filenames, nonbond = None):
        '''Read in itp file, and load mass and charge infomation such as
        molar mass, atomic charge, bonds, angles, dihedrals, etc.'''
        fftypes = {}
#        sessions = ['atoms', 'bondtypes', 'angletypes', 'nonbond']
        for filename in filenames:
            with open(filename) as myfile:
                for line in myfile:
                    if '[' in line and ']' in line:
                        if 'atoms' in line:
                            fftypes['mq'] = filename
                        elif 'bondtypes' in line:
                            fftypes['bondtypes'] = filename
                        elif 'angletypes' in line:
                            fftypes['angletypes'] = filename
                        elif 'dihedraltypes' in line:
                            fftypes['dihedral'] = filename
                        elif 'atomtypes' in line:
                            fftypes['nonbond'] = filename
                        elif 'bonds' in line:
                            fftypes['bonds'] = filename
                        elif 'angles' in line:
                            fftypes['angles'] = filename
        self.load_mq(fftypes['mq'])
        self.load_bondtypes(fftypes['bondtypes'])
        self.load_bonds(fftypes['bonds'])
        self.load_angletypes(fftypes['angletypes'])
        self.load_angles(fftypes['angles'])
        self.load_nonbonds(fftypes['nonbond'])
        self.check_itp_gro()
        self.nmols = self.natoms//len(self.mq)
    def load_mq(self, filename):
        '''Read corresponding *.itp file and save charge and molar mass
        into self.mq and self.mass.
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
                    self.mq.append([float(_) for _ in temp[-2:]])
    def load_bondtypes(self, filename, kcal = True):
        with open(filename) as myfile:
            cnt = 1
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'bondtypes' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    temp = line.split()
                    k = float(temp[-1])
                    b = float(temp[-2]) #angstrom
                    if not kcal:
                        k /= 100*4.184*2 # kJ to kCal/A^2
                        b *= 10 # nm to angstrom
                    self.bondtypes[(temp[0], temp[1])] = [cnt, k, b]
                    self.bondtypes[(temp[1], temp[0])] = self.bondtypes[(temp[0], temp[1])]
                    cnt += 1

    def load_bonds(self, filename):
        with open(filename) as myfile:
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'bonds' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    self.bonds.append([int(_) for _ in line.split()[:2]])
                    
    def load_angletypes(self, filename, kcal = True):
        with open(filename) as myfile:
            cnt = 1
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'angletypes' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    temp = line.split()
                    theta = float(temp[-2]) 
                    k = float(temp[-1])
                    if not kcal:
                        k /= 4.184*2 # kJ/rad^2 to kCal/rad^2
                    self.angletypes[(temp[0], temp[1], temp[2])] = [cnt, k, theta]
                    self.angletypes[(temp[2], temp[0], temp[0])] = self.angletypes[(temp[0], temp[1], temp[2])]
                    cnt += 1
    def load_angles(self, filename):
        with open(filename) as myfile:
            flag = False
            for line in myfile:
                if line.strip() == '' or line.strip()[0] in [';', '#']:
                    continue
                if '[' in line and ']' in line:
                    if 'angles' in line:
                        flag = True
                        continue
                    elif flag:
                        break
                if flag:
                    self.angles.append([int(_) for _ in line.split()[:3]])
    def load_dihedrals(self):
        pass
    def load_nonbonds(self, filename, kcal = True):
        '''Read gromacs topology file and find the Lennard-Jones parameter
        and write it with atomtypes'''
        # if the input nonbonded file is an itp file, need to change the units
        if filename.endswith('.itp') or filename.endswith('.top'):
            kcal = False
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
                        sigma = float(temp[-2])
                        epsilon = float(temp[-1])
                        if not kcal:
                            sigma *= 10
                            epsilon /= 4.184
                        self.lj[temp[0]] = [epsilon, sigma]
    def get_atomtypes(self):
        # Note if we want to deal with multiple input gro
        # file, probably need to define a unique atomtype
        # number under Gro, something like self.typenum
        cnt = 1
        for item in self.info:
            if item[2] not in self.atomtypes:
                self.atomtypes[item[2]] = cnt
                cnt += 1
    def check_itp_gro(self):
        '''Check if the number of atoms in gro file is a multiple of
        the number of atoms in itp file'''
        if self.natoms and self.mq and self.natoms%len(self.mq):
            raise Exception('gro file does not match itp file!')
    @staticmethod
    def read_line(string):
        '''gro file format '%5d%5s%5s%5d%8.3f%8.3f%8.3f'
        '''
        results = []
        for i in range(4):
            results.append(string[i*5:(i+1)*5].strip())
        results.extend([float(_) for _ in string[20:44].split()])
        return results
