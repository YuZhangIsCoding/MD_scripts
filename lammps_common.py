#!/Users/yuzhang/anaconda3/bin/python
# Filename:     lammps_common.py
# Description:  This is a python script that has a class for lammps data file
#               that stores info for atom, bond, angle, etc.
# Date: 10-02-2018  Created
#       10-03-2018  Recursively read in data files using _read_currentline
#                   Add check_consistency

class LAMMPS(object):
    def __init__(self):
        numnames = ['atoms', 'bonds', 'angles', 'dihedrals',
                'atom types', 'bond types', 'angle types', 'dihedral types']
        self.nums = {item: 0 for item in numnames}
        fftypenames = ['Masses', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs']
        self.fftypes = {item: [] for item in fftypenames}
        self.atoms= []
        conn_names = ['Bonds', 'Angles', 'Dihedrals']
        self.connections = {item: [] for item in conn_names}
        self.box = []
        self.lj = []
    def read_data(self, filenames):
        for filename in filenames:
            with open(filename) as myfile:
                filegen = (line for line in myfile if not (line.strip() == '' or
                                line.strip()[0] is '#'))
                for line in filegen:
                    self._read_currentline(line, filegen)
        self.check_consistency()
                
    def _read_currentline(self, currentline, filegen):
        temp = currentline.split()
        if temp[-1] in self.nums:
            self.nums[temp[-1]] = int(temp[0])
        elif temp[-1] == 'types':
            self.nums[' '.join(temp[-2:])] = int(temp[0])
        elif temp[-1] in ['Masses', 'Coeffs']:
            typename = ' '.join(temp)
            self._read_types(typename, filegen)
        elif temp[-1] == 'Atoms':
            self._read_atoms(filegen)
        elif temp[-1] in self.connections:
            self._read_connections(temp[-1], filegen)
        elif temp[0] == 'pair_coeff':
            self.lj.append([int(_) for _ in temp[1:3]])
            self.lj[-1].extend([float(_) for _ in temp[3:]])

    def _read_types(self, typename, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.fftypes[typename].append([int(temp[0])])
                self.fftypes[typename][-1].extend([float(_) for _ in temp[1:]])
            except:
                self._read_currentline(line, filegen)
                return

    def _read_atoms(self, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.atoms.append([int(_) for _ in temp[:3]])
                self.atoms[-1].extend([float(_) for _ in temp[3:]])
            except:
                self._read_currentline(line, filegen)
                return

    def _read_connections(self, typename, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.connections[typename].append([int(_) for _ in temp])
            except:
                self._read_currentline(line, filegen)
                return

    def check_consistency(self):
        if self.nums['atoms'] != len(self.atoms):
            raise ValueError('%d atoms do not match the number of coordinates: %d' 
                    %(self.nums['atom'], len(self.atoms)))
        if self.nums['bonds'] != len(self.connections['Bonds']):
            raise ValueError('%d bonds do not match the number of bonds in connections: %d'
                    %(self.nums['bonds'], len(self.connections['Bonds'])))
        if self.nums['angles'] != len(self.connections['Angles']):
            raise ValueError('%d angles do not match the number of angles in connections: %d'
                    %(self.nums['angles'], len(self.connections['Angles'])))
        if self.nums['dihedrals'] != len(self.connections['Dihedrals']):
            raise ValueError('%d dihedrals do not match the number of dihedrals in connections: %d'
                    %(self.nums['dihedrals'], len(self.connections['Dihedrals'])))
        if self.nums['atom types'] != len(self.fftypes['Masses']):
            raise ValueError('%d atom types do not match the number of Masses: %d'
                    %(self.nums['atom types'], len(self.fftypes['Masses'])))
        if self.nums['bond types'] != len(self.fftypes['Bond Coeffs']):
            raise ValueError('%d bond types do not match the number of Bond Coeffs: %d'
                    %(self.nums['bond types'], len(self.fftypes['Bond Coeffs'])))
        if self.nums['angle types'] != len(self.fftypes['Angle Coeffs']):
            raise ValueError('%d angle types do not match the number of Angle Coeffs: %d'
                    %(self.nums['angle types'], len(self.fftypes['Angle Coeffs'])))
        if self.nums['dihedral types'] != len(self.fftypes['Dihedral Coeffs']):
            raise ValueError('%d dihedral types do not match the number of Dihedral Coeffs: %d'
                    %(self.nums['dihedral types'], len(self.fftypes['Dihedral Coeffs'])))

def test():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs = '*', help = 'input file')
    args = parser.parse_args()
    if args.input:
        lmp = LAMMPS()
        lmp.read_data(args.input)
        print(lmp.nums)
        print(lmp.fftypes)
        print('Atoms:', len(lmp.atoms))
        for item in lmp.connections:
            print(item, len(lmp.connections[item]))
        print('LJ:',lmp.lj)
if __name__ == '__main__':
    test()
