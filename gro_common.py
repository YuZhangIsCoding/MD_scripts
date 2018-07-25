#!/Users/yuzhang/anaconda3/bin/python
# Filename:     gro_common.py
# Description:  This is a python module that includes some common functions
#               to deal with .gro files, including reading coordinates, etc.
# Date: 06-14-2018

class Gro(object):
    def __init__(self):
        self.num = None
        self.info = []
        self.box = None
    def read(self, grofile):
        myfile = open(grofile, 'r')
        for i, line in enumerate(myfile):
            if i == 1:
                self.num = int(line)
            elif i > 1 and i-2 < self.num:
                name = line[:20]
                coords = [float(_) for _ in line[20:44].split()]
                self.info.append((name, coords))
        self.box = [float(_) for _ in line.split()]
        while len(self.box) < 9:
            self.box.append(0)
        myfile.close()

