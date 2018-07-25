#!/Users/yuzhang/anaconda3/bin/python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'Input gro file')
parser.add_argument('-o', '--output', default = 'MXENE.itp', help = 'Output itp file')
args = parser.parse_args()

myfile = open(args.input, 'r')
outfile = open(args.output, 'w')

outfile.write('[ moleculetype ]\nMXN\t 3\n[ atoms ]\n')
count = 0
for i, line in enumerate(myfile):
    if i > 1:
        temp = line.split()
        if len(temp) == 3:
            break
        name = temp[1]
        z = float(temp[5])
        if 'C' in name:
            amass = 12.011
            atype = 'CM'
            charge = -0.74
        count += 1
        if 'Ti' in name:
            atype = 'TiM'
            amass = 47.867
            if z > 0.2 and z < 0.4:
                charge = 0.68
            else:
                charge = 1.04
        elif 'O' in name:
            atype = 'OM'
            amass = 15.99
            charge = -0.64
            if z > 0.5:
                pass
                #charge += 1/20
        outfile.write('%5d%5s%5d%5s%5s%5d%12.7f%12.7f\n' %(count, atype, 1, 'MXN', atype, count, charge, amass))
outfile.close()
myfile.close()
