#!/Users/yuzhang/anaconda3/bin/python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'Input gro file')
parser.add_argument('-o', '--output', default = 'MXENE.itp', help = 'Output itp file')
parser.add_argument('-sc', '--surfacecharge', dest = 'sc', 
                    type = float, help = 'Surface charge density')
parser.add_argument('-c', '--charge', dest = 'charge', 
                    default = 0, type = float, help = 'Atomic charge on the surface atoms')
parser.add_argument('-n', '--name', dest = 'name', default = 'MXN', help = 'Residue name')
args = parser.parse_args()

myfile = open(args.input, 'r')
outfile = open(args.output, 'w')

outfile.write('[ moleculetype ]\n%s\t 3\n[ atoms ]\n' %args.name)
count = 0

#c_dict = dict(C = -0.74,
#        Ti2 = 0.68,
#        Ti1 = 1.04,
#        O = -0.64
#        )
c_dict = dict(
        O = -1.137605,
        Ti1 = 1.929465,
        C = -1.6701,
        Ti2 = 1.70338
        )

if args.sc:
    delta = args.sc*0.30390*0.26319*10/1.602
else:
    delta = args.charge

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
            charge = c_dict['C']
        count += 1
        if 'Ti' in name:
            atype = 'TiM'
            amass = 47.867
            if z > 0.2 and z < 0.4:
                charge = c_dict['Ti2']
            else:
                charge = c_dict['Ti1']
        elif 'O' in name:
            atype = 'OM'
            amass = 15.99
            charge = c_dict['O']
            charge += delta
        outfile.write('%5d%5s%5d%5s%5s%5d%12.7f%12.5f\n' %(count, atype, 1, args.name, atype, count, charge, amass))
outfile.close()
myfile.close()
