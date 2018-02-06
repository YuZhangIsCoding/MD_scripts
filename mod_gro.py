#!/Users/yuzhang/anaconda/bin/python

import argparse, sys

parser = argparse.ArgumentParser(description = 'User specified boundary')
parser.add_argument('-i', '--input', dest = 'filename', default = 'begin.gro', help = 'filename')
parser.add_argument('-s', '--size', type = float, nargs = '*', help = 'amount reduced')
parser.add_argument('-o', '--output', dest = 'outname', default = 'shrink.gro', help = 'output name')
args = parser.parse_args()

myfile = open(args.filename, 'r')
outfile = open(args.outname, 'w')

print 'Input:', args.filename, 'Output:', args.outname,

mylines = myfile.readlines()
for line in mylines[:-1]:
    outfile.write(line)
bound = [float(i) for i in mylines[-1].split()]
if len(args.size) == 0:
    sys.exit('Exit: please specify the length you want to modify')
elif len(args.size) == 1:
    print 'Modify X direction:', args.size
    size = args.size+[0, 0]
elif len(args.size) == 2:
    print 'Modify X, Y direction:', args.size
    size = args.size+[0]
elif len(args.size) == 3:
    print 'Modify all 3 directions:', args.size
    size = args.size
else:
    sys.exit('Exit: too many input lengths')
outfile.write('%10.5f%10.5f%10.5f\n' %(bound[0]-size[0], bound[1]-size[1], bound[2]-size[2]))
print 'Previous box:', bound
print 'Current box:', [bound[i]-size[i] for i in range(3)]

myfile.close()
outfile.close()

