#!/Users/yuzhang/anaconda/bin/python

import argparse, os, sys

parser = argparse.ArgumentParser(description = 'User specified charge density, size, etc.')
parser.add_argument('-sc', '--surfcharge', dest = 'sc', type = float, help = 'surface charge')
#parser.add_argument('-d', '--dim', dest = 'dim', nargs = '*', type = float, help = 'dimension of the system')
parser.add_argument('-i', '--input', dest = 'filename', default = 'neutral.itp', help = 'input file')
parser.add_argument('-c', '--config', dest = 'config', choices = ['Q', 'H2Q', 'Q-0.25', 'Q-0.5', 'Q-0.75', 'graphene'], help = 'Quinone type: Q, H2Q, Q-0.5, graphene')
args = parser.parse_args()

if os.path.isfile(args.filename):
    print 'Reading itp file: %s ...' %args.filename
    infile = open(args.filename)
else:
    sys.exit('Exit: no file named %s found in current directory!' %args.filename)

if args.filename == 'graphene.itp':
    print 'graphene.itp duplicates with the output name, please rename the input file'
    sys.exit()

print 'Appending surface charge density of', args.sc,'C/cm2'

if args.config == 'Q':
    delta = args.sc*5.11302*4.92*10/1.602/12/30
elif args.config == 'H2Q':
    delta = args.sc*5.11302*4.92*10/1.602/14/30
elif args.config == 'Q-0.5':
    delta = args.sc*5.90400*5.11302*10/1.602/26/18
elif args.config == 'Q-0.25':
    delta = args.sc*5.90400*5.11302*10/1.602/54/9
elif args.config == 'Q-0.75':
    delta = args.sc*5.90400*5.11302*10/1.602/50/9
elif args.config == 'graphene':
    delta = args.sc*0.142**2*3*3**0.5/4*10/1.602

outfile = open('graphene.itp', 'w')
for line in infile:
    if 'GPO' == line[:3]:
        flag = 1
    elif 'GNE' == line[:3]:
        flag = -1
    elif 'GPH' == line[:3]:
        flag = 0
    temp = line.split()
    if args.config == 'graphene':
        if len(temp) == 8 and temp[1] == 'CG':
            outfile.write('%s%12.8f%10.3f\n' %(line[:31], flag*delta, float(temp[-1])))
        else:
            outfile.write(line)
    else:
        if len(temp) == 8 and temp[1] != 'CG':
            outfile.write('%s%12.8f%10.3f\n' %(line[:31], float(temp[-2])+flag*delta, float(temp[-1])))
        else:
            outfile.write(line)
outfile.close()
infile.close()
print 'Output file: graphene.itp'
