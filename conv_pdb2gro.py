#!/Users/yuzhang/anaconda/bin/python
# Filename: conv_pdb2gro.py
# Description: This is python script to convert a pdb file to gro file
# Date: 02-15-2016 Make this script executable

import sys
########### Main ############
if len(sys.argv) == 1:
    sys.exit('Exit: Please specify at lest one file you want to modify!')
for filename in sys.argv[1:]:
    myfile = open(filename, 'r')
    outfile = open( filename[:-4]+'_conv'+'.gro', 'w')
    outfile.write('Convert from '+filename+' to '+filename[:-4]+'.pdb')
    outfile.write('\n0\n')
    for line in myfile:
        if line[:4] == 'ATOM':
            myline = '%5s%5s%5s%5s%8.3f%8.3f%8.3f\n' %(line[22:26], line[17:21], line[12:16], line[6:11], float(line[30:38])/10, float(line[38:46])/10, float(line[46:54])/10)
            outfile.write(myline)
        elif line[:3] == 'END':
            outfile.write("%10.5f%10.5f%10.5f" %(0, 0, 0))
    outfile.close()
    myfile.close()
