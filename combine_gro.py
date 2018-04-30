#!/Users/yuzhang/anaconda3/bin/python
# Filename: combine_gro.py
# This is a python script to combine .gro files and assort residue
# Date: 05-20-2015

import sys

############ Main ###########
if len(sys.argv) == 1:
    sys.exit('Exit: please specify at lest one file you want to modify!')
    
res_count = [1]
res_name = []
atom_count = 0
mydict = {}
dict_count = 0
temp = []
mylist = []
for filename in sys.argv[1:]:
    '''Read file and assort each line according to their residue name and index'''
    myfile = open(filename, 'r')
    count = 0
    for line in myfile:
        count += 1
        if count > 2:
            entry = line.split()
            if len(entry) != 3:
                res = line[5:10]
                if res not in res_name:
                    res_name.append(res)
                    mydict[res] = len(res_name)-1
                    temp.append(None)
                    res_count.append(1)
                    mylist.append([])
                if temp[mydict[res]] == None:
                    temp[mydict[res]] = line[:10]
                if line[:10] != temp[mydict[res]]:
                    res_count[mydict[res]] += 1
                    temp[mydict[res]] = line[:10]
                mylist[mydict[res]].append([res_count[mydict[res]],line[5:]])
                atom_count += 1
    myfile.close()
########### Output ##########
outfile = open('myout.gro','w')
outfile.write('This is a combined file from %s\n' %(', '.join(sys.argv[1:])))
outfile.write(str(atom_count)+'\n')
temp = 0
for item in mylist:
    for subitem in item:
        outfile.write('%5d' %(temp+subitem[0])+subitem[1])
    temp += subitem[0]
outfile.write("%10.5f%10.5f%10.5f" %(0, 0, 0))
outfile.close()
