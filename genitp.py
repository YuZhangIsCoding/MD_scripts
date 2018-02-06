print("[ moleculetype ]\n; molname\tnrexcl\nGPH\t\t\t1\n")
print("[ atoms ]")
print("; id type    res residu  at  cg  charge   mass")
# dimensions
l = 16
h = 28
for i in range(1,l*h+1):
    print "%5d%5s%5d%5s%5s%5d%8.4f%8.3f" %(i,'CGP',1,'GPH','CG',i,0,12.011)

#bond
print("[    bonds   ]")
for i in range(1,l+1):
    for j in range(1,h+1):
        num=h*(i-1)+j
        if j == 1:
            below = num+h-1
        else:
            below = num-1
            print "%8d%8d%5d" %(num, below,1)
        if j == h:
            above = num-h+1
        else:
            above = num+1
            print "%8d%8d%5d" %(num, above,1)
        if (i+num)%2 == 0:
            if i == 1:
                left = j+(l-1)*h
            else:
                left = num-h
                print "%8d%8d%5d" %(num, left,1)
#            print "%8d%8d%5d" %(num, left,1)
#            print "%8d%8d%5d" %(num, above,1)
#            print "%8d%8d%5d" %(num, below,1)
#file_in = open('mbuild.top')
#po=0
#for line in file_in:
#    entry = line.split()
#    if len(entry) >1 and entry[0]=="[":
#        if entry[1] == "angles":
#            po=1
#        elif entry[1] == "dihedrals":
#            po=2
#        else:
#            po=0
#    elif len(entry)==0:
#        po=0
#    if po==1:
#        if entry[0]=="[":
#            print entry[0],entry[1],entry[2]
#        else:
#            print("\t"+entry[0]+"\t"+entry[1]+"\t"+entry[2]+"\t1")
#    elif po==2:
#        if len(entry)>3:
#            print("\t"+entry[0]+"\t"+entry[1]+"\t"+entry[2]+"\t"+entry[3]+"\t1")
#        else:
#            print entry[0],entry[1],entry[2]
#print(";improper dihedrals")
#for i in range(1,l+1):
#    for j in range(1,h+1):
#        num=h*(i-1)+j
#        if j == 1:
#            below = num+h-1
#        else:
#            below = num-1
#        if j == h:
#            above = num-h+1
#        else:
#            above = num+1
#        if (i+num)%2 == 0:
#            if i==1:
#                left = j+(l-1)*h
#            else:
#                left = num-h
#            print "%8d%6d%6d%6d%5d" %(num,left, below, above, 2)
#        else:
#            if i==l:
#                right=j
#            else:
#                right=num+h
#            print "%8d%6d%6d%6d%5d" %( num,below, right, above, 2)
