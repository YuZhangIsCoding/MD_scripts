# itp file for pristine graphene with atom types and bonds, but no angles and dihedrals

print("[ moleculetype ]\n; molname\tnrexcl\nGPH\t\t\t1\n")
print("[ atoms ]")
print("; id type    res residu  at  cg  charge   mass")
# dimensions
l = 16
h = 28
for i in range(1,l*h+1):
    print("%5d%5s%5d%5s%5s%5d%8.4f%8.3f" %(i,'CGP',1,'GPH','CG',i,0,12.011))

#bond
print("[    bonds   ]")
for i in range(1,l+1):
    for j in range(1,h+1):
        num=h*(i-1)+j
        if j == 1:
            below = num+h-1
        else:
            below = num-1
            print("%8d%8d%5d" %(num, below,1))
        if j == h:
            above = num-h+1
        else:
            above = num+1
            print("%8d%8d%5d" %(num, above,1))
        if (i+num)%2 == 0:
            if i == 1:
                left = j+(l-1)*h
            else:
                left = num-h
                print("%8d%8d%5d" %(num, left,1))
