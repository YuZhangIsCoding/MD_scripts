#!/Users/yuzhang/anaconda/envs/py3/bin/python
import numpy as np

############ Functions ###########

def gen_frame(l,h,stype):
    '''General frame to generate a graphene sheet'''
    gframe = [[] for row in range(l*h)]
    for i in range(0,l):
        for j in range(0,h):
            seq = i*h+j
            if stype == 1:
                if i%2 == 0:
                    gframe[seq] = np.array([1.5*i+(j%2)/2.0,j*3**0.5/2])
                else:
                    gframe[seq] = np.array([0.5+1.5*i-(j%2)/2.0,j*3**0.5/2])
            elif stype == 2:
                if i%2 == 0:
                    gframe[seq] = np.array([1.5*i-(j%2)/2.0,j*3**0.5/2])
                else:
                    gframe[seq] = np.array([-0.5+1.5*i+(j%2)/2.0,j*3**0.5/2])
    return gframe

def gen_circle(pstart, lx, ly, lz, myfile):
    '''Generate a circle of graphene sheets'''

    print('Dimensions:',lx,ly,lz)
    nx = int(np.floor(lx/bondlen/3*2))+1
    ny = int(np.ceil(ly/bondlen/3**0.5))*2
    nz = int(np.floor(lz/bondlen/3*2))
    
    # Button and top sheets
    while True:
        lframe = gen_frame(nx-2, ny, 1)
        if nx%2 == 0:
            x_real = (1+max(lframe[-1][0],lframe[-2][0])+1)*bondlen
        else:
            x_real = (1+min(lframe[-1][0],lframe[-2][0])+2)*bondlen
        dx = lx-x_real
        print('x_real',x_real)
        print('dx',dx)
        if dx > 2*bondlen:
            nx += 1
        else:
            break
    print('Real length of y direction:', (lframe[-1][1]+3**0.5/2)*bondlen)
    
    # Left and right sheets
    while True:
        vframe= gen_frame(nz, ny, 1)
        lz_real = max(vframe[-1][0],vframe[-2][0])*bondlen
        dz = lz-lz_real
        print('lz_real',lz_real)
        print('dz',dz)
        if dz > 2*bondlen:
            nz += 1
        else:
            break
    
    print('Real dimension numbers:', nx-2, ny, nz)
    for i, item in enumerate(lframe):
        myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',i+1,pstart[0]+(item[0]+1)*bondlen+dx/2,item[1]*bondlen,pstart[1]))
    for i, item in enumerate(lframe):
        myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',i+1,pstart[0]+(item[0]+1)*bondlen+dx/2,item[1]*bondlen,lz+pstart[1]))
    for i, item in enumerate(vframe):
        myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',i+1, pstart[0], item[1]*bondlen, pstart[1]+dz/2+item[0]*bondlen))
    for i, item in enumerate(vframe):
        myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',i+1, lx+pstart[0], item[1]*bondlen, pstart[1]+dz/2+item[0]*bondlen))
    return 2*(nx-2+nz)*ny

#################### Main #####################

#### Inputs ####
global bondlen
bondlen = 0.142 # Specify the bond length for graphene

# Lengths for the out circle
lx = 12
ly = 4
lz = 3.5
myfile = open('slit_out.gro','w')
myfile.write('Undefined Name\n')
myfile.write('%5d\n' %(00000))

n_total = 0
n_total += gen_circle([0, 0], lx, ly, lz, myfile)
n_total += gen_circle([0.341, 0.341], lx-0.682, ly, lz-0.682, myfile)
n_total += gen_circle([0.682, 0.682], lx-0.682*2, ly, lz-0.682*2, myfile)
print('Total number of atoms:', n_total)

myfile.write("%10.5f%10.5f%10.5f" %(lx ,10 ,lz))
myfile.close()

