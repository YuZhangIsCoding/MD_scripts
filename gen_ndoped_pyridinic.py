# python script that generates pyridinic N doped graphene
# dimensions
import numpy as np

unit = [4, 6]
repeat = [4, 4]

bcc = 0.142 # bond length of C-C bond
bcn = 0.1339 # bond length of C-N bond
db = bcc*np.cos(np.radians(60))-(bcc**2*(np.cos(np.radians(60)))**2-(bcc**2-bcn**2))**0.5
dbx = db/2
dby = dbx*3**0.5

l = unit[0]*repeat[0]
h = unit[1]*repeat[1]

area = l*3/2*h/2*3**0.5*bcc**2
e0 = 0.160217657

molname = 'GPH'
#molname = raw_input('Please input the molecule name:\n')
zgap = 5
#zgap = input('Please input the channel height (nm):\n')

sigma = input('Please input the surface charge density (C/m^2):\n')

delta = sigma*area/e0/((unit[0]*unit[1]-1)*repeat[0]*repeat[1])

Nlist = [[1, 1], [1, 3], [2, 2]]
CNlist = [[0, 1], [0, 3], [1, 0], [1, 4], [2, 1], [2, 3]]

outfile = open('strcture.gro','w')
itpfile = open('graphene.itp','w')

outfile.write('N3V graphene\n')
outfile.write(str((unit[0]*unit[1]-1)*repeat[0]*repeat[1]*2)+'\n')

itpfile.write('[ moleculetype ]\n; molname\tnrexcl\n%s\t\t\t1\n' %(molname))
itpfile.write('[ atoms ]\n')
itpfile.write('; id type    res residu  at  cg  charge   mass\n')
natom = 0
for seq in [0, 1]:
    for i in range(0,l):
        for j in range(0,h):
            if i%2 == 0:
                coord = [1.5*i*bcc+bcc/2*(j%2), bcc/2*3**0.5*j, seq*zgap]  
            else:
                coord = [0.5*bcc+1.5*i*bcc-bcc/2*(j%2), bcc/2*3**0.5*j, seq*zgap]
            atomname = 'CG'
            atomtype = 'grC'
            molmass = 12.011
            charge = delta*(-1)**seq
            if i%unit[0] == 1 and j%unit[1] == 1:
                coord[0] = coord[0]-dbx
                coord[1] = coord[1]-dby
                atomname = 'NG'
                atomtype = 'grN'
                charge = delta*(-1)**seq-0.38
                molmass = 14.0067
            elif i%unit[0] == 1 and j%unit[1] == 3:
                coord[0] = coord[0]-dbx
                coord[1] = coord[1]+dby
                atomname = 'NG'
                atomtype = 'grN'
                charge = delta*(-1)**seq-0.38
                molmass = 14.0067
            elif i%unit[0] == 2 and j%unit[1] == 2:
                coord[0] = coord[0]+db
                atomname = 'NG'
                atomtype = 'grN'
                charge = delta*(-1)**seq-0.38
                molmass = 14.0067
            if [i%unit[0], j%unit[1]] in CNlist:
                atomname = 'CG'
                atomtype = 'grC'
                charge = delta*(-1)**seq+0.19
            if i%unit[0] != 1 or j%unit[1] != 2:
                natom += 1
                outfile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(1, 'GPH', atomname, natom,coord[0], coord[1], coord[2]))
                itpfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(natom, atomtype, 1, 'GPH', atomname, natom, charge, molmass))

x=l*3/2*bcc
y=h/2*3**0.5*bcc
z=(x+y)/2
outfile.write( "%10.5f%10.5f%10.5f" %(x, y , z))
itpfile.close()
outfile.close()
