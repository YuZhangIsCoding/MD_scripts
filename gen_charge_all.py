## This is a script to add partial charges to functionalized graphene sheets
## This is to be added to the [atoms] section in .itp file

# Input dimensions and bond length for 2 pristine graphene layers
# And some constants to convert units
l=22
h=38
nlayer=2
bondlen=0.139
area=l*3/2*h/2*3**0.5*bondlen**2
e0=0.160217657


# nper: Number of atoms in each repeated cycle -- to make functional position less linear.
# Only for 5% functionalized, other percentage please change the numbers below.
# p1, p2: Position of dividing points.
# nf: Number of functional sites
nper=60
p1=26
p2=50
nf=((l*h-2)/nper+1)*3 # Number of functional sites

choice = input('Please select your system:\n\
    1. Hydrogenated\n\
    2. Epoxy functionalized\n\
    3. Hydroxyl functionalized\n\
    4. Pristine graphene\n')
sigma = input('Please input the surface delta density(C/m^2):\n')

if choice==1:
    delta=sigma*area/e0/(l*h+nf)
    myfile=open('hydrogen.charge','w')
    for i in range(1,l*h+1):
        if (i%nper==2 or i%nper==p1+1 or i%nper==p2+1):
            myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'CHA',1,'GPH','CHF',i,-0.053+delta,12.011))
        else:
            myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'CGP',1,'GPH','CG',i,0+delta,12.011))
    for i in range(1,nf+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i+l*h,'HGF',1,'GPH','HG',i+l*h,0.053+delta,1.008))
    for i in range(1,nlayer*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i+l*h+nf,'C0',1,'GPH','CG',i+l*h+nf,0,12.011))
    myfile.close()
elif choice==2:
    delta=sigma*area/e0/(l*h+nf)
    myfile=open('epoxy.charge','w')
    for i in range(1,l*h+1):
        if (i%nper==1 or i%nper == 2 or i%nper == p1+1 or i%nper == p1+2 or i%nper == p2+1 or i%nper == p2+2):
            myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i,'CHA',1,'GPH','CHF',i,0.1330+delta,12.011))
        else:
            myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i,'CGP',1,'GPH','CG',i,0+delta,12.011))
    for i in range(1,nf+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i+l*h,'OGF',1,'GPH','OG',i+l*h,-0.2660+delta,15.9994))
    for i in range(1,nlayer*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i+l*h+nf,'C0',1,'GPH','CG',i+l*h+nf,0,12.011))
    myfile.close()
elif choice==3:
    l=16
    h=28
    bondlen=-0.139
    area=l*3/2*h/2*3**0.5*bondlen**2
    nper=20
    nf=(l*h-2)/nper+1
    delta=sigma*area/e0/(l*h+2*nf)
    myfile=open('hydroxyl.charge','w')
    for i in range(1,l*h+1):
        if (i%nper==2):
            myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i,'CHA',1,'GPH','CHF',i,0.1330+delta,12.011))
        else:
            myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i,'CGP',1,'GPH','CG',i,delta,12.011))
    for i in range(1,nf+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i+l*h,'OGF',1,'GPH','OG',i+l*h,delta-0.5321,15.9994))
    for i in range(1,nf+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n" %(i+l*h+nf,'HGF',1,'GPH','HG',i+l*h+nf,0.3991+delta,1.008))
    for i in range(1,nlayer*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i+l*h+2*nf,'C0',1,'GPH','CG',i+l*h+2*nf,0,12.011))
    myfile.close()
elif choice==4:
    [l,h]=input('Please input the dimension numbers [l,h]:\n')
    bondlen=0.142
    area=l*3/2*h/2*3**0.5*bondlen**2
    charge=sigma*area/e0/l/h
    myfile=open('charge.dat','w')
    for i in range(1,l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'grC',1,'GHS','CG',i,charge,12.011))
    for i in range(l*h+1,(1+nlayer)*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'grC',1,'GHS','CG',i,0,12.011))
    for i in range((nlayer+1)*l*h+1,(nlayer+2)*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'grC',1,'GHS','CG',i,-charge,12.011))
    for i in range((nlayer+2)*l*h+1,2*(nlayer+1)*l*h+1):
        myfile.write("%5d%5s%5d%5s%5s%5d%11.7f%8.3f\n" %(i,'grC',1,'GHS','CG',i,0,12.011))
    myfile.close()
else:
    print 'System not recoganized\n'
