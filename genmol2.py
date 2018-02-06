#!/Users/yuzhang/anaconda/bin/python
print("@<TRIPOS>MOLECULE\n*****")

# dimensions
l=16
h=28

print l*h, l*h*3/2, 0, 0,0
print("SMALL\nNO CHARGES\n")
print("@<TRIPOS>ATOM")

wholepair = [[0]*4 for row in range(l*h)]

wholepair[0][0]=0
wholepair[0][1]=0

# general frame

for i in range(0,l):
    if (i%2 == 0):
        wholepair[h*i][0]=1.5*i
        wholepair[h*i+1][0]=wholepair[h*i][0]+0.5
    else:
        wholepair[h*i][0]=wholepair[h*(i-1)][0]+2
        wholepair[h*i+1][0]=wholepair[h*i][0]-0.5
    wholepair[h*i][1]=0
    wholepair[h*i+1][1]=0.5*3**0.5
    for j in range(2,h):
        wholepair[h*i+j][0]=wholepair[h*i+j-2][0]
        wholepair[h*i+j][1]=wholepair[h*i+j-2][1]+3**0.5

for i in range(0,l*h):
    wholepair[i][0]=wholepair[i][0]*0.139
    wholepair[i][1]=wholepair[i][1]*0.139
    wholepair[i][3]='CG'
    print "%5d%5s%8.3f%8.3f%8.3f%5s%5d%5s%5d" %(i+1,wholepair[i][3],wholepair[i][0],wholepair[i][1],wholepair[i][2],wholepair[i][3],1,'GHS',0)

print("@<TRIPOS>BOND")

count=0
for i in range(1,l+1):
    for j in range(1,h+1):
        num=h*(i-1)+j
        if j == 1:
            below = num+h-1
        else:
            below = num-1
        if j == h:
            above = num-h+1
        else:
            above = num+1
        if (i+num)%2 == 0:
            if i == 1:
                left = j+(l-1)*h
            else:
                left = num-h
            count = count + 1
            print "%5d%8d%5d%5d" %(count,num, left,1)
            count=count+1
            print "%5d%8d%5d%5d" %(count,num, above,1)
            count = count +1
            print "%5d%8d%5d%5d" %(count,num, below,1)
