#!/bin/bash
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import datetime

print 'This is a script to generate a quick plot for your data, use -h to see more details'
if '-h' in sys.argv:
    print '-i: input file name\n\
-o: output file name\n\
-h: show help'

if len(sys.argv) == 1:
    print 'Please input at least one command or one file!'
else:
    if '-i' in sys.argv:
        inname = sys.argv[(sys.argv.index('-i')+1)]
    if '-o' in sys.argv:
        outname = sys.argv[(sys.argv.index('-o')+1)]
    else:
        now=datetime.datetime.now()
        outname = 'out_%s-%s-%s.pdf' %(now.month, now.day, now.year)
        print 'No filename specified, gonna save the file as %s' % outname
    while(True):
        if os.path.isfile(outname):
            print 'Opps!!! The filename %s already taken, ' %outname,
            outname = 'r_'+outname
            print 'gonna chage the filename as %s' %outname
        else:
            break
    datainput = np.loadtxt(inname)#,dtype=[('x',float),('y',float)])
    fig1 = plt.figure()
    legendname = ['data'+str(i) for i in range(1, len(datainput.T))]
    for column in datainput.T[1:,]:
        plt.plot(datainput.T[0, :], column)
    plt.legend(legendname)
    fig1.savefig(outname)
