#!/Users/yuzhang/anaconda/bin/python
# Filename:     calc_autoco.py
# Description:  This is a python script to calculate the auto correlation 
#               functions. The auto correlation function could be descriped as
#               f = <c(0)c(t)>/<c(0)c(t)>, where c(t) is a binary function that
#               equals to 1 if the value was continously found in t frame, and
#               equals to 0 otherwise.
#               This algorithm uses (frame+1)/2 origins and frame-(frame+1)/2 
#               max length.
# Date:         10-18-2017 Created
#               10-19-2017 Added counter variable to count the average 
#               appearance of values.
import numpy as np
import matplotlib.pyplot as plt
import argparse, collections

parser = argparse.ArgumentParser(description = 'user specified options')
parser.add_argument('-i', '--input', dest = 'input', help = 'input file')
parser.add_argument('-o', '--output', dest = 'output', help = 'output name')

args = parser.parse_args()

if args.input:
    filename = args.input
else:
    filename = input('Please enter the input file:\n')

if args.output:
    outname = args.output
else:
    outname = 'Autoco'

counters = collections.Counter()

myfile = open(filename, 'r')
data = []
for line in myfile:
    data.append({int(i) for i in line.split()})
    counters.update(data[-1])

orig = (len(data)+1)/2
length = len(data)-orig
temps = data[:orig]
counts = [[] for _ in range(length+1)]
counts[0] = [len(item) for item in temps]
for i, item in enumerate(data):
    for j in range(max(0, i-length), min(i, orig)):
        temps[j] = temps[j].intersection(item)
        counts[i-j].append(len(temps[j]))
out = [sum(item) for item in counts]
out = [float(i)/out[0] for i in out]
plt.plot(range(len(out)), out)
plt.xlabel('Frame')
plt.ylabel('Auto correlation function')
np.savetxt(outname+'.txt', out)
plt.savefig(outname+'.pdf')
#import pdb
#pdb.set_trace()
fig = plt.figure()
vals = np.array(counters.values(), dtype = float)/len(data)
plt.hist(vals, np.arange(0, 1, 0.01), density = True)
plt.xlabel('Fraction of total frame')
plt.ylabel('Probability density')
plt.savefig(outname+'_avg.pdf')
