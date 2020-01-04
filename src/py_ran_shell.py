from __future__ import division
import sys
import random
import numpy as np
#from __future__ import division
print 'Usage: \n\t input_filename output_filename times bins weight(optional, ranw or 1w or now) nosuffix(T,F) 1-eighth-shell(Y or N)'

print '\ntimes is the number of random particles generated with a real particle, bins is the size of the shell in radius.'
print '\nExample py_ran_shell inputf outputf 2 100 1w \n'
print '## the last argument tells whether generating 1/8 shell(xyz>0)'
args = sys.argv

if len(args)<4:
        print 'ERROR (ran_shell.py)! Must have at least 4 arguments.'
        sys.exit()

weight = 'now'
nosuffix=True
noshell_1_8=True
inputf, outputf, nran, bins = args[1],args[2],int(float(args[3])),int(float(args[4]))
if len(args) > 5:
        if args[5] in ['ranw','1w','now']:
                weight = args[5]
        else:
                print 'ERROR(ran_shell.py)! weight type must be ranw, 1w, now: ', weight
                sys.exit()

if len(args) > 6:
        if args[6][0] in ['T','F','Y','N']:
                if args[6][0] == 'T':
                        nosuffix = True
                elif args[6][0] == 'F':
                        nosuffix = False
                elif args[6][0] == 'Y':
                        noshell_1_8 = False
                elif args[6][0] == 'N':
                        noshell_1_8 = True

        else:
                print 'ERROR(ran_shell.py)! nosuffxi must be True/False or 1-eighth-shell must be Y/N'
                sys.exit()
if len(args) > 7:
        if args[7][0] in ['Y','N']:
                if args[7][0] == 'Y':
                        noshell_1_8 = False
                else:
                        noshell_1_8 = True
        else:
                print 'ERROR(ran_shell.py)!  1-eighth-shell must be Y/N', shll_1_8
                sys.exit()


print 'Reading the data: ', inputf

data = np.loadtxt(inputf)
X = data[:,0];Y = data[:,1];Z = data[:,2]
R = [np.sqrt(X[row]**2.0 + Y[row]**2.0 + Z[row]**2.0) for row in range(len(X))]

mR = max(R)
mRrange = mR/bins
precord = np.zeros(int(bins))

for i in range(len(R)):
        j = int(R[i]/mRrange)
        if j == bins:  ### maxR
            j = j - 1
        
        precord[int(j)] = precord[int(j)] + 1

if nosuffix:
        ofilename = outputf
else:
        ofilename = outputf +'.bins.'+str(bins)+'.times.'+str(nran)

print 'write to', ofilename, '...'
outf = open(ofilename, 'w')

for i in range(int(bins)):
        rrange = i*mRrange     
        if noshell_1_8:
            cost = np.random.random(size=int(nran*precord[i])) * 2 - 1
            t = np.arccos(cost)
            p = np.random.random(size=int(nran*precord[i])) * 2 * np.pi 
        else:
            cost = np.random.random(size=int(nran*precord[i])) 
            t = np.arccos(cost)
            p = np.random.random(size=int(nran*precord[i])) * np.pi / 2

        x = np.sin(t)*np.cos(p)
        y = np.sin(t)*np.sin(p)
        z = np.cos(t)
        k_set = np.arange(0,int(nran*precord[i]),1)
        for k in k_set:
#                if i == int(bins-1):
#                        rlen = (mR-rrange)*np.random.random() + rrange
#                else:
                rlen = mRrange*(np.random.random()) + rrange
                x[k] = x[k] * rlen
                y[k] = y[k] * rlen
                z[k] = z[k] * rlen
                if weight == 'now':
                        wstr = ''
                elif weight == '1w':
                        wstr = '1'
                elif weight == 'ranw':
                        wstr = str(random.uniform(0,1))
                outstr = str(x[k])+' '+str(y[k])+' '+str(z[k])+' '+wstr+'\n'
                outf.write(outstr)
outf.close()

print 'Finishing writing. '

#mR = max(R)
#sect = mR//bins + 1

#precord = np.zeros(int(sect))

#for i in range(len(R)):
#        j = R[i]//bins
#        precord[int(j)] = precord[int(j)] + 1

#if nosuffix:
#        ofilename = outputf
#else:
#        ofilename = outputf +'.bins.'+str(bins)+'.random.'+str(len(R)*nran)

#print 'write to', ofilename, '...'
#outf = open(ofilename, 'w')

#for i in range(int(sect)):
#       rlen = bins*np.sqrt(np.random.random())       
#        rbin = i*bins     
#        t = np.random.random(size=int(nran*precord[i])) * np.pi
#        cost = np.random.random(size=int(nran*precord[i])) * 2 - 1
#        t = np.random.random(size=int(nran*precord[i])) * np.pi - np.pi
#        t = np.arccos(cost)
#        p = np.random.random(size=int(nran*precord[i])) * 2 * np.pi 
#        x = np.sin(t)*np.cos(p)
#        y = np.sin(t)*np.sin(p)
#        z = np.cos(t)
#        k_set = np.arange(0,int(nran*precord[i]),1)
#        for k in k_set:
#                if i == int(sect-1):
#                        rlen = (mR-rbin)*np.random.random() + rbin
#                else:
#                        rlen = bins*np.random.random() + rbin
#                x[k] = x[k] * rlen
#                y[k] = y[k] * rlen
#                z[k] = z[k] * rlen
#                if weight == 'now':
#                        wstr = ''
#                elif weight == '1w':
#                        wstr = '1'
#                elif weight == 'ranw':
#                        wstr = str(random.uniform(0,1))
#                outstr = str(x[k])+' '+str(y[k])+' '+str(z[k])+' '+wstr+'\n'
#                outf.write(outstr)
#outf.close()
#
#print 'Finishing writing. '

