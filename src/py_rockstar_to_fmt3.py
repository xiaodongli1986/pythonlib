import sys, os
import numpy as np

printstr = 'Usage: py_rockstar_to_fmt3 inputfilename'



args = sys.argv

if len(args) <=1:
    print(printstr)
    sys.exit()

filename = args[1]
headfile = filename + '.head'
fm3file = filename + '_fmt3'
ascifile = filename + '_ascii'


nowf = open(filename, 'r')
nowf2 = open(ascifile, 'w')

iline = 0
while True:
    nowstr = nowf.readline()
    if nowstr == '':
        break
    else:
        nowstrs = nowstr.split()
        if nowstrs[0] == '#a':
            a = float(nowstrs[2]); redshift = 1/a -1.
        elif nowstrs[0] == '#Om':
            vals = [nowstrs[row] for row in [2,5,8]]
            vals = [xx[:len(xx)-1] for xx in vals]
            om, ol, h = [float(val) for val in vals]
        elif nowstrs[0] == '#Box':
            boxsize = float(nowstrs[2])
        elif nowstrs[0][0] != '#':
            nowf2.write(' '.join([nowstrs[row] for row in [8,9,10,11,12,13,2]])+'\n')
            iline += 1
nowf.close(); nowf2.close()
print ('boxsize, redshift, om, ol, h = ', boxsize, redshift, om, ol, h)
print ('In total read-in ', iline, 'halos.')
print ('ascii file write to :\n\t', ascifile)
print ('Headfile write to :\n\t', headfile)
noutput = iline

head = [noutput, boxsize, 1., redshift, om, h]
nowf = open(headfile,'w')
nowf.write(' '.join([str(xx) for xx in head]))
nowf.close()

cmd = 'py_binary_ascii_conv '+ascifile+' '+fm3file+' ab 3 -weightcol 7 -fmt3_headfile '+headfile
print('Run command: ',cmd, '..')
print(os.popen(cmd).readline())
