import numpy as np
import sys, struct


printstr = 'Usage: py_binary_ascii_conv   inputfile    outputfile    mode(=ba/ab)   format(=1/2) \n\n\tmode: 0/ba means binary-to-ascii; 1/ab means ascii-to-binary (not supported yet); \n\tformat: 1=CUTE_subsample format'+\
        '\n\tformat: 2/simpler_xyzw/xyzw: simple x,y,z,w. all in float32'+\
        '\n\t -weightcol 4 (column of the weight in the ascii file; used in case of mode=1'

print printstr+'\n'

###
args = sys.argv
if len(args) < 5:
    sys.exit()

###
inputfile, outputfile, mode, binaryformat = args[1:5]
infofile = outputfile+'.info'
print 'inputfile, outputfile, mode, binaryformat = ', inputfile, outputfile, mode, binaryformat 

###
weightcol = 4

nopt = len(args)-5
for iopt in range(int(nopt/2)):
    str1, str2 = args[5+2*iopt], args[5+2*iopt+1]
    if str1 in ['-weightcol', '-wcol']:
        weightcol = int(str2)
    else:
        print 'Unknown args! ', str1, str2
        sys.exit()
if binaryformat in ['1', 'CUTE_subsample', 'CUTEsubsample', 'CUTE-subsample']:
    if mode in ['0','ba']: # binary to ascii
        nowf1 = open(inputfile,'rb')
        nowf2 = open(outputfile, 'w')
        nowf3 = open(infofile, 'w')
        head = nowf1.read(7*4)
        head = struct.unpack("6f1i",head)
        print '\tread-in head:',head
        nowf3.write("head = "+" ".join([str(xx) for xx in head]))
        npar = head[6]
        for row in range(npar):
            data = np.array(struct.unpack("6f",nowf1.read(6*4))).astype('float32').reshape(-1)
            nowf2.write(' '.join(str(xx) for xx in data)+'\n')
        npar2 = struct.unpack("i",nowf1.read(4))[0]
        print '\tcheck consistency: npar, npar2 = ', npar, npar2
        if npar != npar2:
            print 'WARNING!!!! Inconsistency!'
            print 'WARNING!!!! Inconsistency!'
            print 'WARNING!!!! Inconsistency!'
        else:
            print '\t\tconfirm consistency!'
        nowf1.close(); nowf2.close(); nowf3.close()
elif binaryformat in ['2', 'simple_xyzw', 'xyzw']:
    if mode in ['1','ab']: # ascii to binary
        nowf1 = open(inputfile,'r')
        nowf2 = open(outputfile, 'wb')
        nowf3 = open(infofile, 'w')
        iline = 0
        while True:
            nowstr = nowf1.readline(); 
            if nowstr == '': break;
            iline+=1;
        print '\tcount line: ', str(iline), 'lines in input ascii file.'
        nowf3.write('count line: '+str(iline)+' lines in input ascii file.\n')
        npar = iline
        if iline > 4294967296:
            print 'ERROR! input_line shall not excess ', 4294967296
            sys.exit()
        nowf1.close(); nowf1= open(inputfile,'r')
        nowf2.write(struct.pack('i',npar))
        while True:
            nowstr = nowf1.readline()
            if nowstr == '': break
            nowstrs = nowstr.split()
            x, y, z, w = float(nowstrs[0]), float(nowstrs[1]), float(nowstrs[2]), float(nowstrs[weightcol-1])
            nowf2.write(struct.pack('4f',x,y,z,w))
        nowf2.write(struct.pack('i',npar))
        nowf2.close(); nowf3.close()
    elif mode in ['2','ba']: # ascii to binary
        nowf1 = open(inputfile,'rb')
        nowf2 = open(outputfile, 'w')
        nowf3 = open(infofile, 'w')
        npar = struct.unpack('i',nowf1.read(4))[0]
        print '\tget iline: ', str(npar), 'lines in input binary file.'
        nowf3.write('get iline: '+str(npar)+' lines in input binary file.\n')
        for iline in range(npar):
            x,y,z,w = struct.unpack('4f',nowf1.read(4*4))
            nowf2.write(str(x)+' '+str(y)+' '+str(z)+' '+str(w)+'\n')
        npar2 = struct.unpack('i',nowf1.read(4))[0]
        if npar != npar2:
            print 'WARNING!!!! Inconsistency!'
            print 'WARNING!!!! Inconsistency!'
            print 'WARNING!!!! Inconsistency!'
        else:
            print '\t\tconfirm consistency!'
        nowf2.close(); nowf3.close()
print 'Finishing writing ', npar, 'particles to ', outputfile

