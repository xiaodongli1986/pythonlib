import numpy as np
import sys, struct


printstr = 'Usage: py_binary_ascii_conv   inputfile    outputfile    mode(=ba/ab)   format(=1/2)/3 '+\
        '\n\tmode: 0/ba means binary-to-ascii; 1/ab means ascii-to-binary (not supported yet); '+\
        '\n\tformats: '+\
        '\n\t        1/CUTE_subsample format: \n\t\t\tboxsize/par-mass/omegam/h/a/redshift/nran/[x,y,z,vx,vy,vz,w]*nran/nran'+\
        '\n\t        2/simpler_xyzw/xyzw: \n\t\t\tinteger-npar + [x,y,z,w]*npar + integer-npar'+\
        '\n\t        3/xyzvxvyvzw/pos_vel_w: \n\t\t\tint-block+64*double+int-block + int-block+[x,y,z,vx,vy,vz,w]*nran+int-block'+\
        '\n\t -weightcol (column of the weight in the ascii file; can use in case of mode=1'+\
        '\n\t -fmt3_headfile headfile for the fmt3 output; should be an ascii with 1 row: (noutput, boxsize, mass, redshfit, Omega0, HubbleParam, De_w)'+\
        '\n\t -fmt3_set_v_as_0   for fmt3 file: input file only provides x,y,z,(weight);  automatically skip read-in of v by setting vx,vy,vz as 0'

print printstr+'\n'

###
args = sys.argv
if len(args) < 5:
    sys.exit()

###
inputfile, outputfile, mode, binaryformat = args[1:5]
headfile = None
infofile = outputfile+'.info'
fmt3_set_v_as_0 = False
print 'inputfile, outputfile, mode, binaryformat = ', inputfile, outputfile, mode, binaryformat 

###
weightcol = -1

nopt = len(args)-5
for iopt in range(int(nopt/2)):
    str1, str2 = args[5+2*iopt], args[5+2*iopt+1]
    if str1 in ['-weightcol', '-wcol']:
        weightcol = int(str2)
        print 'set weightcol = ', weightcol
    elif str1 == '-fmt3_headfile':
        headfile = str2
        print 'set headfile = ', headfile
    elif str1 == '-fmt3_set_v_as_0':
        if str2[0] in 'tT':
            fmt3_set_v_as_0 = True
        elif str2[0] in 'fF':
            fmt3_set_v_as_0 = False
        else:
            print 'Wrong value of fmt3_set_v_as_0: must begin with t/T/f/F!'
            sys.exit()
        print 'set fmt3_set_v_as_0 = ', fmt3_set_v_as_0

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
elif binaryformat in ['3', 'xyzvxvyvzw', 'pos_vel_w']:
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
        nowf2.write(struct.pack('i',64*8))
        try:
            #nowf4 = open(headfile,'r')
            #for nowstr in nowf4.readlines():
            #    if nowstr == '': break
            #    nowstrs = nowstr.split()
            #    if nowstrs[0] == '#': continue
            #    #headinfo = [float(xx) for xx in nowf4.readline().split()]; headinfo[0] = npar
            #    headinfo = [float(xx) for xx in nowstrs]; headinfo[0] = npar; break
            headinfo = np.loadtxt(headfile); headinfo = list(headinfo); headinfo[0] = npar
            print '\tread in headinfo: ', headinfo
            nowf2.write(struct.pack('64d',*(headinfo+[0 for row in range(64-len(headinfo))])))
            #print '\twrite head according to content of file: ',inputfile+'.head'
            #nowf4.close()
        except:
            print '\theadfile read-write error : ',headfile
            print '\theadfile read-write error : ',headfile
            print '\theadfile read-write error : ',headfile
            sys.exit()
            nowf2.write(struct.pack('64d',*[-1.0 for row in range(64)]))

        nowf2.write(struct.pack('i',64*8))
        nowf2.write(struct.pack('i',npar*7*4))
        while True:
            nowstr = nowf1.readline()
            if nowstr == '': break
            nowstrs = nowstr.split()
            if fmt3_set_v_as_0:
                x, y, z, vx, vy, vz = float(nowstrs[0]), float(nowstrs[1]), float(nowstrs[2]), 0, 0, 0
            else:
                x, y, z, vx, vy, vz = float(nowstrs[0]), float(nowstrs[1]), float(nowstrs[2]), float(nowstrs[3]), float(nowstrs[4]), float(nowstrs[5])
            if weightcol == -1:
                w = 1
            else:
                w = float(nowstrs[weightcol-1])
            nowf2.write(struct.pack('7f',x,y,z,vx,vy,vz,w))
        nowf2.write(struct.pack('i',npar*7*4))
        nowf2.close(); nowf3.close()
    elif mode in ['0','ba']: # binary to ascii
        nowf1 = open(inputfile,'rb')
        nowf2 = open(outputfile, 'w')
        nowf3 = open(infofile, 'w')
        nowf4 = open(outputfile+'.head', 'w')
        block1 = struct.unpack('i',nowf1.read(4))[0]
        head  = struct.unpack('64d',nowf1.read(8*64))
        block2 = struct.unpack('i',nowf1.read(4))[0]
#        n-particle/boxsize/part-mass/redshift/omegam/h/
        noutput, boxsize, parmass, redshift, omegam, h, w = head[:7]; noutput = int(noutput+0.1)
        print '\tread-in head: npar-total, boxsize, parmass, redshift, omegam, h, w (w=0 means w=-1) = ', head[:7]
        nowf3.write('read-in head: npar-total, boxsize, parmass, redshift, omegam, h, w = '+' '.join([str(xx) for xx in head[:7]])+'\n')
        nowf4.write(' '.join([str(xx) for xx in head]))
        print '\tread-in noutput = ', noutput
        print '\tget iline: ', str(noutput), 'lines in input binary file.'
        block1 = struct.unpack('i',nowf1.read(4))[0]; npar = int(block1/7/4 + 0.1)
        nowf3.write('get iline: '+str(noutput)+' lines in input binary file.\n')
        for iline in range(noutput):
            contents = struct.unpack('7f',nowf1.read(7*4))
            nowf2.write(' '.join([str(xx) for xx in contents])+'\n')
        block2 = struct.unpack('i',nowf1.read(4))[0]; npar2 = int(block2/7/4 + 0.1)
        if npar != npar2:
            print 'WARNING!!!! Inconsistency!', npar, npar2
            print 'WARNING!!!! Inconsistency!', npar, npar2
            print 'WARNING!!!! Inconsistency!', npar, npar2
        else:
            print '\t\tconfirm consistency!'
        npar = noutput
        nowf2.close(); nowf3.close()
elif binaryformat in ['2', 'simple_xyzw', 'xyzw']:
    if weightcol == -1: weightcol = 4
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
    elif mode in ['0','ba']: # ascii to binary
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

print 'Finishing writing ', npar, 'particles to ', outputfile

