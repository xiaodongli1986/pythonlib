import numpy as np
import sys, struct


printstr = 'Usage: py_binary_ascii_conv inputfile outputfile mode format \n\tmode: 0 means binary-to-ascii; 1 means ascii-to-binary (not supported yet); \n\tformat: 0=simple float32;(not supported yet) 1=CUTE_subsample format'



print printstr+'\n'

# default settings
args = sys.argv

if len(args) != 5:
    sys.exit()

inputfile, outputfile, mode, binaryformat = args[1:]
infofile = outputfile+'.info'
print 'inputfile, outputfile, mode, binaryformat = ', inputfile, outputfile, mode, binaryformat 

if binaryformat in ['1', 'CUTE_subsample', 'CUTEsubsample', 'CUTE-subsample']:
    if mode =='0': # 
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

    

print 'Finishing writing ', npar, 'particles to ', outputfile

