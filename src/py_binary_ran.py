
import numpy as np
import struct, sys

printstr = 'Usage: py_binary_ran filename -format your_format -nran nran -shape box -size boxsize'+\
            '\n\tformat=1/CUTE_subsample: boxsize/par-mass/omegam/h/a/redshift/nran/[x,y,z,vx,vy,vz]*nran/nran'+\
            '\n\tformat=2/xyzw: nran/[x,y,z,w]*nran/nran'

print printstr+'\n'

# default settings
binaryformat = 'CUTE_subsample'
nran = 10000
shape = 'box'
size = 0

# default parameters for head
parmass = 1
omegam = 0.3071
h = 0.7 
a = 1
redshift = 0

args = sys.argv

if len(args) <= 1:
    sys.exit()

filename = args[1]
print '\toutput file-name = ',filename
nopt =int( (len(args)-2)/2)

for iopt in range(nopt):
    str1, str2 = args[2+iopt*2], args[2+iopt*2+1]
    print str1, str2
    if str1 == '-format':
        if str2 in ['1','CUTE_subsample', 'CUTEsubsample', 'CUTE-subsample']:
            binaryformat = 'CUTE_subsample'
        elif str2 in ['2','xyzw', ]:
            binaryformat = 'xyzw'
        else:
            print 'Unknown format: ', str1, str2
            sys.exit()
    elif str1 == '-nran':
        nran = int(str2)
    elif str1 == '-shape':
        if str2 not in ['box']:
            print 'Unknown shape: ', str1, str2
            sys.exit()
        shape = str2
    elif str1 == '-size':
        size = float(str2)
    else:
        print 'Unknown options: ', str1, str2
        sys.exit()

print '\t', nran, 'randoms in', shape, 'with size',size,' | binary-file-format =', binaryformat

nowf = open(filename,'wb')
if binaryformat == 'CUTE_subsample':
    head = [size,parmass,omegam,h,a,redshift,nran]
    print '\tcreate a head for the data:',head
    nowf.write(struct.pack("6f1i",*head))
    for iran in range(nran):
        if shape == 'box':
            nowf.write(struct.pack("3f",*np.random.uniform(0,size,3)))
        nowf.write(struct.pack("3f",0,0,0))
    nowf.write(struct.pack("i",nran))
elif binaryformat == 'xyzw':
    head = [nran]
    print '\tcreate a head for the data:',head
    nowf.write(struct.pack("1i",*head))
    for iran in range(nran):
        if shape == 'box':
            nowf.write(struct.pack("3f",*np.random.uniform(0,size,3)))
        nowf.write(struct.pack("1f",1))
    nowf.write(struct.pack("i",nran))
nowf.close()
print 'Finishing writing ', nran, 'particles to ', filename

