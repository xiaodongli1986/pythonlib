
import numpy as np
import struct, sys

printstr = 'Usage: py_binary_ran filename -format your_format -nran nran -shape box -size boxsize -fmt3_headfile headfile'+\
                '\n\tformats: '+\
                '\n\t        1/CUTE_subsample format: \n\t\t\tboxsize/par-mass/omegam/h/a/redshift/nran/[x,y,z,vx,vy,vz,w]*nran/nran'+\
                '\n\t        2/simpler_xyzw/xyzw: \n\t\t\tinteger-npar + [x,y,z,w]*npar + integer-npar'+\
                '\n\t        3/xyzvxvyvzw/pos_vel_w: \n\t\t\tint-block+64*double+int-block + int-block+[x,y,z,vx,vy,vz,w]*nran+int-block'+\
                '\n\tshape: '+\
                '\n\t        box / sphere / shell (1-eighth of a sphere): uniform random distributed in cavity of given shape'+\
                '\n\t        input_file_spherical / input_file_shell \n\t\t\tGenerating a spherical / shell (1-eighth of sphere) shape sample, whose N(r) varis according to the N(r) of the inputfile \n\t\t\t   Must give: -inputfile -binsize ;  \n\t\t\t   Optional: \n\t\t\t\t-rmin/-rmax (if not given, use the sample\'s rmin and rmax) \n\t\t\t\t-drop_at_rmax (drop some sample near the maximal R boundary of the sample, to avoid the RSD shifting-out) \n\t\t\t\t-input_file_format  by default ascii; also support fmt3\n\t\t\t\t-rat (e.g., 10 means 10 times number in the input_file) '+\
                '\n\t        Example:     py_binary_ran FILENAME   -format 3   -shape input_file_shell   -rmin 0   -binsize 30   -rat 10   -input_file_format 3   -inputfile YOURFILE   '

print printstr+'\n'

# default settings
binaryformat = 'CUTE_subsample'
nran = 10000
shape = 'box'
size = 0
headfile = None

# default parameters for head
parmass = 1
omegam = 0.3071
h = 0.7 
a = 1
redshift = 0
weos = -1

# for input_file_spherical and input_file_shell
inputfile = None
binsize = None
rmin = None
rmax = None
drop_at_rmax = None
rat = None
inputfile_format = 'ascii'

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
        elif str2 in ['3', 'xyzvxvyvzw', 'pos_vel_w']:
            binaryformat = 'pos_vel_w'
        elif str2 in ['4','Nbody', 'nbody', 'gadget' ]:
            binaryformat = 'gadget'
        else:
            print 'Unknown format: ', str1, str2
            sys.exit()
    elif str1 == '-nran':
        nran = int(str2)
    elif str1 == '-fmt3_headfile':
        headfile = str2
    elif str1 == '-shape':
        if str2 not in ['box', 'sphere', 'shell', 'input_file_spherical', 'input_file_shell']:
            print 'Unknown shape: ', str1, str2
            sys.exit()
        shape = str2
    elif str1 == '-size':
        size = float(str2)
    elif str1 == '-inputfile':
        inputfile = str(str2)
    elif str1 == '-rmin':
        rmin = float(str2)
    elif str1 == '-rmax':
        rmax = float(str2)
    elif str1 == '-binsize':
        binsize = float(str2)
    elif str1 == '-drop_at_rmax':
        drop_at_rmax = float(str2)
    elif str1 in ['-inputfile_format', '-input_file_format', '-input_file_fmt', '-inputfile_fmt']:
        inputfile_format = str2
    elif str1 == '-rat':
        rat = float(str2)
    else:
        print 'Unknown options: ', str1, str2
        sys.exit()

print '\t', nran, 'randoms in', shape, 'with size',size,' | binary-file-format =', binaryformat

def gadget_head(nran,size,omegam,h):
    head = [nran,0,0,0,0,0] + [1.]*6 + [0.,0.] + [0,0,nran,0,0,0,0,0,0,1] + [size,omegam,1.-omegam,h] + [0]*24
    return head, struct.pack("6i6d2d10i4d24i",*head)

### tools for input_file_spherical and input_file_shell

if shape in ['input_file_spherical', 'input_file_shell']:
    if binsize == None:
        print ' ERROR! found binsize None'; sys.exit()
    else:
        print ' read-in binsize as ', binsize
    if rat == None:
        print ' Warning! found rat == None. '
    else:
        print ' read-in rat as ', rat
    if inputfile_format == 'ascii':
        print ' read in ascii : ', inputfile
        input_data = np.loadtxt(inputfile)[:,:3]
    elif inputfile_format in ['fmt3', 'format3', '3']:
        print ' read in binary fmt3: ', inputfile
        nowf1 = open(inputfile,'rb')
        block1 = struct.unpack('i',nowf1.read(4))[0]; head  = struct.unpack('64d',nowf1.read(8*64)); block2 = struct.unpack('i',nowf1.read(4))[0]
        noutput, boxsize, parmass, redshift, omegam, h = head[:6]; # need to add weos!
        if head[6] <= -0.1:
            weos = head[6]
        else:
            weos = -1.0
        print '\tread-in head: npar-total, boxsize, parmass, redshift, omegam, h = ', head[:6]
        block1 = struct.unpack('i',nowf1.read(4))[0]; npar = int(block1/7/4 + 0.1)
        input_data = np.array(struct.unpack(npar*'7f',nowf1.read(npar*7*4))).reshape(-1,7)[:,:3]
        block2 = struct.unpack('i',nowf1.read(4))[0]; npar2 = int(block2/7/4 + 0.1)

    r = (input_data[:,0]**2 + input_data[:,1]**2 + input_data[:,2]**2)**0.5
    if rat == None:
        rat = nran / float(len(r))
        print ' set rat = nran / #-line of inputfile = ', rat
    if rmin == None:
        rmin = min(r); print ' set rmin as ', rmin
    if rmax == None:
        rmax = max(r); 
        if drop_at_rmax != None:
            rmax -= drop_at_rmax
            print ' considering drop_at_rmax', drop_at_rmax
        print ' set rmax as ', rmax
    nbin = int((rmax-rmin) / float(binsize)+0.5)
    binsize = (rmax-rmin) / nbin
    print ' set nbin, binsize as ', nbin, binsize

    nums, edges = np.histogram(r, range=(rmin,rmax), bins=nbin)
    nums = nums*rat; nums = nums.astype('int32'); nran = sum(nums)
    print ' set nums as ', nums
    print ' set r-edges as ', edges


nowf = open(filename,'wb')

def yield_xyz():
	if shape =='box':
            for iran in range(nran):
		x,y,z = np.random.uniform(0,size), np.random.uniform(0,size), np.random.uniform(0,size)
                yield [x,y,z]
	elif shape == 'sphere':
            for iran in range(nran):
		x,y,z = np.random.uniform(-size,size), np.random.uniform(-size,size), np.random.uniform(-size,size)
		radius = (x*x+y*y+z*z)**0.5
		while radius > size:
			x,y,z = np.random.uniform(-size,size), np.random.uniform(-size,size), np.random.uniform(-size,size)
			radius = (x*x+y*y+z*z)**0.5
                yield [x,y,z]
	elif shape == 'shell':
            for iran in range(nran):
		x,y,z = np.random.uniform(0,size), np.random.uniform(0,size), np.random.uniform(0,size)
		radius = (x*x+y*y+z*z)**0.5
		while radius > size:
			x,y,z = np.random.uniform(0,size), np.random.uniform(0,size), np.random.uniform(0,size)
			radius = (x*x+y*y+z*z)**0.5
                yield [x,y,z]
        elif shape in ['input_file_spherical', 'input_file_shell']:
            for iedge in range(len(nums)):
              r1, r2 = edges[iedge], edges[iedge+1]; 
              print '   \t  generating %12i'%nums[iedge] ,'randoms in the bin of  %20.3f'%r1, '< r < ', '%.3f'%r2
              for ipar in range(nums[iedge]):
                r = np.random.uniform(r1**3, r2**3)**(1/3.) # uniform distribution with r cube
                cos_theta = np.random.uniform(-1,1); sin_theta = (1-cos_theta**2)**0.5; 
                phi = np.random.uniform(0,2*np.pi); sin_phi = np.sin(phi); cos_phi = np.cos(phi)
                x, y, z = r*sin_theta*cos_phi, r*sin_theta*sin_phi, r*cos_theta
                if shape == 'input_file_shell':
                    x, y, z = abs(x), abs(y), abs(z)
                yield [x,y,z]


if binaryformat == 'CUTE_subsample':
    head = [size,parmass,omegam,h,a,redshift,nran]
    print '\tcreate a head for the data:',head
    nowf.write(struct.pack("6f1i",*head))
    #for iran in range(nran):
    for X in yield_xyz():
        #nowf.write(struct.pack("3f",*np.random.uniform(0,size,3)))
        #nowf.write(struct.pack("3f",*yield_xyz()))
        nowf.write(struct.pack("3f",*X))
        nowf.write(struct.pack("3f",0,0,0))
    nowf.write(struct.pack("i",nran))
elif binaryformat == 'xyzw':
    head = [nran]
    print '\tcreate a head for the data:',head
    nowf.write(struct.pack("1i",*head))
    #for iran in range(nran):
    for X in yield_xyz():
        #nowf.write(struct.pack("3f",*np.random.uniform(0,size,3)))
        #nowf.write(struct.pack("3f",*yield_xyz()))
        nowf.write(struct.pack("3f",*X))
        nowf.write(struct.pack("1f",1))
    nowf.write(struct.pack("i",nran))
elif binaryformat == 'pos_vel_w':
    if headfile !=None:
        for nowstr in open(headfile,'r').readlines():
            if nowstr.split()[0][0] == '#':
                continue
            #head = [float(xx) for xx in open(headfile,'r').readline().split()]
            head = [float(xx) for xx in nowstr.split()]
            head = head + [0 for row in range(64-len(head))]
            head[0] = nran
            break
    else:
        head = [nran, size, parmass, redshift, omegam, h, weos]+[0. for row in range(64-7)]
    nowf.write(struct.pack('i',64*8)); nowf.write(struct.pack('64d',*head)); nowf.write(struct.pack('i',64*8))   
    print 'nran is ', nran
    print '4*7*nran is ', nran*7*4
    nowf.write(struct.pack('i',nran*7*4))
    #for iran in range(nran):
    for iran, X in enumerate(yield_xyz()):
            x,y,z = X
            nowf.write(struct.pack("7f",x,y,z,0,0,0,1))
            if iran ==0:       print('\t begin with x,y,z = '+str(x)+' '+str(y)+' '+str(z))
            if iran ==nran-1:  print('\t end with   x,y,z = '+str(x)+' '+str(y)+' '+str(z))
    nowf.write(struct.pack('i',nran*7*4))
elif binaryformat == 'gadget':
    # head
    nowf.write(struct.pack("1i",666))
    nowf.write(gadget_head(nran,size,omegam,h)[1])
    nowf.write(struct.pack("1i",666))
    # xyz
    nowf.write(struct.pack("1i",6666))
    #for iran in range(nran):
    for X in yield_xyz():
            x,y,z = X
            nowf.write(struct.pack("3f",x,y,z))
            if iran ==0:       print('\t begin with '+str(x)+' '+str(y)+' '+str(z))
            if iran ==nran-1:  print('\t end with   '+str(x)+' '+str(y)+' '+str(z))
            #nowf.write(struct.pack("3f",*np.random.uniform(0,size,3)))
    nowf.write(struct.pack("1i",6666))
    # vxvyvz
    nowf.write(struct.pack("1i",66666))
    for iran in range(nran):
            nowf.write(struct.pack("3f",0.,0.,0.))
    nowf.write(struct.pack("1i",66666))


nowf.close()
print 'Finishing writing ', nran, 'particles to ', filename

