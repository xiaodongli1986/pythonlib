
import sys, random, numpy

print 'Usage: \n\tEXE filename n-particle size shape(box/sphere/shell) weight(optional, ranw or 1w or now) nosuffix(T,F)'
print 'Example:\n\tpy_ran file0 10000 100 box 1w \n'
args = sys.argv

if len(args) <4:
	print 'ERROR (ran.py)! Must have at least 4 arguments.'
	sys.exit()

weight ='now'
nosuffix=True

filestr, nran, size, shape = args[1], int(float(args[2])), float(args[3]), args[4]
if not shape in ['box', 'sphere', 'shell']:
	print 'ERROR (ran.py)! shape must be box, sphere, shell:', shape
	sys.exit()
if len(args) > 5:
	if args[5] in ['ranw', '1w', 'now']:
		weight = args[5]
	else:
		print 'ERROR (ran.py)! weight type must be ranw, 1w, now: ', weight
		sys.exit()
if len(args)>6:
	if args[6][0] in ['T', 'F']:
		if args[6][0] == 'T': 
			nosuffix = True
		else:
			nosuffix = False
	else:
		print 'ERROR (ran.py)! nosuffxi must be True/False ', nosuffix
		sys.exit()
		 
if nosuffix:
	filename = filestr
else:
	filename = filestr+'.'+str(nran)+'random.size%i'%size+'.'+shape
print 'Creating random sample: ', filename

nowf = open(filename, 'w')
iran = 0
while iran < nran:
	if shape =='box':
		x,y,z = random.uniform(0,size), random.uniform(0,size), random.uniform(0,size)
	elif shape == 'sphere':
		x,y,z = random.uniform(-size,size), random.uniform(-size,size), random.uniform(-size,size)
		radius = (x*x+y*y+z*z)**0.5
		while radius > size:
			x,y,z = random.uniform(-size,size), random.uniform(-size,size), random.uniform(-size,size)
			radius = (x*x+y*y+z*z)**0.5
	elif shape == 'shell':
		x,y,z = random.uniform(0,size), random.uniform(0,size), random.uniform(0,size)
		radius = (x*x+y*y+z*z)**0.5
		while radius > size:
			x,y,z = random.uniform(0,size), random.uniform(0,size), random.uniform(0,size)
			radius = (x*x+y*y+z*z)**0.5
	if weight == 'now':
		wstr = ''
	elif weight == '1w':
		wstr = '1'
	elif weight == 'ranw':
		wstr = str(random.uniform(0,1))
	nowstr = str(x)+' '+str(y)+' '+str(z)+' '+wstr+'\n'
	nowf.write(nowstr)
	iran += 1
nowf.close()
print 'Finishing writing ', nran, 'lines.'
	
