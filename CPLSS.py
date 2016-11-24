
import commands
pythonlibPATH=commands.getoutput('echo $pythonlibPATH')
execfile(pythonlibPATH+'/stdA.py')
execfile(pythonlibPATH+'/Tpcftools.py')
pyfile=pythonlibPATH+'/'+'CPLSS.py'

def read_filament(filename, printinfo=False):
	FILAMENTS = []
	FIL_BEGIN = False
	nowf = open(filename, 'r')
	while True:
		nowstr = nowf.readline();
		if nowstr == '':
			break
		if nowstr == '[FILAMENTS]\n':
			break
	nFIL = int(nowf.readline().split()[0])
	if printinfo: print '  ', nFIL, ' filaments identified in the file ', filename
	while True:
		nowstr = nowf.readline();
		if nowstr == '' or nowstr == '[CRITICAL POINTS DATA]\n':
			break
		nsegp = int(nowstr.split()[2])
		FILAMENTS.append ([    [float(xx) for xx in nowf.readline().split()] for row in range(nsegp)  ])
	if printinfo: print len(FILAMENTS), ' FILAMENTS read. FIRST, LAST ONE IS: ', FILAMENTS[0], FILAMENTS[nFIL-1]
	nowf.close()
	return FILAMENTS 
			

