
import commands
import numpy as np
import sys
import stdA



printstr = 'Usage: EXE filename '

cmdargs = sys.argv

if len(cmdargs) < 2:
	print printstr; sys.exit()

filenames = commands.getoutput('ls '+cmdargs[1]).split()
	

#filename = cmdargs[1]


for filename in filenames:
 nowf = open(filename, 'r')

 mins = [1.0e30 for row in range(100000)]
 maxs = [-1.0e30 for row in range(100000)]
 sums = [0 for row in range(100000)]
 sumsqs = [0 for row in range(100000)]
 avgs = [0 for row in range(100000)]
 stds = [0 for row in range(100000)]

 iline = 0
 nowchar = None
 while True:
	nowstr = nowf.readline()
	if nowstr == '': break
        nowstr = nowstr[:len(nowstr)-1]
	iline += 1
	if nowstr[0] in ['#']: continue
	if nowchar == None:
		try:
			nowstrs = nowstr.split()
			ncol = len(nowstrs)
                        xs = [float(x) for x in nowstrs]
		except:
			nowstrs = nowstr.split(',')
			ncol = len(nowstrs)
			xs = [float(x) for x in nowstrs]
			nowchar = ','
	else:
		nowstrs = nowstr.split(nowchar)
		ncol = len(nowstrs)
		xs = [float(x) for x in nowstrs]
		
	for icol in range(ncol):
		mins[icol] = min(mins[icol], xs[icol])
		maxs[icol] = max(maxs[icol], xs[icol])
		sums[icol] += xs[icol]
		sumsqs[icol] += xs[icol]**2.0

 for icol in range(ncol):
	avgs[icol] = sums[icol] / float(iline)
	stds[icol] = (sumsqs[icol]/float(iline) - avgs[icol]**2) ** 0.5

 print 'Finishing processing', iline, 'lines of ', filename, ' / n-col =', ncol
 for icol in range(ncol):
	print ' \ti-col, min, max, avg, std = ', icol, '\t', mins[icol], '\t', maxs[icol], '\t', avgs[icol], '\t',  stds[icol]
	
