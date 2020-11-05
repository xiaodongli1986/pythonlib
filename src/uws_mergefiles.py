
import os #commands
import numpy as np
import sys
import stdA

cmdargs = sys.argv

outputfile = cmdargs[1]

filenames = []
for i in range(2,len(cmdargs)):
	#nowfilenames = commands.getoutput('ls '+cmdargs[i]).split()
	nowfilenames = os.getoutput('ls '+cmdargs[i]).read().split()
	filenames = filenames + nowfilenames

print('Merge ', len(filenames) ,'csv files:')
for nowfile in filenames:
	print('\t', nowfile)
print('Will create new file:\n\t', outputfile)

outputf = open(outputfile, 'w')

totlines = 0
ifile = 0
for nowfile in filenames:
	print('\tOpening ', nowfile, '...')
	nowf = open(nowfile,'r')
	nowlines = 0
	nowstr = nowf.readline()
	if ifile == 0: outputf.write(nowstr)
	while True:
		nowstr = nowf.readline()
		if nowstr == '': break
		outputf.write(nowstr)
		totlines += 1
		nowlines += 1
	print('\t\t', nowlines, 'lines read & write.')
	nowf.close()
	ifile += 1
print('\t', totlines, 'lines read & write in total. Finish.')
outputf.close()
