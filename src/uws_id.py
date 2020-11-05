#!/home/xiaodongli/anaconda/bin/python
import commands
import sys

cmdargs = sys.argv
print(cmdargs)

nowfile=cmdargs[1]
if len(cmdargs) <3:
	outputfile = nowfile+'.ids'
else:
	outputfile = cmdargs[2]
print('Load in result of job submission:   \n\t', nowfile)
print('Extract job ids and output to file: \n\t', outputfile)
numjob = 0
nowf2=open(outputfile,'w')
for nowstr in open(nowfile,'r').readlines():
	if nowstr[0:6]=='Job id':
		nowstr2 = nowstr.split()[2]
		nowf2.write(nowstr2+' ')
		numjob += 1
nowf2.close()
print( numjob, 'job ids. Done.')

