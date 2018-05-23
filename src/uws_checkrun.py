import numpy as np
import sys, commands

#cmdargs = sys.argv
#joblist = cmdargs[1]

rlt = commands.getoutput('uws $cstr  list')

#for nowstr in open(joblist,'r').readlines():
for nowstr in rlt.split('\n'): #open(joblist,'r').readlines():
	if nowstr == '': break
	nowstrs = nowstr.split()
	if len(nowstrs) < 5: continue
	jobid = nowstrs[0]
	nowinfo = [X.split() for X in commands.getoutput('uws $cstr job show '+jobid).split('\n')]
	if len(nowinfo) <5: continue
	print 
	print 'Check job ', jobid, '...'
	if nowinfo[5][1] == 'COMPLETED':
		print '\tJob COMPLETED. Skip.'
	else:
		nowtab = nowinfo[12][2]
		nowquery = '';
		for row in range(2,len(nowinfo[13])):
			nowquery += nowinfo[13][row]; nowquery +=' '
		print '\tResubmit job ',nowquery
		cmd = 'uws $cstr job new queue="long" query="'+nowquery+'" table="'+nowtab+'" --run'
		print '\tDELETE JOB AND EXECUTING CMD: ',cmd
		print commands.getoutput('uws $cstr job delete '+jobid)
		print commands.getoutput(cmd)
	#print nowinfo


