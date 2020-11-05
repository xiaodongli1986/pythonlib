import numpy as np
#import sys, commands
import sys, os

#cmdargs = sys.argv
#joblist = cmdargs[1]

#rlt = commands.getoutput('uws $cstr  list')
rlt = os.popen('uws $cstr  list').read()
print(rlt)

#for nowstr in open(joblist,'r').readlines():
for nowstr in rlt.split('\n'): #open(joblist,'r').readlines():
	if nowstr == '': break
	nowstrs = nowstr.split()
	if len(nowstrs) < 5: continue
	jobid = nowstrs[0]
	nowinfo = [X.split() for X in os.popen('uws $cstr job show '+jobid).read().split('\n')]
	if len(nowinfo) <5: continue
	print(' ')
	print('Check job ', jobid, '...')
	if nowinfo[5][1] == 'COMPLETED':
		print('\tJob COMPLETED. Skip.')
	elif nowinfo[5][1] == 'EXECUTING':
		print('\tJob EXECUTING. Skip.')
	else:
		nowtab = nowinfo[12][2]
		nowquery = '';
		for row in range(2,len(nowinfo[13])):
			nowquery += nowinfo[13][row]; nowquery +=' '
		print('\tResubmit job ',nowquery)
		cmd = 'uws $cstr job new queue="long" query="'+nowquery+'" table="'+nowtab+'" --run'
		print('\tDELETE JOB AND EXECUTING CMD: ',cmd)
		print(os.popen('uws $cstr job delete '+jobid).read())
		print(os.popen(cmd).read())
	#print(nowinfo)


