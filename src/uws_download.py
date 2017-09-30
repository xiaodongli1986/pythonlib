#!/home/xiaodongli/anaconda/bin/python
import commands
import sys
import os
import time


cmdargs = sys.argv
print cmdargs

nowfile=cmdargs[1]
if len(cmdargs) <3:
	outputfile = nowfile+'.downloadinfo'
else:
	outputfile = cmdargs[2]
print 'Load in ids:   \n\t', nowfile
print 'showing rlt:   \n\t', outputfile

os.system('rm $outputfile')

numjob = 0
completed_ids = []
for nowstr in open(nowfile,'r').readlines():
	nowstrs = nowstr.split()
	nowid, nowstatus = nowstrs[0], nowstrs[1]
	if nowstatus == 'COMPLETED':
		completed_ids.append(nowid)
	else:
		print 'Not completed id: ', nowid, nowstatus
print  len(completed_ids), 'completed jobs:\n\t', completed_ids


nowf2=open(outputfile,'w')
time0 = time.time()
numid = 0
for nowid in completed_ids:
	numid += 1
	time1 = time.time()
	print 'Downloding ', nowid, ':       ',numid, 'of', len(completed_ids)
	nowf2.write('Downloding '+str(nowid)+'...\n')
	os.system('uws $cstr job results '+str(nowid)+' csv >> '+outputfile)
	time2 = time.time()
	print '   Time consumed:   this_job/all_jobs = ', time2-time1, time2-time0
	nowf2.write('   Time consumed:   this_job/all_jobs = '+str(time2-time1)+'/'+str(time2-time0))
nowf2.close()
#nowf2.close()
#print 'Results written to ', outputfile



