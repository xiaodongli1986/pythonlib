#!/usr/bin/env python
#import commands
import sys
import os
import time

usagestr = 'Usage: EXE inputfile [options...]\n'+\
	'   Input file format (line by line): jobid, job_status (will be donwloaded if completed), outputfile (optional) \n'+\
	'    [-output/-outputfile  name of output file]\n'+\
	'    [-skip_existed/-skipexisted  True/true/T/t/False/false/F/f] skip existed file and do not download it (requring "outputfile" in the input file)\n'+\
	'    [-parallel  True/true/T/t/False/false/F/f] downloading files in a parallel manner\n'+\
	'    [-timeinterval True/true/T/t/False/false/F/f] when using parallel, gap time among start time of jobs\n'+\
	''

cmdargs = sys.argv
if len(cmdargs) < 2:
	print('Error! Must have at list one arg!')
	print(usagestr)
	sys.exit()
print(cmdargs)

nowfile=cmdargs[1]

outputfile = nowfile+'.downloadinfo'
skipexisted = True
sleepstr = str(100)
parallel = True
if len(cmdargs) > 3:
	i1 = 2
	while i1 <= len(cmdargs) -1:
		i2 = i1+1
		if i2 == len(cmdargs):
			print('Error! overflowing: nothing after ', cmdargs[i1], ', arg -', i1)
			print(usagestr)
			sys.exit
		if cmdargs[i1] in ['-outputfile', '-output']:
			outputfile = cmdargs[i2]
		elif cmdargs[i1] in ['-timeinterval']:
			sleepstr = cmdargs[i2]
		elif cmdargs[i1] in ['-skip_existed', '-skipexisted']:
			nowstr = cmdargs[i2]
			if nowstr[0] in ['t', 'T']:
				skipexisted = True
			elif nowstr[0] in ['f', 'F']:
				skipexisted = False		
			else:
				print('Wrong skipexisted: must start with t/T/f/F, we get ', nowstr)
				sys.exit()
		elif cmdargs[i1] in [ '-parallel']:
			nowstr = cmdargs[i2]
			if nowstr[0] in ['t', 'T']:
				parallel = True
			elif nowstr[0] in ['f', 'F']:
				parallel = False		
			else:
				print('Wrong parallel: must start with t/T/f/F, we get ', nowstr)
				sys.exit()
		else:
			print('Wrong arg!', cmdargs[i1])
			print(usagestr)
		i1 += 2


#if len(cmdargs) <3:
#	outputfile = nowfile+'.downloadinfo'
#else:
#	outputfile = cmdargs[2]
print('Load in ids:   \n\t', nowfile)
print('showing rlt:   \n\t', outputfile)

os.system('rm $outputfile')

numjob = 0
dl_id_files = []
for nowstr in open(nowfile,'r').readlines():
	nowstrs = nowstr.split()
	nowid, nowstatus, nowfile = nowstrs[0], nowstrs[1], nowstrs[2]+'.csv'
	if skipexisted and os.path.isfile(nowfile):
		print('Skip existed file: ', nowfile)
		continue
	if nowstatus != 'COMPLETED':
		print('Skip not completed job: ', nowid, nowstatus)
		continue
	dl_id_files.append([nowid,nowfile])
print(  len(dl_id_files), 'jobs for download:')
for nowidfile in  dl_id_files:
	print('\t', nowidfile[0], '   ', nowidfile[1])


nowf2=open(outputfile,'w')
time0 = time.time()
numid = 0
for nowid_file in dl_id_files:
	numid += 1
	nowid, nowfile = nowid_file
	time1 = time.time()
	if skipexisted and os.path.isfile(nowfile):
		print( 'Skip existed file: ', nowid, nowfile)
		nowf2.write('Skip existed file: '+str(nowid)+' '+nowfile+'\n')
		continue
	print( 'Downloding ', nowid, '   ', nowfile, ':       ',numid, 'of', len(dl_id_files))
	nowf2.write('Downloding '+str(nowid)+'...\n')
	if parallel:
		os.system('uws $cstr job results '+str(nowid)+' csv >> '+outputfile+' &')
		print( 'Start next job after ', sleepstr, 'seconds...')
		os.system('sleep '+sleepstr)
	else:
		os.system('uws $cstr job results '+str(nowid)+' csv >> '+outputfile)
	time2 = time.time()
	print( '   Time consumed:   this_job/all_jobs = ', time2-time1, time2-time0)
	nowf2.write('   Time consumed:   this_job/all_jobs = '+str(time2-time1)+'/'+str(time2-time0))
nowf2.close()
#nowf2.close()
#print 'Results written to ', outputfile



