#!/home/xiaodongli/anaconda/bin/python
import os #commands
import sys

cmdargs = sys.argv
print(cmdargs)

nowfile=cmdargs[1]
if len(cmdargs) <3:
	outputfile = nowfile+'.show'
else:
	outputfile = cmdargs[2]
print('Load in ids:   \n\t', nowfile)
print('showing rlt:   \n\t', outputfile)
numjob = 0
allids = []
for nowstr in open(nowfile,'r').readlines():
		for nowstr2 in nowstr.split():
			allids.append(nowstr2)
			numjob += 1
print( numjob, 'jobs identified:\n\t', allids)

rlts = []

for nowid in allids:
	#nowrlt = commands.getoutput('uws $cstr job show '+str(nowid))
	nowrlt = os.popen('uws $cstr job show '+str(nowid)).read()
	#print(nowrlt

	nowstrs = nowrlt.split()
	for i in range(len(nowstrs)):
		if nowstrs[i] == 'Phase':
			rlts.append([nowid, nowstrs[i+1], nowstrs[i-4]])
			print(nowid, nowstrs[i+1], nowstrs[i-4])
			break
			
nowf2=open(outputfile,'w')
for rlt in rlts:
	nowf2.write(rlt[0]+' '+rlt[1]+' '+rlt[2]+'\n')
nowf2.close()
print('Results written to ', outputfile)



