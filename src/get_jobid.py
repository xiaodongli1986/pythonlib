
nowfile='download.sh.output'
outputfile=nowfile+'.jobids'
nowf2=open(outputfile,'w')
for nowstr in open(nowfile,'r').readlines():
	if nowstr[0:6]=='Job id':
		nowstr2 = nowstr.split()[2]
		nowf2.write(nowstr2+' ')
nowf2.close()

