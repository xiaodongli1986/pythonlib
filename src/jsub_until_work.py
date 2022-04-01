
import os, sys

printstr = '''Usage: jsub_until_work -sleep sleep_time -J job_name -n #-thread -o output -e error your_command
    # keep jsub a job until it suceessfully submitted 
    # sleep gives time sleep (before it make the next trial)
jsub_until_work *** is equivalaent to 
    jsub ***
except that:
    1) keep trails of jsub untial successfully submitted (jobid foundin jjobs)
    2) you can use -sleep sleep_time to specify interval between two trials
'''

args = sys.argv

#print(args)

sleeptime = 1

if len(args) <= 2:
    print(printstr)
    sys.exit()

flag = False
for row in range(len(args)):
    if args[row] == '-J':
        jobname = args[row+1]
        flag = True; 
    if args[row] == '-sleep':
        sleeptime = args[row+1]



cmd = 'jsub '

icmd = 1
while icmd < len(args):
    if args[icmd] == '-sleep':
        icmd += 2; continue
    else:
        cmd = cmd+' '+args[icmd]
    icmd += 1

print('===================================================================================')
print(' (jsub_until_work) Will execute \n\t', cmd)

while True:
    jobinfo = os.popen(cmd).read()
    jobid = jobinfo.split()[1].strip('<>')
    print('   ', jobinfo)

    print(os.popen('sleep '+str(sleeptime)).read())

    jjobs = os.popen('jjobs').read()
    flag = False
    for jjob in jjobs.split('\n'):
        if jjob == '': continue
        jobstrs = jjob.split()
        print('     look for ', jobid, 'in', jobstrs )
        if jobstrs[0] == jobid:
            flag = True
            break
    if flag: 
        print(' (jsub_until_work) check and found job existed: ', jjob, ' Done!')
        break
    else:
        print(' (jsub_unitl_work) check and not found the job, resubmit after ', sleeptime, 'seconds!')



