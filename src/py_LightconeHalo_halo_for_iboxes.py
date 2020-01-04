
import os, numpy, sys

print 'Usage: EXE filename num-of-iboxes quickstart-cfg-file'

#filename = 'test_lightcone.nbox4_overlap10.0_xyz-1800.0to1800.0.ibox'

args = sys.argv

if len(args) < 4:
    sys.exit()

filename, nbox, cfgfile = args[1:]

#nbox = 64

nbox = int(nbox)

nowf0 =open('halo_jsub.sh', 'w')
for ibox in range(1,nbox+1):
    os.popen('mkdir ibox'+str(ibox))
    os.popen('mv '+filename+str(ibox)+'.* ibox'+str(ibox))
    nowshfile = 'halo_ibox'+str(ibox)+'.sh'
    nowf = open(nowshfile, 'w')
    nowcmd = 'cd ibox'+str(ibox)+'; /home/xiaodongli/software/Rockstar/rockstar -c '+cfgfile+' '+filename+str(ibox)+'.*'
    nowf.write(nowcmd)
    nowf0.write('sleep 2 \n')
    nowf0.write('jsub -n 1 -o '+nowshfile+'.out -e '+nowshfile+'.er sh '+nowshfile+'\n')
    nowf.close()
nowf0.close()

