

# This get the *.er files whose size not zero!!!
import os, sys
import numpy  as np

print('Eg: py_erfile er er2\nstandar format: py_erfile suffix1 suffix1\n\tWill find all files with *.suffix1; create a bash file to re-run all jobs having error. \n\tsuffix1 for the previous output; suffix2 for new output ')

argv = sys.argv

suffix1, suffix2 = 'er', 'er2'

if len(argv) > 1:
    suffix1 = argv[1]
if len(argv) > 2:
    suffix2 = argv[2]

file_info = os.popen('ls *.'+suffix1+' -alh ').read();


file_infos = np.array([nowstr.split() for nowstr in file_info.split('\n')])

filenames = []
for X in file_infos:
    try:
        file_size = X[4]
        if file_size != '0':
            print(X[8]); filenames.append(X[8])
            #print(X)
    except:
        pass

finalbashfile = 'jsub_'+suffix1+"_"+suffix2+'.sh'
nowf = open(finalbashfile, 'w')
nowf.write('np=1\n')
for filename in filenames:
    bashfile = filename[:-len(suffix1)-1]
    nowf.write('jsub -n $np  -o '+bashfile+'.output   -e '+bashfile+'.'+suffix2+'    sh '+bashfile+'\nsleep 5\n')
nowf.close()
print('Finish. bash created: ', finalbashfile)


