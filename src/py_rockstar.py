

import os, sys
import numpy as np

printstr ='''Usage: py_rockstar filename 
    -exe/-EXE exe_path 
    -c/-config/-config_file config_file 
    -o/-output/-outputfile outputfile

Example:
    py_rockstar snp04441e.\*

Bash example:
    for ((i=0; i<=8; i++))
    do
	nowfile=Omega_0.3071_binary.nbox2_overlap15.0_xyz0.0to150.0.ibox$i
	echo $nowfile
	py_rockstar $nowfile
        #jsub -n 1 -o ${nowfile}.output -e ${nowfile}.er py_rockstar $nowfile
        #sleep 2
    done
    '''

if len(sys.argv) <=1:
    print(printstr)
    sys.exit()

inputfile = sys.argv[1]
print('\t set inputfile as: ', inputfile)

exe_path = '/home/xiaodongli/software/Rockstar/rockstar '
config_file = '/home/xiaodongli/software/Rockstar/quickstart.cfg' 
outputfile = None

if len(sys.argv) > 2:
    for iarg in range(2,len(sys.argv),2):
        str1, str2 = sys.argv[iarg:iarg+2]
        if str1 in ['-exe', '-EXE']:
            exe_path = str2
            print('\t set exe as: ', exe_path)
        elif str1 in ['-c', '-config', '-config_file']:
            config_file = str2
            print('\t set config_file as: ', config_file)
        elif str1 in ['-o', '-output', '-outputfile']:
            outputfile = str2
            print('\t set outputfile as: ', outputfile)
        else:
            print('Unknown option: ', str1)
            print(printstr)
            sys.exit()

if outputfile == None:
    print(inputfile)
    outputfile = inputfile.replace(".*","").replace(".?","").replace("\\",'').replace("*","").replace("?","")+'_rockstar_halo'
    print('\t automatically set outputfile as: ', outputfile)

iran = np.random.uniform(10000000,99999999); ranpath = 'rs_ranpath_'+str(iran)

cmd1 = 'mkdir -p '+ranpath 
cmd2 = 'cd '+ranpath
cmd3 = 'ln -s ../'+inputfile+' ./'
cmd4 = exe_path+' -c  '+config_file +' ' + inputfile
cmd5 = 'mv halos_0.0.ascii ../'+outputfile+'.ascii' 
cmd6 = 'mv halos_0.0.bin ../'+outputfile+'.bin'
cmd7 = 'mv rockstar.cfg ../'+outputfile+'.cfg'
cmd8 = 'cd ..'
cmd9 = '/usr/bin/rm -rf '+ranpath

cmd = ' && '.join([cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,cmd9])


print('Now execute the command:\n\t\t', cmd, '...')
print(os.popen(cmd).read())

#print('rename the files...\n\t\t', cmd5, '\n\t\t', cmd6); 

print('Convert file to x,y,z,vx,vy,vz,mvir,vmax format...')
cmd = 'LSS_rockstar_select_xyzvxvyvz_mvir_vmax '+outputfile+'.ascii'
print('\t\t', cmd); print(os.popen(cmd).read())
