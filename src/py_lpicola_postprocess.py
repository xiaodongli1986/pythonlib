
import os, sys

printstr='''#Example: \n\tpy_lpicola_postprocess   -basename YOUR_BASE_NAME   -nbox 2   -overlap_distance 15   -xyzmin 0   -xyzmax 600   -output_1eighthLC T   -mv_lightcone_files T   -just_create_bash F   -split_np 4

# Bash Example:

#################################
# This will, 
#  1. Split the binary output of l-picola into subboxes; 
#  2. search for rockstar halos in each subsample 
#  3. convert all halo samples into format of x,y,z,vx,vy,vz,mvir,vmax
#  4. merger these subsamples of halos into one 
#################################

# Split the full LC into nbox*nbox*nbox subboxes, 
nbox=2

# In the range of xyzmin to xyzmax, (if want full-sky, set xyzmin as minus xyzmax)
xyzmin=0
xyzmax=600

# Each subbox has this amount of overlap_distance 
overlap_distance=15

# Run the job using this numer of processors
np=4

output_1eighthLC="T" # only output 1/8 of the LC (e.g., the x>0, y>0, z>0 part)
mv_lightcone_files="F" # mv the many  ${basename}_ligthcone.* files into one directory 
                ## if you want to delete them, set it as 'Delete'

for basename in ............... De_w_-0.7_AnalyticalGrowth   De_w_-1.3_AnalyticalGrowth   Omega_0.3071_AnalyticalGrowth De_w_-0.7    De_w_-1.3  Omega_0.3071 
do
	cmd="py_lpicola_postprocess  -basename $basename  -nbox $nbox -overlap_distance $overlap_distance -xyzmin $xyzmin -xyzmax $xyzmax -output_1eighthLC $output_1eighthLC   -mv_lightcone_files $mv_lightcone_files"
	#$cmd
	jsub -n $np -o pp_${basename}.output  -e pp_${basename}.er $cmd
	sleep 2
done

'''

args = sys.argv

if len(args) <=2:
    print(printstr); sys.exit()

output_1eighthLC = True
mv_lightcone_files = True
just_create_bash = False
split_np = 4

for iarg in range(1,len(args),2):
    str1, str2 = args[iarg], args[iarg+1]
    print(str1, str2)
    if str1 in ['-basename']:
        basename = str2
        print('\t set basename as ', basename)
    elif str1 in ['-nbox']:
        nbox = int(str2)
        print('\t set nbox as ', nbox)
    elif str1 in ['-overlap_distance']:
        overlap_distance = int(str2)
        print('\t set overlap_distance as ', overlap_distance)
    elif str1 in ['-xyzmin']:
        xyzmin  = int(str2)
        print('\t set xyzmin as ', xyzmin)
    elif str1 in ['-xyzmax']:
        xyzmax  = int(str2)
        print('\t set xyzmin as ', xyzmax)
    elif str1 in ['-output_1eighthLC']:
        if str2[0] in ['T', 't']:
            output_1eighthLC = True
        elif str2[0] in ['F', 'f']:
            output_1eighthLC = False
        else:
            print('ERROR! wrong arg for output_1eighthLC: (must start with T or F) ', str2)
            print(printstr); sys.exit()
        print('\t set output_1eighthLC as ', output_1eighthLC)
    elif str1 in ['-mv_lightcone_files']:
        if str2[0] in ['T', 't']:
            mv_lightcone_files= True; print('\t set mv_lightcone_files as ', mv_lightcone_files)
        elif str2[0] in ['F', 'f']:
            mv_lightcone_files= False; print('\t set mv_lightcone_files as ', mv_lightcone_files)
        elif str2[0] in ['D', 'd']:
            mv_lightcone_files= 'Delete'
            print('\t WARNING! set mv_lightcone_files as ', mv_lightcone_files)
            print('\t WARNING! set mv_lightcone_files as ', mv_lightcone_files)
            print('\t WARNING! set mv_lightcone_files as ', mv_lightcone_files)
        else:
            print('ERROR! wrong arg for mv_lightcone_files: (must start with T or F) ', str2)
            print(printstr); sys.exit()
    elif str1 in ['-just_create_bash', '-just_bash', '-just_shfile']:
        if str2[0] in ['T', 't']:
            just_create_bash = True
        elif str2[0] in ['F', 'f']:
            just_create_bash = False
    elif str1 in ['-split_np']:
        split_np = int(str2)
        print('\t set split_np as ', split_np)
    else:
        print('ERROR! unknown arg: ', str1)
        print(printstr); sys.exit()


#basename= 'De_w_-1.3_AnalyticalGrowth'
#nbox=2
#overlap_distance=15.0
#xyzmin = 0.0
#xyzmax = 200.0

filelist = basename+'_filelist'
headfile = basename+'_parameters'
outputname = 'rockstar_halos/' + basename 
cmd1 = 'ls '+basename+'_lightcone.* > '+ filelist
cmd2 = 'mkdir -p rockstar_halos'

output_suffix = '.nbox'+str(nbox)+'_overlap%.1f'%overlap_distance+'_xyz%.1f'%xyzmin+'to%.1f'%xyzmax+'.ibox'

cmd3 = 'LSS_lpicola_lightcone_boxsplit -inputfilelist '+filelist+'   -outputname '+outputname+'   -nbox %.1f'%nbox+'   -overlap_distance %.1f'%overlap_distance+'   -xyzmin %.1f'%xyzmin+'   -xyzmax %.1f'%xyzmax+'  -binary_IO T    -headfile '+headfile
cmd = ' && '.join([cmd1,cmd2,cmd3])

bashfile1 = basename+'_1_run_split.sh'
bashfile2 = basename+'_2_jsub_find_rockstar.sh'
bashfile3 = basename+'_3_run_merge_mv_conv.sh'

if just_create_bash:
    bashf = open(bashfile1, 'w'); bashf.write(cmd); bashf.close()
    

files = [basename+output_suffix+str(ibox) for ibox in range(1,nbox**3+1)]
for ifile, nowfile in enumerate(files):
    cmd += ' &&  cd rockstar_halos && py_rockstar '+nowfile+' &&  cd .. && LSS_rockstar_select_xyzvxvyvz_mvir_vmax rockstar_halos/'+nowfile+'_rockstar_halo.ascii'
    if just_create_bash:
        bashfile = basename+'_2_run'+str(ifile)+'.sh'; bashf = open(bashfile, 'w'); bashf.write('cd rockstar_halos && py_rockstar '+nowfile+' &&  cd .. && LSS_rockstar_select_xyzvxvyvz_mvir_vmax rockstar_halos/'+nowfile+'_rockstar_halo.ascii'); bashf.close()
        if ifile == 0:
            bashf = open(bashfile2, 'w'); 
        else:
            bashf = open(bashfile2, 'a'); 
        bashf.write('jsub -n '+str(split_np)+' -o '+bashfile+'.output  -e '+bashfile+'.er   sh '+bashfile+'\nsleep 2 \n'); bashf.close()

    

cmd4 = 'py_LightconeHalo_halo_for_iboxes_merge '+str(nbox)+' %.1f'%overlap_distance+' %.1f'%xyzmin+' %.1f'%xyzmax+' rockstar_halos '+basename+' '+str(output_1eighthLC)+' T '
if just_create_bash:
    bashf = open(bashfile3,'w'); bashf.write(cmd4+'\n'); bashf.close()

cmd += (' && '+cmd4)

if mv_lightcone_files == True:
    cmd5 = ' mkdir -p  '+basename+'_lightcone_files && mv '+basename+'_lightcone.* '+basename+'_lightcone_files'
    cmd += (' && '+cmd5)
    if just_create_bash:
        bashf = open(bashfile3,'a'); bashf.write(cmd5+'\n'); bashf.close()
elif mv_lightcone_files == "Delete": 
    cmd5 = ' /usr/bin/rm -rf  '+basename+'_lightcone.* '
    cmd += (' && '+cmd5)
    if just_create_bash:
        bashf = open(bashfile3,'a'); bashf.write(cmd5+'\n'); bashf.close()

file1 = 'rockstar_halos/'+basename+output_suffix+'ALL_rockstar_halo.ascii.xyzvxvyvz_mvir_vmax'
file2 = file1+'_fmt3'
cmd6 = 'py_binary_ascii_conv '+file1+' '+file2+' ab 3 -fmt3_headfile '+headfile
cmd += (' && '+cmd6)
if just_create_bash:
    bashf = open(bashfile3,'a'); bashf.write(cmd6+'\n'); bashf.close()

if not just_create_bash:
    print('Executing cmd:\n\t', cmd)
    print(os.popen(cmd).read())


