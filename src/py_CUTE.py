
import numpy as np
import sys
import commands
import os

print 'This is python args version of CUTE.'
#global RRfile_create = True


allkeys = [   'data_filename', 'random_filename', 'data_filename_2', 'random_filename_2', 'input_format', 'output_filename', 'corr_type', 'omega_M', 'omega_L', 'w', 'log_bin', 'dim1_max', 'dim1_min_logbin', 'dim1_nbin', 'dim2_max', 'dim2_nbin', 'dim3_min', 'dim3_max', 'dim3_nbin', 'radial_aperture', 'use_pm', 'n_pix_sph', 'RR_filename', 'xplus', 'yplus', 'zplus', 'weight_pow', 'readbinary_fmt', 'addrsd_shiftz_pbbox', 'addrsd_redshift', 'addrsd_pbboxsize', 'addrsd_om']

output_dict = {
    'data_filename': 'test/shell.dat',
    'random_filename': 'test/random.dat',
    #The next 2 files are only needed if you want to
    #cross-correlate two datasets
    'data_filename_2': 'test/shell.dat',
    'random_filename_2': 'test/random.dat',
    'input_format': '0',
    'output_filename': 'test/corr_full_pm.dat',
    # estimation parameters
    'corr_type': 'angular',
    # cosmological parameters
    'omega_M': '0.315',
    'omega_L': '0.685',
    'w': '-1',
    # binning
    'log_bin': '0',
    'dim1_max': '10.',
    'dim1_min_logbin': '0.01',
    'dim1_nbin': '30',
    'dim2_max': '0.1',
    'dim2_nbin': '30',
    'dim3_min': '0.4',
    'dim3_max': '0.7',
    'dim3_nbin': '1',
    # pixels for radial correlation
    'radial_aperture': '1.',
    # pm parameters
    'use_pm': '0',
    'n_pix_sph': '2048',
    'RR_filename':  None, 
    'xplus': '0',
    'yplus': '0',
    'zplus': '0',
    'weight_pow': '1',
    'readbinary_fmt': '0',
    'addrsd_shiftz_pbbox': '0',
    'addrsd_redshift': '-1',
    'addrsd_pbboxsize': '-1',
    'addrsd_om': '-1',
 #   'have_weight': '1',
            }

printstr = 'Usage:\n\tpy_CUTE -cute_exe /home/xiaodongli/software/CUTE/CUTE/CUTE -cute_ini_filename ./tmp_cute_ini ...\nDefault values of optional options:\n\t'

example_str = '''bashf = open(???.sh, 'w') # bash script of py_CUTE commands
bashf.write('export OMP_NUM_THREADS=48\\n\\n') # running CUTE in 48-cores 
for ...:

        datafile = ... 
        ranfile = .... 

        wei=1;  ### weight -> weight ^wei ; convenient for mark weighted CF
        zplus=0 #zplus = 1000000000
        input_format = 4 # 4 means directly read-in x,y,z;

        ### binary format sample
        #  0 means a usual ascii file; 
        #  1 means CUTE subsample binary file (weight be forece to be 1);  
        #  2 means xyzw binary (head and end contains a 4-bytes integer specifying n-sample)
        readbinary_fmt = 0 

        ### automatically add rsd to the z-direction
        # Use them only if readbinary_fmt != 0 and velocites provided in the format !!!
        addrsd_shiftz_pbbox = 0  ### 0: don't do it; 1: add it.
        addrsd_pbboxsize = -1  ### boxsize; in case of readbinary_fmt==CUTE-subsample, you do not need to give it (set as -1)
        addrsd_redshift = -1  ### redshift of the sample; in case of readbinary_fmt==CUTE-subsample, you do not need to give it (set as -1)
        addrsd_om = -1 ### omegam of the sample; in case of readbinary_fmt==CUTE-subsample, you do not need to give it (set as -1)

        smax = 150; sbin = 150; mubin = 120
        suffixstr_rr = '.'+str(sbin)+'s0to'+str(smax)+'.'+str(mubin)+'mu'
        suffixstr = '.weipow'+str(wei)+suffixstr_rr
        if addrsd_shiftz_pbbox != 0:
            suffixstr = '.shiftz'+suffixstr
        if zplus!= 0: suffixstr = '.zplus'+str(zplus)+suffixstr

        inifile = datafile+suffixstr+'.ini'
        Tpcffile = datafile+suffixstr+'.2pcf'
        rrfile = ranfile+suffixstr_rr+'.rr'

        print('datafile, ranfile = \\n\\t', datafile, '\\n\\t',ranfile)
        print('We will generate: inifile, 2pcffile, rrfile: \\n\\t',inifile,'\\n\\t',Tpcffile,'\\n\\t',rrfile)
        print('Start running CUTE...')

        py_CUTE_cmd = 'py_CUTE    -cute_exe /home/xiaodongli/software/CUTE/CUTE/CUTE    -cute_ini_filename '+inifile+'    -corr_type 3D_rm  -input_format '+str(input_format)+'    -log_bin 0 -dim1_max '+str(smax)+'   -dim1_nbin '+str(sbin)+' -dim2_max 1   -dim2_nbin '+str(mubin)+'  -omega_M 0.3071  -omega_L 0.6929 -w -1    -data_filename '+str(datafile)+'   -random_filename '+str(ranfile)+' -output_filename '+str(Tpcffile)+'   -RR_filename '+str(rrfile)+' -weight_pow '+str(wei)+' -zplus '+str(zplus)+' -readbinary_fmt '+str(readbinary_fmt)+' -addrsd_shiftz_pbbox '+str(addrsd_shiftz_pbbox)+' -addrsd_pbboxsize '+str(addrsd_pbboxsize)+' -addrsd_redshift '+str(addrsd_redshift)+' -addrsd_om '+str(addrsd_om)
        print(py_CUTE_cmd)

        bashf.write('echo \\'datafile, ranfile = '+str(datafile)+' '+str(ranfile)+'\\'\\n')
        bashf.write('echo \\'Will generate:\\'\\n')
        bashf.write('echo \\'  inifile:   '+str(inifile)+'\\'\\n')
        bashf.write('echo \\'  2pcffile:  '+str(Tpcffile)+'\\'\\n')
        bashf.write('echo \\'  rrfile:    '+str(rrfile)+'\\'\\n')
        bashf.write('echo \\'Start running CUTE...\\'\\n')
        bashf.write(py_CUTE_cmd+'\\n')
'''

for nowkey in allkeys:
    printstr += ('-'+nowkey+' '+str(output_dict[nowkey])+'   ')

def cute_ini(ini_file_name = None, **kws):
    global RRfile, RRfile_create
    RRfile_create = True
    for nowkey in output_dict.keys():
        output_dict[nowkey] = None
    output_dict['dim3_nbin'] = 1; output_dict['dim3_min'] = 0.4; output_dict['dim3_max'] = 0.7
    if ini_file_name == None:
        ini_file_name = './tmp_cute_ini.ini'
    for key in kws.keys():
        if key not in output_dict.keys():
            print '(cute_ini) ERROR! unkwon key!: key = ', key
            return  'Nothing!!!!'
        else:
            output_dict[key] = kws[key]
    RRfile = output_dict['RR_filename']
    if RRfile != None:
        if not os.path.exists(RRfile):
            print '(cute_ini) Can not find given RRfile... will be created.'
            #if RRfile_create:
            RRfile_create = True
            output_dict['RR_filename'] = None
    print 'create ini file: ', ini_file_name, '...\n'
    nowf = open(ini_file_name, 'w')
    for key in allkeys:
        if output_dict[key] != None:
            nowstr = str(key)+'= '+str(output_dict[key])
            nowf.write(nowstr+'\n')
            print nowstr
    nowf.close()
    return ini_file_name

cute_ini_filename = None
cute_exe = '/home/xiaodongli/software/CUTE/CUTE/CUTE'
arg_dict = {}
cmdargs = sys.argv

if (len(cmdargs) -1)%2 != 0 or len(cmdargs) < 2:
    print 'ERROR! Number of options + values must be a even number: we get len(cmdargs)-1 = ', len(cmdargs)-1
    print printstr; print '\n###################\nexample of calling py_CUTE and write commands (to bash) in python:\n\n\n', example_str; sys.exit()

for iarg in range(1, len(cmdargs), 2):
    key = cmdargs[iarg]
    key = key[1:len(key)]
    value = cmdargs[iarg+1]
    if key == 'cute_ini_filename':
        cute_ini_filename = value
    elif key == 'cute_exe':
        cute_exe = value
    else:
        arg_dict[key] = value

cute_ini_filename = cute_ini(cute_ini_filename, **arg_dict)
nowcmd = cute_exe + ' ' + cute_ini_filename
print '\nRun command:\n\t'+nowcmd+'\n'
print commands.getoutput(cute_exe + ' ' + cute_ini_filename)
print output_dict
if RRfile_create:
    print 'Create RRfile...'
    nowcmd = 'cp '+output_dict['output_filename'] + ' ' + RRfile
    print '\t'+nowcmd
    print commands.getoutput(nowcmd)
