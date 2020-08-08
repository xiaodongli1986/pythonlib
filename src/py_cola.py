
import numpy as np
import sys, os

default_keys = {
        'random_seed': 4000, 

        'nc': 256,               'boxsize': 512, 
        'nrealization': 1,       # multiple realisations for random_seed, random_seed+1, ...

        'ntimestep': 30,        'a_final': 1,     
        'output_redshifts': '{5.5, 5.0, 4.5, 4.0, 3.5, 3.0, 2.75, 2.5, 2.25, 2.0, 1.75, 1.5, 1.25,  1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0}' ,

        'omega_m': 0.3071,      'h': 0.6787,     'sigma8': 0.8228, 
        'omega_l': -1.,         'de_w_add10': 9.,

        #'De_w': None,

        # pm_nc_factor= 3            -- Particle Mesh grid pm_nc_factor*nc per dimension
        # np_alloc_factor= 5.        -- Amount of memory allocated for particle
        # loglevel=2                 -- 0=debug, 1=verbose, 2=normal, increase the value to reduce output msgs
        #'pm_nc_factor': 3,      'np_alloc_factor':  1.5,    
        'pm_nc_factor': 5,      'np_alloc_factor':  18,    
        'loglevel':  2,
        'powerspectrum':        '"BigMDPL_matterpower.dat"', 

        #'fof':   '"fof"',   'linking_factor' : 0.2,
        'snapshot': '"snap"',
        'subsample': '"subsample"', 'subsample_factor' : 0.01,

        #'coarse_grid' : None, 'coarse_grid_nc' : 32,
        'init': None, 
        
        ## for lightcone
        }

append_param_dict = {
        'omegab': 0.048206,   'ns': 0.96  
        }

basename  = None; outputdir=None
EXE_path = "/home/xiaodongli/software/cola_halo/cola_alo"
nproc = 48
de_w = -1.
simple_basename = True

args = sys.argv

printstr = '''# Example:
   py_cola   -np 4   -basename None   -nc 256   -boxsize 256   -random_seed 4000   -ntimestep 20   -powerspectrum "BigMDPL_matterpower.dat"   -EXE_path "/home/xiaodongli/software/cola_halo/cola_halo"    -omegam 0.3071   -h 0.6787    -omegal 0.6929  -sigma8 0.8228    -omegab 0.048206    -ns 0.96  -de_w -1.0        -pm_nc_factor 5   -np_alloc_factor 1.5  # nc means nmesh 
   ## !!! Still need to add: De_w, zmax, basename, omgeal!!!  use_solve_growth!  # firstly, add them to cola_halo!
   #  recommend to set use_solve_growth=1 (if de_w != -1)
        baename: if set as None, automatic generating a name based on parameters'''


if len(args) <=1 or len(args)%2 != 1:
    print(printstr)
    sys.exit()

for iarg in range(1,len(args),2):
    arg1, arg2 = args[iarg], args[iarg+1]; 
    if arg1[0] != '-':
        print('Unknown arg1 :', arg1); print(printstr); sys.exit()
    arg1=arg1[1:]
    #print(arg1, arg2)
    if arg1 == 'EXE_path':
        EXE_path= arg2; print('\t EXE_path= ', EXE_path)
    elif arg1 == 'np':
        nproc = arg2; print('\tnproc= ', nproc)
    elif arg1 == 'basename':
        basename = arg2; print('\t set basename (filename)= ', basename)
    elif arg1 in ['omegab', 'ns', ]:
        append_param_dict[arg1] = arg2;  print('\t set append_param_dict[', arg1, '] = ', arg2)
    elif arg1 in ['omegam', 'Omegam', 'omega_m']:
        default_keys['omega_m'] = arg2;  print('\t set default_keys[ omega_m ] = ', arg2)
    elif arg1 in ['omegal', 'Omegal', 'omega_l']:
        default_keys['omega_l'] = arg2;  print('\t set default_keys[ omega_l ] = ', arg2)
    elif arg1 in ['w', 'de_w']:
        default_keys['de_w_add10'] = float(arg2) + 10.; print('\t set default_keys[ de_w_add10 ] = ', default_keys['de_w_add10'])
    elif arg1 in ['simple_basename', 'simplebasename']:
        if arg2[0] in ['t', 'T']:
            simple_basename = True
        elif arg2[0] in ['f', 'F']:
            simple_basename = False
        else:
            print('ERROR! value for simple_basename must be logic, we get ', arg2); sys.exit()
    else:
        if arg1 in default_keys.keys():
            default_keys[arg1] = arg2;  print('\t set default_keys[', arg1, '] = ', arg2)
        else: 
            print('Unknown arg1 :', arg1); print(printstr); sys.exit()

default_keys['powerspectrum'] = '"'+default_keys['powerspectrum'] + '"'
boxsize=default_keys['boxsize'];  nc=default_keys['nc'];  ntimestep=default_keys['ntimestep'];  random_seed=default_keys['random_seed']
runname=''.join([ str(boxsize),'box_', str(nc), 'npar_', str(default_keys['pm_nc_factor']), 'pmfactor_',   str(ntimestep), 'steps_ranseed', str(random_seed), ])

# ----------------------------------------
#  * set names

if outputdir in [None, 'None']:
    outputdir = runname
else:
    outputdir = outputdir # +'__'+runname

omegam=default_keys['omega_m']; sigma8=default_keys['sigma8']; h=default_keys['h']; w=float(default_keys['de_w_add10']) - 10.
omegab=append_param_dict['omegab']; omegal=default_keys['omega_l']; ns=append_param_dict['ns'] 
cosmostr = 'om%.4f'%(float(omegam))+'_omb%.4f'%(float(omegab))+'_oml%.4f'%(float(omegal))+'_w%.4f'%(float(w))+'_Hubble%.2f'%(float(h))+'_sig8%.4f'%(float(sigma8))+'_ns%.4f'%(float(ns))


if basename in [None, 'None']:
    if not simple_basename:
        basename = cosmostr
    else:
        outputdir = runname + '/' + cosmostr
        basename = 'sim'
else:
    basename = basename #+'__'+cosmostr

os.popen('mkdir -p '+outputdir); print(os.popen('sleep 2').read())

modelname = outputdir+'/'+basename
default_keys['lightcone_basename'] = '"'+modelname+'"'


# ----------------------------------------
#  * parameter file
nowfile = 'param_'+runname+'_'+cosmostr+'.lua'

print(default_keys)
nowf = open(nowfile,'w')
for nowkey in default_keys.keys():
    if default_keys[nowkey] not in ['None', None]:
        nowf.write(str(nowkey) + '=  '+str(default_keys[nowkey])+'\n')
nowf.close()


# ----------------------------------------
# * bash file
bashfile =  'run_'+runname+'_'+cosmostr+'.sh'
nowf = open(bashfile, 'w')
## libarary
nowf.write('export LD_LIBRARY_PATH=/home/xiaodongli/software/gsl2.5/lib/:$LD_LIBRARY_PATH\n\n' )
nowf.write('export OMP_NUM_THREADS=1\n\n')
## cmd
runcmd = 'mpirun -np '+str(nproc)+' '+EXE_path+' '+nowfile
print('# generate command:\n\t',runcmd)
nowf.write(runcmd+'\n')
nowf.close()


# ----------------------------------------
# * jsub bash file
jsubcmd = 'jsub -n '+str(nproc)+' -o '+(bashfile+'.output')+' -e '+(bashfile+'.er')+' sh '+bashfile
print('\t',jsubcmd)
jsubfile = 'jsub_'+runname+'_'+cosmostr+'.sh'
nowf = open(jsubfile,'w')
nowf.write(jsubcmd+'\n')
nowf.close()


# ----------------------------------------
# * particle mass
G = 6.67428*(10**-8) #cm3*g-1*s-2
M_sun = 1.98892*(10**33) #g
Mpc = 30.85677*(10**23) #cm
h = float(h); omegam = Omega_m = float(omegam); boxsize=float(boxsize); nparticle=float(nc)

rho = (3*((100*h)**2)/(8*np.pi*G)) * ((10**10)*Mpc/M_sun) #M_sun*h2/Mpc-3
rho_h = rho/(h**2)
m_p1 = (((boxsize**3)*3*omegam*((100)**2))/(nparticle**3*8*np.pi*G)) * ((10**10) *Mpc /M_sun) #M_sun/h
m_p2 = (rho_h*omegam*(boxsize**3))/(nparticle**3)
#print ('m_p=',m_p2)

# summarizing the values of parameters
parfile = outputdir+'/'+basename+'_parameters'
nowf = open(parfile,'w')
nowf.write('# nran, size, parmass, redshift, omegam, h, w\n')
nowf.write(' '.join([str(xx) for xx in [0, boxsize, m_p2, -1.0, omegam, h, w]]))
nowf.close()


print('\n# parameterfile: \n\t', nowfile)
print('# command to run it\n\t',bashfile)
print('# command to submit the job\n\t',jsubfile)
print('# values of parameters\n\t',parfile)




