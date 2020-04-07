
import numpy as np
import sys, os

Buffer = 1.3
NumFilesWrittenInParallel = None
Omega		 = 0.31        #% Total matter density (CDM + Baryons at z=0).
OmegaBaryon	 = 0.048       # % Baryon density (at z=0).
OmegaLambda	 = 0.69        #% Dark Energy density (at z=0)
HubbleParam	 = 0.69        # % Hubble parameter, 'little h' (only used for power spectrum parameterization).
Sigma8		 = 0.83        #% Power spectrum normalization (power spectrum may already be normalized correctly).
PrimordialIndex	 = 0.96        #% Used to tilt the power spectrum for non-tabulated power spectra (if != 0.4.0 and nongaussian, generic flag required)
De_w = -1

nproc = 0

FileWithInputSpectrum  = '/home/xiaodongli/software/l-picola/files/input_spectrum.dat'
EXE_path = '/home/xiaodongli/software/l-picola-wcdm/L-PICOLA'
outputdir= None
basename = None

args = sys.argv

#-outputdir outputdir   -basename basename   -omegam/-om omegam   -w w   -omegab omegab   -omegal omegal   -H0 H0   -sigma8 sigma8   -ns ns   -Pkfile Pkfile    -np np   - start_redshfit start_redshift  -end_redshift end_redshift   -nstep1 nstep1   -nstep2 nstep2   -boxsize boxsize   -nmesh nmesh   -nparticle nparticle   -zinit zinit   -buffer buffer   -seed seed

printstr = '''# Example:
   py_lpicola  -np 48   -nparticle 256   -nmesh 256   -boxsize 512   -seed 5000   -buffer 2.0   -EXE_path "/home/xiaodongli/software/l-picola-wcdm/L-PICOLA"      -NumFilesWrittenInParallel 0   -zinit 9.0   -nstep1 30   -start_redshift 0.05   -end_redshift 0.0   -nstep2 30  -omegam 0.3000   -w -1.0000   -omegab 0.0480   -omegal 0.7000   -h 69.00   -sigma8 0.8300   -ns 0.9600   -basename None   -outputdir None   -Pkfile "/home/xiaodongli/software/l-picola/files/input_spectrum.dat"   
        NumFilesWrittenInParallel: if set as 0, then will be set to np 
        baename/outputdir: if set as None, automatic generating a name based on parameters'''



if len(args) <=1 or len(args)%2 != 1:
    print(printstr)
    sys.exit()

for iarg in range(1,len(args),2):
    arg1, arg2 = args[iarg], args[iarg+1]
    #print(arg1, arg2)
    if arg1 in ['-outputdir']:
        outputdir = arg2
        print('\tOutpudir = ', outputdir)
    elif arg1 in ['-basename']:
        basename = arg2
        print('\t set basename (filename)= ', basename)
    elif arg1 in ['-EXE_path']:
        EXE_path = arg2
        print('\t set EXE_path = ', EXE_path)
    elif arg1 in ['-np']:
        nproc = arg2
        print('\t set np (num-of-processes) as ',nproc)
    elif arg1 in ['-zinit', '-init_redshift']:
        Init_redshift = arg2
        print('\t set zinit as ',Init_redshift)
    elif arg1 in ['-start_redshift']:
        start_redshift = arg2
        print('\t set start_redshift as ',start_redshift)
    elif arg1 in ['-end_redshift']:
        end_redshift = arg2
        print('\t set end_redshift as ',end_redshift)
    elif arg1 in ['-nstep1']:
        nstep1 = arg2
        print('\t set nstep1 as ', nstep1) 
    elif arg1 in ['-nstep2']:
        nstep2 = arg2
        print('\t set nstep2 as ', nstep2) 
    elif arg1 in ['-boxsize']:
        boxsize = arg2
        print('\t set boxsize as ', boxsize) 
    elif arg1 in ['-nparticle']:
        nparticle = arg2
        print('\t set nparticle as ', nparticle) 
    elif arg1 in ['-NumFilesWrittenInParallel' ]:
        NumFilesWrittenInParallel = arg2
        print('\t set NumFilesWrittenInParallel as ',NumFilesWrittenInParallel)
    elif arg1 in ['-nmesh']:
        nmesh = arg2
        print('\t set nmesh as ', nmesh) 
    elif arg1 in ['-seed']:
        seed = arg2
        print('\t set random-seed as ', seed) 
    elif arg1 in ['-Omega', '-omegam', '-om']:
        Omega = arg2
        print('\t set omegam as ', Omega)
    elif arg1 in ['-omb', '-OmegaBaryon', '-omegab']:
        OmegaBaryon = arg2
        print('\t set omegab as ', OmegaBaryon)
    elif arg1 in ['-OmegaLambda', '-oml', '-omegalabmda', '-omegade', '-omegal']:
        OmegaLambda = arg2
        print('\t set omega_de as ', OmegaLambda)
    elif arg1 in ['-h', '-Hubble', '-HubbleParam']:
        HubbleParam = str(float(arg2))
        print('\t set h as ', HubbleParam)
    elif arg1 in ['-sigma8', '-sig8', '-Sigma8']:
        Sigma8 = arg2
        print('\t set sigma8 as ', Sigma8)
    elif arg1 in ['-ns', '-PrimordialIndex']:
        PrimordialIndex = arg2
        print('\t set ns as ', PrimordialIndex)
    elif arg1 in ['-w', '-De_w', '-de_w']:
        De_w = arg2
        print('\t set w as ', De_w)
    elif arg1 in ['-Pkfile', '-pkfile', '-']:
        FileWithInputSpectrum  = arg2
        print('\t set Pkfile as ', FileWithInputSpectrum)
    elif arg1 in ['-buffer', '-Buffer']:
        Buffer = arg2
        print('\t set Buffer as ', Buffer)
    else:
        print('Unknown arg1 :', arg1)
        print(printstr)
        sys.exit()

#print('NumFilesWrittenInParallel = ')
#print(NumFilesWrittenInParallel)

if NumFilesWrittenInParallel in [None, 0, '0']:
    NumFilesWrittenInParallel = nproc
    print('\t reset NumFilesWrittenInParalle as ', NumFilesWrittenInParallel)

#sys.exit()
runname=''.join([ boxsize,'box_', nmesh, 'mesh_', nparticle, 'par_z', start_redshift, '_startstep', nstep1, '_', end_redshift, '_endstep', nstep2])

if outputdir in [None, 'None']:
    outputdir = runname
else:
    outputdir = outputdir # +'__'+runname

os.popen('mkdir -p '+outputdir)
os.popen('sleep 1')

cosmostr = 'om%.4f'%(float(Omega))+'_omb%.4f'%(float(OmegaBaryon))+'_oml%.4f'%(float(OmegaLambda))+'_w%.4f'%(float(De_w))+'_Hubble%.2f'%(float(HubbleParam))+'_sig8%.4f'%(float(Sigma8))+'_ns%.4f'%(float(PrimordialIndex))

if basename in [None, 'None']:
    basename = cosmostr
else:
    basename = basename #+'__'+cosmostr

## redshift file
modelname = outputdir+'/'+basename
redshiftfile = modelname+'_output_redshifts.dat'

setstr = '''% =============================== %
% This is the run parameters file % 
% =============================== %

% Simulation outputs
% ==================
OutputDir      OutputDir_OutputDir             %files/                              % Directory for output.
FileBase                    test                     % Base-filename of output files (appropriate additions are appended on at runtime) 
OutputRedshiftFile          Outred_Outred     % The file containing the redshifts that we want snapshots for
NumFilesWrittenInParallel   64                                    % limits the number of files that are written in parallel when outputting.

% Simulation Specifications
% =========================
UseCOLA		1           % Whether or not to use the COLA method (1=true, 0=false).
Buffer		1.3         % The amount of extra memory to reserve for particles moving between tasks during runtime.
Nmesh		Nmesh_Nmesh         % This is the size of the FFT grid used to compute the displacement field and gravitational forces.

Nsample		Nsample_Nsample         % This sets the total number of particles in the simulation, such that Ntot = Nsample^3.
Box		Box_Box       % The Periodic box size of simulation.
Init_Redshift   9.0         % The redshift to begin timestepping from (redshift = 9 works well for COLA)
Seed		5001        % Seed for IC-generator
SphereMode	1 %0           % If "1" only modes with |k| < k_Nyquist are used to generate initial conditions (i.e. a sphere in k-space), 
                             % otherwise modes with |k_x|,|k_y|,|k_z| < k_Nyquist are used (i.e. a cube in k-space).

WhichSpectrum	1           % "0" - Use transfer function, not power spectrum
                             % "1" - Use a tabulated power spectrum in the file 'FileWithInputSpectrum'
                             % otherwise, Eisenstein and Hu (1998) parametrization is used
                             % Non-Gaussian case requires "0" and that we use the transfer function

WhichTransfer	0           % "0" - Use power spectrum, not transfer function
                             % "1" - Use a tabulated transfer function in the file 'FileWithInputTransfer' 
                             % otherwise, Eisenstein and Hu (1998) parameterization used 
                             % For Non-Gaussian models this is required (rather than the power spectrum) 

% FileWithInputSpectrum  files/input_spectrum.dat    % filename of tabulated input spectrum (if used)
FileWithInputSpectrum  /home/xiaodongli/software/l-picola/files/input_spectrum.dat % files/input_spectrum.dat    % filename of tabulated input spectrum (if used)
                                                   % expecting k and Pk 

FileWithInputTransfer  files/input_transfer.dat    % filename of tabulated transfer function (if used)
                                                   % expecting k and T (unnormalized)

% Cosmological Parameters
% =======================
Omega		0.31        % Total matter density (CDM + Baryons at z=0).
OmegaBaryon	0.048        % Baryon density (at z=0).
OmegaLambda	0.69        % Dark Energy density (at z=0)
HubbleParam	0.69         % Hubble parameter, 'little h' (only used for power spectrum parameterization).
Sigma8		0.83        % Power spectrum normalization (power spectrum may already be normalized correctly).
PrimordialIndex	0.96        % Used to tilt the power spectrum for non-tabulated power spectra (if != 0.4.0 and nongaussian, generic flag required)
De_w            -1          % Dark energy Equation-of-state // added by xiaodong

% Timestepping Options
% ====================
StepDist	0           % The timestep spacing (0 for linear in a, 0.4 for logarithmic in a)
DeltaA		0           % The type of timestepping: "0" - Use modified COLA timestepping for Kick and Drift. Please choose a value for nLPT.
                             % The type of timestepping: "1" - Use modified COLA timestepping for Kick and standard Quinn timestepping for Drift. Please choose a value for nLPT.
                             % The type of timestepping: "2" - Use standard Quinn timestepping for Kick and Drift
                             % The type of timestepping: "3" - Use non-integral timestepping for Kick and Drift
nLPT		-2.5        % The value of nLPT to use for modified COLA timestepping


% Units
% =====                                                                                             
UnitLength_in_cm		3.085678e24       % defines length unit of output (in cm/h) 
UnitMass_in_g			1.989e43          % defines mass unit of output (in g/h)
UnitVelocity_in_cm_per_s	1e5               % defines velocity unit of output (in cm/sec)
InputSpectrum_UnitLength_in_cm	3.085678e24       % defines length unit of tabulated input spectrum in cm/h. 
                                                  % Note: This can be chosen different from UnitLength_in_cm

% ================================================== %
% Optional extras (must comment out if not required) %
% ================================================== %

% Non-Gaussianity
% ===============
%Fnl			-400                              % The value of Fnl.
%Fnl_Redshift		49.0                              % The redshift to apply the nongaussian potential
%FileWithInputKernel  /files/input_kernel_ortog.txt     % the input kernel for generic non-gaussianity (only needed for GENERIC_FNL)


% Lightcone simulations
% =====================
Origin_x	0.0                % The position of the lightcone origin in the x-axis
Origin_y	0.0                % The position of the lightcone origin in the y-axis
Origin_z	0.0                % The position of the lightcone origin in the z-axis
Nrep_neg_x	100000 %0                  % The maximum number of box replicates in the negative x-direction
Nrep_pos_x	100000 %0                  % The maximum number of box replicates in the positive x-direction
Nrep_neg_y	100000 %0                  % The maximum number of box replicates in the negative y-direction
Nrep_pos_y	100000 %0                  % The maximum number of box replicates in the positive y-direction
Nrep_neg_z	100000 %0                  % The maximum number of box replicates in the negative z-direction
Nrep_pos_z	100000 %0                  % The maximum number of box replicates in the positive z-direction'''

#print(nowstr)


setstrs = setstr.split('\n')
for istr, nowstr in enumerate(setstrs):
    if nowstr == '': continue
    nowstrs = nowstr.split()
    if nowstrs[0] == 'Omega':
        setstrs[istr] = nowstrs[0] + '     ' + Omega
    elif nowstrs[0] == 'OmegaBaryon':
        setstrs[istr] = nowstrs[0] + '     ' + OmegaBaryon
    elif nowstrs[0] == 'OmegaLambda':
        setstrs[istr] = nowstrs[0] + '     ' + OmegaLambda
    elif nowstrs[0] == 'HubbleParam':
        setstrs[istr] = nowstrs[0] + '     ' + HubbleParam
    elif nowstrs[0] == 'Sigma8':
        setstrs[istr] = nowstrs[0] + '     ' + Sigma8
    elif nowstrs[0] == 'PrimordialIndex':
        setstrs[istr] = nowstrs[0] + '     ' + PrimordialIndex
    elif nowstrs[0] == 'De_w':
        setstrs[istr] = nowstrs[0] + '     ' + De_w
    elif nowstrs[0] == 'FileWithInputSpectrum':
        setstrs[istr] = nowstrs[0] + '     ' + FileWithInputSpectrum
    elif nowstrs[0] == 'Box':
        setstrs[istr] = nowstrs[0] + '     ' + boxsize
    elif nowstrs[0] == 'Nmesh':
        setstrs[istr] = nowstrs[0] + '     ' + nmesh
    elif nowstrs[0] == 'Nsample':
        setstrs[istr] = nowstrs[0] + '     ' + nparticle
    elif nowstrs[0] == 'Init_redshift':
        setstrs[istr] = nowstrs[0] + '     ' + Init_redshift
    elif nowstrs[0] == 'Seed':
        setstrs[istr] = nowstrs[0] + '     ' + seed
    elif nowstrs[0] == 'Buffer':
        setstrs[istr] = nowstrs[0] + '     ' + Buffer
    elif nowstrs[0] == 'Box':
        setstrs[istr] = nowstrs[0] + '     ' + boxsize
    elif nowstrs[0] == 'OutputDir':
        setstrs[istr] = nowstrs[0] + '     ' + outputdir
    elif nowstrs[0] == 'FileBase':
        setstrs[istr] = nowstrs[0] + '     ' + basename
    elif nowstrs[0] == 'OutputRedshiftFile':
        setstrs[istr] = nowstrs[0] + '     ' + redshiftfile
    elif nowstrs[0] == 'NumFilesWrittenInParallel':
        setstrs[istr] = nowstrs[0] + '     ' + str(NumFilesWrittenInParallel)
setstr = '\n'.join(setstrs)
setstr = '\n'.join(setstrs)

# parameter file
nowfile = outputdir+'_'+basename+'.dat'
nowf = open(nowfile,'w') 
nowf.write(setstr)
nowf.close()

# bash file
bashfile =  'run_'+outputdir+'_'+basename+'.sh'
nowf = open(bashfile, 'w')
## libarary
nowf.write('export LD_LIBRARY_PATH=/home/xiaodongli/software/gsl2.5/lib/:$LD_LIBRARY_PATH\n\n' )
## cmd
runcmd = 'mpirun -np '+str(nproc)+' '+EXE_path+' '+nowfile
print('# generate command:\n\t',runcmd)
nowf.write(runcmd+'\n')
nowf.close()

# jsub bash file
jsubcmd = 'jsub -n '+str(nproc)+' -o '+(bashfile+'.output')+' -e '+(bashfile+'.er')+' sh '+bashfile
print('\t',jsubcmd)
jsubfile = 'jsub_'+outputdir+'_'+basename+'.sh'
nowf = open(jsubfile,'w')
nowf.write(jsubcmd+'\n')
nowf.close()

### particle mass
G = 6.67428*(10**-8) #cm3*g-1*s-2
M_sun = 1.98892*(10**33) #g
Mpc = 30.85677*(10**23) #cm
h = float(HubbleParam); Omega_m = float(Omega); boxsize=float(boxsize); nparticle=float(nparticle)

rho = (3*((100*h)**2)/(8*np.pi*G)) * ((10**10)*Mpc/M_sun) #M_sun*h2/Mpc-3
rho_h = rho/(h**2)
m_p1 = (((boxsize**3)*3*Omega_m*((100)**2))/(nparticle**3*8*np.pi*G)) * ((10**10) *Mpc /M_sun) #M_sun/h
m_p2 = (rho_h*Omega_m*(boxsize**3))/(nparticle**3)
#print ('m_p=',m_p2)

# summarizing the values of parameters
parfile = outputdir+'/'+basename+'_parameters'
nowf = open(parfile,'w')
nowf.write('# nran, size, parmass, redshift, omegam, h, w\n')
nowf.write(' '.join([str(xx) for xx in [0, boxsize, m_p2, -1.0, Omega, HubbleParam, De_w]]))
nowf.close()


print('\n# parameterfile: \n\t', nowfile)
print('# command to run it\n\t',bashfile)
print('# command to submit the job\n\t',jsubfile)
print('# values of parameters\n\t',parfile)


nowf=open(redshiftfile, 'w')
nowf.write('''
% Here we place a list of comma-delimited output redshifts and the number of steps to take between
% that redshift and the previous output redshift (or the initial redshift specified in the run
% parameters file if a given output redshift is the first one to be output)
%
% We can't put an output redshift greater than the initial redshift, but we can ask for output
% at the initial redshift, i.e. 2LPT initial conditions only, by putting one of the output redshifts
% as the initial redshift. In fact, as comparing floating point numbers is hard, for any output redshift
% within 0.0001 percent of the initial redshift we will assume you meant to put the initial redshift
% or that rounding errors have occured and output at the initial redshift anyway.
%
% A final point is that all the number of steps must be greater than 0 except for if the
% corresponding redshift satisfies the condition above and is equivalent to the initial redshift.
% In this case it doesn't matter what value you put.
%
% Example: output a simulation at the initial redshift of 9.0 and then take 6 steps to redshift 1.0, output,
% and take 4 steps to redshift 0.0 then output:
%   9.0, 0
%   1.0, 6
%   0.0, 4
'''+str(start_redshift)+', '+str(nstep1)+'\n'+str(end_redshift)+', '+str(nstep2)+'\n')
nowf.close()
print('\b# create redshfit file:\n\t',redshiftfile)


