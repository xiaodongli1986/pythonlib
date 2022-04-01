
#import commands
import numpy as np
import sys, os
import stdA_py3


printstr = '''test!test!test!
        Postprocess a Snapshot sample
        1. Convert background to another cosmology
        2. Do a mass-cut (to reach certain number density)
        3. Measure numberdensity (using SPH, using given numNB)
        
        Assuming z-direction is the Line-of-Sight
            i.e. rescale  z  according to Hz(Cosmo_out)/Hz(Cosmo_in), rescale  x, y   according to DA(Cosmo_out)/DA(Cosmo_in)

#################################
####  Usage: 
            py_SnapCosmoConv_MassCut_MeasureRho   -inputfile inputfile    -omin omin -win win -omout omout -wout wout    -boxsize -redshift redshift   -nbar nbar    -numNBs numNBs  -margins margins   -rm_tmpfiles F -rmcmd /us/bin/rm  -xcol xcol -ycol ycol -zcol zcol -masscol masscol  -create_random T -random_rat 10  (col start with 1)

#################################
##### Example: 
            py_SnapCosmoConv_MassCut_MeasureRho    -inputfile om0.3471_omb0.0482_oml0.6529_w-1.0000_Hubble0.68_sig80.8228_ns0.9600_simsnap01810r_rockstar_halo_xyzvxvyvz_mvir_vmax.shiftz    -omin 0.3471 -win -1 -omout 0.26 -wout -1    -boxsize 680    -numNBs 30_100 -margins 30_100     -xcol 1 -ycol 2 -zcol 3 -masscol 7     -redshift  0.6 -nbar 3e-4     -create_random T -random_rat 10   -rm_tmpfiles T -rmcmd /usr/bin/rm

#################################
#### Example (a bash file):

import os
bashfile = '1_mksample.sh'

### options

files = os.popen('ls *shiftz').read()

snap_red_dict = {'r': 0.6, 's': 0.5, 't': 0.4,  }

omout, wout = 0.26, -1
numNBs = '30_100'
margins = '30_100'
nbar = '3e-4'
xcol, ycol, zcol, masscol = 1, 2, 3, 7
boxsize = '680'

bashf = open(bashfile, 'w')
for nowfile in files.split():
    paramstrs = nowfile.split('_')
    ## get omegam, w
    omstr, wstr = paramstrs[0], paramstrs[3]
    om, w = float(omstr.strip('om')), float(wstr.strip('w'))
    ## get redshift
    for paramstr in paramstrs:
        if 'simsnap' in paramstr:
            snapstr = paramstr[-1]
    if snapstr in snap_red_dict.keys():
        redshift = snap_red_dict[snapstr]

        cmd = 'py_SnapCosmoConv_MassCut_MeasureRho    -inputfile  '+nowfile+'   -omin %.4f'%(om)+' -win %.4f'%w+' -omout %.4f'%omout+' -wout %.4f'%wout+'    -boxsize '+str(boxsize)+'    -numNBs '+str(numNBs)+' -margins '+str(margins)+'    -xcol '+str(xcol)+' -ycol '+str(ycol)+' -zcol '+str(zcol)+' -masscol '+str(masscol)+'     -redshift  '+str(redshift)+' -nbar '+str(nbar)+'     -create_random T -random_rat 10   -rm_tmpfiles T -rmcmd /usr/bin/rm'
        print(cmd)
        bashf.write('echo \necho '+cmd+'\n')
        bashf.write(cmd+'\necho \n\n')
bashf.close()
print('commands written to ', bashfile)


        '''

cmdargs = sys.argv

if len(cmdargs) < 2:
    print (printstr); sys.exit()


omin = 0.3071
win = -1
omout = omin; wout=win
boxsize = None
nbar = None
numNBs = []
margins = []
redshift = None

rm_tmpfiles = False
rmcmd = '/usr/bin/rm'

xcol, ycol, zcol = 1, 2, 3
masscol  = None

create_random = True
random_rat = 10

for iarg in range(len(cmdargs)):
        if np.mod(iarg, 2) == 0:
            continue
        else:
            opt1, opt2 = cmdargs[iarg:iarg+2]
            #print('opt1, opt2 = ', opt1, opt2)
            if opt1 == '-inputfile': 
                inputfile = opt2
            elif opt1 == '-omin': 
                omin = float(opt2)
            elif opt1 == '-omout':
                omout = float(opt2)
            elif opt1 == '-win': 
                win = float(opt2)
            elif opt1 == '-wout':
                wout = float(opt2)
            elif opt1 == '-boxsize':
                boxsize = float(opt2)
            elif opt1 == '-nbar':
                nbar = opt2
            elif opt1 == '-redshift':
                redshift = float(opt2)
            elif opt1 == '-numNBs':
                numNBs = [int(xx) for xx in opt2.split('_')]
            elif opt1 == '-margins':
                margins = [int(xx) for xx in opt2.split('_')]
            elif opt1 == '-xcol':
                xcol = int(opt2)
            elif opt1 == '-ycol':
                ycol = int(opt2)
            elif opt1 == '-zcol':
                zcol = int(opt2)
            elif opt1 == '-masscol':
                masscol = int(opt2)
            elif opt1 == '-rmcmd':
                rmcmd = opt2
            elif opt1 == '-rm_tmpfiles':
                if opt2[0]  in 'tT':
                    rm_tmpfiles = True
                elif opt2[0]  in 'fF':
                    rm_tmpfiles = False
                else:
                    print('wrong -rm_tmpfiles: (must start with t, T, f, F!) ', opt2); sys.exit()
            elif opt1 == '-create_random':
                if opt2[0]  in 'tT':
                    create_random = True
                elif opt2[0]  in 'fF':
                    create_random = False
                else:
                    print('wrong -create_random: (must start with t, T, f, F!) ', opt2); sys.exit()
            elif opt1 == '-random_rat':
                random_rat = float(opt2)

print(' (py_Makerandom_MeasureRho) Start.')
print('  Read in parameters:')
print('\tinputfile = ', inputfile)
print('\tomin,  win  = ', omin, win)
print('\tomout, wout = ', omout, wout)
print('\tboxsize  = ', boxsize)
print('\tredshift = ', redshift)
print('\tnumNBs  = ', numNBs)
print('\tmargins = ', margins)
print('\txcol, ycol, zcol, masscol = ', xcol, ycol, zcol, masscol)
print('\trm_tmpfiles, rmcmd= ', rm_tmpfiles, rmcmd)
print('\tcreate_random, random_rat= ', create_random, random_rat)

option_list = [inputfile, omin, win, omout, wout, boxsize, redshift, numNBs, margins, xcol, ycol, zcol, masscol]
if None in option_list:
    print(' (py_Makerandom_MeasureRho) Error! One option is None!')
    sys.exit()

tmpfiles = []


DAin, DAout = stdA_py3.DA(omin, win, 0.7, redshift), stdA_py3.DA(omout, wout, 0.7, redshift)
Hin, Hout = stdA_py3.Hz(omin, win, 0.7, redshift), stdA_py3.Hz(omout, wout, 0.7, redshift)
DArat, Hrat = DAout/DAin, Hout/Hin

Volrat = DArat**2 / Hrat

### Step 1. numdegrade

numdegrade = int(boxsize**3 * Volrat * float(nbar))
'''
file_numdegrade = inputfile.strip('.tmpfile') + '.numdegrade_nbar'+nbar+'.tmpfile'; tmpfiles.append(file_numdegrade)
print('\n ----------------------------------------------------------------------------------------------------------------')
print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 1) numdegrade.\n\t DArat, Hrat, Volrat = ', DArat, Hrat, Volrat )
print('\t Do a mass-cut to reach nbar = ',nbar,'. target #-object = ', numdegrade)
print('    input       = ', inputfile)
print('    output      = ', file_numdegrade)
cmdstr = 'LSS_rmass_stat    -input '+inputfile+ ' -rmin 0 -rmax 1e15 -numrbin 1 -nummassbin 500 -xcol '+str(xcol)+' -ycol '+str(ycol)+' -zcol '+str(zcol)+' -masscol '+str(masscol)+' -dodegrade T -degradedfile  '+file_numdegrade+' -numdegrade '+str(numdegrade)+' -const_nbar_degrade T'
print(' (Executing command...)  \n\t', cmdstr)
os.popen('sleep 1')
print(os.popen(cmdstr).read())
print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 2) numdegrade: Done.')

### Step 2.  Cosmo Conv

print('\n ----------------------------------------------------------------------------------------------------------------')
if abs(DArat - 1.0) < 1e-6 and abs(Hrat - 1) < 1e-6:
    file_cosmoconv = file_numdegrade
    print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 2) cosmo conv: Skip since DArat, Hrat close to 1')
else:
    file_cosmoconv = file_numdegrade+'.CosmoConv.omout%.4f'%omout+'wout%.4f'%wout+'.tmpfile'; tmpfiles.append(file_cosmoconv)
    print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 2) cosmo conv: DArat, Hrat = ', DArat, Hrat )
    print('    input   = ', file_numdegrade)
    print('    output  = ', file_cosmoconv)
    f1, f2 = open(file_numdegrade, 'r'), open(file_cosmoconv, 'w')
    iline = 0
    while True:
        nowstr = f1.readline()
        if nowstr == '': break
        nowstrs = nowstr.split()
        x, y, z, mass = nowstrs[xcol-1], nowstrs[ycol-1], nowstrs[zcol-1], nowstrs[masscol-1]
        x = str(float(x)*DArat);  y = str(float(y)*DArat); z = str(float(z)/Hrat)
        f2.write(' '.join([x,y,z,mass])+'\n')
        iline += 1
    f1.close(); f2.close()
    print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 2) cosmo conv: Done. ',iline,' lines processed.')

### Step 3.  Measure density

print('\n ----------------------------------------------------------------------------------------------------------------')
print(' (py_SnapCosmoConv_Masscut_MeasureRho) (Step 3) Measure rho: numNBs, margins = ', numNBs, margins)

outputfiles = []

for numNB, margin in zip(numNBs, margins):
    file_rho_drho = inputfile +'.SnapCosmoConv.omout%.4f'%omout+'wout%.4f'%wout+'.nbar'+nbar+'.numNB'+str(numNB)+'.margin'+str(margin)
    outputfiles.append(file_rho_drho)
    xsize = ysize = boxsize * DArat; zsize = boxsize / Hrat
    print('    input   = ', file_cosmoconv)
    print('    output  = ', file_rho_drho)
    print('    x, y, z size = ', xsize, ysize, zsize)
    cmdstr = 'LSS_rho_knn_rectangular_box -inputfile '+str(file_cosmoconv)+' -xcol 1 -ycol 2 -zcol 3 -margin '+str(margin)+' -numNB '+str(margin)+' -outputfile '+file_rho_drho+' -xsize '+str(xsize)+' -ysize '+str(ysize)+' -zsize '+str(zsize)
    print(' (Executing command...)  \n\t', cmdstr)
    os.popen('sleep 1')
    print(os.popen(cmdstr).read())
'''

if create_random:
    file_random = inputfile +'.SnapCosmoConv.omout%.4f'%omout+'wout%.4f'%wout+'.nbar'+nbar+'.x'+str(random_rat)+'ran'
    nran = int(numdegrade * float(random_rat))
    f_ran = open(file_random, 'w')
    print(' (py_Makerandom_MeasureRho) Create random file: nran = ', nran)
    print('    output = ', file_random)
    for row in range(nran  ):
        x = np.random.uniform(0, xsize)
        y = np.random.uniform(0, ysize)
        z = np.random.uniform(0, zsize)
        f_ran.write(' '.join([str(x) for x in [x, y, z, 1]])+'\n')
    f_ran.close()
    for numNB, margin in zip(numNBs, margins):
      file_rho_drho = file_random + '.numNB'+str(numNB*random_rat)+'.margin'+str(margin)
      outputfiles.append(file_rho_drho)
      xsize = ysize = boxsize * DArat; zsize = boxsize / Hrat
      print('    input   = ', file_random)
      print('    output  = ', file_rho_drho)
      print('    x, y, z size = ', xsize, ysize, zsize)
      cmdstr = 'LSS_rho_knn_rectangular_box -inputfile '+str(file_random)+' -xcol 1 -ycol 2 -zcol 3 -margin '+str(margin)+' -numNB '+str(int(numNB*random_rat))+' -outputfile '+file_rho_drho+' -xsize '+str(xsize)+' -ysize '+str(ysize)+' -zsize '+str(zsize)
      print(' (Executing command...)  \n\t', cmdstr)
      os.popen('sleep 1')
      print(os.popen(cmdstr).read())



if rm_tmpfiles:
    print('\t Clearing up tmp files... ')
    for nowfile in tmpfiles:
        print('\t', rmcmd,  nowfile)
        os.popen(rmcmd+' '+nowfile).read()
print(' (py_Makerandom_MeasureRho) (Step 3) Measure rho: Done. Outputs = ')
for nowfile in outputfiles:
    print('\t', nowfile)

