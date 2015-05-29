
import commands
pythonlibPATH=commands.getoutput('echo $pythonlibPATH')
execfile(pythonlibPATH+'/stdA.py')
pyfile=pythonlibPATH+'/'+'bossdatamock.py'

#########################
## Information

# catname means DR12-CMASS-N, DR12-CMASS-S, ...
# catname2 means HR3, L93, B08, ...

catnameliststr = 'DR12v4-CMASS-N DR12v4-CMASS-S DR12v4-LOWZ-N DR12v4-LOWZ-S'
catnamelist = ['DR12v4-CMASS-N', 'DR12v4-CMASS-S', 'DR12v4-LOWZ-N', 'DR12v4-LOWZ-S']

catinfo_npatch = {
			'DR12v4-CMASS-N' : 4,
			'DR12v4-CMASS-S' : 8,
			'DR12v4-LOWZ-N' : 4,
			'DR12v4-LOWZ-S' : 8
			}

catinfo_redshiftrange = {
			'DR12v4-CMASS-N' : [0.43, 0.6929],
			'DR12v4-CMASS-S' : [0.43, 0.6929],
			'DR12v4-LOWZ-N' : [0.15, 0.43],
			'DR12v4-LOWZ-S' : [0.15, 0.43]	
			}

catinfo_rrange = {
			'DR12v4-CMASS-N' : [1172.74, 1772.34],
			'DR12v4-CMASS-S' : [1172.74, 1772.34],
			'DR12v4-LOWZ-N' : [436.07, 1172.74],
			'DR12v4-LOWZ-S' : [436.07, 1172.74]	
			}

catinfo_3binsplit = {
			'DR12v4-CMASS-N' : [1172.74,    1367.47,    1507.22,    1772.34],
			'DR12v4-CMASS-S' : [1172.76,    1365.31,    1510.98,    1772.34],
			'DR12v4-LOWZ-N' : [436.07,     775.25,     975.38,    1172.74],
			'DR12v4-LOWZ-S' : [436.07,     767.21,     972.80,    1172.73]	
			}

catinfo_4binsplit = {
			'DR12v4-CMASS-N':   [1172.74,    1333.84,    1434.71,    1550.67,    1772.34],
			'DR12v4-CMASS-S':   [1172.76,    1330.04,    1437.35,    1552.80,    1772.34],
			'DR12v4-LOWZ-N':    [436.07,     702.40,     892.16,    1017.54,    1172.74],
			'DR12v4-LOWZ-S':    [436.07,     696.79,     884.27,    1016.24,    1172.73],
			}

catinfo_5binsplit = {
			'DR12v4-CMASS-N':    [1172.74,    1312.32,    1394.34,    1477.22,    1579.24,    1772.34],
			'DR12v4-CMASS-S':   [1172.76,    1308.79,    1393.88,    1480.64,    1580.34,    1772.34],

			'DR12v4-LOWZ-N':   [436.07,     657.42,     827.88,     943.81,    1044.16,    1172.74],
			'DR12v4-LOWZ-S':   [436.07,     653.53,     820.32,     939.13,    1042.28,    1172.73]
			}

HR4catname2list = ['LC93', 'J08.dat.z_0', 'J08.dat.z_0.5']
catname2list = ['HR3'] + HR4catname2list

miniscan_omwlist = [[0.26, -1.0], [1.0, -1.0], [0.0, -1.0], [0.26, -3.0], [0.26, -0.2] ]
miniscan_scanname = 'miniscan'

#########################################
### Files (name, check)

execfile(pythonlibPATH+'/bossdatamock_files.py')

######################################
#### commands

### These are filenames of 2pcf but we think it is more related with Sec. commands rather than files, so put them in this section

def Tpcf_suffix(rmax=150, nbins=150, mubins=120):
	return '.rmax'+str(rmax)+'.'+str(nbins)+'rbins.'+str(mubins)+'mubins'

def Tpcfrltfilename(galfile, rmax=150, nbins=150, mubins=120):
	return galfile+Tpcf_suffix(rmax, nbins, mubins)+'.2pcf'

def Tpcfrrfilename(ranfile, rmax=150, nbins=150, mubins=120):
	return ranfile+Tpcf_suffix(rmax, nbins, mubins)+'.rr'

def Tpcf_create_shfile(galfile, ranfile, rmax=150, nbins=150, mubins=120, ncpu=160, 
	shfilename='', jobname='', printinfo=True,
	exedir='/home/xiaodongli/software/csabiu-kstat-8b4343c51687/bin/2pcf',
	exit_if_filenotexit=True,
	checkexit=True,
	maxlenjobname=30):
	if checkexit or exit_if_filenotexit:
		if not isfile(galfile):
			print '\n (Tpcf_creat_shfile) WARNING! galfile not found: ', galfile, '\n'
			if exit_if_filenotexit:
				sys.exit()
		if not isfile(ranfile):
			print ' (Tpcf_creat_shfile) WARNING! ranfile not found: ', ranfile, '\n'
			if exit_if_filenotexit:
				sys.exit()
	cmd = cmd_Tpcf(galfile, ranfile, rmax, nbins, mubins)
	path, filename = separate_path_file(galfile)
	if jobname == '':
		jobname = filename+Tpcf_suffix(rmax, nbins, mubins)
		if len(jobname) > maxlenjobname:
			jobname = jobname[0:maxlenjobname]
	if shfilename == '':
		shfilename = galfile + Tpcf_suffix(rmax, nbins, mubins) + '.run_sge_baekdu.sh'
	cmdstr = '\n\n## Tpcf job file created by Tpcf_create_shfile()\n\n# shfile dir: '+shfilename+'\n# galfile : '+galfile+'\n# ranfile : '+ranfile+'\n# rmax, nbins, mubins = '+str((rmax, nbins, mubins))+'\n# exedir : '+exedir+'\n\n'+cmd_Tpcf(galfile, ranfile, rmax, nbins, mubins, exedir)
	shfile_baekdu(shfilename, jobname, ncpu, cmdstr, printinfo = printinfo )
	#if printinfo:
	#	print ' (Tpcf_create_shfile) shfile created: ', shfilename	
	return shfilename

def Tpcf_create_shfiles(galfilelist, ranfilelist, rmax=150, nbins=150, mubins=120, ncpu=160, 
	shfilename='', jobname='', printinfo=True,
	exedir='/home/xiaodongli/software/csabiu-kstat-8b4343c51687/bin/2pcf',
	exit_if_filenotexit=True,
	checkexit=True,
	maxlenjobname=1000):
	cmdstr = ''
	for ifile in range(len(galfilelist)):
		galfile, ranfile = galfilelist[ifile], ranfilelist[ifile]
		if checkexit or exit_if_filenotexit:
			if not isfile(galfile):
				print '\n (Tpcf_creat_shfile) WARNING! galfile not found: ', galfile, '\n'
				if exit_if_filenotexit:
					sys.exit()
			if not isfile(ranfile):
				print '\n (Tpcf_creat_shfile) WARNING! ranfile not found: ', ranfile, '\n'
				if exit_if_filenotexit:
					sys.exit()
		cmd = cmd_Tpcf(galfile, ranfile, rmax, nbins, mubins)
		path, filename = separate_path_file(galfile)
		if ifile == 0:
			if jobname == '':
				jobname = filename+Tpcf_suffix(rmax, nbins, mubins)+'.'+str(len(galfilelist))+'jobs'
				if len(jobname) > maxlenjobname:
					jobname = jobname[0:maxlenjobname]
			if shfilename == '':
				shfilename = galfile + Tpcf_suffix(rmax, nbins, mubins) +'.'+str(len(galfilelist))+'jobs'+ '.run_sge_baekdu.sh'
		cmdstr += '\n\n## Tpcf job file created by Tpcf_create_shfile()\n\n# shfile dir: '+shfilename+'\n# galfile : '+galfile+'\n# ranfile : '+ranfile+'\n# rmax, nbins, mubins = '+str((rmax, nbins, mubins))+'\n# exedir : '+exedir+'\n\n'+cmd_Tpcf(galfile, ranfile, rmax, nbins, mubins, exedir)
	shfile_baekdu(shfilename, jobname, ncpu, cmdstr, printinfo = printinfo )
	#if printinfo:
	#	print ' (Tpcf_create_shfile) shfile created: ', shfilename	
	return shfilename


def cmd_Tpcf(galfile, ranfile, rmax, nbins, mubins, exedir='/home/xiaodongli/software/csabiu-kstat-8b4343c51687/bin/2pcf'):
	cmd='mpirun '+exedir+' -gal '+galfile+' -ran '+ranfile+' -rmax '+str(rmax)+' -nbins '+str(nbins)+' -tbins '+str(mubins)+' -out '+Tpcfrltfilename(galfile, rmax, nbins, mubins)+' -iso ANISO -decp SMU -RR '+Tpcfrrfilename(ranfile, rmax, nbins, mubins)+' -wgt .true.'
	return cmd

def cmd_filesplit(filename, redges):
    cmd = 'LSS_Rbinned-Split-Data '+filename+' '+str(len(redges)-1)+' '+fmtstrlist(redges, '%.2f', ' ')
    for i in range(len(redges)-1):
	cmd += ' T '
    filelist = []
    for nowi in range(len(redges)-1):
        r1, r2 = redges[nowi:nowi+2]
        filelist.append(filename+'.r%.2f'%r1+'to%.2f'%r2)
    return cmd, filelist

def boss2pcf_creatbashfile_split(catname, catinfo_binsplit, givenfilelist=''):
        numbinsplit = len(catinfo_binsplit[catname])-1

        shfileA = datamockdir+'/Step4.'+catname+'.finesplit.sh'
        shfileB = datamockdir+'/Step4.'+catname+'.'+str(numbinsplit)+'binsplit.sh'

        fA = open(shfileA, 'w'); fA.write('#!/bin/bash\n\n');
        fB = open(shfileB, 'w'); fB.write('#!/bin/bash\n\n');
        rrange = catinfo_rrange[catname]
        splitinfo = catinfo_binsplit[catname]
        print 'rrange, splitinfo = ', rrange, splitinfo
	if givenfilelist == '':
		givenfilelist = allfilelist_xyzw(catname)
        for nowfile in givenfilelist:

            ## cmd, filelist
            cmdA, filelistA = cmd_filesplit(nowfile, rrange)
            cmdB, filelistB = cmd_filesplit(nowfile, splitinfo)

            ## finesplit
            fileA = filelistA[0]; fileAmv = nowfile+'.finesplitted'
            fA.write(cmdA+'\n')
            fA.write('mv '+fileA+' '+fileAmv+'\n')
            path, filename = separate_path_file(fileAmv)
            fA.write('mv '+path+'*.finesplitted '+path+'../xyzw.finesplitted/ \n')
            fA.write('mv '+path+'*.rbinned-split-data.info '+path+'../xyzw.finesplitted/ \n\n')

            ### bin split

            fB.write(cmdB+'\n')
            for ifile in range(len(filelistB)):
                path, fileB = separate_path_file(filelistB[ifile])
                fileBmv = nowfile+'.'+str(ifile+1)+'of'+str(len(splitinfo)-1)
                fB.write('mv '+filelistB[ifile]+' '+fileBmv+' \n')
            fB.write('mv '+path+'*.rbinned-split-data.info '+path+'../xyzw.binsplitted/ \n')
            fB.write('mv '+path+'*.*of* '+path+'../xyzw.binsplitted/ \n\n')
        fA.close(); fB.close()
        print 'file created: \n\t', shfileA, '\n\t', shfileB
	print
	os.system('cd '+datamockdir+'; chmod a+x *.sh ')

def boss2pcf_creatbashfile_split_nocpfile(catname):

        shfileA = datamockdir+'/Step4.'+catname+'.finesplit.cpignored.sh'

        fA = open(shfileA, 'w'); fA.write('#!/bin/bash\n\n');
        rrange = catinfo_rrange[catname]
        splitinfo = catinfo_binsplit[catname]
        print 'rrange, splitinfo = ', rrange, splitinfo
        for nowfile in allfilelist(catname):
	    nowfile = xyzwfilename(nowfile+'.cpignored')
	    if not isfile(nowfile):
		continue

            ## cmd, filelist
            cmdA, filelistA = cmd_filesplit(nowfile, rrange)

            ## finesplit
            fileA = filelistA[0]; fileAmv = nowfile+'.finesplitted'
            fA.write(cmdA+'\n')
            fA.write('mv '+fileA+' '+fileAmv+'\n')
            path, filename = separate_path_file(fileAmv)
            fA.write('mv '+path+'*.finesplitted '+path+'../xyzw.finesplitted/ \n')
            fA.write('mv '+path+'*.rbinned-split-data.info '+path+'../xyzw.finesplitted/ \n\n')

        fA.close()
        print 'file created: \n\t', shfileA
	print
	os.system('cd '+datamockdir+'; chmod a+x *.sh ')

### Columns in data, mock
icol_data_x = 0
icol_data_y = 1
icol_data_z = 2
icol_data_red = 3
icol_data_index = 4
icol_data_w_nofkp = 5
icol_data_wfkp = 6
icol_data_completeness = 7
icol_data_wfkp2 = 8
icol_data_wcp = 9
icol_data_wnoz =10
icol_data_wstar =11
icol_data_wsee =12

icol_mock_x = 0
icol_mock_y = 1
icol_mock_z = 2
icol_mock_weight = 3
icol_mock_mass = 4
icol_mock_ID = 5
icol_mock_AW = 6
icol_mock_VW = 7
icol_mock_CPW = 8
icol_mock_vx = 9
icol_mock_vy = 10
icol_mock_vz =11
icol_mock_nbar =12
icol_mock_wfkp =13

### Useful quantities: distance at r = 0.43; r = 0.7
r0p69 = comov_r_dft(0.69)
r0p7 = comov_r_dft(0.7)
r0p71 = comov_r_dft(0.71)
r0p43 = comov_r_dft(0.43)
r0p44 = comov_r_dft(0.44)
r0p42 = comov_r_dft(0.42)
r0p15 = comov_r_dft(0.15)
r0p14 = comov_r_dft(0.14)

def print_mock_fmt():
    print 'fmt of mock: \n\t x,y,z, weight [~1/nbar],mass, ID,AndersonWeight,VetoWeight,CPWeight, vx,vy,vz, nbar, wfkp'
    print ' \n#Name\ti-column \n0\t x \n1\t y \n2\t z \n3\t weight [~1/nbar] \n4\t mass \n5\t  ID', \
    ' \n6\t AndersonWeight \n7\t VetoWeight \n8\t CPWeight \n9\t vx \n10\t vy \n11\t vz \n12\t  nbar \n13\t  wfkp'


def print_dataxyzwei_fmt():
    print 'fmt of file named as *Data-xyz-red-wei.txt*: '
    print 'x y z redshift Index W_NOFKP WFKP Completeness WFKP2 WCP WNOZ WSTAR WSEE'
    print ' \n#Name\ti-column \n0\t x \n1\t y \n2\t z \n3\t redshift  \n4\t Index  \n5\t W_NOFKP  \n6\t WFKP',\
    ' \n7\t Completeness \n8\t WFKP2  \n9\t WCP  \n10\t WNOZ  \n11\t WSTAR  \n12\t WSEE'

### Reduce a loaded mock [dropping gals with zero weight]
def reduce_mock(mockdata, 
                Andersonweicol=6,Vetoweicol=7,CPweicol=8,
                Andersonweireduce=True,Vetoweireduce=True,CPweireduce=True):
    reduceddata = []; 
    reducedindex = [];
    ndat = len(mockdata)
    nCP0 = 0
    nCP2 = 0
    nCPlarge = 0
    print '##################################'
    print 'Reducing loaded mock data...'
    print '\tAndersonweireduce =',Andersonweireduce
    print '\tVetoweireduce =',Vetoweireduce
    print '\tCPweireduce =',CPweireduce
    for i in range(ndat):
        AW, VW, CPW = mockdata[i][Andersonweicol], mockdata[i][Vetoweicol], mockdata[i][CPweicol]
        if CPW ==0:
            nCP0+=1
        elif CPW ==2:
            nCP2+=1
        elif CPW >2:
            nCPlarge+=1
        if Andersonweireduce and AW <= 0.1:
            continue
        if Vetoweireduce and VW <= 0.1:
            continue
        if CPweireduce and CPW <= 0.1:
            continue
        reduceddata.append([x for x in mockdata[i]])
        reducedindex.append(i)
    print 'Finishing processing ', ndat, 'lines; #-left = ', len(reduceddata)
    print '# of CP=0:', nCP0
    print '# of CP=2:', nCP2
    print '# of CP>2:', nCPlarge
    return reduceddata, reducedindex

## Reducing the redshift range of data/mock. 
## Although we use the name "BOSSdata", it is applicatable for both data and mock.
def redshiftreduce_data(BOSSdata, zmin=0.43, zmax=0.7, ix=0, iy=1, iz=2):
    BOSSX, BOSSY, BOSSZ = XYZfromdata(BOSSdata, ix,iy,iz)
    BOSSRA, BOSSDEC, BOSSR = list_xyz_to_radecr(BOSSX, BOSSY, BOSSZ)
    ### Reduce BOSS Data
    if catname[7] == 'C':
       indice = Xinrange(BOSSR, comov_r_dft(0.43), comov_r_dft(0.7), onlyreturnindex=True)
    else:
       indice = Xinrange(BOSSR, comov_r_dft(0.15), comov_r_dft(0.43), onlyreturnindex=True)
    return Xs(BOSSdata, indice), indice
