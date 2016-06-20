
import commands
pythonlibPATH=commands.getoutput('echo $pythonlibPATH')
execfile(pythonlibPATH+'/stdA.py')
execfile(pythonlibPATH+'/Tpcftools.py')
pyfile=pythonlibPATH+'/'+'bossdatamock.py'


#########################
## Information

# catname means DR12-CMASS-N, DR12-CMASS-S, ...
# catname2 means HR3, L93, B08, ...

catnameliststr = 'DR12v4-CMASS-N DR12v4-CMASS-S DR12v4-LOWZ-N DR12v4-LOWZ-S DR12v4-CMASS DR12v4-LOWZ'
catnamelist = ['DR12v4-CMASS-N', 'DR12v4-CMASS-S', 'DR12v4-LOWZ-N', 'DR12v4-LOWZ-S', 'DR12v4-CMASS', 'DR12v4-LOWZ']

catinfo_npatch = {
			'DR12v4-CMASS-N' : 4,
			'DR12v4-CMASS-S' : 8,
			'DR12v4-LOWZ-N' : 4,
			'DR12v4-LOWZ-S' : 8,
			'DR12v4-CMASS' : 4,
			'DR12v4-LOWZ' : 4
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

SixBinRedshift_NS = [0.2154098242, 0.3156545036, 0.3869159269, 0.4785988874, 0.5408209467, 0.6186561000]


SixBinRedshift_N = [0.21600642244568982, 0.31646577619767519, 0.38721946672537055, 0.47897499723114539, 0.54059313570594303, 0.61850917062883259]

SixBinRedshift_S = [0.21410525733734107, 0.31388039435706311, 0.38625224442260986, 0.47757238122842527, 0.54144257236794946, 0.61905709858835356]

Seven_Edges = [436.07,     775.25,     975.38,    1172.74,    1367.47,    1507.22,    1772.34 ]

Sky_catname = {
		'N': ['DR12v4-CMASS-N', 'DR12v4-LOWZ-N'],
		'S': ['DR12v4-CMASS-S', 'DR12v4-LOWZ-S'],
		'NS': ['DR12v4-CMASS', 'DR12v4-LOWZ']
}


BinSplittedSamplesInfo = BSSInfo ={'DR12v4-CMASS-N.1of3': [188227,
  [1172.741406785507, 1367.469046283597],
  1291.0535627205716,
  0.47897499723114539],
 'DR12v4-CMASS-N.2of3': [188206,
  [1367.4732224621805, 1507.2199077616378],
  1435.5345281492446,
  0.54059313570594303],
 'DR12v4-CMASS-N.3of3': [188227,
  [1507.2201048559464, 1772.3398525139473],
  1611.3844045605133,
  0.61850917062883259],
 'DR12v4-CMASS-S.1of3': [68966,
  [1172.7636691935306, 1365.3090581727679],
  1287.7082096923421,
  0.47757238122842527],
 'DR12v4-CMASS-S.2of3': [68973,
  [1365.3116190881847, 1510.9783351604913],
  1437.4925575154489,
  0.54144257236794946],
 'DR12v4-CMASS-S.3of3': [68968,
  [1510.98262251195, 1772.3354903962963],
  1612.5944239099224,
  0.61905709858835356],
 'DR12v4-LOWZ-N.1of3': [82746,
  [436.07187411187164, 775.24922582434101],
  618.95992640133159,
  0.21600642244568982],
 'DR12v4-LOWZ-N.2of3': [82745,
  [775.25446374475644, 975.37966607041835],
  886.42475118915456,
  0.31646577619767519],
 'DR12v4-LOWZ-N.3of3': [82745,
  [975.38246906830341, 1172.7353404558635],
  1066.8380569241624,
  0.38721946672537055],
 'DR12v4-LOWZ-S.1of3': [37841,
  [436.08566381217133, 767.20797202023391],
  613.77080078879521,
  0.21410525733734107],
 'DR12v4-LOWZ-S.2of3': [37838,
  [767.2228201274321, 972.79080390696379],
  879.70764728807524,
  0.31388039435706311],
 'DR12v4-LOWZ-S.3of3': [37844,
  [972.80004343547478, 1172.7239164963337],
  1064.4160778934527,
  0.38625224442260986]}
SixBinGalNums = [ 82746 +37841,82745+37838 ,82745+37844 , 188227+68966,188206+68973,188227+68968, ]

HR4LCcatname2list = ['J08', 'LC93', 'M12', 'V13', 'B08',]
HR4catname2list = ['HR4PSB',  'J08.dat.z_0', 'J08.dat.z_0.5', ] + HR4LCcatname2list
catname2list = ['HR3', 'MD'] + HR4catname2list
RSDstrlist = ['noRSD', 'RSD']

def numHR3mock(catname):
	if catname in ['DR12v4-CMASS-N', 'DR12v4-LOWZ-N' ]:
		return 108
	elif catname in ['DR12v4-CMASS-S', 'DR12v4-LOWZ-S' ]:
		return 216
	elif catname in ['DR12v4-CMASS']:
		return 72
	elif catname in ['DR12v4-LOWZ']:
		return 72

def catinfo_nummock(catname, catname2):
	if catname2 == 'HR3':
		return numHR3mock(catname)
	elif catname2 == 'MD':
		return 10
	elif catname2 in catname2list:
		return catinfo_npatch[catname]
	else:
		print 'ERROR (catinfo_nummock): Unknown catname2! ', catname2
		return

miniscan_omwlist = [[0.26, -1.0], [1.0, -1.0], [0.0, -1.0], [0.26, -3.0], [0.26, -0.2] ]
miniscan_scanname = 'miniscan'

#########################################
### NS scan

NSscan_omlist = np.linspace(0.06, 0.86, 9)
NSscan_wlist = np.linspace(-3.0, 0.0, 7)
NSscan_omwlist = sumlist([[[om,w] for om in NSscan_omlist] for w in NSscan_wlist])
NSscan_scanname = 'NS-97Scan'


#########################################
#########################################
### Dense1 scan

#Dense1subscan_omlist = np.linspace(0.06, 0.56, 21)
#Dense1subscan_wlist = np.linspace(-2.0, 0.0, 41)

Dense1subscan_omlist = np.linspace(0.06, 0.66, 25)
Dense1subscan_wlist = np.linspace(-2.5, 0.0, 51)

Dense1subscan_omwlist = sumlist([[[om,w] for om in Dense1subscan_omlist] for w in Dense1subscan_wlist])
Dense1subscan_scanname = 'Dense1sub--'+str(len(Dense1subscan_omlist))+'-'+str(len(Dense1subscan_wlist))+'-Scan'

## Total scan: all om, w in the previous NS-97Scan are included into the omlist, wlist

Dense1scan_omlist = [x for x in Dense1subscan_omlist]
Dense1scan_wlist = [x for x in Dense1subscan_wlist]

Dense1scan_omlist = merge_two_list_finitedigits(Dense1scan_omlist, NSscan_omlist, sortlist=True)
Dense1scan_wlist = merge_two_list_finitedigits(Dense1scan_wlist, NSscan_wlist, sortlist=True)
Dense1scan_omwlist = sumlist([[[om,w] for om in Dense1scan_omlist] for w in Dense1scan_wlist])

Dense1scan_scanname = 'Dense1--'+str(len(Dense1scan_omlist))+'-'+str(len(Dense1scan_wlist))+'-Scan'


#########################################
#########################################
### Dense2 scan
#########################################

Dense2subscan_omlist = np.linspace(0.06, 0.41, 14*5+1)

Dense2subscan_wlist = np.linspace(-1.5, -0.4, 22*2+1)

Dense2subscan_omwlist = sumlist([[[om,w] for om in Dense2subscan_omlist] for w in Dense2subscan_wlist])
Dense2subscan_scanname = 'Dense2sub--'+str(len(Dense2subscan_omlist))+'-'+str(len(Dense2subscan_wlist))+'-Scan'

## Total scan: all om, w in the previous NS-97Scan are included into the omlist, wlist

Dense2scan_omlist = [x for x in Dense2subscan_omlist]
Dense2scan_wlist = [x for x in Dense2subscan_wlist]

Dense2scan_omlist = merge_two_list_finitedigits(Dense2scan_omlist, Dense1scan_omlist, sortlist=True)
Dense2scan_wlist = merge_two_list_finitedigits(Dense2scan_wlist, Dense1scan_wlist, sortlist=True)
Dense2scan_omwlist = sumlist([[[om,w] for om in Dense2scan_omlist] for w in Dense2scan_wlist])

Dense2scan_scanname = 'Dense2--'+str(len(Dense2scan_omlist))+'-'+str(len(Dense2scan_wlist))+'-Scan'

#########################################
#########################################
#########################################
### Dense2sub + Dense1
Dense2sub_Dense1_scanomwlist = [[omw[0], omw[1]] for omw in Dense2subscan_omwlist]

omw_er = 0.000001
#Dense2sub_Dense1_num_omw_added = 0
for i in range(len(Dense1scan_omwlist)):
    om, w = Dense1scan_omwlist[i]
    if not ((om > 0.060 - omw_er) and (om < 0.41 + omw_er) and (w > -1.5-omw_er) and (w < -0.4 + omw_er)):
        Dense2sub_Dense1_scanomwlist.append([om, w])
        #Dense2sub_Dense1_num_omw_added += 1

#########################################
#########################################


#########################################
### Files (name, check)

execfile(pythonlibPATH+'/bossdatamock_files.py')

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

execfile(pythonlibPATH+'/bossdatamock_smu.py')
execfile(pythonlibPATH+'/bossdatamock_ChisqInfo.py')
