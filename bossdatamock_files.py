
### Names of 'original files'

datadir_orig = '/home/xiaodongli/SparseFilaments/data/input/DR12/'
mockdir_orig = '/home/xiaodongli/SparseFilaments/data/input/HR4/DR12_mock/'
datamockdir = '/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/'

### files
def datafile(catname, suffix = ''):
	return datamockdir+catname+'/data'+suffix
def dataranfile(catname, suffix = ''):
	return datamockdir+catname+'/dataran'+suffix
def mockfile(catname, catname2, imock, RSDstr='RSD', suffix = ''):
	return datamockdir+catname+'/'+catname2+'.'+RSDstr+'.%03i'%imock+suffix
def mockranfile(catname, suffix = ''):
	return datamockdir+catname+'/mockran'+suffix

def Tpcf_finesplitte_keyfilename(catname, catname2, imock=0, RSDstr='RSD', rmax=150, nbins=150, mubins=120, shortkeyname=False):

	if catname2 != 'data':
		nowfile = separate_path_file(mockfile('', catname2, imock, RSDstr))[1]
		if not shortkeyname:
			keyname = catname+'--'+catname2+'--'+RSDstr+'--mock'+str(imock+1)
		else:
			keyname = catname2+'--'+RSDstr+'--mock'+str(imock+1)
	else:
		nowfile = keyname = 'data'
	filename = datamockdir+'/'+catname+'/xyzw.finesplitted/'+nowfile+'.xyzw.finesplitted'+Tpcf_suffix(rmax, nbins, mubins)+'.2pcf'
	return keyname, filename

def Tpcf_binsplitted_keyfilename(catname, catname2, imock=0, RSDstr='RSD', ibin=1, totbin=3, rmax=150, nbins=150, mubins=120, shortkeyname=False ):

	if catname2 != 'data':
		nowfile = separate_path_file(mockfile('', catname2, imock, RSDstr))[1]
		if not shortkeyname:
			keyname = catname+'--'+catname2+'--'+RSDstr+'--mock'+str(imock+1)
		else:
			keyname = catname2+'--'+RSDstr+'--mock'+str(imock+1)
	else:
		nowfile = keyname = 'data'
	binstr = str(ibin)+'of'+str(totbin)
	keyname += ('--binsplitted-'+binstr)
	filename = datamockdir+'/'+catname+'/xyzw.binsplitted/'+nowfile+'.xyzw.'+binstr+Tpcf_suffix(rmax, nbins, mubins)+'.2pcf'
	return keyname, filename


### orig files
def datafile_orig(catname):
	return '/home/xiaodongli/SparseFilaments/data/input/DR12/OurRandoms/'+catname+'/Om0.260w-1.000/Data-xyz-red-wei.txt'

def dataranfile_orig(catname):
	if catname == 'DR12v4-CMASS-N':
		return datadir_orig+'/cmass-dr12v4-N-Reid.ran.ra_dec_z_wfkp_nz.txt'
	elif catname == 'DR12v4-CMASS-S':
		return datadir_orig+'/cmass-dr12v4-S-Reid.ran.ra_dec_z_wfkp_nz.txt'
	elif catname == 'DR12v4-LOWZ-N':
		return datadir_orig+'/lowz-dr12v4-N-Reid.ran.ra_dec_z_wfkp_nz.txt'
	elif catname == 'DR12v4-LOWZ-S':
		return datadir_orig+'/lowz-dr12v4-S-Reid.ran.ra_dec_z_wfkp_nz.txt'

def mockfile_orig(catname, catname2, imock, RSDstr='RSD'):
	isky, ipatch = imock_to_iskypatch(imock, npatch=catinfo_npatch[catname])
	if catname2 == 'HR3':
		return mockdir_orig+catname+'/Xiao-dong.all.pv.000%02i'%isky+\
			'.dat.compact.patch%1i'%ipatch+'.'+RSDstr+'-radial-selected'
	elif catname2 == 'J08.dat.z_0':
		return mockdir_orig+catname+'/'+catname2+'.mge8e12.%02i'%isky+'.compact.patch%1i'%ipatch+'.'+RSDstr+'-radial-selected'
	elif catname2 == 'J08.dat.z_0.5':
		return mockdir_orig+catname+'/'+catname2+'.mge4.05e12.%02i'%isky+'.compact.patch%1i'%ipatch+'.'+RSDstr+'-radial-selected'

# J08.dat.z_0.mge8e12.00.compact.patch4.RSD-radial-selected
	elif catname2 == 'LC93':
		return mockdir_orig+catname+'/'+catname2+'_massge2.0e12_rle1817.%02i'%isky+'.compact.patch%1i'%ipatch+'.'+RSDstr+'-radial-selected'

def mockranfile_sep_orig(catname, imock):
	return mockdir_orig+catname+'/random.%02i'%imock+'.0vxvyvz.log10massfrom10to14.compact.patch2.noRSD-radial-selected'
def mockranfile_orig(catname):
	return mockdir_orig+catname+'/mockran'

### Useful file list
def HR3noRSDmocklist(catname):
	return [mockfile(catname, 'HR3', imock, 'noRSD') for imock in range(27*catinfo_npatch[catname])]
def HR3RSDmocklist(catname):
	return [mockfile(catname, 'HR3', imock, 'RSD') for imock in range(27*catinfo_npatch[catname])]
def HR4noRSDmocklist(catname,catname2):
	return [mockfile(catname, catname2, imock, 'noRSD') for imock in range(catinfo_npatch[catname])]
def HR4RSDmocklist(catname,catname2):
	return [mockfile(catname, catname2, imock, 'RSD') for imock in range(catinfo_npatch[catname])]
def HR3noRSDmocklist_orig(catname):
	return [mockfile_orig(catname, 'HR3', imock, 'noRSD') for imock in range(27*catinfo_npatch[catname])]
def HR3RSDmocklist_orig(catname):
	return [mockfile_orig(catname, 'HR3', imock, 'RSD') for imock in range(27*catinfo_npatch[catname])]
def HR4noRSDmocklist_orig(catname,catname2):
	return [mockfile_orig(catname, catname2, imock, 'noRSD') for imock in range(catinfo_npatch[catname])]
def HR4RSDmocklist_orig(catname,catname2):
	return [mockfile_orig(catname, catname2, imock, 'RSD') for imock in range(catinfo_npatch[catname])]

### all file list: Weill be used to creating shortcuts of files; original and shortcut files must have 1-to-1 mapping in right order; Very important!

def allfilelist(catname ):
	return [datafile(catname), dataranfile(catname), mockranfile(catname)] + \
		HR3noRSDmocklist(catname) + \
		HR3RSDmocklist(catname) + \
		sumlist([HR4noRSDmocklist(catname,catname2) for catname2 in HR4catname2list]) + \
		sumlist([HR4RSDmocklist(catname,catname2) for catname2 in HR4catname2list])
def mockfilelist(catname):
	return	HR3noRSDmocklist(catname) + \
		HR3RSDmocklist(catname) + \
		sumlist([HR4noRSDmocklist(catname,catname2) for catname2 in HR4catname2list]) + \
		sumlist([HR4RSDmocklist(catname,catname2) for catname2 in HR4catname2list])

def HR4mockfilelist(catname):
	return sumlist([HR4noRSDmocklist(catname,catname2) for catname2 in HR4catname2list]) + \
		sumlist([HR4RSDmocklist(catname,catname2) for catname2 in HR4catname2list])

def allfilelist_orig(catname):
	return [datafile_orig(catname), dataranfile_orig(catname), mockranfile_orig(catname)] + \
		HR3noRSDmocklist_orig(catname) + \
		HR3RSDmocklist_orig(catname) + \
		sumlist([HR4noRSDmocklist_orig(catname,catname2) for catname2 in HR4catname2list]) + \
		sumlist([HR4RSDmocklist_orig(catname,catname2) for catname2 in HR4catname2list])

def allfilelist_xyzw(catname):
	return [xyzwfilename(nowfile) for nowfile in allfilelist(catname)]

def imock_to_iskypatch(imock, npatch):
    isky = imock / npatch
    ipatch = imock - isky*npatch+1
    return isky, ipatch

def xyzwfilename(datafilename):
	path, filename = separate_path_file(datafilename)
	return path+'/xyzw/'+filename+'.xyzw'

def finesplittedfilename(datafilename):
	path, filename = separate_path_file(datafilename)
	return path+'/xyzw.finesplitted/'+filename+'.xyzw.finesplitted'

def binsplittedfilename(datafilename, ibin, totbin=3):
	path, filename = separate_path_file(datafilename)
	return path+'/xyzw.binsplitted/'+filename+'.xyzw.'+str(ibin)+'of'+str(totbin)

def binsplittedfilename_list(datafilelist, totbin=3, ibin=''):
	rlt = []
	for nowfile in datafilelist:
		if ibin == '':
			for ibin in range(totbin):
				rlt.append(binsplittedfilename(nowfile, ibin+1, totbin))
		else:
			rlt.append(binsplittedfilename(nowfile, ibin, totbin))
	return rlt

	
### Useful routines checking files...

def boss2pcf_lnfiles(catname, HR4catname2list=HR4catname2list, do_link = False):
    ####################################################
    ### 1. Check whether original files exist?
    
    filelist_orig = allfilelist_orig(catname)
    print 'Check ', catname, '...'
    print '\tcheck data, dataran, mockran;'
    print '\tcheck HR3 mocks: ',27*catinfo_npatch[catname],' noRSD, ',27*catinfo_npatch[catname],' RSD '
    print '\tcheck HR4 mocks: ',catinfo_npatch[catname],' noRSD, ',\
		catinfo_npatch[catname],' RSD; HR4cat2namelist =  ', HR4catname2list
    nfile = 1 + 1 + 1 + 27*catinfo_npatch[catname]*2 + catinfo_npatch[catname]*2*len(HR4catname2list)
    print '\tTotal #-files = 1 + 1 + 1 + 27*catinfo_npatch[catname]*2 + 27*catinfo_npatch[catname]*2*len(HR4catname2list) = ', nfile
    if nfile != len(filelist_orig):
        print 'ERROR! #-file mismatching: ', nfile, len(filelist_orig)
        return
    nmissfile = 0
    isfilelist = []
    missfilelist = []
    for nowfile in filelist_orig:
        if not isfile(nowfile):
            nmissfile += 1
            missfilelist.append(nowfile)
        else:
            isfilelist.append(nowfile)
    print '\n\tFinishing checking ', nfile, 'files. ', nmissfile, 'files missing.'
    if nmissfile != 0:
        print '\tlist of missing files: ', missfilelist

    ####################################################
    ### 2. ln original files to files
    
    filelist = allfilelist(catname)
    if len(filelist) != len(filelist_orig):
        print 'Error! #-file mismatches: ', len(filelist), len(filelist_orig)
        return
        
    if do_link:
        print '\n\tCreating shortcuts of files...'
        for i in range(nfile):
            print '#### ', i
            print '\t', filelist_orig[i]
            print '\t', filelist[i]
            cmd = 'ln -s '+filelist_orig[i]+' '+filelist[i]
            os.system(cmd)
    else:
        print '\n\tWill not creating shortcuts of files.'


def boss2pcf_check_xyzwfiles(catname, HR4catname2list=HR4catname2list):
    ####################################################
    ### 1. Check whether original files exist?
    
    filelist = allfilelist(catname)
    print 'Check ', catname, '...'
    print '\tcheck data, dataran, mockran;'
    print '\tcheck HR3 mocks: ',27*catinfo_npatch[catname],' noRSD, ',27*catinfo_npatch[catname],' RSD '
    print '\tcheck HR4 mocks: ',catinfo_npatch[catname],' noRSD, ',\
		catinfo_npatch[catname],' RSD; HR4cat2namelist =  ', HR4catname2list
    nfile = 1 + 1 + 1 + 27*catinfo_npatch[catname]*2 + catinfo_npatch[catname]*2*len(HR4catname2list)
    print '\tTotal #-files = 1 + 1 + 1 + 27*catinfo_npatch[catname]*2 + 27*catinfo_npatch[catname]*2*len(HR4catname2list) = ', nfile
    if nfile != len(filelist):
        print 'ERROR! #-file mismatching: ', nfile, len(filelist)
        return
    nmissfile = 0
    isfilelist = []
    missfilelist = []
    for nowfile in filelist:
        nowfile = xyzwfilename(nowfile)
        if not isfile(nowfile):
            nmissfile += 1
            missfilelist.append(nowfile)
        else:
            isfilelist.append(nowfile)
    print '\n\tFinishing checking ', nfile, 'files. ', nmissfile, 'files missing.'
    if nmissfile != 0:
        print '\tlist of missing files: ', missfilelist

def boss2pcf_check_splitfiles(catname, totbin = 3, printinfo=True):
	splitfilelist = []
	for nowfile in allfilelist(catname):
	    splitfilelist.append(finesplittedfilename(nowfile))
	    for ibin in range(totbin):
	        splitfilelist.append(binsplittedfilename(nowfile, ibin+1,totbin))
	isfiles(splitfilelist,printinfo)
