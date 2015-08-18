######################################
#### commands

### These are filenames of 2pcf but we think it is more related with Sec. commands rather than files, so put them in this section

def Tpcf_keyname(catname, catname2, RSDstr='noRSD', imock=0, ibin=0, totbin=0):
    keyname = catname+'-'+catname2
    if catname2 != 'data':
        keyname += ('-'+RSDstr+'-mock'+str(imock))
    if totbin >= 1:
        keyname += ('--'+str(ibin)+'bin-of-'+str(totbin))
    return keyname

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
