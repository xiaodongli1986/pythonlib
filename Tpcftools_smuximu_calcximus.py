
### Calculating a list of chisq contours using a list of options
def smu_ximu_calcchisqs2(
	omlist, wlist,### 1. List of omegam, w
	baseoutputfile, ### 2. Basic name of outputfile
	mumins = [0.01, 0.02, 0.03], ### 3. Range of mumin
	#smu__intxi__settings = smu__intxi__settings_std,
	smusettings_data = smusettings_smax51, 
	smusettings_mock = { 'smin':0.0, 'smax':150.0, 'numsbin':150, 'mumin':0.0, 'mumax':1.0, 'nummubin':120, 'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[] } , 
	s1=6.0, s2=40, NumMubins= [16,21,26],   Smax = 51, Smax_covmock = 150, Smax_sysmock=150, ### 4. Things about s
	cov_catname2='HR3', cov_RSDstr = 'noRSD', cov_nummock=72, ### 5. settings for covmat and systematic
	xifunction = intxi_FractionalS,
	sys_catname2= 'J08', sys_RSDstr ='RSD', sys_catnamesuffix = '', syscor_imocks = range(4), polyfitdeg=3,
	make_plot=False, showfig=False, savefig=True, 
	totbin = 3, delta_time = 100, nummubins = 120, 	normfun=norm_OneSkip, ### Other miscellous
	outputdir = '/home/xiaodongli/software/', suffix = '',
        use_omw_as_sys = False,
	cosmo_rescaled_s = False, z_of_each_bin = [0.2154098242, 0.3156545036, 0.3869159269, 0.4785988874, 0.5408209467, 0.6186561000],
	consider_tilt_in_chisq = False,
	covmat_nondiag_rescale = False,
        usingmapping_for_nonstd_omw = False, basepar_for_smu_mapping = [0.26,-1],
	rebinxi=False, # In early settings we re-bin the values of DD/DR/RR
	polyfit_dxi = None,
	simple_replacement=False,
	use_DenseMapping=False, DM_nbins=750, DM_nummubin=600, DM_smax=150, method2str='divided_pixel',
		):
	ommin, ommax, wmin, wmax = min(omlist), max(omlist), min(wlist), max(wlist)
	#omwlist = sumlist([[[om,w] for om in omlist] for w in wlist])
	omwlist = [[0.31,-1.5]]+sumlist([[[om,w] for om in omlist] for w in wlist]) #DEBUG
	#sarray = smu_smids(s1, s2)
	#sarray = [ x**sfact for x in sarray ]
	
	mumin_mubins = [[[mumin, NumMubin] for mumin in mumins] for NumMubin in NumMubins]
	smu__intxi__settings = [[{ 'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,## In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
				 } for mumin in mumins] for NumMubin in NumMubins]
	smu__intxi__settings_orig = [[{ 'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,## In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
				 } for mumin in mumins] for NumMubin in NumMubins]
	#for mumin in mumins:
	if True:
		nowchisqstr = baseoutputfile
		nowchisqstr += ('.CovMock_'+str(cov_nummock)+'_'+cov_catname2+'_'+cov_RSDstr)
		if cosmo_rescaled_s != False:
			CRbaseom, CRbasew = cosmo_rescaled_s
			nowchisqstr += ('.CosmoRescaledS_'+omwstr(CRbaseom,CRbasew))
		if use_omw_as_sys !=False:
			omsys, wsys = use_omw_as_sys
			nowchisqstr += '.force_bestfit_as_'+omwstr(omsys,wsys)
		if sys_catname2 != 'J08' and use_omw_as_sys == False:
			nowchisqstr += ('.SysCor' + sys_catname2)
		if sys_catnamesuffix != '':
			if sys_catnamesuffix == '-N':
				nowchisqstr += '.SysCor_NGC'
			elif sys_catnamesuffix == '-S':
				nowchisqstr += '.SysCor_SGC'
			else: 
				nowchisqstr += ('.SysCor_'+sys_catnamesuffix)
		if consider_tilt_in_chisq != False:
			nowchisqstr += '.Chisq_ConsiderTilt'
		if polyfitdeg != 3:
			nowchisqstr += ('.SysCor_poly'+str(polyfitdeg))
		if covmat_nondiag_rescale != False:
			nowchisqstr += ('.CovNDRescal%.3f'%covmat_nondiag_rescale)
			def covmat_adjust(covmat):
				for rowa in range(len(covmat)):
					for rowb in range(len(covmat)):
						if rowa != rowb:
							covmat[rowa][rowb] *= covmat_nondiag_rescale
				return covmat
		if usingmapping_for_nonstd_omw == True:
			nowchisqstr += '.SMU_Mapping'
			if simple_replacement:
				nowchisqstr += '.simple_replacement'
			if use_DenseMapping: 
				nowchisqstr += ('.DenseMapping.'+method2str)
				DM_smusettings = {'deltamu': 0,'deltas': 0,'mulist': [],'mumax': 1.0,'mumin': 0.0,\
					'nummubin': DM_nummubin,'numsbin': DM_nbins,'slist': [], 'smax': DM_smax,'smin': 0.0};
				smu__initsmusettings(DM_smusettings)
				smu__initsmusettings(smusettings_data)
				#print DM_smusettings, smusettings_data
		if rebinxi == True:
			nowchisqstr += '.rebinxi'
		if polyfit_dxi!= None:
			nowchisqstr += '.polyfit_dxi_deg'+str(polyfit_dxi)

	        ximudir = smu_ximu_covchisqdir
		#nowfile=outputdir+'/'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
		nowfile = [[outputdir+'/'+nowchisqstr+'_'+smu__intxi_str(smu__intxi__settings_orig[row1][row2])+'.txt' for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
		nowf = [[open(nowfile[row1][row2],'w') for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
		
		print '    Output chisq to file: ' 
		for row2 in range(len(mumins)): 
			for row1 in range(len(NumMubins)):
				print '\t\t', nowfile[row1][row2]
				nowf[row1][row2].write('### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor\n')
		### 6. Compute
		#### 6.1 Load in 
		##### 6.1.1 Mock
		#xis_mock = smu_xis_loadmockrlt(mumin, 
		#				      print_nummock=False)  #fmt: catname, catname2, RSDstr, ibin
		##### 6.1.2 Data
		#xis_data_file = smu_xis_ximudir+'/'+baseoutputfile+'.smax'+str(smax)+'.mumin%03i'%mumin+'.xis'
		#xis_data = smu_xis_loaddatarlt(xis_data_file) # fmt: catname, ibin, omwstr
		#print xis_data['DR12v4-CMASS'][1][omwstr(0.31,-1)]
		
		#### 6.2 Compuate chisqs for all ...
		chisqs_nosyscor, chisqs_syscor = [[{} for row2 in range(len(mumins))] for row1 in range(len(NumMubins))], \
			[[{} for row2 in range(len(mumins))] for row1 in range(len(NumMubins))] # fmt: omwstr, ibin
		iomw = 0
		covmats = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
		xisys_list = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
		dxisys_list = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]

		t0 = time.clock(); t1 = t0; 
		if usingmapping_for_nonstd_omw:
			omstd, wstd = basepar_for_smu_mapping
		smutabstd_list = []
		for omw in omwlist:
			#print omw ## BOSSLi
			t2 = time.clock();
			if t2 - t1 > 10:
				print t2-t0, 'sec ellapses.   iomw/numomw = ', iomw,'/', len(omwlist), \
					'     rat = %.4f'%(float(iomw)/float(len(omwlist)))
				t1 = t2
			om,w = omw; nowomwstr=omwstr(om,w); 
			for row2 in range(len(mumins)): 
				for row1 in range(len(NumMubins)):
					chisqs_nosyscor[row1][row2][nowomwstr], chisqs_syscor[row1][row2][nowomwstr] = [], []
			i_redshiftbin = 0
			xihatdata_now = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
			
			for catname in ['DR12v4-LOWZ', 'DR12v4-CMASS', ]:
				for ibin in range(totbin):
					if usingmapping_for_nonstd_omw:
				          if catname == 'DR12v4-LOWZ':
					    zmedian = SixBinRedshift_NS[ibin]
				          elif  catname == 'DR12v4-CMASS':
					    zmedian = SixBinRedshift_NS[ibin+3]
					  #print 'ibin / zmedian =   ', ibin, zmedian
    				          DAstd, Hstd = DA(omstd, wstd, 0.7, zmedian), Hz(omstd, wstd, 0.7, zmedian)
    				          DAnew, Hnew = DA(om,    w,    0.7, zmedian), Hz(om,    w,    0.7, zmedian)
					if iomw == 0 and usingmapping_for_nonstd_omw:
						if not use_DenseMapping:
							smufile = Tpcfrltfilename(cosmoconvertedfilename(\
	 						 	binsplittedfilename(datafile(catname), ibin+1, totbin),omstd,wstd), \
							 	mubins=nummubins,nbins=Smax,rmax=Smax );
							smutabstd_list.append(smu__loadin(smufile,smusettings_data ))
						else:
							smufile = Tpcfrltfilename(cosmoconvertedfilename(\
	 						 	binsplittedfilename(datafile(catname), ibin+1, totbin),omstd,wstd), \
							 	mubins=DM_nummubin,nbins=DM_nbins,rmax=DM_smax );
							smutabstd_list.append(smu__loadin(smufile,DM_smusettings ))
							
						#print 'smax of the smutabstd_list: ', len(smu__loadin(smufile,smusettings_data ))
						#print smusettings_data
					if cosmo_rescaled_s != False:
						nowz = z_of_each_bin[i_redshiftbin]
						VolBase = DA(CRbaseom,CRbasew,0.7,nowz)**2.0 / Hz(CRbaseom,CRbasew,0.7,nowz)
						VolNow  = DA(om,w,0.7,nowz)**2.0 / Hz(om,w,0.7,nowz)
						ReScale = (VolNow/VolBase)**(1.0/3.0)
						for row2 in range(len(mumins)): 
							for row1 in range(len(NumMubins)):
								s1,s2=smu__intxi__settings[row1][row2]['smin'],smu__intxi__settings[row1][row2]['smax']
								CRs1, CRs2 = s1*ReScale, s2*ReScale
								smu__intxi__settings[row1][row2]['smin'] = CRs1
								smu__intxi__settings[row1][row2]['smax'] = CRs2
						#print om, w, nowz, ': ', CRs1, CRs2
					if i_redshiftbin in [0]:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	
						smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))] 
						#print smu__intxi_file
						if False: #isfile(smu__intxi_file):
							xidata_base = [[smu__intxi_quickload(smu__intxi_file[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						else:
						  if not usingmapping_for_nonstd_omw:
							xidata_base = [[smu__intxi_calcwrite(smufile, smusettings_data,\
								writetofile=False, smu__intxi__settings=smu__intxi__settings[row1][row2],\
								rebinxi=rebinxi)[2] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						  elif not use_DenseMapping:
							smudata = mapping_smudata_to_another_cosmology(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								deltamu=1.0/120.0, smin_mapping=1, smax_mapping=51,\
								simple_replacement=simple_replacement )         
							xidata_base = [[smu__intxi_calcwrite(smufile, smusettings_data, \
								writetofile=False,smudata=smudata,smu__intxi__settings=smu__intxi__settings[row1][row2],\
								rebinxi=rebinxi)[2] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						  else:
							smudata = mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								method=method2str,\
								smin_mapping=1, smax_mapping=51,\
                    						deltas1=DM_smusettings['deltas'],  deltamu1=DM_smusettings['deltamu'],\
					                    	deltas2=smusettings_data['deltas'],deltamu2=smusettings_data['deltamu'])
							xidata_base = [[smu__intxi_calcwrite(smufile, smusettings_data, \
								writetofile=False, smudata=smudata, \
								smu__intxi__settings=smu__intxi__settings[row1][row2],rebinxi=rebinxi)[2]\
								for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						## B. normalize the amplitude
						if polyfit_dxi== None:
						  xihatdata_base =  [[ normfun(xidata_base[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						else:
						  xihatdata_base =  [[ normfun(polyfitY(range(len(xidata_base[row1][row2])),xidata_base[row1][row2],polyfit_dxi)) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						for row2 in range(len(mumins)):
							for row1 in range(len(NumMubins)):
								xihatdata_now[row1][row2].append([x for x in xihatdata_base[row1][row2]])
						if iomw == 0:
							xicov_base = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							for imock in range(cov_nummock):
								smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_covmock,rmax=Smax_covmock ); 
								smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
								for row2 in range(len(mumins)):
								 for row1 in range(len(NumMubins)):
								  if isfile(smu__intxi_file[row1][row2]):
									xicov_base[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
								  else:
									xicov_base[row1][row2].append(smu__intxi_calcwrite(smufile, smusettings_mock,\
										writetofile=False, smu__intxi__settings=\
										smu__intxi__settings_orig[row1][row2],rebinxi=rebinxi)[2])
							xihatcov_base  = [[ [ normfun(X) for X in xicov_base[row1][row2]] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							#if polyfit_dxi== None:
							# xihatcov_base  = [ normfun(X) for X in xicov_base]
							#else:
							# xihatcov_base  = [ normfun(polyfitY(range(len(X)),X,polyfit_dxi)) for X in xicov_base]
							for row2 in range(len(mumins)):
                                                         for row1 in range(len(NumMubins)):
							  covmats[row1][row2].append([])
							xisys_base = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_sysmock,rmax=Smax_sysmock ); 
									smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=\
									  smu__intxi__settings_orig[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
									#print 'smu__intxi_file', smu__intxi_file # HB
									for row2 in range(len(mumins)):
									 for row1 in range(len(NumMubins)):
									  if isfile(smu__intxi_file[row1][row2]):
										xisys_base[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
									  else:
										xisys_base[row1][row2].append(smu__intxi_calcwrite(smufile, \
										  smusettings_mock, writetofile=False, \
										  smu__intxi__settings=smu__intxi__settings_orig[row1][row2],\
										  rebinxi=rebinxi)[2])
							else:
								omsys, wsys = use_omw_as_sys
								smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); 
								smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
								for row2 in range(len(mumins)):
								 for row1 in range(len(NumMubins)):
                                                	 	   if isfile(smu__intxi_file[row1][row2]):
		                                                        xisys_base[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
		                                                   else:
                		                                        xisys_base[row1][row2].append(smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings[row1][row2],rebinxi=rebinxi)[2])
							for row2 in range(len(mumins)): 
							 for row1 in range(len(NumMubins)):
							  xisys_list[row1][row2].append(get_avg_array([normfun(X) for X in xisys_base[row1][row2]]))
							  dxisys_list[row1][row2].append([])
					else:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	
						smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings[row1][row2]) \
							for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						#print smu__intxi_file; sys.exit()
						xidata = [[ 0 for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						smudata2_done = False
						for row2 in range(len(mumins)):
						 for row1 in range(len(NumMubins)):
						  if isfile(smu__intxi_file[row1][row2]):
							#print 'load in smu__intxi_file: ', smu__intxi_file
							xidata[row1][row2] = smu__intxi_quickload(smu__intxi_file[row1][row2])
						  else:
                                                    if not usingmapping_for_nonstd_omw:
                                                        xidata[row1][row2] = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
                                                    elif not use_DenseMapping:
							if not smudata2_done:
                                                         smudata2 = mapping_smudata_to_another_cosmology(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
                                                                deltamu=1.0/120.0, smin_mapping=1, smax_mapping=51,\
								simple_replacement=simple_replacement )
							 smudata2_done = True
							#print  'smutabstd_list[i_redshiftbin][10][10][9], smudata[10][10][9] = ', smutabstd_list[i_redshiftbin][10][10][9], smudata[10][10][9]
                                                        xidata[row1][row2] = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, \
                                                           smudata=smudata2, smu__intxi__settings=smu__intxi__settings[row1][row2],rebinxi=rebinxi)[2]
						    else: 
							if not smudata2_done:
                                                         smudata2 = mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								method=method2str,\
                                                                smin_mapping=1, smax_mapping=51,\
                    						deltas1=DM_smusettings['deltas'],  deltamu1=DM_smusettings['deltamu'],\
					                    	deltas2=smusettings_data['deltas'],deltamu2=smusettings_data['deltamu'])
							 smudata2_done = True
                                                        xidata[row1][row2] = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, \
                                                                smudata=smudata2, smu__intxi__settings=smu__intxi__settings[row1][row2],rebinxi=rebinxi)[2]
						  #print 'om, w, i_redshiftbin, normto1(xidata) = ', om, w, i_redshiftbin, normto1(xidata) #DEBUG
						  ## B. normalize the amplitude
						  #xihatdata =   normfun(xidata)
						if polyfit_dxi== None:
						   xihatdata =  [[ normfun(xidata[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						else:
						   xihatdata =  [[ normfun(polyfitY(range(len(xidata[row1][row2])),xidata[row1][row2],polyfit_dxi)) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
						xihatdata_now[row1][row2].append([x for x in xihatdata])
						dxidata = [[ get_diffarray(xihatdata_base[row1][row2], xihatdata[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
					#if om == 0.26 and w == -1.0:
						#	print om,w, ':', ibin, xidata, xihatdata
						if iomw == 0:
							xisys = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
	                                                if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_sysmock,rmax=Smax_sysmock ); 
									smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
									#print 'smu__intxi_file', smu__intxi_file # HB
									for row2 in range(len(mumins)): 
									 for row1 in range(len(NumMubins)):
									  if isfile(smu__intxi_file[row1][row2]):
										xisys[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
									  else:
										xisys[row1][row2].append(smu__intxi_calcwrite(smufile, smusettings_mock, writetofile=False, smu__intxi__settings=smu__intxi__settings_orig[row1][row2],rebinxi=rebinxi)[2])
	                                                else:
        	                                                 omsys, wsys = use_omw_as_sys
                	                                         smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); 
								 smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=\
									smu__intxi__settings[row1][row2]) for row2 in range(len(mumins))] \
									for row1 in range(len(NumMubins))]
								 for row2 in range(len(mumins)):
								  for row1 in range(len(NumMubins)):
                        	                                   if isfile(smu__intxi_file[row1][row2]):
                                	                                xisys[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
                                        	                   else:
                                                	                xisys[row1][row2].append(smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings[row1][row2],rebinxi=rebinxi)[2])

							polyfitdeg = polyfitdeg
							xihatsys_base  = [[ [normfun(X) for X in xisys_base[row1][row2]] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							xihatsys  = [[ [normfun(X) for X in xisys[row1][row2]] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							#print get_avg_array(xisys)
							#print 'i_redshiftbin, xihatsys, xihatdata', i_redshiftbin, get_avg_array(xihatsys), xihatdata  # HB
							dxisys = [[ get_diffarray(get_avg_array(xihatsys_base[row1][row2]),get_avg_array(xihatsys[row1][row2])) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							dxisys = [[ polyfitY(range(len(dxisys[row1][row2])), dxisys[row1][row2], deg=polyfitdeg) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
                                                        for row2 in range(len(mumins)):
                                                         for row1 in range(len(NumMubins)):
                                                          xisys_list[row1][row2].append(get_avg_array(xihatsys[row1][row2]))
						   	  dxisys_list[row1][row2].append([x for x in dxisys[row1][row2]])

                                                        xicov = [[[] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
                                                        for imock in range(cov_nummock):
                                                                smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_covmock,rmax=Smax_covmock ); 
								smu__intxi_file = [[smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
                                                                for row2 in range(len(mumins)):
                                                                 for row1 in range(len(NumMubins)):
								  if isfile(smu__intxi_file[row1][row2]):
                                                                	xicov[row1][row2].append(smu__intxi_quickload(smu__intxi_file[row1][row2]))
								  else:
									xicov[row1][row2].append(smu__intxi_calcwrite(smufile, smusettings_mock, writetofile=False, smu__intxi__settings=smu__intxi__settings_orig[row1][row2],rebinxi=rebinxi)[2])

                                                        xihatcov  = [[ [ normfun(X) for X in xicov[row1][row2]] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							#if polyfit_dxi== None:
                                                        # xihatcov  = [ normfun(X) for X in xicov]
							#else:
                                                        # xihatcov  = [ normfun(polyfitY(range(len(X)),X,polyfit_dxi)) for X in xicov]
                                                        #X0=range(len(xihatcov[row1][row2][0]))
                                                        dxicov  = [[ [ get_diffarray(xihatcov_base[row1][row2][row], xihatcov[row1][row2][row]) \
                                                                for row in range(len(xihatcov[row1][row2]))] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
							#print 'cov_nummock = ', cov_nummock
							#print 'xicov[0][0] = ', xicov[0][0]
							#print 'dxicov[0][0] = ', dxicov[0][0]
							for row2 in range(len(mumins)):
							 for row1 in range(len(NumMubins)):
                                                          if covmat_nondiag_rescale == False:
                                                                covmats[row1][row2].append(get_covmat(transpose(dxicov[row1][row2])))
                                                          else:
                                                                covmats[row1][row2].append(covmat_adjust(get_covmat(transpose(dxicov[row1][row2]))))
                                                        #print '\n\n ibin, covmat = ', ibin, get_covmat(transpose(dxicov))
                                                covmat = [[covmats[row1][row2][i_redshiftbin] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
                                                #if iomw == 0: print np.mat(covmat).I # BOSSLi
                                                chisq_nosyscor = [[ chisq_like_cov_xbar(dxidata[row1][row2], covmat[row1][row2])[0] for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]
                                                chisq_syscor   = [[ chisq_like_cov_xbar(XplusY(dxidata[row1][row2],dxisys_list[row1][row2][i_redshiftbin],b=-1), covmat[row1][row2]) for row2 in range(len(mumins))] for row1 in range(len(NumMubins))]

                                                #print nowomwstr, '/ redshiftbin=', i_redshiftbin, '/ dxi = ', dxidata, chisq_nosyscor, chisq_syscor
						for row2 in range(len(mumins)):
						 for row1 in range(len(NumMubins)):
                                                  chisqs_nosyscor[row1][row2][nowomwstr].append(chisq_nosyscor[row1][row2])
                                                  chisqs_syscor[row1][row2][nowomwstr].append(chisq_syscor[row1][row2])
                                        i_redshiftbin += 1
			for row2 in range(len(mumins)):
			 for row1 in range(len(NumMubins)):
                          X1, X2 = chisqs_nosyscor[row1][row2][nowomwstr], chisqs_syscor[row1][row2][nowomwstr]
                          str1, str2, strA, strB = \
                                array_to_str(X1), array_to_str(X2), str(sum(X1)), str(sum(X2))
                          #print str(mumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2
                          nowf[row1][row2].write(str(mumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2+'\n')
                        iomw += 1
                        #if iomw == 2:
                        #       sys.exit()
                        #print om, w, chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
                #print chisqs_nosyscor.keys()
                #sys.exit(' ')
		for row2 in range(len(mumins)):
			 for row1 in range(len(NumMubins)):
		                nowf[row1][row2].close()
                if make_plot:
			print 'Error! Make plot not supported!' 
                	chisqlist_nosyscor = [ [  sum(chisqs_nosyscor[omwstr(om,w)])    for om in omlist] for w in wlist ]
                	chisqlist_syscor =   [ [  sum(chisqs_syscor[omwstr(om,w)])      for om in omlist] for w in wlist ]
                        fig = plt.figure(figsize=(16,6))
                        ax1 = fig.add_subplot(121)
                        ax2 = fig.add_subplot(122)
                        sigA, sigB, sigC = sig1, sig2, sig3
                        plot_contour(ax1, omlist, wlist, chisqlist_nosyscor, label='NO Sys Cor',
                            ommin = ommin, ommax = ommax, wmin = wmin, wmax = wmax,
                            do_smooth=False, smsigma=0.5, sigA = sigA, sigB = sigB, sigC = sigC)
                        plot_contour(ax2, omlist, wlist, chisqlist_syscor, label='After Sys Cor',
                            ommin = ommin, ommax = ommax, wmin = wmin, wmax = wmax,
                            do_smooth=False, smsigma=0.5, sigA = sigA, sigB = sigB, sigC = sigC)
                        fig.suptitle(nowchisqstr+'.mumin'+str(mumin))
                        if savefig:
                                figname = outputdir+'/'+nowchisqstr+'.mumin'+str(mumin)+'.png';
                                print '   Figure saved: ', figname
                                plt.savefig(figname, fmt='png')
                        if showfig: plt.show()
