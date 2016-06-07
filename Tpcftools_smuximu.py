
def norm_OneSkip(X, deg=1.0):
                        avg = sum(X) / float(len(X))
                        imid = len(X)/2 + 1
                        Y = [X[row]/avg for row in range(imid)] + [X[row]/avg for row in range(imid+1,len(X))]
                        return Y
def norm_powerC(X, deg=1.0):
                        avg = sum([x**deg for x in X])**(1.0/deg) / float(len(X))
                        return [X[row]/avg for row in range(len(X)-2)]


def smu_ximu_calcchisqs(
	omlist, wlist,### 1. List of omegam, w
	baseoutputfile, ### 2. Basic name of outputfile
	mumins = [0.01, 0.02, 0.03], ### 3. Range of mumin
	#smu__intxi__settings = smu__intxi__settings_std,
	smusettings_data = smusettings_smax51, 
	smusettings_mock = { 'smin':0.0, 'smax':150.0, 'numsbin':150, 'mumin':0.0, 'mumax':1.0, 'nummubin':120, 'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[] } , 
	s1=6.0, s2=40, NumMubin=21,   Smax = 51, Smax_mock = 150, ### 4. Things about s
	cov_catname2='HR3', cov_RSDstr = 'noRSD', cov_nummock=72, ### 5. settings for covmat and systematic
	xifunction = intxi_FractionalS,
	sys_catname2= 'J08', sys_RSDstr ='RSD', sys_catnamesuffix = '', syscor_imocks = range(4), polyfitdeg=3,
	make_plot=True, showfig=False, savefig=True, 
	totbin = 3, delta_time = 100, nummubins = 120, 	normfun=norm_OneSkip, ### Other miscellous
	outputdir = '/home/xiaodongli/software/', suffix = '',
        use_omw_as_sys = False,
	cosmo_rescaled_s = False, z_of_each_bin = [0.2154098242, 0.3156545036, 0.3869159269, 0.4785988874, 0.5408209467, 0.6186561000],
		):
	ommin, ommax, wmin, wmax = min(omlist), max(omlist), min(wlist), max(wlist)
	omwlist = sumlist([[[om,w] for om in omlist] for w in wlist])
	#sarray = smu_smids(s1, s2)
	#sarray = [ x**sfact for x in sarray ]
	for mumin in mumins:
		smu__intxi__settings = {
            	'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,# In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
         	      }
		smu__intxi__settings_orig = {
            	'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,## In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
         	      }
		nowchisqstr = baseoutputfile+'_'+smu__intxi_str(smu__intxi__settings)

		nowchisqstr += ('.NumCovMock'+str(cov_nummock))
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
		if polyfitdeg != 3:
			nowchisqstr += ('.SysCor_poly'+str(polyfitdeg))

	        ximudir = smu_ximu_covchisqdir
		#nowfile=outputdir+'/'+baseoutputfile+'_'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
		nowfile=outputdir+'/'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
		print '    Output chisq to file: ', nowfile
		nowf.write('### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor\n')
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
		chisqs_nosyscor, chisqs_syscor = {}, {} # fmt: omwstr, ibin
		iomw = 0
		covmats = []
		dxisys_list = []

		t0 = time.clock(); t1 = t0; 
		for omw in omwlist:
			t2 = time.clock();
			if t2 - t1 > 10:
				print t2-t0, 'sec ellapses.   iomw/numomw = ', iomw,'/', len(omwlist), \
					'     rat = %.4f'%(float(iomw)/float(len(omwlist)))
				t1 = t2
			om,w = omw; nowomwstr=omwstr(om,w); chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr] = [], []
			i_redshiftbin = 0
			for catname in ['DR12v4-LOWZ', 'DR12v4-CMASS', ]:
				for ibin in range(totbin):
					if cosmo_rescaled_s != False:
						nowz = z_of_each_bin[i_redshiftbin]
						VolBase = DA(CRbaseom,CRbasew,0.7,nowz)**2.0 / Hz(CRbaseom,CRbasew,0.7,nowz)
						VolNow  = DA(om,w,0.7,nowz)**2.0 / Hz(om,w,0.7,nowz)
						ReScale = (VolNow/VolBase)**(1.0/3.0)
						CRs1, CRs2 = s1*ReScale, s2*ReScale
						smu__intxi__settings['smin'] = CRs1
						smu__intxi__settings['smax'] = CRs2
						#print om, w, nowz, ': ', CRs1, CRs2
					if i_redshiftbin in [0]:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
						#print smu__intxi_file
						if isfile(smu__intxi_file):
							xidata_base = smu__intxi_quickload(smu__intxi_file)
						else:
							xidata_base = smu__intxi_calcwrite(smufile, smusettings_data, smu__intxi__settings=smu__intxi__settings)[2]
						## B. normalize the amplitude
						xihatdata_base =   normfun(xidata_base)
						if iomw == 0:
							xicov_base = []
							for imock in range(cov_nummock):
								smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1), mubins=nummubins,nbins=Smax_mock,rmax=Smax_mock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
								if isfile(smu__intxi_file):
									xicov_base.append(smu__intxi_quickload(smu__intxi_file))
								else:
									xicov_base.append(smu__intxi_calcwrite(smufile, smusettings_mock, smu__intxi__settings=smu__intxi__settings_orig)[2])
								xihatcov_base  = [ normfun(X) for X in xicov_base]
							covmats.append([])
							xisys_base = []
							if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1), mubins=nummubins,nbins=Smax_mock,rmax=Smax_mock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
									if isfile(smu__intxi_file):
										xisys_base.append(smu__intxi_quickload(smu__intxi_file))
									else:
										xisys_base.append(smu__intxi_calcwrite(smufile, smusettings_mock, smu__intxi__settings=smu__intxi__settings_orig)[2])
							else:
								omsys, wsys = use_omw_as_sys
								smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                                                		if isfile(smu__intxi_file):
		                                                        xisys_base.append(smu__intxi_quickload(smu__intxi_file))
		                                                else:
                		                                        xisys_base.append(smu__intxi_calcwrite(smufile, smusettings_data, smu__intxi__settings=smu__intxi__settings)[2])
							dxisys_list.append([])
					else:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
						if isfile(smu__intxi_file):
							xidata = smu__intxi_quickload(smu__intxi_file)
						else:
							xidata = smu__intxi_calcwrite(smufile, smusettings_data, smu__intxi__settings=smu__intxi__settings)[2]
						## B. normalize the amplitude
						xihatdata =   normfun(xidata)
						dxidata = get_diffarray(xihatdata_base, xihatdata)

						#if om == 0.26 and w == -1.0:
						#	print om,w, ':', ibin, xidata, xihatdata
						if iomw == 0:
							xisys = []
	                                                if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1), mubins=nummubins,nbins=Smax_mock,rmax=Smax_mock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
									if isfile(smu__intxi_file):
										xisys.append(smu__intxi_quickload(smu__intxi_file))
									else:
										xisys.append(smu__intxi_calcwrite(smufile, smusettings_mock, smu__intxi__settings=smu__intxi__settings_orig)[2])
	                                                else:
        	                                                 omsys, wsys = use_omw_as_sys
                	                                         smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                        	                                 if isfile(smu__intxi_file):
                                	                                xisys.append(smu__intxi_quickload(smu__intxi_file))
                                        	                 else:
                                                	                xisys.append(smu__intxi_calcwrite(smufile, smusettings_data, smu__intxi__settings=smu__intxi__settings)[2])

							polyfitdeg = polyfitdeg
							xihatsys_base  = [ normfun(X) for X in xisys_base]
							xihatsys  = [ normfun(X) for X in xisys]
							dxisys = get_diffarray(get_avg_array(xihatsys_base),get_avg_array(xihatsys))
							dxisys = polyfitY(range(len(dxisys)), dxisys, deg=polyfitdeg)
							dxisys_list.append([x for x in dxisys])

                                                        xicov = []
                                                        for imock in range(cov_nummock):
                                                                smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1), mubins=nummubins,nbins=Smax_mock,rmax=Smax_mock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
								if isfile(smu__intxi_file):
                                                                	xicov.append(smu__intxi_quickload(smu__intxi_file))
								else:
									xicov.append(smu__intxi_calcwrite(smufile, smusettings_mock, smu__intxi__settings=smu__intxi__settings_orig)[2])

                                                                xihatcov  = [ normfun(X) for X in xicov]
		
							dxicov  = [ get_diffarray(xihatcov_base[row], xihatcov[row]) \
								for row in range(len(xihatcov))]
							covmats.append(get_covmat(transpose(dxicov)))
							#print '\n\n ibin, covmat = ', ibin, get_covmat(transpose(dxicov))
						covmat = covmats[i_redshiftbin]
						chisq_nosyscor, like = chisq_like_cov_xbar(dxidata, covmat)
						#chisq_syscor, like   = chisq_like_cov_xbar(get_diffarray(dxidata,dxisys_list[i_redshiftbin]), covmat)
						chisq_syscor, like   = chisq_like_cov_xbar(XplusY(dxidata,dxisys_list[i_redshiftbin],b=-1), covmat)
						#if i_redshiftbin in [5]:
						#	chisq_nosyscor, chisq_syscor = 0, 0
						chisqs_nosyscor[nowomwstr].append(chisq_nosyscor)  
						chisqs_syscor[nowomwstr].append(chisq_syscor) 
					i_redshiftbin += 1
			X1, X2 = chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
			str1, str2, strA, strB = \
				array_to_str(X1), array_to_str(X2), str(sum(chisq_nosyscor)), str(sum(chisq_syscor))
			nowf.write(str(mumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2+'\n')
			iomw += 1	
			#print om, w, chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
		#print chisqs_nosyscor.keys()
		#sys.exit(' ')
		chisqlist_nosyscor = [ [  sum(chisqs_nosyscor[omwstr(om,w)])    for om in omlist] for w in wlist ]
		chisqlist_syscor =   [ [  sum(chisqs_syscor[omwstr(om,w)])      for om in omlist] for w in wlist ]
		nowf.close()
		if make_plot:
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
