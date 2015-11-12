import copy

execfile(pythonlibPATH+'/Tpcftools_smuintxi_list.py')

smu__intxi__settings_std = {
            'xifunction': intxi_FractionalS,
            'smin':6.0, 'smax':50, 'mumin':0.01, 'mumax':1.0, 'nummubin':41, # In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
               }

smu__intxi__settings__list_all = [smu__intxi__settings_std]

### Different smin, with 19, 20, 21 bins
smu__intxi__settings_std_21bin_mumin0p01 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_21bin_mumin0p01['mumin'] = 0.01
smu__intxi__settings_std_21bin_mumin0p01['nummubin'] = 21
smu__intxi__settings__list_all.append(smu__intxi__settings_std_21bin_mumin0p01)

smu__intxi__settings_std_20bin_mumin0p05 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_20bin_mumin0p05['mumin'] = 0.05
smu__intxi__settings_std_20bin_mumin0p05['nummubin'] = 21
smu__intxi__settings__list_all.append(smu__intxi__settings_std_20bin_mumin0p05)

smu__intxi__settings_std_19bin_mumin0p10 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_19bin_mumin0p10['mumin'] = 0.10
smu__intxi__settings_std_19bin_mumin0p10['nummubin'] = 19
smu__intxi__settings__list_all.append(smu__intxi__settings_std_19bin_mumin0p10)

smu__intxi__settings_std_18bin_mumin0p15 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_18bin_mumin0p15['mumin'] = 0.15
smu__intxi__settings_std_18bin_mumin0p15['nummubin'] = 18
smu__intxi__settings__list_all.append(smu__intxi__settings_std_18bin_mumin0p15)

### Different scales
smu__intxi__settings_std_s5to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s5to50['smin'] = 5.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s5to50)

smu__intxi__settings_std_s4to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s4to50['smin'] = 4.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s4to50)

smu__intxi__settings_std_s3to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s3to50['smin'] = 3.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s3to50)

smu__intxi__settings_std_s2to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s2to50['smin'] = 2.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s2to50)

smu__intxi__settings_std_s1to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s1to50['smin'] = 1.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s1to50)

smu__intxi__settings_std_s0to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s0to50['smin'] = 0.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s0to50)

smu__intxi__settings_std_s10to50 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s10to50['smin'] = 10.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s10to50)

smu__intxi__settings_std_s6to40 = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_s6to40['smax'] = 40.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_s6to40)


### different nummubins
smu__intxi__settings_std_1bin = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_1bin['nummubin'] = 2
smu__intxi__settings__list_all.append(smu__intxi__settings_std_1bin)

smu__intxi__settings_std_1bin_mumin0 = copy.deepcopy(smu__intxi__settings_std_1bin)
smu__intxi__settings_std_1bin_mumin0['mumin'] = 0.0
smu__intxi__settings__list_all.append(smu__intxi__settings_std_1bin_mumin0)

smu__intxi__settings_std_120bin = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_std_120bin['nummubin'] = 120
smu__intxi__settings__list_all.append(smu__intxi__settings_std_120bin)


### Different intxi function
smu__intxi__settings_s2xi = copy.deepcopy(smu__intxi__settings_std)
smu__intxi__settings_s2xi['xifunction'] = int_s_square_xi
smu__intxi__settings__list_all.append(smu__intxi__settings_s2xi)

smu__intxi__settings_s2xi_1bin = copy.deepcopy(smu__intxi__settings_s2xi)
smu__intxi__settings_s2xi_1bin['nummubin'] = 2
smu__intxi__settings__list_all.append(smu__intxi__settings_s2xi_1bin)


### filenames

def smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_std,
                        mufmt='automatic', sfmt='%.2f'):
    ### for mumax we use %.2f; for mumin we check how many digits
    mumin = smu__intxi__settings['mumin']
    if mufmt == 'automatic':
	if round(mumin,5) == round(mumin,2):
	   mustr = '%.2f'%mumin
	elif round(mumin,5) == round(mumin,3):
	   mustr = '%.3f'%mumin
	elif round(mumin,5) == round(mumin,4):
	   mustr = '%.4f'%mumin
	else:
	   mustr = '%.5f'%mumin
    else:
	mustr = '%.2f'%mumin
    return smufile+'__intxis/'+smu__xifunname(smu__intxi__settings['xifunction'])+\
        '.'+str(smu__intxi__settings['nummubin'])+\
        'mu'+mustr+'to'+'%.2f'%smu__intxi__settings['mumax']+'.s'+\
        sfmt%smu__intxi__settings['smin']+'to'+sfmt%smu__intxi__settings['smax']

def smu__intxi_get_muedgemids(smu__intxi__settings = smu__intxi__settings_std):
    mumin = smu__intxi__settings['mumin']
    mumax = smu__intxi__settings['mumax']
    nummubin = smu__intxi__settings['nummubin']
    muedges = np.linspace(mumin, mumax, nummubin)
    return get_mid_array1d(muedges)
def smu__intxi_str(smu__intxi__settings=smu__intxi__settings_std):
	return separate_path_file(smu__intxi_filename('', smu__intxi__settings))[1]

### Something about getting information from smu chisq file name
def smu_chisq__getSky(chisqfile):
    if chisqfile.find('-NGC-') >= 0:
        return 'N'
    elif chisqfile.find('-NSGC-') >= 0:
        return 'NS'
    elif chisqfile.find('-SGC-') >= 0:
        return 'S'
    else:
        return 'NS'
def smu_chisq__getSixBinRedshift(Sky):
    if Sky == 'N':
        SixBinRedshift = SixBinRedshift_N
    elif Sky == 'S':
        SixBinRedshift = SixBinRedshift_S
    else:
        SixBinRedshift = SixBinRedshift_NS
    return [x for x in SixBinRedshift]
def smu_chisq__getrefibin(chisqfile):
    refibin = int(chisqfile[chisqfile.find('refibin')+7])
    return refibin
def smu_chisq__getzref_zcomps(refibin, SixBinRedshift):
    zref = SixBinRedshift[refibin]
    zcomps = []
    izcomps = []
    for i in range(len(SixBinRedshift)):
        if i == refibin: continue
        zcomps.append(SixBinRedshift[i])
        izcomps.append(i)
    return zref, [x for x in zcomps], [x for x in izcomps]

### Compute intxi

def smu__intxi_calcwrite(smufile, smusettings, smu__intxi__settings=smu__intxi__settings_std,
                         mufmt='automatic', sfmt='%.2f', smudata=None):
    if smudata == None:  smudata = smu__loadin(smufile, smusettings)
    DDlist, DRlist, RRlist = Xsfrom2ddata(smudata, [4,5,6])
    smu__initsmusettings(smusettings)
    muedges = np.linspace(smu__intxi__settings['mumin'], smu__intxi__settings['mumax'], smu__intxi__settings['nummubin'])
    i_muedges = smu__get_imuedges(muedges,smusettings);
    muedgemids = smu__get_muedgemids(i_muedges, smusettings);
    intxis = [smu__intxi__settings['xifunction'](DDlist, DRlist, RRlist,  
                        smu__intxi__settings['smin'], smu__intxi__settings['smax'], 
                         i_muedges[imubin][0], i_muedges[imubin][1]) 
                          for imubin in range(len(i_muedges))]
    
    smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings, 
                                          mufmt=mufmt, sfmt=sfmt)
    
    callsys('mkdir -p '+separate_path_file(smu__intxi_file)[0])
    
    f0 = open(smu__intxi_file, 'w')
    for intxi in intxis:
        f0.write(str(intxi)+'\n')
    f0.close()
    return smu__intxi_file, muedgemids, intxis

def smu__intxi_calcwrite_list(smufilelist, smusettings = smusettings_sample, 
				smu__intxi__settings__list = [smu__intxi__settings_std], 
                              checkexit = True, just_check_intxifileexit = False, just_get_intxifilelist=False,
				skip_exited_smufile=True, do_not_preserve_files = False):
    print exms100rt + ' Computing intxi for ',len(smufilelist),' smufiles' + exms100rt
    ### Check existence of smu files
    if checkexit:
        if (not just_check_intxifileexit) and (not just_get_intxifilelist):
            print exmsrt + ' Checking whether all files existed...'
            checkrlt = isfiles(smufilelist)
            if not checkrlt:
                print exmsrt + ' Some files missing: Re-create missing files!'
                return
    ifile = 0
    ijob_skip = 0
    ijob_tot = 0
    smu__intxi_file_list = []
    timebegin = time.time()
    time1 = timebegin
    print exmsrt
    for smufile in smufilelist:
        smufileloaded = False
        ifile += 1
        time2 = time.time()
        if time2 - time1 > 600:
            print '\t', time2 - timebegin, 'seconds elapsed; files processed = ', ifile, \
                '; rat = ', ifile / float(len(smufilelist)), '; skipped job = ', ijob_skip
            print '\t\t', smu__intxi_file
            time1 = time2
        for smu__intxi__settings in smu__intxi__settings__list:
            smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
            if not do_not_preserve_files:
            	smu__intxi_file_list.append(smu__intxi_file)
            ijob_tot += 1
            if just_check_intxifileexit or just_get_intxifilelist:
                ijob_skip += 1
                continue
            if skip_exited_smufile and isfile(smu__intxi_file):
                ijob_skip += 1
                continue
            else:
                if smufileloaded == False:
                    smudata = smu__loadin(smufile, smusettings)
                    smufileloaded = True
            smu__intxi_file, muedgemids, intxis = \
                smu__intxi_calcwrite(smufile, smusettings, smu__intxi__settings=smu__intxi__settings, smudata=smudata)
    if just_check_intxifileexit:
        print exmsrt + ' Checking existence of intxi files...'
        isfiles(smu__intxi_file_list)
        return
    timeend = time.time()
    print exmsrt + ' Done.\n   #-files = ', ifile, '    #-job=', ijob_tot, '    #-job skipped = ', ijob_skip
    print '   Total time taken: ', timeend - timebegin
    print exmsrt
    return smu__intxi_file_list

def time_used_for_calc_intxi(timeload = 0.3, timecalc = 0.02, nsmufile = 5052, nsmusetting = 1):
    '''Roughly estimate time requred to compute intxi'''
    time_in_s = nsmufile * (timeload + timecalc * nsmusetting)
    print time_in_s / 3600.0,  ' hours'
    print time_in_s / 60.0,  ' minutes'


### Load files; plot files;...

def smu__intxi_quickload(intxifile):
    return [float(x) for x in open(intxifile,'r').readlines()]
def smu__intxi__plot_files( intxifilelist = [], mus = None, prefix='smuplot-', savefig=False):
    time1 = time.time()
    intxis_list = [smu__intxi_quickload(nowfile) for nowfile in intxifilelist]
    time2 = time.time()
    if mus == None:
	mus = np.linspace(0,1,len(intxis_list[0]))
    figname = str(len(intxifilelist))+'-files'
    fig, ax = figax(title = figname)
    for i in range(len(intxis_list)):
        ax.plot(mus, intxis_list[i])
    ax.grid()
    plt.show()
    if savefig:
	fig.savefig(prefix+figname+'.png', format='png')

### RSD Correction
import copy
totbin = 3

refibin = 0

def smu__intxi_RSDCovbasename(Sky, catname2, RSDstr, smu__intxi__settings, refibin = 0, 
	totbin=3, normedintxi=True):
	rlt =  Sky+'GC-'+str(totbin*2)+'bins-'+catname2+'-'+RSDstr+'-refibin'+str(refibin)+\
        '--'+smu__intxi_str(smu__intxi__settings)
	if normedintxi == False:
		rlt += '.unnormed'
	return rlt

def smu__intxi_covmatfilename(icovmat, Sky, catname2, RSDstr, smu__intxi__settings, refibin = 0, 
	totbin=3, normedintxi=True, refine_intxi=None):
    covmatname = datamockdir+'/covmats'+'/'+smu__intxi_RSDCovbasename(Sky, catname2, RSDstr, smu__intxi__settings, refibin, totbin, normedintxi)+'.cov'+str(icovmat)
    if refine_intxi != None:
	covmatname += ('.'+refine_intxi)
    return covmatname

def smu__intxi_chisqfilename(Sky, catname2, RSDstr, smu__intxi__settings, refibin = 0, 
	totbin=3, numchisq=10, prefix='', normedintxi=True, refine_intxi=None):
    chisqfilename = datamockdir+'/covmats'+'/'+prefix+smu__intxi_RSDCovbasename(Sky, catname2, RSDstr, smu__intxi__settings, refibin, totbin, normedintxi)+'.'+str(numchisq)+'chisqs'
    if refine_intxi != None:
	chisqfilename += ('.'+refine_intxi)
    return chisqfilename


def smu__inxi_RSDCor_Covmat(Sky = 'N', 
                            catname2='data', 
                            smu__intxi__settings = smu__intxi__settings_std, 
                            totbin=3, 
                            cosmoconvertomw = None,
                            smusettings = smusettings_sample,
                            RSDstr = 'RSD', 
                            makeplot=True, savefig=False, figfmt='png', 
                            docovmat = False, savecovmat=True, 
                            refibin = 0, 
                            normedintxi=True,
                            use_omw_testfun=False, ## using a test function
                            refine_intxi = None,     
                            ):
    
    if Sky == 'N':
        catnames = ['DR12v4-LOWZ-N', 'DR12v4-CMASS-N']
    elif Sky == 'S':
        catnames = ['DR12v4-LOWZ-S', 'DR12v4-CMASS-S']
    elif Sky == 'NS':
        catnames = ['DR12v4-LOWZ', 'DR12v4-CMASS']
    else:
        print 'Error (smu__inxi_RSDCor): Sky must be N, S: ', Sky
        return
    mus = smu__intxi_get_muedgemids(smu__intxi__settings)
    
    if True:
        if makeplot: fig, ax = figax(figxsize=16)
        Tpcfs = []
        labels = []
        
        #AllTpcfs[isambin][imock][iquan]
        
        isambin = 0
        AllTpcfs = []
            
        for catname in catnames:
          for ibin in range(totbin):
            if catname2 == 'data':
                nummock = 1
            else:
                nummock = catinfo_nummock(catname, catname2)
            AllTpcfs.append([])
            labels.append(catname+'--'+str(ibin+1) +'-th bin')
            NowTpcfs = []
            
            for imock in range(nummock):
                if catname2 == 'data':
                    file0 = datafile(catname)
                else:
                    file0 = mockfile(catname, catname2, imock, RSDstr)
                    
                rmax, nbins, mubins = smusettings_get_rmax_nbins_mubins(smusettings)
                if cosmoconvertomw == None:
                    smufile = Tpcfrltfilename(binsplittedfilename(file0, ibin+1), 
                                          mubins=mubins,nbins=nbins,rmax=rmax )
                else:
                    om, w = cosmoconvertomw
                    smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(file0, ibin+1), om, w),
                                          mubins=mubins,nbins=nbins,rmax=rmax )
                smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                #print isfile(smufile)
                #print isfile()
                
                if use_omw_testfun:
                    if cosmoconvertomw == None:
                        om, w = 0.26, -1
                    else:
                        om, w = cosmoconvertomw
                    sigma = 2.0
                    noisesig = 0.1
                    X = [ (((om-0.6)/sigma*1.5)**2.0 + ((w+2.0)/sigma)**2.0) 
                                 * ibin * x * random.uniform(1-noisesig, 1+noisesig) \
                                 + random.uniform(1-noisesig/3, 1+noisesig/3)
                            for x in mus]
                    #X = polyfitY(mus, X, 3)
                    savecovmat = False
                else:
                    X = quickload_1col(smu__intxi_file)
		if normedintxi: X = normto1(X)
		if refine_intxi != None:
			if refine_intxi[0:8] == 'polyfit_':
			    polyfitorder = int(refine_intxi[8:len(refine_intxi)])
			    X = polyfitY(mus, X, polyfitorder)
			elif refine_intxi in ['muge0p15', 'muge0p15_polyfit3']:
				mus, X = smu__intxi_refine_intxi(mus, X, None, refine_intxi)
			elif refine_intxi[0:13] == 'normedpolyfit':
				if refine_intxi[len(refine_intxi)-2] == '_':
					polyfitorder = int(refine_intxi[len(refine_intxi)-1])
				else:
					polyfitorder = int(refine_intxi[len(refine_intxi)-2:len(refine_intxi)])
			        X = polyfitY(mus, X, polyfitorder)
				X = normto1(X)
                NowTpcfs.append([x for x in X])
                if docovmat: 
                    AllTpcfs[isambin].append([x for x in X])
            isambin += 1
            Tpcfs.append(avgarray_2d(NowTpcfs))
        #print NowTpcfs 
        #print exmsrt
        #print avgarray_2d(NowTpcfs)
        #break
        Ys = []
        for i in range(len(Tpcfs)):
            #print Tpcfs[i]
            if i == refibin: continue
            X = mus# = np.linspace(0,1,len(Tpcfs[0]))
            Y = [Tpcfs[i][row] -  Tpcfs[refibin][row] for row in range(len(X))]
            Ys.append([y for y in Y])
            if makeplot: ax.plot(X, Y, label = labels[i] + ' - ' + labels[0], lw=2)
        if makeplot:     
            mocktypestr = 'mock-type = $\\rm '+catname2+',\\ '+RSDstr+\
		'$ (averaged over $'+str(nummock)+'$ realizations)'
            catloguestr = 'Difference among the $'+str(len(Tpcfs))+'$  bins from $\\rm'+\
                array_to_str(catnames,',\\ ')+',$ reference = $'+str(refibin+1)+'$-th bin'
            smuintxistr = separate_path_file(smu__intxi_filename('', smu__intxi__settings))[1]
            integrationstr = 'integration = $\\rm ' + smuintxistr + '$'
            figtilt = catloguestr + '\n' + mocktypestr + '\n' + integrationstr
            ax.set_title(figtilt, fontsize=16)
            ax.grid()
            ax.set_ylim(-0.15, 0.15)
            ax.legend(loc='best', frameon=False)
            fig.subplots_adjust(top=0.85);
            plt.show()
            if savefig:
                figname = Sky+'-'+catname2+'-'+RSDstr+'-'+smuintxistr+'.'+figfmt
                fig.savefig(figname, figfmt=figfmt)

        avgdiffs = Ys
        avgintxis = Tpcfs
        
        if not docovmat:
            return mus, avgdiffs, avgintxis
        else:
            
            #AllTpcfs[isambin][imock][iquan]
            #AllQuans[iquan][imock]
            #AvgDiffs[isambin][iquan]
            
            numsambin = len(AllTpcfs) # sample bin
            numquan = len(AllTpcfs[0][0])
            
            covmats = []
            AvgDiffs = [[0.0 for row2 in range(numquan)] for row1 in range(numsambin)]
            
            for isambin in range(numsambin):
                if isambin == refibin: continue
                AllQuans =  [[0 for row2 in range(nummock)] for row1 in range(numquan)]
                for imock in range(nummock):
                    for iquan in range(numquan):
                        nowdiff = AllTpcfs[isambin][imock][iquan] - AllTpcfs[refibin][imock][iquan]
                        AllQuans[iquan][imock] = nowdiff
                        AvgDiffs[isambin][iquan] += nowdiff / float(nummock)
                covmats.append(get_covmat(AllQuans))
            if savecovmat:
            	for icovmat in range(len(covmats)):
		    nowfile = smu__intxi_covmatfilename(icovmat, Sky, catname2, RSDstr, smu__intxi__settings, refibin = refibin,
                totbin=totbin, normedintxi=normedintxi, refine_intxi=refine_intxi)
		    np.savetxt(nowfile, covmats[icovmat])
		    print '\tfilesaved: ', nowfile
            return mus, avgdiffs, avgintxis, covmats
### Plot intxi files (convenient for checking...)

###########################################3
#### smu__intxi_ChisqContour
###########################################3

### TO run me, you must run these settings first
## 0. sigmas

#sigs = [sig3,  sig5,  sig7]

# 1. Parameters of om,w

#omws = NSscan_omwlist
#scanname = NSscan_scanname
#oms = NSscan_omlist
#ws = NSscan_wlist
#numw=len(ws); numom = len(oms)

# 2. The scheme for 2pcf


## scheme used in covmat, simulation mock
#import copy;
#smusettings_simulation_mock = smusettings_covmat = copy.deepcopy(smusettings_sample)

## scheme used in chisq, real observational data
#smusettings_observational_data = smusettings_chisq = copy.deepcopy(smusettings_smax60)

#if True:
#    print exmsrt
#    print 'Scheme of smu for simulation mock (covmat): \n\t', smusettings_covmat
#    print 'Scheme of smu for observational data (chisq): \n\t', smusettings_smax60
#smu__initsmusettings(smusettings_simulation_mock)
#smu__initsmusettings(smusettings_observational_data)
if not ('oms' in locals() or 'oms' in globals()): 
	oms, ws, omws, scanname =  None, None, None, None 
def smu__intxi_ChisqContour(Skylist=['N'], 
                            smu__intxi__settings_list = [smu__intxi__settings_std],
                            normedintxilist = [True],
                            refibinlist = [0], 
                            compute_covmat = True, covmat_RSDstr = 'noRSD',
                            catname2_RSDCor = 'J08', plot_RSDcor = False,
                            printchisqinfo = False, printchisqroughinfo = False,
                            savecontourfig = False,
                            use_omw_testfun=False,
                            sigs = [0.683, 0.954, 0.997],
                            skip_calc_if_chisqfile_found = True,
                            refine_intxi = None,
			    #oms=oms, ws=ws, omws=omws, scanname=scanname,
                            ):

    if use_omw_testfun:
        scanname = 'omw-Test-Scan'
        oms = np.linspace(0, 1.0, 20)
        ws = np.linspace(-2.5, -0.5, 20)
        numw=len(ws); numom = len(oms)
        omws = sumlist( [[[om,w] for om in oms] for w in ws])
    else:
	global oms, ws, omws, scanname
        numw=len(ws); numom = len(oms)

    for refibin in refibinlist:
     print '* refibin = ', refibin

     for smu__intxi__settings in smu__intxi__settings_list:
      print '** smu__intxi = ', smu__intxi_str(smu__intxi__settings)

      for Sky in Skylist:
       print '*** Sky = ', Sky

       for normedintxi in normedintxilist:
               print '**** normedintxi = ', normedintxi

               ###########################################
               ### Step 3.1 Compute covmat

               print exms
               print '******** Step 3.1 Covmat'

                        ###!!!!!  RSDCor_covmat : 4
	       numchisq = len(omws)
               chisqfilename = smu__intxi_chisqfilename(Sky, catname2_RSDCor, 'RSD', smu__intxi__settings, 
                                                         refibin = refibin, totbin=totbin, 
                                                         numchisq=numchisq, prefix=scanname+'--',
                                                         normedintxi=normedintxi, refine_intxi=refine_intxi)

	       if skip_calc_if_chisqfile_found and isfile(chisqfilename):
			print '\t\tSkip computation: using chisqs loaded from file: ', chisqfilename
			pass
	       else:
                if compute_covmat:
                    print '\t\tcompute covmats...'

                    ###!!!!!  RSDCor_covmat : 1 !!!Very Important!!!
                    ### Key product: covmats

                    mus, avgdiffs, avgintxis, covmats = \
                            smu__inxi_RSDCor_Covmat(Sky, 'HR3', smu__intxi__settings, 
                                        totbin=totbin,
                                        cosmoconvertomw=None,
                                        smusettings=smusettings_simulation_mock, ## smusettings = simulation_mock
                                        RSDstr=covmat_RSDstr, ## RSDstr = covmat_RSDstr
                                        makeplot=False, 
                                        docovmat=True, savecovmat=True,
                                        refibin = refibin,
                                        normedintxi=normedintxi,
                                        use_omw_testfun=use_omw_testfun,
					refine_intxi=refine_intxi)

                else: # load in rather than do computation
                    covmats = []
                    for icovmat in range(totbin*2-1):
                        ###!!!!!  RSDCor_covmat : 2
                        covmatfilename = smu__intxi_covmatfilename(icovmat, Sky, 'HR3', 
                                    RSDstr = covmat_RSDstr, 
                                    smu__intxi__settings = smu__intxi__settings, 
                                    refibin = refibin,
                                    totbin = totbin,
                                    normedintxi = normedintxi,
                                    refine_intxi=refine_intxi)


                        print '\t\tload in covmat from:    ', separate_path_file(covmatfilename)[1]
                        covmats.append(loadtxt_rand(covmatfilename, rat=1.1))

                ###########################################
                ### Step 3.2 Calc RSD Cor

                print exms
                print '******** Step 3.2 RSD Correction'


                        ###!!!!!  RSDCor_covmat : 3 !!!Very Important!!!
                        ### Key product: avgdiffs_simulation_mock

                mus, avgdiffs_simulation_mock, avgintxis_simulation_mock = \
                    smu__inxi_RSDCor_Covmat(Sky, catname2_RSDCor, smu__intxi__settings, 
                                        totbin=totbin,
                                        cosmoconvertomw=None,
                                        smusettings=smusettings_simulation_mock, ## smusettings = simulation_mock
                                        RSDstr='RSD', 
                                        makeplot=plot_RSDcor, savefig=False, 
                                        docovmat = False,
                                        refibin = refibin,
                                        normedintxi = normedintxi,
                                        use_omw_testfun=use_omw_testfun,
                                        refine_intxi=refine_intxi)

                ## Use 3-rd order polynomial to fit the RSD cor???

                polyfitdeg = 3

                avgdiffs_simulation_mock_polyfit = [polyfitY(mus, Y, polyfitdeg) for Y in avgdiffs_simulation_mock]
                if plot_RSDcor:
                    fig, ax = figax(16,8, title='RSD correction (solid) and its '+str(polyfitdeg)+\
                                    '-rd polynomial fitted curve (dashed)')
                    for iY in range(len(avgdiffs_simulation_mock)):
                        nowcolor = PLOT_COLOR_ARRAY [iY]
                        ax.plot(mus, avgdiffs_simulation_mock[iY], lw=1, c=nowcolor)
                        ax.plot(mus, avgdiffs_simulation_mock_polyfit[iY], lw=3, c=nowcolor, ls='--', label=str(iY)+'th curve')
                    ax.legend(loc='best', fontsize=16, frameon=False)
                    ax.grid()
                    plt.show()

                ###########################################
                ### Step 3.3 Compute Chisq

                print exms
                print '******** Step 3.3 Compute Chisq'

                noRSDCorchisqlist = []
                RSDCorchisqlist = []
                noRSDCorchisq_sep_list = []
                RSDCorchisq_sep_list = []


                ### Begin of omw loop
                for omw in omws:
                    om, w = omw

                            ###!!!!!  RSDCor_covmat : 5 !!!Very Important!!!
                            ### Key product: avgdiffs_observational_data
                    mus, avgdiffs_observational_data, avgintxis_observational_data = \
                        smu__inxi_RSDCor_Covmat(Sky, 'data', smu__intxi__settings, 
                                        totbin = totbin,
                                        cosmoconvertomw=[om, w],
                                        smusettings=smusettings_observational_data, ## smusettings = observational_data
                                        RSDstr='RSD', ## in fact it does not matter
                                        makeplot=False,  savefig=False, 
                                        docovmat=False,
                                        refibin = refibin, 
                                        normedintxi=normedintxi,
                                        use_omw_testfun=use_omw_testfun,
                                        refine_intxi=refine_intxi)  


                    noRSDCorchisqs = []
                    RSDCorchisqs = []

                    if printchisqinfo:
                        print '\t', omwstr(om, w), ':'

                    ## chisq in many bins
                    for iY in range(len(avgdiffs_simulation_mock)):
                        ## diff
                        noRSDCorX = avgdiffs_observational_data[iY]
                        RSDCorX = XplusY(avgdiffs_observational_data[iY], avgdiffs_simulation_mock_polyfit[iY], b=-1)
                        ## cov
                        covmat = covmats[iY]
                        ## compute  chisq
			if refine_intxi in ['Wei1']:
				noRSDCorX = smu__intxi_refine_intxi(mus, noRSDCorX, covmat, refine_intxi)
				RSDCorX = smu__intxi_refine_intxi(mus, RSDCorX, covmat, refine_intxi)
			#elif refine_intxi in ['muge0p15', 'muge0p15_polyfit3']:
			#	mus, noRSDCorX = smu__intxi_refine_intxi(mus, noRSDCorX, covmat, refine_intxi)
			#	mus, RSDCorX = smu__intxi_refine_intxi(mus, RSDCorX, covmat, refine_intxi)
			if False:
                        	noRSDCorchisq, noRSDCorlike = chisq_like_cov_xbar(noRSDCorX, covmat)
                        	RSDCorchisq, RSDCorlike = chisq_like_cov_xbar(RSDCorX, covmat)
			else:
				npcov = np.mat(covmat); invcov = npcov.I;
				noRSDCorchisq = chisq_invcov_xbar(noRSDCorX, invcov)
				RSDCorchisq = chisq_invcov_xbar(RSDCorX, invcov)
                        noRSDCorchisqs.append(noRSDCorchisq); RSDCorchisqs.append(RSDCorchisq); 
                    ## summation of chisq
                    noRSDCorchisq = sum(noRSDCorchisqs)
                    RSDCorchisq = sum(RSDCorchisqs)
                    if printchisqinfo:
                        print '\t\t  noRSDCor-chisqs = %.4f'%noRSDCorchisq, '(',\
                            array_to_str_fmtted(noRSDCorchisqs, fmt='%.4f'), ')'
                        print '\t\t  RSDCor-chisqs   = %.4f'%RSDCorchisq, '(', \
                            array_to_str_fmtted(RSDCorchisqs, fmt='%.4f'), ')'

                    noRSDCorchisqlist.append(noRSDCorchisq)
                    RSDCorchisqlist.append(RSDCorchisq)

                    noRSDCorchisq_sep_list.append([x for x in noRSDCorchisqs])
                    RSDCorchisq_sep_list.append([x for x in RSDCorchisqs])

                ### End of of omw loop

                chisqrlt = [ [omws[iomw][0], omws[iomw][1], noRSDCorchisqlist[iomw], RSDCorchisqlist[iomw]] \
				+ noRSDCorchisq_sep_list[iomw] + RSDCorchisq_sep_list[iomw] \
                           	 for iomw in range(len(omws))]
                print '\t\tResult saving to file: ', chisqfilename
                np.savetxt(chisqfilename, chisqrlt, fmt='%.7f')
                if printchisqroughinfo:
                    print_file(chisqfilename)

               ###########################################
               ## Step 3.4 Plot Contour
	       if True:
                print exms
                print '******** Step 3.3 Plot Contour'

                # Load in data
                X1, X2, X3, X4 = Xsfromdata(loadtxt_rand(chisqfilename, rat=1.1), [0,1,2,3])
                X3 = [x for x in X3]
                X4 = [x for x in X4]
                noRSDCor_2dchisqlist = get_2darray_from_1d(X3, numw, numom)
                RSDCor_2dchisqlist = get_2darray_from_1d(X4, numw, numom)

                # Plot Contour
                fig = plt.figure(figsize=(24,6))
                ax1 = fig.add_subplot(121); ax2 = fig.add_subplot(122)

                contourlabel = ''
                for isig in range(len(sigs)-1):
                    nsig = nsigofsigCL(sigs[isig])
                    contourlabel += ('%g'%(round(nsig, 1))+', ')
                contourlabel += ('%g'%(round(nsigofsigCL(sigs[len(sigs)-1]), 1))+'')
                contourlabel += ' $\\sigma$ contours'

                plot_contour(ax1, oms, ws, noRSDCor_2dchisqlist, do_smooth=False, sigs=sigs,
                             ommin=min(oms),  ommax=max(oms), wmin=min(ws), wmax=max(ws), plotformat=1,  
                                 sigA=1,sigB=2,sigC=3,
                                 label=contourlabel+'\nBefore RSD Correction', legftsize=20)

                plot_contour(ax2, oms, ws, RSDCor_2dchisqlist, do_smooth=False, sigs=sigs, 
                             ommin=min(oms),  ommax=max(oms), wmin=min(ws), wmax=max(ws), plotformat=1,  
                                 sigA=1,sigB=2,sigC=3,
                                 label=contourlabel+'\nAfter RSD Correction', legftsize=20)
                ax1.grid(); ax2.grid()
                fig.suptitle(separate_path_file(chisqfilename)[1], fontsize=24)
                if savecontourfig:
                    fig.savefig(separate_path_file(chisqfilename)[1]+'.contour.png', format='png')
                plt.show()

    if use_omw_testfun and savecontourfig:            
        fig.suptitle('To test the code, I use some mathematical curves as xi(mu) in six redshift bins\n'+
                     'They are designed to have no redshift evolution at om/w = 0.6/-2.0\n'+
                     'The code can reproduce the properties of these functions (shown in left)\n'+
                     'And, after correction I found it can return to 0.26, -1 (shown in right)', fontsize=18)
        fig.subplots_adjust(top=0.75)
        fig.savefig('Check_Correctness_of_Code.png', format='png')


###########################################3
#### Distance Prior FItting
###########################################3

############################################
### List of Distance Priors

### Now only have DA * H; 

def DistancePrior__DA_pro_H__rat(om, w, zcomp, zref):
    fcomp = DA_pro_H(om, w, 0.6, zcomp)
    fref = DA_pro_H(om, w, 0.6, zref)
    return fcomp / fref

def DistancePrior__DA_pro_H__sqrat(om, w, zcomp, zref):
    fcomp = DA_pro_H(om, w, 0.6, zcomp)
    fref = DA_pro_H(om, w, 0.6, zref)
    return (fcomp / fref)**2.0

def DistancePrior__DA_pro_H__diff(om, w, zcomp, zref):
    fcomp = DA_pro_H(om, w, 0.6, zcomp)
    fref = DA_pro_H(om, w, 0.6, zref)
    return abs(fcomp - fref)

def DistancePrior__rz__rat(om, w, zcomp, zref):
    fcomp = comov_r(om, w, 0.6, zcomp)
    fref = comov_r(om, w, 0.6, zref)
    return fcomp / fref

def DistancePrior__Hz__rat(om, w, zcomp, zref):
    fcomp = Hz(om, w, 0.6, zcomp)
    fref = Hz(om, w, 0.6, zref)
    return fcomp / fref

def DistancePrior_name(DP, zcomp, zref, zfmt = '%.3f', nametype='plot'):
    if DP == DistancePrior__DA_pro_H__rat:
        if nametype == 'plot':
            fname = 'D_A H'
            return '${D_A H('+zfmt%zcomp+')}\\ /\\ {D_A H('+zfmt%zref+')}$'
        elif nametype == 'txt':
            return 'DA_pro_H__rat'
    elif DP == DistancePrior__DA_pro_H__sqrat:
        if nametype == 'plot':
            fname = 'D_A H'
            return '$({D_A H('+zfmt%zcomp+')}\\ /\\ {D_A H('+zfmt%zref+')})^2$'
        elif nametype == 'txt':
            return 'DA_pro_H__sqrat'
    elif DP == DistancePrior__DA_pro_H__diff:
        if nametype == 'plot':
            fname = 'D_A H'
            return '${D_A H('+zfmt%zcomp+')}\\ - \\ {D_A H('+zfmt%zref+')}$'
        elif nametype == 'txt':
            return 'DA_pro_H__diff'
    elif DP == DistancePrior__rz__rat:
        if nametype == 'plot':
            fname = 'r'
            return '${r('+zfmt%zcomp+')}\\ /\\ {r('+zfmt%zref+')}$'
        elif nametype == 'txt':
            return 'r__rat'
    elif DP == DistancePrior__Hz__rat:
        if nametype == 'plot':
            fname = 'H'
            return '${H('+zfmt%zcomp+')}\\ /\\ {H('+zfmt%zref+')}$'
        elif nametype == 'txt':
            return 'H__rat'
    else:
        return '$Distance Prior (zcomp='+zfmt%zcomp+', zref = '+zfmt%zref+')}$'
    
def smu_chisq__auto_omws(chisqdata):
    omws = [[x[0], x[1]] for x in chisqdata]
    oms, ws = XYfromdata(chisqdata)
    oms = merge_two_list_finitedigits(oms, [])
    ws = merge_two_list_finitedigits(ws, [])
    return oms, ws, omws

def smu_chisq__DPfitted(om, w, DistancePrior, zref, zcomps, Polys):
                    chisq = 0
                    num_sep_chisq = len(Polys)
                    for ichisq in range(num_sep_chisq):
                        zcomp = zcomps[ichisq]
                        DP = DistancePrior(om, w, zcomp, zref)
                        #f0, sig, chisq0 = noRSDCor_f0sigchisq0_all[ichisq]; 
                        #noRSDCorchisq += chisq__f0_sig_chisq0(DP, f0, sig, chisq0)
                        #f0, sig, chisq0 = RSDCor_f0sigchisq0_all[ichisq]; 
                        #RSDCorchisq += chisq__f0_sig_chisq0(DP, f0, sig, chisq0)
                        chisq += polyval(Polys[ichisq], DP)
                    return chisq
def smu_chisq__DPfitted_list(oms, ws, DistancePrior, zref, zcomps, Polys):
                    chisqs = []
                    for iomw in range(len(oms)):
                        om, w = oms[iomw], ws[iomw]
                        chisqs.append( smu_chisq__DPfitted(om, w, DistancePrior, zref, zcomps, Polys) )
                    return chisqs

import stdA; execfile(stdA.pyfile)

def smu__intxi_DPfit(chisqfile, multifiles = False, print_all_filename=False,
		     nowdir = '/home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/covmats/', 
                     DistancePrior = DistancePrior__DA_pro_H__rat,
                     polyfitdegs = [6], truncated_chisq_fit=None,
                     sigs = [sig1, sig2, sig3],
                     chisq__calibrating_for_diff_omw = None, 
			chisq__calibrating_for_diff_omw__dft_omw = [0.26,-1.0],
                     Plot_DP_chisq = True, labelfs = 8, legfs = 16, save_DP_chisq_curve = False, figfmt='pdf',
                         DP_chisq_figxsize=16, DP_chisq_figysize=3.5, figtiltfs=18,
		     do_DPfit_RandomWalk_search = None, given_SixBinRedshift=None,
                     plot_tot_fitted_chisq_contour = True, save_tot_contour = False, plot_fitted_rlt=True,
                     print_fittingrlt = False, only_print_f0sigchisq0 = True, 
                     print_convenient_sum = False,
                     check_tot_sep_match = False, 
                     saverlt = True,
                     EnhancedOmWGrid=[50,50],
                     ):
     DPname_txt = DistancePrior_name(DistancePrior, 0, 0, nametype='txt' )
     if (not multifiles) or print_all_filename:
	     print '* chisqfile = ', chisqfile
     else:
	     print '* chisqfile = ', len(chisqfile), ' files:  \n\t\t\t', chisqfile[0], \
			'\n\t\t\t\t...\n\t\t\t', chisqfile[len(chisqfile)-1]
     if not multifiles:
         Sky = smu_chisq__getSky(chisqfile)
         refibin = smu_chisq__getrefibin(chisqfile)
     else:
         Sky = smu_chisq__getSky(chisqfile[0])
         refibin = smu_chisq__getrefibin(chisqfile[0])
                
     if given_SixBinRedshift == None:
	     SixBinRedshift = smu_chisq__getSixBinRedshift(Sky)
     else:
	     SixBinRedshift = given_SixBinRedshift
     SixBinRedshift_orig = SixBinRedshift

     zref, zcomps, izcomps = smu_chisq__getzref_zcomps(refibin, SixBinRedshift_orig)
     DPname_general = DistancePrior_name(DistancePrior, 'z', '%.3f'%zref, '%s')
     print '\t\t ### Sky = ', Sky, '  ### refibin = ', refibin, '  ### zref = %.3f'%zref, \
        '  zcomps = ', array_to_str_fmtted(zcomps, div=' ', fmt='%.3f')
        
     if not multifiles:
                chisqfilename = nowdir + chisqfile
                chisqdata = loadtxt_rand(chisqfilename)
     else:
                chisqdatas = []
                iskip = 0
                for nowfile in chisqfile:
                    try:
                        chisqfilename = nowdir + nowfile
                        chisqdatas.append(loadtxt_rand(chisqfilename))
                    except:
                        iskip += 1
                        pass
                chisqdata = get_avg_array2D(chisqdatas)
                chisqfilename += ('.'+str(len(chisqfile))+'files')
                if iskip != 0:
                    chisqfilename += ('.'+str(iskip)+'skiped')

     if do_DPfit_RandomWalk_search == None:
	RW_max_step = 1
	RW_min_BC = 1
     else:
	RW_max_step, RW_min_BC = do_DPfit_RandomWalk_search

     for polyfitdeg in polyfitdegs:
      print '** polyfit degree = ', polyfitdeg

      ### RW
      RW_BC = 0.01
      SixBinRedshift_bf = SixBinRedshift_orig
      RW_chisq_min = 1.0e10
      if do_DPfit_RandomWalk_search == None:
	RW_break = True
      else:
	RW_break = False

      for RW_step in range(1,RW_max_step+1):
            ### RW
	    if do_DPfit_RandomWalk_search != None and RW_step > 1:
		SixBinRedshift = [x + random.normal(0, RW_BC) for x in SixBinRedshift_bf]
     		zref, zcomps, izcomps = smu_chisq__getzref_zcomps(refibin, SixBinRedshift)
		if RW_step == RW_max_step or RW_BC < RW_min_BC:
			RW_break = True

            ################################################
            ################################################
            ################################################
            ### Step 1. chisqs
            oms, ws, omws = smu_chisq__auto_omws(chisqdata)
            numw, numom = len(ws), len(oms)

            num_sep_chisq = len(chisqdata[0]) / 2 - 2
            oms_all, ws_all = Xsfromdata(chisqdata, [0,1])
            numomw_all = len(oms_all)

	    ### Step 1.0 calibrate the chisq for different om, w 
		### simply because we use the covariance matrix compuated at default om, w for all sets of om, w
		###  but that's not true
		### so we estimate "how much the error bar at a different set of om, w is wrongly estimated"
		###  and use this to correct the chisq values
            if chisq__calibrating_for_diff_omw == None:
		pass
            elif chisq__calibrating_for_diff_omw == 'simple_volume':
		omdft, wdft = chisq__calibrating_for_diff_omw__dft_omw = [0.26,-1.0]
		for ichisq in range(num_sep_chisq):
			zcomp = zcomps[ichisq];
#			noRSDCor_sep_chisqs, RSDCor_sep_chisqs = Xsfromdata(chisqdata, [col1, col2])
                	col1 = 4+ichisq;   col2 = col1 + num_sep_chisq
			for iomw in range(len(chisqdata)):
				om = oms_all[iomw]; w = ws_all[iomw];
				volcomp = comov_r(om, w, 0.7, zcomp) ** 2.0 / Hz_Mpc_to_h(om, w, 0.7, zcomp)
				volref = comov_r(om, w, 0.7, zref) ** 2.0 / Hz_Mpc_to_h(om, w, 0.7, zref)

				volcomp_dft = comov_r(omdft, wdft, 0.7, zcomp) ** 2.0 / Hz_Mpc_to_h(omdft, wdft, 0.7, zcomp)
				volref_dft = comov_r(omdft, wdft, 0.7, zref) ** 2.0 / Hz_Mpc_to_h(omdft, wdft, 0.7, zref)
				## We assume sig^2 = sigzcomp^2 + sigzref^2 ~ 1/Npair_zcomp + 1/Npair_zref
				### we assume Npair ~ 1/vol (larger volume means lower density, smaller number of galaxy pairs)
				## chisq ~ 1/sig^2
				## e.g. small omegam -> large volume -> ersq is actually larger than the ersq we use (which is adopted as the ersq in the default cosmology) -> chisq shall be smaller
				ersq = volcomp + volref
				ersq_dft = volcomp_dft + volref_dft
				ersq_rat = ersq / ersq_dft
				chisq_rat = ersq_dft / ersq 
				#chisq_rat = volcomp_dft / volcomp + volref_dft / volref
				#chisq_rat = chisq_rat ** 0.3
				if ichisq == 0:
					#print '\t\tichisq = ', ichisq, '; zref, zcomp = ', zref, zcomp, \
					#	'; om, w =', om, w, \
					#	'; chisq rat = ', chisq_rat
					pass
				chisqdata[iomw][col1] *= chisq_rat
				chisqdata[iomw][col2] *= chisq_rat
            for iomw in range(len(chisqdata)):
		chisqdata[iomw][2] = sum([chisqdata[iomw][4+ichisq] for ichisq in range(num_sep_chisq)])
		chisqdata[iomw][3] = sum([chisqdata[iomw][4+ichisq+num_sep_chisq] for ichisq in range(num_sep_chisq)])



            noRSDCor_tot_chisqs, RSDCor_tot_chisqs = Xsfromdata(chisqdata, [2,3])

            noRSDCor_2dchisqlist = get_2darray_from_1d(noRSDCor_tot_chisqs, numw, numom)
            RSDCor_2dchisqlist = get_2darray_from_1d(RSDCor_tot_chisqs, numw, numom)

            print '\t\t #-omw  = ', numomw_all, '; #-sep_chisq = ', num_sep_chisq

            ################################################
            ################################################
            ################################################
            ### Step 2. Get sep_chisqs, do polynomial fit

            DPs_all = []; DPname_all = []
            noRSDCor_sep_chisqs_all = []; RSDCor_sep_chisqs_all = []
            noRSDCorpoly_all = []; RSDCorpoly_all = []
            noRSDCor_f0sigchisq0_all = []; RSDCor_f0sigchisq0_all = []

            if Plot_DP_chisq and RW_break: 
                print exmsrt
                fig = plt.figure(figsize=(DP_chisq_figxsize,DP_chisq_figysize*num_sep_chisq));
                figtilt = separate_path_file(chisqfilename)[1]+'\nFitting Result of '+DPname_general
                fig.suptitle(figtilt, fontsize=figtiltfs)

	    RW_chisq = 0.0
            for ichisq in range(num_sep_chisq):

                ################################################
                ## Step 2.1: Load in sep_chisq

                col1 = 4+ichisq;   col2 = col1 + num_sep_chisq
                noRSDCor_sep_chisqs, RSDCor_sep_chisqs = Xsfromdata(chisqdata, [col1, col2])
                noRSDCor_sep_chisqs_all.append([x for x in noRSDCor_sep_chisqs]); 
                RSDCor_sep_chisqs_all.append([x for x in RSDCor_sep_chisqs])

                ################################################
                ## Step 2.2: redshift, DP
                zcomp = zcomps[ichisq]; izcomp = izcomps[ichisq]
                DPname = DistancePrior_name(DistancePrior, zcomp, zref)
                DPs = []
                for iomw in range(numomw_all):
                    om = oms_all[iomw]; w = ws_all[iomw];
                    DP = DistancePrior(om, w, zcomp, zref); DPs.append(DP)
                DPname_all.append(DPname)
                DPs_all.append([DP for DP in DPs])

                ################################################
                ## Step 2.3: Polynomial fitted : chisq as a function of DP
                ################################################
                ### Step 2.3.1: f0, sig, chisq0
                noRSDCorpoly = polyfit(DPs, noRSDCor_sep_chisqs, deg=2);
                RSDCorpoly = polyfit(DPs, RSDCor_sep_chisqs, deg=2)

		if truncated_chisq_fit != None:
	                a,b,c = noRSDCorpoly; 
	                [f0, sig, chisq2 ] = abc_to_f0sigchisq0(a, b, c)
	                a,b,c = RSDCorpoly; 
	                [f0, sig, chisq3 ] = abc_to_f0sigchisq0(a, b, c)
			x2s, y2s, x3s, y3s = [], [], [], []
			for ixy in range(len(DPs)):
				x2, y2, x3, y3 = DPs[ixy], noRSDCor_sep_chisqs[ixy], \
					DPs[ixy], RSDCor_sep_chisqs[ixy]
				if y2 <= chisq2 + truncated_chisq_fit:
					x2s.append(x2); y2s.append(y2)
				if y3 <= chisq3 + truncated_chisq_fit:
					x3s.append(x3); y3s.append(y3)
			if len(x2s) <= 3:
				x2s, y2s = DPs, noRSDCor_sep_chisqs
			if len(x3s) <= 3:
				x3s, y3s = DPs, RSDCor_sep_chisqs
	                noRSDCorpoly = polyfit(x2s, y2s, deg=2);
	                RSDCorpoly = polyfit(x3s, y3s, deg=2)

                a,b,c = noRSDCorpoly; 
                noRSDCor_f0sigchisq0 = [f0, sig, chisq0 ] = abc_to_f0sigchisq0(a, b, c)
                noRSDCor_labadd = '\n$\\'+chisq__f0_sig_chisq0_str(f0, sig, chisq0)+'$'
                noRSDCor_f0sigchisq0_all.append([f0, sig, chisq0 ])
                a,b,c = RSDCorpoly; 
                RSDCor_f0sigchisq0 = [f0, sig, chisq0 ] = abc_to_f0sigchisq0(a, b, c)
                RSDCor_labadd = '\n$\\'+chisq__f0_sig_chisq0_str(f0, sig, chisq0)+'$'
                RSDCor_f0sigchisq0_all.append([f0, sig, chisq0 ])

                ################################################
                ### Step 2.3.2: Polynomial fit
		if truncated_chisq_fit == None:
	                noRSDCorpoly = polyfit(DPs, noRSDCor_sep_chisqs, deg=polyfitdeg)
        	        RSDCorpoly = polyfit(DPs, RSDCor_sep_chisqs, deg=polyfitdeg)
		else:
	                noRSDCorpoly = polyfit(x2s, y2s, deg=polyfitdeg);
	                RSDCorpoly = polyfit(x3s, y3s, deg=polyfitdeg)
                noRSDCorpoly_all.append([x for x in noRSDCorpoly])
                RSDCorpoly_all.append([x for x in RSDCorpoly])

		### RW
		for iDP in range(len(DPs)):
			RW_chisq += ((RSDCor_sep_chisqs[iDP] - polyval(RSDCorpoly, DPs[iDP])))**2.0

                ################################################
                ### Step 2.3.3: Plot DP .vs. chisq and the polynomial fitted result
                if Plot_DP_chisq and RW_break:
                    ### Step 2.3.3.1. axs
                    ax1 = fig.add_subplot(num_sep_chisq,2,ichisq*2+1); 
                    ax2 = fig.add_subplot(num_sep_chisq,2,ichisq*2+2)
                    ### Step 2.3.3.2. Scatter Plot of DP, chisq
                    ax1.scatter(DPs, noRSDCor_sep_chisqs);
                    ax2.scatter(DPs, RSDCor_sep_chisqs);
                    ### Step 2.3.3.3. Line Plot of DP, poly-fitted-chisq
                    X = np.linspace(min(DPs), max(DPs), 100)
                    Y1 = [polyval(noRSDCorpoly, x) for x in X]
                    Y2 = [polyval(RSDCorpoly, x) for x in X]
                    ### label with polyfit str
                    #label1 = 'No RSD Cor \n$'+polystr(noRSDCorpoly, polyfitdeg)+'$'+noRSDCor_labadd
                    #label2 = 'After RSD Cor \n$'+polystr(RSDCorpoly, polyfitdeg)+'$'+RSDCor_labadd
                    label1 = '   No RSD Cor, $z = %.3f'%zcomp+'$'+noRSDCor_labadd
                    label2 = 'After RSD Cor, $z = %.3f'%zcomp+'$'+RSDCor_labadd
                    ax1.plot(X, Y1, lw=2, label=label1)
                    ax2.plot(X, Y2, lw=2, label=label2)
                    ### Step 2.3.3.4. extra settings and elements
                    DPWMAP5 = DistancePrior(0.26, -1.0, zcomp, zref)
                    DPPlanck = DistancePrior(0.31, -1.0, zcomp, zref)
                    DPw1 = DistancePrior(0.26, -0.9, zcomp, zref)
                    DPw2 = DistancePrior(0.26, -1.1, zcomp, zref)
                    for ax in [ax1, ax2]:
                        if ax == ax1:
                            [f0, sig, chisq0 ] = noRSDCor_f0sigchisq0
                        else:
                            [f0, sig, chisq0 ] = RSDCor_f0sigchisq0
                        ax.set_xlabel(DPname, fontsize=labelfs)
                        ax.set_ylabel('$\\chi^2$', fontsize=labelfs)
                        nowylim = ax.get_ylim();
                        if ichisq == 0:
                            ax.plot([DPWMAP5,DPWMAP5], nowylim, lw=3, c='g', label='WMAP5')
                            ax.plot([DPPlanck,DPPlanck], nowylim, lw=3, c='k',label='Planck')
                            ax.plot([f0-sig, f0-sig], nowylim, lw=3, ls='--', c='b', label='1$\\sigma$ range')
                            ax.plot([f0+sig, f0+sig], nowylim, lw=3, ls='--', c='b')
                            ax.legend(loc='best', frameon=False, fontsize=legfs)
                        else:
                            ax.plot([DPWMAP5,DPWMAP5], nowylim, lw=3, c='g', )
                            ax.plot([DPPlanck,DPPlanck], nowylim, lw=3, c='k')
                            ax.plot([f0-sig, f0-sig], nowylim, lw=3, ls='--', c='b')
                            ax.plot([f0+sig, f0+sig], nowylim, lw=3, ls='--', c='b')
                            ax.legend(loc='best', frameon=False, fontsize=legfs+5)
                        ax.set_ylim(nowylim)
                        ax.grid()
            if Plot_DP_chisq and RW_break:
                fig.subplots_adjust(top=0.95)
                #fig.tight_layout()
                plt.show()
                if save_DP_chisq_curve:
                    figname = chisqfilename+'.'+DPname_txt+'.scatter_fit_curve.'+figfmt
                    print 'Figure saved: ', figname
                    fig.savefig(figname, fmt=figfmt)
                    
                    
            ################################################
            ################################################      
            ################################################
            ## Step 3: check matching between tot - sep chisqs
            if check_tot_sep_match:
                diffchisqs = []
                for iomw in range(numomw_all):
                    chisq1 = noRSDCor_tot_chisqs[iomw]
                    chisq2 = RSDCor_tot_chisqs[iomw]
                    chisq1B = sum(noRSDCor_sep_chisqs_all[ichisq][iomw] for ichisq in range(num_sep_chisq))
                    chisq2B = sum(RSDCor_sep_chisqs_all[ichisq][iomw] for ichisq in range(num_sep_chisq))
                    diffchisqs.append(abs(chisq1 - chisq1B))
                    diffchisqs.append(abs(chisq2 - chisq2B))
                print '\t\tCheck matching betwee tot, sep chisq: max diff = ', max(diffchisqs)

                
            ################################################
            ################################################    
            ################################################
            ## Step 4: Print the final result
            if print_fittingrlt:
                print exmsrt
                print 'Fitting result of quantity: '
                print '\t\t f \t= ', DPname_general
                print '\t fitting order \t= ', polyfitdeg
                for ichisq in range(num_sep_chisq):
                    print
                    print '\tzcomp, zref = ', zcomps[ichisq], zref
                    if not only_print_f0sigchisq0:
                        print '\tnoRSDCor-polynomial: '
                        print '\t\tcoef = ', noRSDCorpoly_all[ichisq]
                        print '\t\tchisq = ', polystr(noRSDCorpoly_all[ichisq], polyfitdeg, valstr='f')
                    f0, sig, chisq0 = noRSDCor_f0sigchisq0_all[ichisq]
                    print '\tnoRSDCor-f0sigchisq0: '
                    print '\t\t', chisq__f0_sig_chisq0_str(f0, sig, chisq0)
                    if not only_print_f0sigchisq0:
                        print '\tRSDCor-polynomial: '
                        print '\t\tcoef = ', RSDCorpoly_all[ichisq]
                        print '\t\tchisq = ', polystr(RSDCorpoly_all[ichisq], polyfitdeg, valstr='f')
                    f0, sig, chisq0 = RSDCor_f0sigchisq0_all[ichisq]
                    print '\tRSDCor-f0sigchisq0: '
                    print '\t\t', chisq__f0_sig_chisq0_str(f0, sig, chisq0)
            if print_convenient_sum:
                    print '###################################'
                    print 'Convenient summary: '
                    print 'DP = ', DPname_general, '; zref = ',zref,  ';   zcomps = ', zcomps
                    print ' noRSDCor_f0sigchisq0_all = ', noRSDCor_f0sigchisq0_all
                    print ' RSDCor_f0sigchisq0_all = ', RSDCor_f0sigchisq0_all
                    print '###################################'
                    print ' noRSDCorpoly_all = ', noRSDCorpoly_all
                    print ' RSDCorpoly_all = ', RSDCorpoly_all
                    print '###################################'
            if saverlt:
                rltfilename = chisqfilename+'.'+str(polyfitdeg)+\
                    'rdpolyfit_by_'+DistancePrior_name(DistancePrior, 0.0, 0.0, nametype='txt')
                nowf = open(rltfilename, 'w')
                nowf.write('### zref = '+str(zref)+'; Distance Prior = '+\
                           DPname_general+'; zcomps = '+array_to_str(zcomps)+'\n')
                for ichisq in range(num_sep_chisq):
                    nowarray = [zcomps[ichisq], '   ']
                    nowarray = nowarray+noRSDCor_f0sigchisq0_all[ichisq]+['       ']
                    nowarray = nowarray+RSDCor_f0sigchisq0_all[ichisq]+['               ']
                    nowarray = nowarray+noRSDCorpoly_all[ichisq]+['       ']
                    nowarray = nowarray+RSDCorpoly_all[ichisq]
                    nowf.write(array_to_str(nowarray)+'\n')
                nowf.close()
		if do_DPfit_RandomWalk_search == None:
	                print 'Fitting result saving to :'
        	        print '\t\t', separate_path_file(rltfilename)[1]

            ################################################
            ################################################        
            ################################################
            ## Step 5: Make a contour plot
            if plot_tot_fitted_chisq_contour:
                ################################################
                ## Step 5.1: poly-fitted chisq vals
                if EnhancedOmWGrid == None:
                    noRSDCor_tot_chisqs_polyfitted = \
                        smu_chisq__DPfitted_list(oms_all, ws_all, DistancePrior, zref, zcomps, noRSDCorpoly_all)
                    RSDCor_tot_chisqs_polyfitted = \
                        smu_chisq__DPfitted_list(oms_all, ws_all, DistancePrior, zref, zcomps, RSDCorpoly_all)
                    oms2, ws2, numom2, numw2 = oms, ws, numom, numw
                else:
                    #global oms2, ws2, omws2, oms_all2, ws_all2, oms, ws, omws, oms_all, ws_all
                    numom2, numw2 = EnhancedOmWGrid
                    oms2 = np.linspace(min(oms), max(oms), numom2)
                    ws2 = np.linspace(min(ws), max(ws), numw2)
                    omws2 = omws_from_om_w(oms2, ws2, TwoDlist=False)
                    oms_all2, ws_all2 = XYfromdata(omws2)
                    noRSDCor_tot_chisqs_polyfitted = \
                        smu_chisq__DPfitted_list(oms_all2, ws_all2, DistancePrior, zref, zcomps, noRSDCorpoly_all)
                    RSDCor_tot_chisqs_polyfitted = \
                        smu_chisq__DPfitted_list(oms_all2, ws_all2, DistancePrior, zref, zcomps, RSDCorpoly_all)
                noRSDCor_2dchisqlist_polyfitted = get_2darray_from_1d(noRSDCor_tot_chisqs_polyfitted, numw2, numom2)
                RSDCor_2dchisqlist_polyfitted = get_2darray_from_1d(RSDCor_tot_chisqs_polyfitted, numw2, numom2)
                ################################################
                ## Step 5.2: plot contour
                fig = plt.figure(figsize=(18,6))
                ax1 = fig.add_subplot(121); ax2 = fig.add_subplot(122)
                contourlabel = ''
                for isig in range(len(sigs)-1):
                    nsig = nsigofsigCL(sigs[isig])
                    contourlabel += ('%g'%(round(nsig, 1))+', ')
                contourlabel += ('%g'%(round(nsigofsigCL(sigs[len(sigs)-1]), 1))+'')
                contourlabel += ' $\\sigma$ contours'
                plot_contour(ax1, oms, ws, noRSDCor_2dchisqlist, do_smooth=False, sigs=sigs, show_marg_rlt=False,
                    ommin=min(oms),  ommax=max(oms), wmin=min(ws), wmax=max(ws), plotformat=1,  
                                             label=contourlabel+'\nBefore RSD Correction', legftsize=20)
		if plot_fitted_rlt:
                 plot_contour(ax1, oms2, ws2, noRSDCor_2dchisqlist_polyfitted, do_smooth=False, sigs=sigs, 
                    ommin=min(oms2),  ommax=max(oms2), wmin=min(ws2), wmax=max(ws2), plotformat=2,  
                                             label='', legftsize=20)
                plot_contour(ax2, oms, ws, RSDCor_2dchisqlist, do_smooth=False, sigs=sigs, show_marg_rlt=False,
                    ommin=min(oms),  ommax=max(oms), wmin=min(ws), wmax=max(ws), plotformat=1,  
                                            label=contourlabel+'\nAfter RSD Correction', legftsize=20)
		if plot_fitted_rlt:
                 plot_contour(ax2, oms2, ws2, RSDCor_2dchisqlist_polyfitted, do_smooth=False, sigs=sigs, 
                    ommin=min(oms2),  ommax=max(oms2), wmin=min(ws2), wmax=max(ws2), plotformat=2,  
                                             label='', legftsize=20)
                
                figtilt = separate_path_file(chisqfilename)[1]
                figtilt += '\n Fitted by '+DPname_general+'; order = '+str(polyfitdeg)
                fig.suptitle(figtilt, fontsize=figtiltfs)
                ax1.grid(); ax2.grid()
                fig.subplots_adjust(top=0.85)
                print ; plt.show()
                if save_tot_contour:
                    figname = chisqfilename+'.'+DPname_txt+'.fitted_contour.'+figfmt
                    print 'Figure saved: ', figname
                    fig.savefig(figname, fmt=figfmt)

            if RW_chisq < RW_chisq_min:
		RW_chisq_min = RW_chisq
		SixBinRedshift_bf = SixBinRedshift
		RW_BC  *= 1.06
	    else:
		RW_BC *= 0.99
	    print '\tRW_step, SixBinRedshift_bf = ', RW_step, SixBinRedshift_bf
	    print '\tRW_BC, RW_chisq_min = ', RW_BC, RW_chisq_min
	    if RW_break:
		break
     return zref, zcomps, noRSDCor_f0sigchisq0_all, noRSDCorpoly_all, RSDCor_f0sigchisq0_all, RSDCorpoly_all
