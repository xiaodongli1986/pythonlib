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
                        mufmt='%.2f', sfmt='%.2f'):
    return smufile+'__intxis/'+smu__xifunname(smu__intxi__settings['xifunction'])+\
        '.'+str(smu__intxi__settings['nummubin'])+\
        'mu'+mufmt%smu__intxi__settings['mumin']+'to'+mufmt%smu__intxi__settings['mumax']+'.s'+\
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
    if chisqfile.find('NGC') >= 0:
        return 'N'
    elif chisqfile.find('SGC') >= 0:
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
                         mufmt='%.2f', sfmt='%.2f', smudata=None):
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
				skip_exited_smufile=True):
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
	totbin=3, normedintxi=True):
    return datamockdir+'/covmats'+'/'+smu__intxi_RSDCovbasename(Sky, catname2, RSDstr, smu__intxi__settings, refibin, totbin, normedintxi)+'.cov'+str(icovmat)

def smu__intxi_chisqfilename(Sky, catname2, RSDstr, smu__intxi__settings, refibin = 0, 
	totbin=3, numchisq=10, prefix='', normedintxi=True):
    return datamockdir+'/covmats'+'/'+prefix+smu__intxi_RSDCovbasename(Sky, catname2, RSDstr, smu__intxi__settings, refibin, totbin, normedintxi)+'.'+str(numchisq)+'chisqs'


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
                            ):
    
    if Sky == 'N':
        catnames = ['DR12v4-LOWZ-N', 'DR12v4-CMASS-N']
    elif Sky == 'S':
        catnames = ['DR12v4-LOWZ-S', 'DR12v4-CMASS-S']
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
                totbin=totbin, normedintxi=normedintxi)
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
                            skip_calc_if_chisqfile_found = True
                            ):

    if use_omw_testfun:
        scanname = 'omw-Test-Scan'
        oms = np.linspace(0, 1.0, 20)
        ws = np.linspace(-2.5, -0.5, 20)
        numw=len(ws); numom = len(oms)
        omws = sumlist( [[[om,w] for om in oms] for w in ws])
    else:
	global oms, ws, omws, scanname
        scanname = scanname
        oms = oms
        ws = ws
        numw=len(ws); numom = len(oms)
        omws = omws

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
                                                         normedintxi=normedintxi)

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
                                        use_omw_testfun=use_omw_testfun)

                else: # load in rather than do computation
                    covmats = []
                    for icovmat in range(totbin*2-1):
                        ###!!!!!  RSDCor_covmat : 2
                        covmatfilename = smu__intxi_covmatfilename(icovmat, Sky, 'HR3', 
                                    RSDstr = covmat_RSDstr, 
                                    smu__intxi__settings = smu__intxi__settings, 
                                    refibin = refibin,
                                    totbin = totbin,
                                    normedintxi = normedintxi)


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
                                        use_omw_testfun=use_omw_testfun)

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
                                        use_omw_testfun=use_omw_testfun)  


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
                        noRSDCorchisq, noRSDCorlike = chisq_like_cov_xbar(noRSDCorX, covmat)
                        RSDCorchisq, RSDCorlike = chisq_like_cov_xbar(RSDCorX, covmat)
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



