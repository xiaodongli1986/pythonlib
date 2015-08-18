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
        if time2 - time1 > 10:
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


def smu__inxi_RSDCor_Covmat(Sky='N', catname2='data', smu__intxi__settings = smu__intxi__settings_std, totbin=3, 
	RSDstr = 'RSD', savefig=False, figfmt='png', makeplot=True, docovmat = False, savecovmat=True,refibin = 0, normedintxi=True):
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
                smufile = Tpcfrltfilename(binsplittedfilename(file0, ibin+1))
                smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                #print isfile(smufile)
                #print isfile()
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
		    nowfile = smu__intxi_covmatfilename(icovmat, Sky, catname2, RSDstr, smu__intxi__settings, refibin = refibin)
		    np.savetxt(nowfile, covmats[icovmat])
		    print 'filesaved: ', nowfile
            return mus, avgdiffs, avgintxis, covmats


### Plot intxi files (convenient for checking...)
