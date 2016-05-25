def smu__plot(filename, smax=50, numsbin=50, nummubin=120, 
		isMAX=None, 
		is_sig_pi=False, 
		deltais = 2,
		no_s_sq_in_y = False,
		savefig=True,
		):

	smusettings = {
            'smin':0.0, 'smax':smax, 'numsbin':numsbin,
            'mumin':0.0, 'mumax':1.0, 'nummubin':nummubin,
#            'mumin':0.0, 'mumax':1.0, 'nummubin':120,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }

	smu__initsmusettings(smusettings)

	if True:

		DDlist, DRlist, RRlist = Xsfrom2ddata(smu__loadin(filename, smusettings), [4,5,6])
		if isMAX == None:
			isMAX = numsbin
		fig, ax1 = figax()

                ### packed count of xi as a function of s
                ismin = 0;
                imumin = 0; imumax = nummubin;
                sasx = []; packedxiasy = [];

                while ismin < isMAX:
                    ismax = min(ismin+deltais,numsbin-1);
                    nows=((slist[ismin]**3.0+slist[ismax]**3.0)/2.0)**(1.0/3.0);   sasx.append(nows)
                    if not no_s_sq_in_y:
			    if is_sig_pi:
			    	Y = []
			    	for nowsig in range(ismax-1):
				 nowpi = int(np.sqrt(nows**2.0 - nowsig**2.0)) 
				 nowpi = max(nowpi,0)
				 Y.append(packedxi(DDlist,DRlist,RRlist,nowsig,nowsig+1,nowpi,nowpi+1))
			    	packedxiasy.append(nows*nows*sum(Y) / (len(Y)+0.0))
			    else:
                                packedxiasy.append(nows*nows*packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
                    else:
                            packedxiasy.append(packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
                    ismin += deltais;
                if True:
                        ax1.plot(sasx, packedxiasy, marker='o', markersize=1)
                        ax1.set_xlabel('$s\ [\\rm Mpc/h]$', fontsize=25)
                        if not no_s_sq_in_y:
                                ax1.set_ylabel('$s^2\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                        else:
                                ax1.set_ylabel('$\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                        ax1.set_xlim(min(sasx),max(sasx))
			ax1.grid()
                        #ax1.set_xlim(0,50)
		if savefig:
			fig.savefig(filename+'.png', format = 'png')

def smu__plot_Tpcfdatas(Tpcfdatas,  smusettings, key_leg_list=None, muedges = np.linspace(0,1,11), 
	ismin_of_intxi = 6, ismax_of_intxi = 50, deltais = 2, 
	xifunction = intxi_FractionalS, lsfun=None, lcfun=None, 
	figxsize=16, figysize=8, markersize=1, lw=3, showleg=True, plot_theta=True, legfs=13, givenfigaxs=None,
	leftpanel_xrange = [], figtilt = None,
	only_plot_shape=False, # do not display the left panel; only show xi as a function of mu
	autoleg=True, # automatic create label
	legiskey=False,  # use key as label in legend
	ignore_catname_in_leg=True, # when autoleg, do not display catname in label
	normedintxi=False, # normalize amplitude of intxi
	no_s_sq_in_y=False, #
	return_XY = False,
	 ):
	""" Plot xi(s,mu) as a funtion of scale and angle"""
	itcolors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray'])
	if givenfigaxs == None:
		fig=plt.figure(figsize=(figxsize,figysize));  
		if not only_plot_shape:
			ax1=fig.add_subplot(121); ax2=fig.add_subplot(122)
		else:
			ax2 = fig.add_subplot(111); ax1=ax2
	else:
		fig, ax1, ax2 = givenfigaxs
	ylabel2 = '$\\int_{'+str(ismin_of_intxi)+'}^{'+str(ismax_of_intxi)+'} \\xi(s,\\mu) ds$ '
	if xifunction != intxi_FractionalS:
		ylabel2 = smu__xifunname(xifunction)+'$(s_{\\rm min}='+str(ismin_of_intxi)+',\\ s_{\\rm max}='+str(ismax_of_intxi)+')$'

	if key_leg_list == None:
		key_leg_list = list(Tpcfdatas.keys())
		autoleg=True

	Xs1, Ys1, Xs2, Ys2 = [], [], [], []
	for key_leg in key_leg_list:                
            	if autoleg:
			nowleg = keyname = key_leg
			if ignore_catname_in_leg:
				if nowleg[0:4] != 'data':
					if nowleg.find('--') >= 1:
						nowleg = nowleg[nowleg.find('--')+2:len(nowleg)]
		elif legiskey:
			nowleg = keyname = key_leg
		else:
	                keyname, nowleg = key_leg
		
                DDlist, DRlist, RRlist = Xsfrom2ddata(Tpcfdatas[keyname], [4,5,6])

		if lsfun == None:
	                ls = '-' 
		else:
			ls = lsfun(keyname)

		if lcfun == None:
			lc = itcolors.next()
		else:
			lc = lcfun(keyname)


                ### packed count of xi as a function of s
                ismin = 0;
                imumin = 0; imumax = nummubin;
                sasx = []; packedxiasy = [];

                while ismin < numsbin:
                    ismax = min(ismin+deltais,numsbin-1);
                    nows=(slist[ismin]+slist[ismax])/2.0;   sasx.append(nows)
		    if not no_s_sq_in_y:
	                    packedxiasy.append(nows*nows*packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
		    else:
	                    packedxiasy.append(packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
                    ismin += deltais;
		if not only_plot_shape:
			Xs1.append([x for x in sasx])
			Ys1.append([x for x in packedxiasy])
	                ax1.plot(sasx, packedxiasy, marker='o', markersize=1, ls=ls, lw=lw, label = nowleg, c=lc)
	                ax1.set_xlabel('$s\ [\\rm Mpc/h]$', fontsize=25)
			if not no_s_sq_in_y:
		                ax1.set_ylabel('$s^2\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
			else:
		                ax1.set_ylabel('$\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
			if leftpanel_xrange != []:
				ax1.set_xlim(leftpanel_xrange)
			else:
		                ax1.set_xlim(min(sasx),max(sasx))

                ### \int xi ds as a function of mu
                i_muedges = smu__get_imuedges(muedges,smusettings);
                muedgemids = smu__get_muedgemids(i_muedges, smusettings);
                intxis = [xifunction(DDlist, DRlist, RRlist, 
                    ismin_of_intxi, ismax_of_intxi, i_muedges[imubin][0], i_muedges[imubin][1]) 
                          for imubin in range(len(i_muedges))];
		
		if normedintxi:
			intxis = normto1(intxis);

		Xs2.append([x for x in muedgemids])
		Ys2.append([x for x in intxis])
                ax2.plot(muedgemids, intxis, marker='o', markersize=markersize, ls=ls, lw=lw, label = nowleg, c=lc)
                ax2.set_xlabel('$1-\\mu$', fontsize=25)
                ax2.set_ylabel(ylabel2, fontsize=25)
		
		if plot_theta:
			ax3 = ax2.twiny()
			ax3.set_xlim(0,1.0)
			thetas = [0,20, 30, 40, 50, 60, 70, 80, 90]
			xticks2 = [1-np.cos(theta/(180.0/np.pi)) for theta in thetas]
			#xtickslabels2 = [frommutotheta(xtick) for xtick in xticks]
			ax3.set_xticks(xticks2)
			xtick2labels = [str(theta)+'$^{\\circ}$' for theta in thetas]
			ax3.set_xticklabels(xtick2labels, color='r')
			ax3.set_xlabel('$\\theta$', fontsize=27, color='r')

	for ax in [ax1, ax2]:
	    ax.grid(ls='--'); 
	    if showleg: ax.legend(loc='best', frameon=False, fontsize=legfs)

	fig.tight_layout()
	if figtilt != None:
		fig.set_title(figtilt, fontsize=24)
	if not return_XY: 
		return fig, ax1, ax2
	else:
		return fig, ax1, ax2, [Xs1, Ys1, Xs2, Ys2]





#### Very old programme...

def make_2pcf_plottings(nowsmins=[5],nowsmaxs=[100],iomws=[0,1,2,3,4],irbinlist=[0,1,2,3,4], noplot=False, plotRSDeff=False,
                   nowxrange=[],nowyrange=[],plotnoRSD=True,plotRSD=True,printsetttings=False, xifunction=intxi,
		   chisqmethod = 'minus',
                   iimocklist=range(len(imocks)),muedges=np.linspace(0,1,11), minmucut = -1, RSDcor=True, RSDcortozero=False,
		   consistsrange=False,
                   maxmucut=1000000,useourerrorbar=False,savfig=False,shifttozero=False,divto1=False, usecovmat=False,
                   reprocessed_shifttozero_covmat=False,reprocessed_divto1_covmat=False,
                   set_last_mubindiff_tozero = False, separate_RSDCor=False,
                   diffchisqmethod='use_lastrbin_as_ref', weight_for_avgdiff=[],
                   ConsiderVarVar=False, CorrectedCovmat=False, usegivencovmats=[], chisqprintinfo=True,
		   notitle = False):

    ### Range of s
    # Minimal/Maximal indice of s
    if xifunction == intxi or xifunction == packed_count:
	    nowismins = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmins]
	    nowismaxs = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmaxs]
    else:
	    nowismins = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmins]
	    nowismaxs = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmaxs]
#    else:
#    	print 'ERROR (make_2pcf_plottings)! Unknown xifunction!'
#    	return
    	
    ### Binning of mu
    i_muedges = [[min(max(0,  int((muedges[row]-mumin)/deltamu)  ),nummubin), 
                 min(max(0,  int((muedges[row+1]-mumin)/deltamu)-1  ),nummubin)]
                 for row in range(len(muedges)-1)]
    # Minimal/Maximal cut of mu
    if minmucut > 0:
        mini_mu = get_i_mu(minmucut)
        i_muedges2 = copy.deepcopy(i_muedges); i_muedges=[]
        for i_muedge in i_muedges2:
            if i_muedge[1] < mini_mu:
                continue
            else:
                i_muedges.append([max(i_muedge[0],mini_mu),i_muedge[1]])
        if printsetttings:
            print 'Re-define i_muedge due to minimal mucut ', minmucut, '.\n\torig:', i_muedges2, '\n\tnow:', i_muedges
    if maxmucut < 1:
        maxi_mu = get_i_mu(maxmucut)
        i_muedges2 = copy.deepcopy(i_muedges); i_muedges=[]
        for i_muedge in i_muedges2:
            if i_muedge[0] > maxi_mu:
                continue
            else:
                i_muedges.append([i_muedge[0],min(maxi_mu,i_muedge[1])])
        if printsetttings:
            print 'Re-define i_muedge due to maximal mucut ', maxmucut, '.\n\torig:', i_muedges2, '\n\tnow:', i_muedges
    # indice of mu for each binning; middle value at each bin
    i_muranges = [range(i_muedge[0],i_muedge[1]+1) for i_muedge in i_muedges]
    muedgemids = [(get_mid_smu(0,i_muedge[0])[1]+get_mid_smu(0,i_muedge[1])[1])*0.5 for i_muedge in i_muedges]
    
    ### print settings of the program
    if printsetttings:
        print 'Settings:'
        print '\tlist of irbins: \n\t\t',irbinlist,'\n\t\t', [rbinstrs[irbin] for irbin in irbinlist]
        print '\tlist of mocks: \n\t\t',iimocklist,'\n\t\t', [imocks[iimock] for iimock in iimocklist]
        print '\tlist of muedges:  \n\t\t',muedges, '\n\t\t', i_muedges
        print '\t\tBinning of i_mu:  \n\t\t\t', i_muranges
        print '\t\t# of i_mu in each bin:  \n\t\t\t', [len(i_muranges[row]) for row in range(len(i_muranges))]
        print '\t\tMiddles of mu: \n\t\t\t', muedgemids
    ### loop: range of s & cosmology (one by one)
    i_of_smin = i_of_smax = -1;
    for nowismin in nowismins:
     i_of_smin +=1;
     for nowismax in nowismaxs:
        i_of_smax += 1;
        nowisrange = range(int(nowismin),int(nowismax))
        chisqrlts = [];
        
        ### use loaded data (fast)
        use_loaded=True; plot_eachmock = False
            
        ### Settings for plot    
        lw=4;
        lslist = ['-', '--']
        lclist = ['b','g','k','r','c','m']
        markerlist = ['o', '*', 'D', 'p', '+', 's']
        ils=-2; iplot = -1; 
        
        ### loop of cosmology
	i_of_omw = -1
        for iomw in iomws:
            tiltstr = ''
	    i_of_omw += 1
            if printsetttings:
                print '\n\n\tCalculating cosmology:   \n\t\t',iomw, '\n\t\t', omws[iomw]
                print '\trange of i_s: \n\t\t', nowismin, '('+str(get_mid_smu(nowismin,0)[0])+') to ', nowismax-1,\
                    '('+str(get_mid_smu(nowismax-1,0)[0])+')\n\t\t', nowisrange
        
            om,w = omws[iomw]; 
            ils+=2;ilc=0; 
           
            ### binned xi
            binnedxis1 = [[0 for irbin in irbinlist] for iimock in iimocklist]
            binnedxis2 = [[0 for irbin in irbinlist] for iimock in iimocklist]
            intxiavgs1 = [[0 for irbin in irbinlist] for iimock in iimocklist]
            intxiavgs2 = [[0 for irbin in irbinlist] for iimock in iimocklist]
            
            ### for each cosmology, loop of different rbin (redshift)
            i_of_bin = -1; 
            xi1_diffrs, xi1std_diffrs, xi2_diffrs, xi2std_diffrs = [],[],[],[]
            
            for irbin in irbinlist:
                i_of_bin +=1;
                zetalist1,zetalist2,zeta_mu_array1s,zeta_mu_array2s,mucompactarray1s,mucompactarray2s = [],[],[],[],[],[]
                ### loop of indice of mock
                i_of_mock = -1;
                if consistsrange:
                    nowsmin = nowsmins[i_of_smin]*rescalfac[iomw][irbin]
                    nowsmax = nowsmaxs[i_of_smin]*rescalfac[iomw][irbin]
                    if xifunction == intxi or xifunction == packed_count:
                    	nowismin = min(max(0,int((nowsmin-smin)/deltas)),numsbin)
	                nowismax = min(max(0,int((nowsmax-smin)/deltas)),numsbin)
	            else:
	                nowismin = min(max(0,(nowsmin-smin)/deltas),numsbin)
	                nowismax = min(max(0,(nowsmax-smin)/deltas),numsbin)
                    #print 'om/w = ', omws[iomw], ', rbin = ', rbinstrs[irbin], ': srange = ', nowsmin, nowsmax, \
                    #    '; isrange = ', nowismin, nowismax
                for iimock in iimocklist:
                    i_of_mock += 1;
                    if not use_loaded:
                        file1 = get_2pcffile(imocks[iimock], 'noRSD', rbinstrs[irbin], om, w)
                        file2 = get_2pcffile(imocks[iimock], 'RSD', rbinstrs[irbin], om, w)
                        zetalist1 = loadin_smu(file1); zetalist2 = loadin_smu(file2)
                        DD1, DR1, RR1 = loadin_counts(file1)
                        DD2, DR2, RR2 = loadin_counts(file2)
                    else:
                        #zetalist1 = zetalist1s[iomw][irbin][iimock]
                        #zetalist2 = zetalist2s[iomw][irbin][iimock]
                        DD1, DR1, RR1 = DDlist1s[iomw][irbin][iimock], DRlist1s[iomw][irbin][iimock], RRlist1s[iomw][irbin][iimock]
                        DD2, DR2, RR2 = DDlist2s[iomw][irbin][iimock], DRlist2s[iomw][irbin][iimock], RRlist2s[iomw][irbin][iimock]
                        
                    ### calculate integral xi,  save result to binnedxis
                    binnedxis1[i_of_mock][i_of_bin] = \
                        [ xifunction(DD1, DR1, RR1, nowismin, nowismax, muedge[0], muedge[1]) for muedge in i_muedges]
                    binnedxis2[i_of_mock][i_of_bin] = \
                        [ xifunction(DD2, DR2, RR2, nowismin, nowismax, muedge[0], muedge[1]) for muedge in i_muedges]

            ### Normalization & Plotting & Chisq calculation 
            # Figure
            if not noplot:
                fig=plt.figure(figsize=(18,6));ax=fig.add_subplot(111); 
                ax.set_xlabel('$\\mu$',fontsize=18); ax.set_ylabel('$\\sum \\xi$',fontsize=18)
            # binned xi and its STD at different r
            mucompactarray1_diffrs, mucompactarray1std_diffrs, mucompactarray2_diffrs, mucompactarray2std_diffrs = [],[],[],[]
            for i_of_bin in range(len(irbinlist)):
                mucompactarray1s,mucompactarray2s = [],[]
                for i_of_mock in range(len(iimocklist)):
                    if divto1:
                        intxiavgs1[i_of_mock][i_of_bin], binnedxis1[i_of_mock][i_of_bin] = \
                        	normto1(binnedxis1[i_of_mock][i_of_bin],returnavg=True)
                        intxiavgs2[i_of_mock][i_of_bin], binnedxis2[i_of_mock][i_of_bin] = \
                        	normto1(binnedxis2[i_of_mock][i_of_bin],returnavg=True)
                        
                    if shifttozero:
                        intxiavgs1[i_of_mock][i_of_bin], binnedxis1[i_of_mock][i_of_bin] = \
                        	shiftozero(binnedxis1[i_of_mock][i_of_bin],returnavg=True)
                        intxiavgs2[i_of_mock][i_of_bin], binnedxis2[i_of_mock][i_of_bin] = \
                        	shiftozero(binnedxis2[i_of_mock][i_of_bin],returnavg=True)
                        
                    intxi1, intxi2 = binnedxis1[i_of_mock][i_of_bin], binnedxis2[i_of_mock][i_of_bin]
                    mucompactarray1s.append(copy.deepcopy(intxi1)); mucompactarray2s.append(copy.deepcopy(intxi2))
                mucompactarray1, mucompactarray1std = get_avgstd_array(mucompactarray1s,ConsiderVarVar=ConsiderVarVar)
                mucompactarray2, mucompactarray2std = get_avgstd_array(mucompactarray2s,ConsiderVarVar=ConsiderVarVar)
                ### Correction for Covmat using arXiv:1312.4841
                if CorrectedCovmat:
                    nownummubin = len(muedges)-1
                    if nownummubin + 1 >= len(iimocklist) -1:
                        print 'Too many mu-bins, not enough # of mocks! Will not correct covmat!'
                    else:
                        CapD = (nownummubin+1)/(len(iimocklist)-1.0)
                        CorcFactor = 1.0 / (1.0 - CapD) ## Eq6 of arXiv:1312.4841
                        CorcFactor = np.sqrt(CorcFactor)
                        for row in range(len(mucompactarray1std)):
                            mucompactarray1std[row] *= CorcFactor
                            mucompactarray2std[row] *= CorcFactor
                mucompactarray1_diffrs.append(copy.deepcopy(mucompactarray1)); 
                mucompactarray1std_diffrs.append(copy.deepcopy(mucompactarray1std));
                mucompactarray2_diffrs.append(copy.deepcopy(mucompactarray2)); 
                mucompactarray2std_diffrs.append(copy.deepcopy(mucompactarray2std));                    
                
                if noplot:
                    continue
                    
                ### label of the noRSD/RSD curve
                labelstr1 = '${\\rm om/w=%.2f'%om+'/%.2f'%w+';\ \ '+rbinstrs[irbin]+\
                    '}$ No RSD: $\\bar\\sum_\\xi=%.2f'%(avgarray(mucompactarray1))+'$'
                labelstr2 = '${\\rm om/w=%.2f'%om+'/%.2f'%w+';\ \ '+rbinstrs[irbin]+\
                    '}$ With RSD: $\\bar\\sum_\\xi=%.2f'%(avgarray(mucompactarray2))+'$'            

                ### plot the curve
                if plotnoRSD:
                    if useourerrorbar:
                        ourerrorbar(ax, muedgemids, mucompactarray1, mucompactarray1std, lw=lw, c=lclist[mod(ilc,len(lclist))],
                                ls=lslist[mod(ils,len(lslist))], capsize=10, label=labelstr1, 
                                marker=markerlist[mod(iplot,len(markerlist))], ms=10, not_dof_div=False, polyorder=1)
                    else:
                        ax.errorbar(muedgemids,mucompactarray1,mucompactarray1std,label=labelstr1,
                            lw=lw,ls=lslist[mod(ils,len(lslist))],capsize=10,c=lclist[mod(ilc,len(lclist))],
                            marker=markerlist[mod(iplot,len(markerlist))],markersize=10)
                if plotRSD:
                    if useourerrorbar:
                        ourerrorbar(ax, muedgemids, mucompactarray2, mucompactarray2std, lw=lw, c=lclist[mod(ilc,len(lclist))], 
                                ls=lslist[mod(ils+1,len(lslist))], capsize=10, label=labelstr1, 
                                marker=markerlist[mod(iplot,len(markerlist))], ms=10, not_dof_div=False, polyorder=1)
                    else:
                        ax.errorbar(muedgemids,mucompactarray2,mucompactarray2std,label=labelstr2,lw=lw,
                                    ls=lslist[mod(ils+1,len(lslist))],capsize=10,c=lclist[mod(ilc,len(lclist))],
                                    marker=markerlist[mod(iplot,len(markerlist))],markersize=10)
                iplot+=1;ilc+=1
            
            ### After the loop of rbins:
            # Calculate the chisq; print it in the title
            if RSDcor == True:
		if not separate_RSDCor:
		    if i_of_omw ==0:
			if RSDcortozero:
	                    	RSDdiff = [ mucompactarray2_diffrs[ir] for ir in range(len(irbinlist))]
			else:
	                    	RSDdiff = [get_diffarray(mucompactarray1_diffrs[ir], mucompactarray2_diffrs[ir]) 
        	                	for ir in range(len(irbinlist))]
		else:
			if RSDcortozero:
	                    	RSDdiff = [ mucompactarray2_diffrs[ir] for ir in range(len(irbinlist))]
			else:
	                    RSDdiff = [get_diffarray(mucompactarray1_diffrs[ir], mucompactarray2_diffrs[ir]) 
        	            	for ir in range(len(irbinlist))]
                if plotRSDeff:
                        fig2 = plt.figure(figsize=(18,6));
                        ax2=fig2.add_subplot(111);
                        for ir in range(len(irbinlist)):
                            ax2.plot(RSDdiff[ir], c=lclist[mod(ir,len(lclist))], label=rbinstrs[ir], lw=4)
                        ax2.plot(mucompactarray2_diffrs[0], c=lclist[mod(0,len(lclist))], label=rbinstrs[0]+' Curve', 
                                 lw=4,ls='--')
                        ax2.legend(prop={'size':13},loc='best',frameon=False);
                        plt.show();
                mucompactarray2_diffrs = [get_diffarray(RSDdiff[ir], mucompactarray2_diffrs[ir]) for ir in range(len(irbinlist))]
            ### Chisqs using error bars
            chisq_norsd = chisq_of_intxiatdiffr(mucompactarray1_diffrs, mucompactarray1std_diffrs, chisqmethod = chisqmethod)
            chisq_rsd = chisq_of_intxiatdiffr(mucompactarray2_diffrs, mucompactarray2std_diffrs, chisqmethod = chisqmethod)
            if not noplot:
	            tiltstr += '${\\rm om/w=%.2f'%om+'/%.2f'%w
        	    tiltstr += ';\ \ NoRSD/RSD\ \\chi^2= %.3f'%chisq_norsd+'/%.3f'%chisq_rsd+'(erbar);}$'
            ### Chisqs using covmats
            if not usecovmat:
                chisqrlts.append([[om,w], chisq_norsd, chisq_rsd])
            else:
                if chisqprintinfo:
                    print '## covmat: norsd '
                if usegivencovmats == []:
                    chisq_norsd, like_norsd, covmats_norsd = chisq_of_intxiatdiffr_cov(binnedxis1, 
                    	printinfo=chisqprintinfo,diffchisqmethod=diffchisqmethod)
                if chisqprintinfo:
                    print '## covmat: rsd'
                if usegivencovmats == []:
                    chisq_rsd, like_rsd, covmats_rsd = chisq_of_intxiatdiffr_cov(binnedxis2, printinfo=chisqprintinfo,
                    diffchisqmethod=diffchisqmethod)
                chisq_norsd2, chisq_rsd2 = 0, 0
                numr = len(irbinlist)
                ### if covmat given, using them rather than our owns; 
                if usegivencovmats != []:
                    covmats_norsd, covmats_rsd = usegivencovmats
                if diffchisqmethod=='use_lastrbin_as_ref':
		        for ir in range(numr-1):
		            diff_norsd = get_diffarray(mucompactarray1_diffrs[ir],mucompactarray1_diffrs[numr-1])
		            diff_rsd = get_diffarray(mucompactarray2_diffrs[ir],mucompactarray2_diffrs[numr-1])
		            if set_last_mubindiff_tozero:
		            	diff_norsd[len(diff_norsd)-1]=0; diff_rsd[len(diff_rsd)-1]=0; 
		            	#print "diff_norsd: ", diff_norsd
		            	#print "diff_rsd: ", diff_rsd
			    chisq_norsd2 += chisq_like_cov_xbar(diff_norsd, covmats_norsd[ir]) [0]
			    chisq_rsd2 += chisq_like_cov_xbar(diff_rsd, covmats_rsd[ir])[0]
		elif diffchisqmethod=='use_weightedavg_as_ref':
			avg1 = avgarray_2d(mucompactarray1_diffrs, wei=weight_for_avgdiff)
			avg2 = avgarray_2d(mucompactarray2_diffrs, wei=weight_for_avgdiff)
			for ir in range(numr):
		            diff_norsd = get_diffarray(mucompactarray1_diffrs[ir],avg1)
		            diff_rsd = get_diffarray(mucompactarray2_diffrs[ir],avg2)
		            if set_last_mubindiff_tozero:
		            	diff_norsd[len(diff_norsd)-1]=0; diff_rsd[len(diff_rsd)-1]=0; 
		            	#print "diff_norsd: ", diff_norsd
		            	#print "diff_rsd: ", diff_rsd
			    chisq_norsd2 += chisq_like_cov_xbar(diff_norsd, covmats_norsd[ir]) [0]
			    chisq_rsd2 += chisq_like_cov_xbar(diff_rsd, covmats_rsd[ir])[0]
		elif diffchisqmethod <= numr-1:
			for ir in range(numr):
		            diff_norsd = get_diffarray(mucompactarray1_diffrs[ir],mucompactarray1_diffrs[diffchisqmethod])
		            diff_rsd = get_diffarray(mucompactarray2_diffrs[ir],mucompactarray2_diffrs[diffchisqmethod])
		            if set_last_mubindiff_tozero:
		            	diff_norsd[len(diff_norsd)-1]=0; diff_rsd[len(diff_rsd)-1]=0; 
		            	#print "diff_norsd: ", diff_norsd
		            	#print "diff_rsd: ", diff_rsd
			    chisq_norsd2 += chisq_like_cov_xbar(diff_norsd, covmats_norsd[ir]) [0]
			    chisq_rsd2 += chisq_like_cov_xbar(diff_rsd, covmats_rsd[ir])[0]
#			mucompactarray1_diffrs
#			for ir in range(numr):
			
                #chisq_like_cov_xbar(diffX, covmat)
                #tiltstr += (' ${\\rm'+str(chisq_norsd)+'/'+str(chisq_rsd)+'}$')
                chisqrlts.append([[om,w], chisq_norsd2, chisq_rsd2])
                if not noplot:
	                tiltstr += (' ${\\rm  %.3f'%chisq_norsd2+'/%.3f'%chisq_rsd2+'(covmat)}$')
            ### Other minor settings of the figure
            if not noplot:
		    if nowxrange != []:
		        ax.set_xlim(nowxrange[0],nowxrange[1]) 
		    if nowyrange != []:
		        ax.set_ylim(nowyrange[0],nowyrange[1]) 
		    fig.subplots_adjust(top=0.85);
		    if not notitle:
		      ax.set_title(tiltstr, fontsize=14);
		    ax.legend(loc='best',frameon=False);
		    plt.show();
		    if savfig:
		        fig.savefig('om%.2f'%om+'_w%.2f'%w+'.png',format='png');
    if not usecovmat:
        return chisqrlts, mucompactarray1_diffrs,mucompactarray2_diffrs
    else:
        return chisqrlts, covmats_norsd, covmats_rsd, mucompactarray1_diffrs,mucompactarray2_diffrs
