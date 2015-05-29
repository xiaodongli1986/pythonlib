def smu__plot_Tpcfdatas(Tpcfdatas,  smusettings, key_leg_list=None,muedges = np.linspace(0,1,11), ismin_of_intxi = 6, ismax_of_intxi = 50, lsfun=None, lcfun=None, xifunction = intxi_FractionalS, figxsize=16, figysize=8, deltais = 2, autoleg=True, ignore_catname_in_leg=True,markersize=1,normedintxi=False,no_s_sq_in_y=False):
	""" Plot xi(s,mu) as a funtion of scale and angle"""
	itcolors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray'])
	fig=plt.figure(figsize=(figxsize,figysize));  
	ax1=fig.add_subplot(121); ax2=fig.add_subplot(122)
	ylabel2 = '$\\int_{'+str(ismin_of_intxi)+'}^{'+str(ismax_of_intxi)+'} \\xi(s,\\mu) ds$ '
	if xifunction != intxi_FractionalS:
		ylabel2 = smu__xifunname(xifunction)+'$(s_{\\rm min}='+str(ismin_of_intxi)+',\\ s_{\\rm max}='+str(ismax_of_intxi)+')$'

	if key_leg_list == None:
		key_leg_list = list(Tpcfdatas.keys())
		autoleg=True

	for key_leg in key_leg_list:                
            	if autoleg:
			nowleg = keyname = key_leg
			if ignore_catname_in_leg:
				if nowleg[0:4] != 'data':
					if nowleg.find('--') >= 1:
						nowleg = nowleg[nowleg.find('--')+2:len(nowleg)]
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

                ax1.plot(sasx, packedxiasy, marker='o', markersize=1, ls=ls, lw=3, label = nowleg, c=lc)
                ax1.set_xlabel('$s\ [\\rm Mpc/h]$', fontsize=25)
		if not no_s_sq_in_y:
	                ax1.set_ylabel('$s^2\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
		else:
	                ax1.set_ylabel('$\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                ax1.set_xlim(0,150)

                ### \int xi ds as a function of mu
                i_muedges = smu__get_imuedges(muedges,smusettings);
                muedgemids = smu__get_muedgemids(i_muedges, smusettings);
                intxis = [xifunction(DDlist, DRlist, RRlist, 
                    ismin_of_intxi, ismax_of_intxi, i_muedges[imubin][0], i_muedges[imubin][1]) 
                          for imubin in range(len(i_muedges))];
		
		if normedintxi:
			intxis = normto1(intxis);

                ax2.plot(muedgemids, intxis, marker='o', markersize=markersize, ls=ls, lw=3, label = nowleg, c=lc)
                ax2.set_xlabel('$1-\\mu$', fontsize=25)
                ax2.set_ylabel(ylabel2, fontsize=25)

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
	    ax.grid(); ax.legend(loc='best', frameon=False, fontsize=13)

	fig.tight_layout()
	return fig, ax1, ax2
