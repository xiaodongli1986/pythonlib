

def ChisqContour_avgfiles(omlist, wlist, chisqfiles, outputchisqfile = None, given_cols = False,
	mult_facts = None, output_corrected_files = True, suffix=''):
	chisqs_dict_NoCor = {}
	chisqs_dict_Cor = {}
	for om in omlist:
		for w in wlist:
			chisqs_dict_NoCor[omwstr(om,w)] = []
			chisqs_dict_Cor[omwstr(om,w)] = []
	ichisqfile = -1
	for chisqfile in chisqfiles:
		ichisqfile += 1
		nowf = open(chisqfile, 'r')
		while True:
			nowstr = nowf.readline();
			if nowstr == '': break
			nowstrs = nowstr.split();
			if nowstrs[0][0] == '#': continue
			if given_cols == False:
				nowomwstr, nowchisq1, nowchisq2 = nowstrs[1], float(nowstrs[2]), float(nowstrs[3])
			else:
				nowomwstr, nowchisq1, nowchisq2 = nowstrs[1], sum([float(nowstrs[i]) for i in [4,5,6,7,8]]), sum([float(nowstrs[i]) for i in given_cols])
			if mult_facts != None:
				nowchisq1 *= mult_facts[ichisqfile]
				nowchisq2 *= mult_facts[ichisqfile]
			chisqs_dict_NoCor[nowomwstr].append(nowchisq1)
			chisqs_dict_Cor[nowomwstr].append(nowchisq2)
		nowf.close()
	chisqavg_dict = {}
	for om in omlist:
		for w in wlist:
			nowomwstr = omwstr(om,w)
			chisqavg_dict[nowomwstr] = [avgarray(chisqs_dict_NoCor[nowomwstr]), avgarray(chisqs_dict_Cor[nowomwstr])]
#			chisqavg_dict[nowomwstr] = [min(chisqs_dict_NoCor[nowomwstr]), min(chisqs_dict_Cor[nowomwstr])]
#	print nowomwstr
#	print chisqs_dict_NoCor[nowomwstr]
#	print chisqs_dict_Cor[nowomwstr]
#	print chisqavg_dict[nowomwstr]
	if outputchisqfile == None: 
		outputchisqfile = chisqfiles[0]+'.'+str(len(chisqfiles))+'files'
		if given_cols != False:
			outputchisqfile+='.givencols'
			for col in given_cols:
				outputchisqfile += ('.'+str(col))
		outputchisqfile += suffix
	nowf = open(outputchisqfile, 'w')
		
	for om in omlist:
		for w in wlist:
			nowomwstr = omwstr(om,w)
			nowf.write(' 0.0000 '+nowomwstr+'  '+str(chisqavg_dict[nowomwstr][0])+' '+str(chisqavg_dict[nowomwstr][1])+' \n')
	print ' (ChisqContour_avgfiles) Averaged result of ', len(chisqfiles),' files outputed to file:\n\t\t', outputchisqfile
	nowf.close()
	if output_corrected_files:
	        nowf = open(outputchisqfile+'.Corrected', 'w')

                for w in wlist:
        		for om in omlist:
                        	nowomwstr = omwstr(om,w)
	                        nowf.write(str(om)+' '+str(w)+'  '+str(chisqavg_dict[nowomwstr][1])+' \n')
        	print ' (ChisqContour_avgfiles) Averaged result of ', len(chisqfiles),' files outputed to file:\n\t\t', outputchisqfile+'.Corrected'

	        nowf.close()
	return outputchisqfile

def ChisqContour_Plot(omlist, wlist, chisqfile,
			make_plot = True,
			savefig = True,
			showfig = True,
			ommin=None,ommax=None,wmin=None,wmax=None,
			return_chisqlist=False):
        chisqs_dict_NoCor = {}
        chisqs_dict_Cor = {}
        nowf = open(chisqfile, 'r')
        while True:
                        nowstr = nowf.readline();
                        if nowstr == '': break
                        nowstrs = nowstr.split();
                        if nowstrs[0][0] == '#': continue
                        nowomwstr, nowchisq1, nowchisq2 = nowstrs[1], float(nowstrs[2]), float(nowstrs[3])
                        chisqs_dict_NoCor[nowomwstr]=nowchisq1
                        chisqs_dict_Cor[nowomwstr]=nowchisq2
        nowf.close()
	chisqlist_nosyscor = [[chisqs_dict_NoCor[omwstr(om,w)] for om in omlist] for w in wlist]
	chisqlist_syscor = [[chisqs_dict_Cor[omwstr(om,w)] for om in omlist] for w in wlist]
	#chisqlist_syscor = [[ ((om-0.26)/0.1)**2.0 + ((w+1)/0.1)**2.0  for om in omlist] for w in wlist]
	if make_plot:
                        fig = plt.figure(figsize=(16,6))
                        ax1 = fig.add_subplot(121)
                        ax2 = fig.add_subplot(122)
                        sigA, sigB, sigC = sig1, sig2, sig3
			if ommin == None: ommin = min(omlist)
			if ommax == None: ommax = max(omlist)
			if wmin == None: wmin = min(wlist)
			if wmax == None: wmax = max(wlist)
                        plot_contour(ax1, omlist, wlist, chisqlist_nosyscor, label='NO Sys Cor',
                            ommin = ommin, ommax = ommax, wmin = wmin, wmax = wmax,
                            do_smooth=False, smsigma=0.5, sigA = sigA, sigB = sigB, sigC = sigC)
                        plot_contour(ax2, omlist, wlist, chisqlist_syscor, label='After Sys Cor',
                            ommin = ommin, ommax = ommax, wmin = wmin, wmax = wmax,
                            do_smooth=False, smsigma=0.5, sigA = sigA, sigB = sigB, sigC = sigC)
                        fig.suptitle(separate_path_file(chisqfile)[1])
                        if savefig:
                                figname = chisqfile+'.png';
                                print '   Figure saved: ', figname
                                plt.savefig(figname, fmt='png')
                        if showfig: plt.show()

	if return_chisqlist:
		return chisqlist_nosyscor, chisqlist_syscor
