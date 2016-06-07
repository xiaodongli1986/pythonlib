### calcchisqs
from numpy import *
import copy 


def smu_xi_str(
	baseoutputfile = 'Dense1subscan', ### 2. Basic name of outputfile
	s1=10, s2=50,sbins=5,sfact=1, 
	cov_catname2='HR3', cov_RSDstr = 'RSD', ### 5. settings for covmat and systematic
	sys_catname='DR12v4-CMASS-N', sys_catname2= 'J08', sys_RSDstr ='RSD',
	suffix = ''
	):
		return baseoutputfile+'.'+str(sbins)+'s'+str(s1)+'to'+str(s2)+'.sfact'+str(1)+'.'\
			+'cov_'+str(cov_catname2)+'_'+cov_RSDstr+'.syscor_'+sys_catname2+'_'+sys_RSDstr+suffix

def smu_xi_calcchisqs(
	omlist, wlist,### 1. List of omegam, w
	baseoutputfile, ### 2. Basic name of outputfile
	imumins = range(20)+[30,40,50,60], ### 3. Range of imumin
	s1=10, s2=50, sbins=5, sfact=1,    smax = 51, smax_mock = 150, ### 4. Things about s
	cov_catname2='HR3', cov_RSDstr = 'RSD', ### 5. settings for covmat and systematic
	sys_catname_suffix='-N', sys_catname2= 'J08', sys_RSDstr ='RSD',
	make_plot=True, showfig=False, savefig=True, 
	totbin = 3, delta_time = 100, nummubins = 120, 	normfun=norm_power, ### Other miscellous
	outputdir = '/home/xiaodongli/software/', suffix = '',
	use_omw_as_sys = False,
	refbins = [0],
	calcbins = [1,2,3,4,5],
		):
	ommin, ommax, wmin, wmax = min(omlist), max(omlist), min(wlist), max(wlist)
	omwlist = sumlist([[[om,w] for om in omlist] for w in wlist])
	sarray = smu_smids(s1, s2)
	sarray = [ x**sfact for x in sarray ]
	nowchisqstr = smu_xi_str(
        	baseoutputfile, s1, s2,sbins,sfact, cov_catname2, cov_RSDstr, sys_catname_suffix, sys_catname2, sys_RSDstr,
		suffix = '.'+str(len(imumins))+'mus.'+str(len(omwlist))+'omws'+suffix )
	if use_omw_as_sys !=False:
		omsys, wsys = use_omw_as_sys
		nowchisqstr += '.force_bestfit_as_'+omwstr(omsys,wsys)
        xisdir = smu_xis_covchisqdir
	nowfile=outputdir+'/'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
	#print '    Check carefully whether systematic correctoni properly done!'; return
	print '    Output chisq to file: ', nowfile
	nowf.write('### imumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor\n')
	### 6. Compute
	for imumin in imumins:
		#### 6.1 Load in 
		##### 6.1.1 Mock
		xis_mock = smu_xis_loadmockrlt(imumin, 
					      print_nummock=False)  #fmt: catname, catname2, RSDstr, ibin
		##### 6.1.2 Data
		xis_data_file = smu_xis_xisdir+'/'+baseoutputfile+'.smax'+str(smax)+'.imumin%03i'%imumin+'.xis'
		xis_data = smu_xis_loaddatarlt(xis_data_file) # fmt: catname, ibin, omwstr
		#print xis_data['DR12v4-CMASS'][1][omwstr(0.31,-1)]
		
		#### 6.2 Compuate chisqs for all ...
		chisqs_nosyscor, chisqs_syscor = {}, {} # fmt: omwstr, ibin
		iomw = 0
		covmats = []

		for omw in omwlist:
			om,w = omw; nowomwstr=omwstr(om,w); chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr] = [], []
			i_redshiftbin = 0
			for catname in ['DR12v4-LOWZ', 'DR12v4-CMASS', ]:
				for ibin in range(totbin):
					if i_redshiftbin in refbins:
						## A. picku up the xis
						xidata_base = xis_data[catname][ibin][nowomwstr]
						if not use_omw_as_sys:
							xisys_base  = xis_mock[catname+sys_catname_suffix][sys_catname2][sys_RSDstr][ibin]
						else:
							omsys, wsys = use_omw_as_sys
							xisys_base = [xis_data[catname][ibin][omwstr(omsys,wsys)]]
						## B. rebin it in s1:s2;   normalize the amplitude
						xihatdata_base =   normfun(array_fracbin(XmultY(sarray, xidata_base[s1:s2]), sbins))
						xihatsys_base  = [ normfun(array_fracbin(XmultY(sarray, X[s1:s2]), sbins)) for X in xisys_base]
						if iomw == 0:
							xicov_base  = xis_mock[catname][cov_catname2][cov_RSDstr][ibin]
							xihatcov_base  = [ normfun(array_fracbin(XmultY(sarray, X[s1:s2]), sbins)) for X in xicov_base]
						covmats.append([])
					else:
						## A. picku up the xis
						xidata = xis_data[catname][ibin][nowomwstr]
						if not use_omw_as_sys:
							xisys  = xis_mock[catname+sys_catname_suffix][sys_catname2][sys_RSDstr][ibin]
						else:
							omsys, wsys = use_omw_as_sys
							xisys = [xis_data[catname][ibin][omwstr(omsys,wsys)]]
						## B. rebin it in s1:s2;   normalize the amplitude
						xihatdata =   normfun(array_fracbin(XmultY(sarray, xidata[s1:s2]), sbins))
						xihatsys  = [ normfun(array_fracbin(XmultY(sarray, X[s1:s2]), sbins)) for X in xisys]
						dxidata = get_diffarray(xihatdata_base, xihatdata)
						dxisys  = [ get_diffarray(xihatsys_base[row], xihatsys[row]) \
							for row in range(len(xihatsys))]
						dxisys = get_avg_array(dxisys)
						dxisys = polyfitY(array_fracbin(sarray,sbins), dxisys, deg=3)

						if iomw == 0:
							xicov  = xis_mock[catname][cov_catname2][cov_RSDstr][ibin]
							xihatcov  = [ normfun(array_fracbin(XmultY(sarray, X[s1:s2]), sbins)) for X in xicov]
							dxicov  = [ get_diffarray(xihatcov_base[row], xihatcov[row]) \
								for row in range(len(xihatcov))]
							covmats.append(get_covmat(transpose(dxicov)))
						covmat = covmats[i_redshiftbin]
						chisq_nosyscor, like = chisq_like_cov_xbar(dxidata, covmat)
						chisq_syscor, like   = chisq_like_cov_xbar(get_diffarray(dxidata,dxisys), covmat)
						if i_redshiftbin not in calcbins:
							chisq_nosyscor = chisq_syscor = 0.0
						#if i_redshiftbin in [2,5]:
						#	chisq_nosyscor, chisq_syscor = 0, 0
						chisqs_nosyscor[nowomwstr].append(chisq_nosyscor)  
						chisqs_syscor[nowomwstr].append(chisq_syscor) 
					i_redshiftbin += 1
			X1, X2 = chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
			str1, str2, strA, strB = \
				array_to_str(X1), array_to_str(X2), str(sum(chisq_nosyscor)), str(sum(chisq_syscor))
			nowf.write(str(imumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2+'\n')
			iomw += 1	
			#print om, w, chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
		#print chisqs_nosyscor.keys()
		#sys.exit(' ')
		chisqlist_nosyscor = [ [  sum(chisqs_nosyscor[omwstr(om,w)])    for om in omlist] for w in wlist ]
		chisqlist_syscor =   [ [  sum(chisqs_syscor[omwstr(om,w)])      for om in omlist] for w in wlist ]
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
			fig.suptitle(nowchisqstr+'.imumin'+str(imumin))
			if savefig: 
				figname = outputdir+'/'+nowchisqstr+'.imumin'+str(imumin)+'.png';
				print '   Figure saved: ', figname
				plt.savefig(figname, fmt='png')
			if showfig: plt.show()

