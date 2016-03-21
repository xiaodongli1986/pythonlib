from scipy import stats

def nsigofsigCL(sig):
	if sig == sig1:
		return 1
	elif sig == sig2:
		return 2
	elif sig == sig3:
		return 3
	elif sig == sig4:
		return 4
	elif sig == sig5:
		return 5
	elif sig == sig6:
		return 6
	elif sig == sig7:
		return 7
	else:
		return stats.norm.interval(sig, loc=0.0, scale=1.0)[1]		

### Convert    a x^2 + b x + c     to    (x-x0)^2/sig^2 + chisq0 
def abc_to_f0sigchisq0(a, b, c):
    f0 = - b / (2.0*a)
    sig = 1.0 / np.sqrt(float(a))
    chisq0 = c - b*b/4.0/a
    return f0, sig, chisq0
def chisq__f0_sig_chisq0(f, f0, sig, chisq0=0):
    return (f-f0)**2.0 / sig**2.0 + chisq0
def chisq__f0_sig_chisq0_str(f0, sig, chisq0, f0sigfmt='%.5f', chisq0fmt='%.1f' ):
    return ('chi^2 = (f-'+f0sigfmt%f0+')^2 / '+f0sigfmt%sig+'^2 + '+chisq0fmt%chisq0)  


def binned_quan(X, Y, nbin, Xrange=''): ## binned quantities
    Xs, Ys = [[] for row in range(nbin)], [[] for row in range(nbin)]
    xmin,xmax = min(X),max(X); deltax = (xmax-xmin)/float(nbin)
    if Xrange != '':
	xmin,xmax = Xrange; deltax = (xmax-xmin)/float(nbin)
    for row in range(len(X)):
        x = X[row]; y = Y[row];
	if x<xmin or x>xmax:
		continue
        ibin = int((x-xmin)/deltax); ibin = max(ibin,0);ibin=min(ibin,nbin-1)
        Xs[ibin].append(x); Ys[ibin].append(y)
    return Xs, Ys
def binned_quan_meaner(X, Y, nbin, Xrange=''): ## binned quantities
    Xs, Ys = binned_quan(X,Y,nbin,Xrange)
    Xavs, Xavers, Yavs, Yavers = [],[],[],[]
    for ibin in range(nbin):
	if len(Xs[ibin]) > 2:
	   xav, xvar, xaver = get_stat_from_list(Xs[ibin])
	   yav, yvar, yaver = get_stat_from_list(Ys[ibin])
	   Xavs.append(xav); Xavers.append(xaver)
	   Yavs.append(yav); Yavers.append(yaver)
    return Xavs, Xavers, Yavs, Yavers

def binned_quan(X, Y, nbin, Xrange=''): ## binned quantities
    Xs, Ys = [[] for row in range(nbin)], [[] for row in range(nbin)]
    xmin,xmax = min(X),max(X); deltax = (xmax-xmin)/float(nbin)
    if Xrange != '':
	xmin,xmax = Xrange; deltax = (xmax-xmin)/float(nbin)
    for row in range(len(X)):
        x = X[row]; y = Y[row];
	if x<xmin or x>xmax:
		continue
        ibin = int((x-xmin)/deltax); ibin = max(ibin,0);ibin=min(ibin,nbin-1)
        Xs[ibin].append(x); Ys[ibin].append(y)
    return Xs, Ys

### covariance matrix
def get_covmat(Xs):
	n = len(Xs)
	Xibars = [0 for row in range(n)]
	XiXjbars = [ [0 for row1 in range(n)] for row2 in range(n)]
	covmat = [ [0 for row1 in range(n)] for row2 in range(n)]
	ndat = len(Xs[0])
	for i in range(n):
		## \bar Xi
		for idat in range(ndat):
			Xibars[i] += Xs[i][idat]
		Xibars[i] /= float(ndat)
		## \bar (XiXj)
		for j in range(i,n):
			for idat in range(ndat):
				XiXjbars[i][j] += (Xs[i][idat] * Xs[j][idat])
			XiXjbars[i][j] /= float(ndat)
	for i in range(n):
		for j in range(i,n):
			XiXjbars[j][i] = XiXjbars[i][j]
	#print 'XiXjbars: '
	#for i in range(n):
	#	print XiXjbars[i]
	for i in range(n):
		for j in range(n):
			covmat[i][j] = (XiXjbars[i][j] - Xibars[i]*Xibars[j]) * ndat / float(ndat-1.0)
#				= E ((Xi -Xibar)(Xj-Xjbar) ) = E(XiXj) - E(Xi)*Xjbar - E(Xj) * Xibar + Xibar * Xjbar
	return covmat

### ps. code to test covmat function: python + Mathematica

#X0 = [random.gauss(0,1) for i in range(1,10001)]
#X1 = [random.gauss(0,1) for i in range(1,10001)]
#X2 = [random.gauss(0,1) for i in range(1,10001)]
#for row in range(len(X0)):
#    #X1[row] += X0[row]
#    X2[row] += X0[row]
#np.savetxt('X0.txt', X0)
#np.savetxt('X1.txt', X1)
#np.savetxt('X2.txt', X2)
#covmat = get_covmat([X0,X1,X2])
#for i in range(len(covmat)):
#    print covmat[i]
    
    #X0 = Import[
    #   "/home/lixiaodong/SparseFilaments/code/A_Standard_Models/B_Study_\
		    #GF/X0.txt", "Data"];
		    #X1 = Import[
		    #   "/home/lixiaodong/SparseFilaments/code/A_Standard_Models/B_Study_\
				    #GF/X1.txt", "Data"];
				    #X2 = Import[
				    #   "/home/lixiaodong/SparseFilaments/code/A_Standard_Models/B_Study_\
						    #GF/X2.txt", "Data"];
						    #M = Table[{X0[[i, 1]], X1[[i, 1]], X2[[i, 1]]}, {i, 1, Length[X0]}];
						    #Covariance[M] // MatrixForm    


def chisq_like_cov_xbar(X, Cov, Xbar = [], chisq0=0):
    p=len(X);
    npcov = np.mat(Cov);
    invcov = npcov.I;
    if Xbar == []:
        Xdiff = X
    else:
        Xdiff = [X[row]-Xbar[row] for row in range(p)]
    #print 'Xdiff:', Xdiff
    #print 'covmat: ', npcov
    #print 'invcov:', invcov
    detcov = np.linalg.det(npcov);
    chisq = 0.0;
    for i in range(p):
        for j in range(p):
            chisq += Xdiff[i]*invcov[i,j]*Xdiff[j]
    like = (2.0*np.pi)**(-p/2.0) * (detcov**(-0.5)) * exp(-0.5 * (chisq-chisq0));
    return chisq, like

def chisq_invcov_xbar(X, InvCov, Xbar = [], chisq0=0):
    p=len(X);
    #invcov = np.mat(InvCov);
    if Xbar == []:
        Xdiff = X
    else:
        Xdiff = [X[row]-Xbar[row] for row in range(p)]
    #print 'Xdiff:', Xdiff
    #print 'covmat: ', npcov
    #print 'invcov:', invcov
    chisq = 0.0;
    for i in range(p):
        for j in range(p):
            chisq += Xdiff[i]*InvCov[i,j]*Xdiff[j]
    return chisq

def get_er_erofer(X):
    x_av, x_var, x_av_er, x_var_er = get_stat_from_list(X)
    erofvar = x_var_er
    var1sup = x_var + erofvar
    var1slow = x_var - erofvar
    sqrtvar1sup = np.sqrt(x_var + erofvar)
    sqrtvar1slow = np.sqrt(x_var - erofvar)
    sqrtvarcenter = (sqrtvar1sup+sqrtvar1slow)/2.0
    sqrtvarcenter2 = np.sqrt(x_var)
    sqrtvarer = (sqrtvar1sup-sqrtvar1slow)/2.0
    print 'x = ', x_av, '+/-', np.sqrt(x_var)
    #print 'Error bar = ', sqrtvarcenter, '+/-', sqrtvarer, '(', sqrtvarcenter-sqrtvarer,sqrtvarcenter+sqrtvarer,')'
    #print 'Error bar = ', sqrtvarcenter2, '+/-', sqrtvarer
    print 'Correct estimation of x: ', x_av, '+/-', np.sqrt(x_var)+sqrtvarer
    return sqrtvarcenter2, sqrtvarer



###############################
### chisq calculation
###############################

### calculate chisq from Y
def get_chisq_from_yyerlist(Y, YEr, diffchisqmethod = 'use_weightedavg_as_ref'):
    if diffchisqmethod == 'use_weightedavg_as_ref':
            Ymean = meanY(Y)
	    chisq = sum( [(Y[row]-Ymean)**2.0/YEr[row]**2.0 for row in range(len(Y))] )
    elif diffchisqmethod == 'use_lastrbin_as_ref':
     	    iref = len(Y)-1
    	    YErsq = YEr[iref]**2.0
	    chisq = sum( [(Y[row]-Y[iref])**2.0/(YEr[row]**2.0+YErsq) for row in range(len(Y))] )
    elif diffchisqmethod < len(Y):
    	    YErsq = YEr[diffchisqmethod]**2.0
	    chisq = sum( [(Y[row]-Y[diffchisqmethod])**2.0/(YEr[row]**2.0+YErsq) for row in range(len(Y))] )
    else:
    	    print 'Error (get_chisq_from_yyerlist)! Wrong diffchisqmethod: ', diffchisqmethod, \
    	    	'; must be use_weightedavg_as_ref or a number smaller than ', len(Y)
    	    return
    return chisq

### calculate chisqs from list of Ys
def get_chisqs_from_yyerlists(Ys, YErs, diffchisqmethod = 'use_weightedavg_as_ref'):
    return [get_chisq_from_yyerlist(Ys[row],YErs[row],diffchisqmethod=diffchisqmethod) for row in range(len(Ys))]
    
##############################################################
## Useful functions used to determine CL, error bar;
##  dealing with likelihood
##############################################################

def get_cut(Xlike, CL=0.683):
    import copy
    copylike = copy.copy(Xlike);
    copylike.sort(reverse=True);
    totlike = sum(copylike);
    goalike = totlike * CL;
    nowlike = 0
    for i in range(len(copylike)):
        nowlike += copylike[i]
        if nowlike > goalike:
            break
    cut1, like1 = copylike[i], nowlike
    cut2, like2 = copylike[i-1], nowlike - copylike[i]
    cut = cut1 + (goalike-like1)/(like2-like1)*(cut2-cut1)
    return cut
    
def get_erbar(X, Xlike, CL=0.683, printinfo=False):
    lablist = [];
    CL = get_cut(Xlike, CL)
    if printinfo:
    	print 'Determine cut of likelihood as ', CL
    for row in range(len(X)):
        if Xlike[row] > CL:
            lablist.append(row)
    minlab = min(lablist);
    maxlab = max(lablist);
    if printinfo:
    	print 'Determine minlab/minlabquan/minlablike as ', minlab, X[minlab], Xlike[minlab]
    	print 'Determine maxlab/minlabquan/maxlablike as ', maxlab, X[maxlab], Xlike[maxlab]
    if minlab > 0:
        xleft = X[minlab-1] + (X[minlab]-X[minlab-1]) * ((CL-Xlike[minlab-1])/(Xlike[minlab]-Xlike[minlab-1]))
    else:
        xleft = X[minlab] + (X[minlab+1]-X[minlab]) * ((CL-Xlike[minlab])/(Xlike[minlab+1]-Xlike[minlab]))
    if maxlab < len(Xlike) -1:
        xright = X[maxlab] + (X[maxlab+1]-X[maxlab]) * ((CL-Xlike[maxlab])/(Xlike[maxlab+1]-Xlike[maxlab]))
    else:
        xright = X[maxlab-1] + (X[maxlab]-X[maxlab-1]) * ((CL-Xlike[maxlab-1])/(Xlike[maxlab]-Xlike[maxlab-1]))
    xbf = X[index_of_max(Xlike)]
    return xbf, xleft, xright

def integral_like(chisqlist, chisqcut):
        like = 0
        numrow1 = len(chisqlist); numrow2 = len(chisqlist[0]);
        for row1 in range(numrow1):
            for row2 in range(numrow2):
                if chisqlist[row1][row2] < chisqcut:
                    like += np.exp(-chisqlist[row1][row2]*0.5)
        return like
    
def find_CL_chisqcut(chisqlist, CL, chisq1 = 0, chisq2 = 30):
        totlike = integral_like(chisqlist, 1.0e30)
        like0 = totlike * CL; chisqA = chisq1; chisqB = chisq2;
        likeA = integral_like(chisqlist, chisqA)
        likeB = integral_like(chisqlist, chisqB)
        if not(likeA < like0 < likeB):
            print 'ERROR!: likeA, like0, likeB = ', likeA, like0, likeB
            return
        while abs(chisqA - chisqB) > 0.01:
            chisqC = chisqA + (like0-likeA)/(likeB-likeA)*(chisqB-chisqA)
            likeC = integral_like(chisqlist, chisqC)
            if likeC > like0:
                likeB = likeC; chisqB = chisqC;
            else:
                likeA = likeC; chisqA = chisqC;
        return chisqC, likeC/totlike  
    
def list_find_CL_chisqcut(chisqlist, CLlist, chisq1 = 0, chisq2 = 30):
    chisq1d = [];
    for row in range(len(chisqlist)):
        chisq1d += chisqlist[row]
    chisqmin = min(chisq1d)
    like1d = [exp(-(chisq1d[row]-chisqmin)/2.0) for row in range(len(chisq1d))]
    chisqcutlist = []
    for CL in CLlist:
        cutlike = get_cut(like1d, CL)
        chisqcut = chisqmin - 2.0 * log(cutlike)
        chisqcutlist.append(chisqcut)
    return chisqcutlist  
    


def get_margconstraint(chisqlist, omlist, wlist, do_smooth=False, smsigma=0.8, smorder=0, printinfo=False):    
    # w, om is the first, second quantity. Anyway just names. not important.
    numw = len(chisqlist); numom = len(chisqlist[0]);
      
    # Marginalized constraints
    margwlike = [0 for row2 in range(numw)] 
    margomlike = [0 for row2 in range(numom)]
    nowlikelist = [[0 for row in range(numom)] for row in range(numw)]
    chisqmin = 1.0e30
    for row1 in range(numw):
        for row2 in range(numom):
            chisqmin= min(chisqmin, chisqlist[row1][row2])
    for row1 in range(numw):
        for row2 in range(numom):
            nowlikelist[row1][row2] = exp(-0.5*(chisqlist[row1][row2]-chisqmin)) 
    # Marginalized constraint on Omegam
    for iom in range(numom):
        margomlike[iom] = sum([nowlikelist[row][iom] for row in range(numw)])
        # Some times there are some too small values, leading to numerical error.
        if not margomlike[iom] > 1.0e-8:
            margomlike[iom] = random.uniform(1.0e-9,1.0e-8) 
    maxlike = max(margomlike[0:numom])
    for iom in range(numom):
        margomlike[iom] /= maxlike
    Z = margomlike[0:numom]
    if do_smooth:
        Z = ndimage.gaussian_filter(Z, sigma=smsigma, order=smorder)
    maxZ = max(Z); Z = [Z[row]/maxZ for row in range(len(Z))]
    ombf,oml,omr = get_erbar(omlist, Z, CL=0.683,printinfo=printinfo)
    if printinfo:
        print 'Constraints on Omegam: ', ombf, '+', omr-ombf, '-', ombf-oml
    
    # Marginalized constraint on w
    for iw in range(numw):
        margwlike[iw] = sum([nowlikelist[iw][row] for row in range(numom)])
        if not margwlike[iw] > 1.0e-8:
            margwlike[iw] = random.uniform(1.0e-9,1.0e-8)
    maxlike = max(margwlike[0:numw])
    for iw in range(numw):
        margwlike[iw] /= maxlike
    Z = margwlike[0:numw]
    if do_smooth:
        Z = ndimage.gaussian_filter(Z, sigma=smsigma, order=smorder)
    maxZ = max(Z); Z = [Z[row]/maxZ for row in range(len(Z))]
    wbf,wl,wr = get_erbar(wlist, Z, CL=0.683,printinfo=printinfo)
    
    if printinfo:
        print 'Constraints on w: ', wbf, '+', wr-wbf, '-', wbf-wl
    
    return margwlike, margomlike,wbf, wl, wr, ombf, oml, omr

### Useful test for plot_contour  and  marg like (seems there is serious problem -- improve it later!!!)

#omlist = np.linspace(-5, 5, 70)
#wlist = np.linspace(-5, 5, 80)
#chisqlist = [[(omlist[iom])**2.0 + (wlist[iw])**2.0  for iom in range(len(omlist))] for iw in range(len(wlist))]
#margwlike, margomlike, wbf, wl, wr, ombf, oml, omr = get_margconstraint(chisqlist, omlist, wlist)
#fig, ax = figax()
#plot_contour(ax, omlist, wlist, chisqlist, ommin=min(omlist), ommax=max(omlist), wmin=min(wlist), wmax=max(wlist), )
#ax.grid()


########################################
### 2d linear interpolations
def LinearInterpolation_2d(x1, y1,  x2, y2,  f_x1y1, f_x1y2, f_x2y1, f_x2y2,   x3, y3):
    '''
        Value of f(x3, y3) from f(x1,y1), f(x1,y2), f(x2,y1), f(x2,y2)
        Method = Linear Interpolation in 2d
        
        ### Testing code:
        
            import stdA; execfile(stdA.pyfile);
            deltaxy = 0.1
            x1, y1, x2, y2 = 1, 1, 1+deltaxy, 1+deltaxy
            x3, y3 = 1+deltaxy*0.17, 1+deltaxy*0.5
            def fun(x,y):
                return x**2+y*100
            print LinearInterpolation_2d(x1, y1, x2, y2, fun(x1,y1), fun(x1,y2), fun(x2,y1), fun(x2,y2), x3, y3), fun(x3,y3)
            print 'fun(x1,y1), fun(x1,y2), fun(x2,y1), fun(x2,y2) = ', fun(x1,y1), fun(x1,y2), fun(x2,y1), fun(x2,y2)
    '''
    fracx = (x3-x1) / (x2-x1)
    fracy = (y3-y1) / (y2-y1)
    #print fracx, fracy
    f_x3y1 = f_x1y1 + (f_x2y1-f_x1y1) * fracx
    f_x3y2 = f_x1y2 + (f_x2y2-f_x1y2) * fracx
    f_x3y3 = f_x3y1 + (f_x3y2-f_x3y1) * fracy
    #print 'f_x3y1, f_x3y2, f_x3y3 = ', f_x3y1, f_x3y2, f_x3y3
    return f_x3y3

def LinearInterpolation_2d_UniformGrid(xmin, xmax, numx,   ymin, ymax, numy,   F_2d,   X3, Y3):
    '''
        Values of f(x3, y3) from F_2d
        F_2d: 2d list with numx rows, numy columns
        X3, Y3: 1d list with same length
        Method = Linear Interpolation in 2d
        
        A useful sentence:
            LinearInterpolation_2d_UniformGrid(min(X), max(X), len(X), min(Y),max(Y),len(Y), 
                chisqlist, X3, Y3)
    '''
    if len(X3) != len(Y3):
        print 'Error (LinearInterpolation_2d_UniformGrid) !!!! Length of X3, Y3 not matching: ', len(X3), len(Y3)
        return
    deltax = (xmax-xmin) / (numx - 1.0)
    deltay = (ymax-ymin) / (numy - 1.0)
    F3s = []
    for row in range(len(X3)):
        x, y = X3[row], Y3[row]
        ### ix and iy
        ix, iy = int((x-xmin)/deltax), int((y-ymin)/deltay)
        ix = min(ix, numx-2)
        iy = min(iy, numy-2)
        ix = max(ix,0); iy = max(iy,0)
        ### x1,x2,y1,y2 and the fs
        x1,y1 = xmin+ix*deltax, ymin+iy*deltay
        x2,y2 = x1+deltax, y1+deltay
        f_x1y1, f_x1y2, f_x2y1, f_x2y2 = F_2d[iy][ix], F_2d[iy+1][ix], F_2d[iy][ix+1], F_2d[iy+1][ix+1]
        ### 
        f_xy = LinearInterpolation_2d(x1, y1,  x2, y2,  f_x1y1, f_x1y2, f_x2y1, f_x2y2,   x,y)         
        F3s.append(f_xy)
    return F3s

########################################
########################################
### Tools for cosmomc

def COSMOMC_load_contours(filename, upper_limit_of_loadin=100, upper_limit_of_loadiny=100):
    ''' 
	Load in the cosmomc contour boundary outputed by matlab
	To have outputed boundary in cosmomc add this in Getdist.f90: (cosmomc version: October 12; delete the three slashes)
		              if (PlotContMATLAB(aunit,rootname,j,j2,do_shading) .and. &
                  num_contours /= 0) then

               write (aunit,'(a)') '[C h] = contour(x1,x2,pts,cnt,lineM{1});'
               write (aunit,'(a)') 'set(h,''LineWidth'',lw1);'
               write (aunit,*) 'hold on; axis manual; '

                plotfile = numcat(trim(numcat(trim(rootname)//'_2D_',colix(j)-2))//'_',colix(j2)-2) !xiaodongli
                write(aunit,'(a)') 'csvwrite(\'\'\'//trim(plotfile)//'.2dcon.txt'',C'');' !xiaodongli 
'''
    data = np.loadtxt(filename, delimiter=',')
    contours = []
    now_index = -1
    for row in range(len(data)):
        if data[row][0] > upper_limit_of_loadin or data[row][1] > upper_limit_of_loadiny:
            now_index += 1
            contours.append([])
        else:
            contours[now_index].append([x for x in data[row]])
    return contours


def MCMC_cosmomc_fmt_convert(MCMCfile, suffix = '.CosmomcFmtConverted', ipt_outputfile = '', fmtstr = '%20.10e'):
	'''Convert a chisq result file into the format of cosmomc.
	    fmt before conversion: par1, par2, ..., chisq
	    fmt after conversion:  weight, chisq/2, par1, par2, ...
	'''
	### determine the outputfile name
	if ipt_outputfile == '':
		outputfile = MCMCfile+suffix
	else:
		outputfile = ipt_outputfile
	print '(MCMC_cosmomc_fmt_convert) Open: \n\t', MCMCfile, '\nOutput to: \n\t', outputfile
	f1 = open(MCMCfile, 'r')
	f2 = open(outputfile, 'w')

	### find out the minimal chisq
	data = np.loadtxt(MCMCfile);
	chisqs = Xfromdata(data, len(data[0])-1);
	chisqmin = min(chisqs)

	### generate the new file
	nlines = 0
	while True:
		nowstr = f1.readline()
		if nowstr == '':
			break
		A = str_to_numbers(nowstr)
		chisq = A[len(A)-1]
		B = [np.exp(-0.5*(chisq-chisqmin)), chisq/2.0] + A[0:len(A)-1]
		f2.write(fmtstrlist(B, fmtstr=fmtstr)+'\n')
		nlines += 1
	print '	Finishing processing ', nlines, 'lines.'
	return

def MCMC_extension(MCMCfile, keyfun, fmtstr='', prefix = '', suffix = '.extended', ipt_outputfile = '', only_new_quan=True):
	nlines = 0
	if ipt_outputfile == '':
		outputfile = prefix+MCMCfile
	else:
		outputfile = ipt_outputfile
	f1 = open(MCMCfile, 'r')
	f2 = open(outputfile, 'w')
	print '(MCMC_extension) Open: \n\t', MCMCfile, '\nOutput to: \n\t', outputfile
	while True:
		nowstr = f1.readline()
		if nowstr == '':
			break
		A = str_to_numbers(nowstr)
		B = keyfun(A)
		if only_new_quan:
			f2.write(fmtstrlist([A[0],A[1]]+B, fmtstr=fmtstr)+'\n')		
		else:
			f2.write(nowstr[0:len(nowstr)-1]+' '+fmtstrlist(B, fmtstr=fmtstr)+'\n')
		nlines += 1
	print 'Finishing processing ', nlines, 'lines.'
	return

def COSMOMC_post_chisq(MCMCfile, 
                       chisqfun, 
                       rejectfun=None,
                       suffix = '_post', CharSkip_of_Suffix=6, 
                       ipt_outputfile = '', 
                       fmtstr = '%20.10e'):
    '''
        Add an additional chisq value of the MCMC chain.
        
        ### This tests COSMOMC_post_chisq: post a x^2+y^2 to constant chisq...
        import stdA; execfile(stdA.pyfile)
        import Tpcftools; execfile(Tpcftools.pyfile)
        import bossdatamock; execfile(bossdatamock.pyfile)
        import copy
        
        if False:
            MCMCfile = bossdatamock_cosmomcchaindir+'../const_chisq_1.txt'
            nowf = open(MCMCfile, 'w')
            for i in range(100000):
                x, y = random.uniform(-3,3), random.uniform(-3,3)
                nowf.write('1.0 1.0 '+str(x)+' '+str(y)+'\n')
            nowf.close()

            def chisqfun(Par):
                return Par[0]**2.0 + Par[1]**2.0

            COSMOMC_post_chisq(MCMCfile, chisqfun, suffix='_post_xsqysq')
            print 'Please use COSMOMC to analyze the outputed file to check whether it is 2d standard gaussian!'

    '''
    if ipt_outputfile == '':
        outputfile = MCMCfile[0:len(MCMCfile)-CharSkip_of_Suffix]+suffix \
            +MCMCfile[len(MCMCfile)-CharSkip_of_Suffix:len(MCMCfile)]
    else:
        outputfile = ipt_outputfile
        
    print 'Open: \n\t', MCMCfile, '\nOutput to: \n\t', outputfile
    f1 = open(MCMCfile, 'r')
    f2 = open(outputfile, 'w')

    ### Main loop
    nlines = 0
    nrejected = 0
    while True:
        
        nowstr = f1.readline()
        if nowstr == '':
            break
        nlines += 1
        
        ### Read in the numbers
        A = str_to_numbers(nowstr)
        Weight = A[0]
        Par = A[2:len(A)]
        Loglike = A[1]
        
        if rejectfun != None:
            if rejectfun(Par):
                nrejected += 1
                continue
        
        ### Add the addtional chisq
        AddChisq = chisqfun(Par)
        Loglike += AddChisq / 2.0
        Weight *= exp(  - AddChisq / 2.0)
        B = [Weight, Loglike] + Par
        f2.write(fmtstrlist(B, fmtstr=fmtstr)+'\n')
        
    print '	Finishing processing ', nlines, 'lines;     ', nrejected, ' rejected.'
    return


