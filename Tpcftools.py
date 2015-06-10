#!/home/xiaodongli/software/anaconda/bin/python
# Filename: Tpcftools.py 

import commands
pythonlibPATH = commands.getoutput('echo $pythonlibPATH')
pyfile=pythonlibPATH+'/'+'Tpcftools.py'

### Settings for smu

def smu__rowofsmu(i_s, i_mu, smusettings):
    """		Row number of given i_s, i_mu; start from 0 [0 means first row]"""
    return i_s*smusettings['nummubin'] + i_mu

def smu__get_srange(i_s, smusettings):
    """		Range of s with given i_s"""
    return smusettings['smin'] + smusettings['deltas'] * (i_s), \
		smusettings['smin'] + smusettings['deltas'] * (i_s+1)
def smu__get_murange(i_mu, smusettings):
    """		Range of mu with given i_mu"""
    return smusettings['mumin'] + smusettings['deltamu'] * (i_mu), \
		smusettings['mumin'] + smusettings['deltamu'] * (i_mu+1)


def smu__get_mid_smu(i_s, i_mu, smusettings):
    """		Middle value of s, mu with given i_s, i_mu"""
    return smusettings['smin'] + smusettings['deltas'] * (i_s+0.5), smusettings['mumin'] + smusettings['deltamu']*(i_mu+0.5)
def smu__get_mids(i_s, smusettings, avgdim=1):
    """		Middle value of s in given dimension, given i_s"""
    s1, s2 = smu__get_srange(i_s, smusettings)
    return dimavg(s1, s2, avgdim)
def smu__get_midmu(i_mu, smusettings, avgdim=1):
    """		Middle value of s in given dimension, given i_s"""
    mu1, mu2 = smu__get_murange(i_mu, smusettings)
    return dimavg(mu1, mu2, avgdim)
def smu__get_imuedges(muedges,smusettings):
    """		Index of mubin given edges of mu"""
    i_muedges = [[min(max(0,  int((muedges[row]-smusettings['mumin'])/smusettings['deltamu'])  ),smusettings['nummubin']), min(max(0,  int((muedges[row+1]-smusettings['mumin'])/smusettings['deltamu'])-1  ),smusettings['nummubin'])] for row in range(len(muedges)-1)]
    return i_muedges
def smu__get_muedgemids(i_muedges, smusettings):
    """		Middle value of mu in bins"""
    return [(smu__get_midmu(i_muedge[0], smusettings)+smu__get_midmu(i_muedge[1], smusettings))*0.5 
                      for i_muedge in i_muedges]


def smu__initsmusettings(smusettings):
	"""	Initializing deltas, deltamu, slist, mulist for smusettings"""
	global mumin, mumax, smin, smax, mulist, slist, deltamu, deltas, nummubin, numsbin
	smusettings['deltas'] = (smusettings['smax']-smusettings['smin'])/float(smusettings['numsbin']); 
	smusettings['deltamu'] = (smusettings['mumax']-smusettings['mumin'])/float(smusettings['nummubin']);
	smusettings['slist'] = [smu__get_mid_smu(i_s,0, smusettings)[0] for i_s in range(smusettings['numsbin'])]
	smusettings['mulist'] = [smu__get_mid_smu(0,i_mu, smusettings)[1] for i_mu in range(smusettings['nummubin'])]
	mumin = smusettings['mumin']
	mumax = smusettings['mumax']
	smin = smusettings['smin']
	smax = smusettings['smax']
	mulist = smusettings['mulist']
	slist = smusettings['slist']
	deltamu = smusettings['deltamu']
	deltas = smusettings['deltas']
	nummubin = smusettings['nummubin']
	numsbin = smusettings['numsbin']

def smu__initsmusettings_quick(smin_ipt=0, smax_ipt=100, numsbin_ipt=100, mumin_ipt=0, mumax_ipt=1, nummubin_ipt=100):
	"""	Initializing deltas, deltamu, slist, mulist for smusettings; quick """
	global mumin, mumax, smin, smax, mulist, slist, deltamu, deltas, nummubin, numsbin
	smusettings = {}
	smusettings['smin']=smin=smin_ipt
	smusettings['smax']=smax=smax_ipt
	smusettings['numsbin']=numsbin=numsbin_ipt
	smusettings['mumin']=mumin=mumin_ipt
	smusettings['mumax']=mumax=mumax_ipt
	smusettings['nummubin']=nummubin=nummubin_ipt
	smusettings['deltas'] = (smusettings['smax']-smusettings['smin'])/float(smusettings['numsbin']); 
	smusettings['deltamu'] = (smusettings['mumax']-smusettings['mumin'])/float(smusettings['nummubin']);
	smusettings['slist'] = [smu__get_mid_smu(i_s,0, smusettings)[0] for i_s in range(smusettings['numsbin'])]
	smusettings['mulist'] = [smu__get_mid_smu(0,i_mu, smusettings)[1] for i_mu in range(smusettings['nummubin'])]
	mumin = smusettings['mumin']
	mumax = smusettings['mumax']
	smin = smusettings['smin']
	smax = smusettings['smax']
	mulist = smusettings['mulist']
	slist = smusettings['slist']
	deltamu = smusettings['deltamu']
	deltas = smusettings['deltas']
	nummubin = smusettings['nummubin']
	numsbin = smusettings['numsbin']
	return smusettings

def smu__briefprint(smusettings):
    """		Very briefly print the info of smusettings"""
    print smusettings['numsbin'], 's bins within', (smusettings['smin'], smusettings['smax']),\
        ',',smusettings['numsbin'], 'mu bins within', (smusettings['mumin'], smusettings['mumax'])
    print '\tfirst sbin:', smu__get_srange(0, smusettings)
    print '\tlast mubin:', smu__get_murange(smusettings['nummubin']-1, smusettings)

### Load in 2pcf result of smu

def smu__loadin(smufile, smusettings, icol=''):
    """		Load in 2pcf file in s, mu space
			if icol=='', load in all columns; else only load in icol column
    """
    data = np.loadtxt(smufile)

    if icol == '':
	Tpcfrlt = [[data[smu__rowofsmu(i_s,i_mu,smusettings)] for i_mu in range(smusettings['nummubin'])] for i_s in range(smusettings['numsbin'])]
    else:
	Tpcfrlt = [[data[smu__rowofsmu(i_s,i_mu,smusettings)][icol] for i_mu in range(smusettings['nummubin'])] for i_s in range(smusettings['numsbin'])]

    return Tpcfrlt

def smu__contour(ax,mulist,slist,xilist,xlabel='$\\mu$',ylabel='$s$',labelfs=18,srange=(5,30),title='',tiltfs=15,levels=[], lw=2, ls='-'):
    """		contour plot for xi, in space of s, mu    """
    ax.set_xlabel(xlabel,fontsize=labelfs);ax.set_ylabel(ylabel,fontsize=labelfs);
    if levels==[]:
        CS = ax.contour(mulist,slist, xilist, cmap=plt.cm.Accent, lw=lw, ls=ls )
    else:
        CS = ax.contour(mulist,slist, xilist, levels=levels, cmap=plt.cm.gray, linewidths=lw, linestyles=ls  )
    #CB = plt.colorbar(CS, shrink=0.8, extend='both' )
    #im = ax.imshow(zetalist, interpolation='bilinear', origin='lower',cmap=cm.gray_r,extent=(0,1,srange[0],srange[1]), aspect='auto')
    #CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
    ax.set_ylim(srange)
    ax.clabel(CS, fontsize=18,inline=3)
    if title!='':
        ax.set_title(title,fontsize=tiltfs)

def smu__Tpcfplot_xis(ax, DDlist, DRlist, RRlist, deltais=2, nowleg=''):
        ismin = 0; 
        imumin = 0; 
        imumax = nummubin;
	sasx=[]
	packedxiasy=[]
        while ismin < numsbin:
            ismax = ismin + deltais;
            ismax = min(ismax, numsbin-1)
            nows = (slist[ismin]+slist[ismax])/2.0
            sasx.append(nows)
            packedxiasy.append(nows*nows*packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
            ismin += deltais;

        #ax.set_title(nowfile, fontsize=12)
        ax.plot(sasx, packedxiasy, marker='o', markersize=1, lw=3, label = nowleg)
        ax.set_xlabel('$r\ [\\rm Mpc/h]$', fontsize=25)
        ax.set_ylabel('$r^2\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
        #ax.set_xlim(0,150)
        #ax.set_xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150])
        #ax.set_yticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150])
        #ax.set_ylim(0, 120)
        ax.grid()
        ax.legend(frameon=False)
	return ax

### key quantities representing properties of 2pcf (integral of \xi over a range of s;  and so on)

def get_funame(function):
	totstr = str(function)
	nowstr = ''
	for i in range(10, len(totstr)):
		if totstr[i] == ' ':
			break
		else:
			nowstr += totstr[i]
	return nowstr
	
def packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax):
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    dd = packed_count(DDlist, ismin, ismax, imumin, imumax)
    dr = packed_count(DRlist, ismin, ismax, imumin, imumax)
    rr = packed_count(RRlist, ismin, ismax, imumin, imumax)
    return calc_xils(dd, dr, rr, dim=0)

def calc_xils(DD, DR, RR,dim=0):
    """	calculate xiLS from DD, RR..."""
    if dim == 0:
	return (DD-2*DR)/RR+1
    elif dim == 1:
	return [(DD[row]-2*DR[row])/RR[row]+1 for row in range(len(DD))]
    elif dim == 2:
	return [[(DD[row1][row2]-2*DR[row1][row2])/RR[row1][row2]+1 for row2 in range(len(DD[row1]))] for row1 in range(len(DD))]

def packed_count(DR, ismin, ismax, imumin, imumax):
    """	count # of pairs within certain range of s, mu"""
    #ismin, ismax, imumin, imumax = get_is(smin), get_is(smax), get_imu(mumin), get_imu(mumax)
    Sumed_DR = 0
    for i_s in range(ismin, ismax+1):
        for i_mu in range(imumin, imumax+1):
            Sumed_DR += DR[i_s][i_mu]
    #print ismin, ismax, imumin, imumax
    return Sumed_DR

def smu__xifunname(xifunction):
    """	Name of the function which calculate the  """
    if xifunction in [intxi, intxi_FractionalS]:
        name =  'int--xi'
    elif xifunction == int_s_square_xi:
        name =  'int--s-square-xi'
    elif xifunction == int_s_xi:
        name =  'int--s-xi'
    elif xifunction == int_s3_xi:
        name =  'int--s3-xi'
    elif xifunction == int_1overs_xi:
        name =  'int--1overs-xi'
    elif xifunction == intxi_rat:
        name =  'int--xi_rat'
    elif xifunction == int_s_square_xi_rat:
        name =  'int--s-square-xi_rat'
    elif xifunction == int_weightedmu_s2xi:
        name =  'int--weightedmu_s2xi'
    elif xifunction == int_weighteds_s2xi:
        name =  'int--weighteds_s2xi'
    elif xifunction == int_weightedmu_xi:
        name =  'int--weightedmu_xi'
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range:
	name =  'int--xi_normed'
    elif xifunction == int_s_squre_xi_normed__by_full_mu_range:
	name =  'int--s-squre-xi_normed'
    elif xifunction == int_s_xi_normed__by_full_mu_range:
	name =  'int--s-xi_normed'
    return name


def intxi(DDlist, DRlist, RRlist, ismin, ismax, i_mumin, i_mumax):
    sumxi = 0
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    if ismin > ismax or i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        sumxi += calc_xils(dd, dr, rr)
        #print s
    return sumxi

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def intxi_mu_over_intxi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        #print 'nowis = ', nowis
        sumxi += calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    #print 'nowis = ', ismin-1
    sumxi += (calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        #print 'nowis = ', ismax+1
        sumxi += (calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        #print 'nowis = ', nowis
        sumxi += calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    #print 'nowis = ', ismin-1
    sumxi += (calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        #print 'nowis = ', ismax+1
        sumxi += (calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 


#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        sumxi += now_s_square*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s_square*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s_square*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi    

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_s3_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        sumxi += now_s_square**1.5 *calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s_square**1.5 *calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s_square**1.5 *calc_xils(dd, dr, rr) * smax_tail)
    return sumxi    

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_weightedmu_s2xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square*nowmu
		sumwei += nowxi*now_s_square
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square*nowmu * smin_tail
		sumwei += nowxi*now_s_square  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square*nowmu * smax_tail
	        sumwei += nowxi*now_s_square  * smax_tail
	        
    return sumxi / sumwei


#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_weightedmu_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*nowmu
		sumwei += nowxi
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*nowmu * smin_tail
		sumwei += nowxi  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*nowmu * smax_tail
	        sumwei += nowxi  * smax_tail
	        
    return sumxi / sumwei

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_weighteds_s2xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square**1.5
		sumwei += nowxi*now_s_square
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square**1.5 * smin_tail
		sumwei += nowxi*now_s_square  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowi_mu in range(i_mumin,i_mumax+1):
        	nowmu = mumin + (nowi_mu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowi_mu], DRlist[nowis][nowi_mu], RRlist[nowis][nowi_mu])
	        sumxi += nowxi*now_s_square**1.5 * smax_tail
	        sumwei += nowxi*now_s_square  * smax_tail
	        
    return sumxi / sumwei

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        sumxi += now_s*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi  

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_1overs_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    i_mumin = max(i_mumin,0); i_mumax = min(i_mumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if i_mumin > i_mumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, i_mumin, i_mumax)
        dr = packed_count(DRlist, nowis, nowis, i_mumin, i_mumax)
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        sumxi += 1.0/now_s*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (1.0/now_s*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (1.0/now_s*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 

def intxi_rat(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax): 
	ismid = (ismin_float+ismax_float)/2.0
	return intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismid,  i_mumin, i_mumax) / intxi_FractionalS(DDlist, DRlist, RRlist, ismid,ismax_float, i_mumin, i_mumax)
	
def int_s_square_xi_rat(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax): 
	ismid = (ismin_float+ismax_float)/2.0
	return int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismid,  i_mumin, i_mumax) / int_s_square_xi(DDlist, DRlist, RRlist, ismid,ismax_float, i_mumin, i_mumax)

def intxi_FractionalS_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax, 
	minimali_mumin=0,maximali_mumax=119):
	return intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax) /\
		intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimali_mumin, maximali_mumax)

def intxi_FractionalS_normed__by_full_mu_range_minimalmuge10(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax, 
	minimali_mumin=10,maximali_mumax=119):
	return intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax) /\
		intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimali_mumin, maximali_mumax)

def int_s_squre_xi_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax, 
	minimali_mumin=0,maximali_mumax=119):
	return int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax) /\
		int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimali_mumin, maximali_mumax)

def int_s_xi_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax, 
	minimali_mumin=0,maximali_mumax=119):
	return int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax) /\
		int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimali_mumin, maximali_mumax)


def settingstr_2pcf(nowsmin=5, nowsmax=100, nownummubin=10, nowminmucut=0.1, RSDcor=True, flexiblesmin=False, flexiblesmax=False, 	
	diffchisqmethod = 'use_lastrbin_as_ref', xifunction = 0):
    nowstr = 's'+str(nowsmin)+'to'+str(nowsmax)+'--'+str(nownummubin)+'mubins--muge%.3f'%nowminmucut+'--RSDCor'+str(RSDcor)
    if flexiblesmin:
        nowstr += '--flexiblesmin'
    if flexiblesmax:
        nowstr += '--flexiblesmax'
    if diffchisqmethod == 'use_weightedavg_as_ref':
        nowstr += '.use_weightedavg_as_ref'
    if xifunction == int_s_square_xi:
        nowstr += '.xifunction--s-square-xi'
    elif xifunction == int_s_xi:
        nowstr += '.xifunction--s-xi'
    elif xifunction == int_s3_xi:
        nowstr += '.xifunction--s3-xi'
    elif xifunction == int_1overs_xi:
        nowstr += '.xifunction--1overs-xi'
    elif xifunction == intxi_rat:
        nowstr += '.xifunction--intxi_rat'
    elif xifunction == int_s_square_xi_rat:
        nowstr += '.xifunction--int_s_square_xi_rat'
    elif xifunction == int_weightedmu_s2xi:
        nowstr += '.int_weightedmu_s2xi'
    elif xifunction == int_weighteds_s2xi:
        nowstr += '.int_weighteds_s2xi'
    elif xifunction == int_weightedmu_xi:
        nowstr += '.int_weightedmu_xi'
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range'
    elif xifunction == int_s_squre_xi_normed__by_full_mu_range:
	nowstr += '.int_s_squre_xi_normed__by_full_mu_range'
    elif xifunction == int_s_xi_normed__by_full_mu_range:
	nowstr += '.int_s_xi_normed__by_full_mu_range'
    return nowstr

#def filename_2pcfcovmat(nowsmin=5, nowsmax=100, nownummubin=2, nowminmucut=0.1, rsdstr = 'norsd', nowdir='./2pcfdata/covmatchisqs-Oct31/', 
def filename_2pcfcovmat(nowsmin=5, nowsmax=100, nownummubin=2, nowminmucut=0.1, rsdstr = 'norsd', nowdir='./2pcfdata/covmatchisqs-Dec15/', 
	extraprefix='', diffchisqmethod = 'use_lastrbin_as_ref', xifunction = 0):
    nowstr = nowdir+'/'+extraprefix+'covmats--s'+str(nowsmin)+'to'+str(nowsmax)+'--'+str(nownummubin)\
           +'mubins--muge%.3f'%nowminmucut+'--'+rsdstr+'.covmat'
    if diffchisqmethod == 'use_weightedavg_as_ref':
        nowstr += '.use_weightedavg_as_ref'
    if xifunction == int_s_square_xi:
        nowstr += '.xifunction--s-square-xi'
    elif xifunction == int_s_xi:
        nowstr += '.xifunction--s-xi'
    elif xifunction == int_s3_xi:
        nowstr += '.xifunction--s3-xi'
    elif xifunction == int_1overs_xi:
        nowstr += '.xifunction--1overs-xi'
    elif xifunction == intxi_rat:
        nowstr += '.xifunction--intxi_rat'
    elif xifunction == int_s_square_xi_rat:
        nowstr += '.xifunction--int_s_square_xi_rat'
    elif xifunction == int_weightedmu_s2xi:
        nowstr += '.int_weightedmu_s2xi'
    elif xifunction == int_weighteds_s2xi:
        nowstr += '.int_weighteds_s2xi'
    elif xifunction == int_weightedmu_xi:
        nowstr += '.int_weightedmu_xi'
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range'
    elif xifunction == int_s_squre_xi_normed__by_full_mu_range:
	nowstr += '.int_s_squre_xi_normed__by_full_mu_range'
    elif xifunction == int_s_xi_normed__by_full_mu_range:
	nowstr += '.int_s_xi_normed__by_full_mu_range'
    elif xifunction != intxi and xifunction != intxi_FractionalS:
	print 'Error!!!! Unkonwn quantity!'
	nowstr += '.unkown-quantity'
    return nowstr

def fmtprint_chisqrlt(chisqrlts,addstr='\t\t'):
    [om,w], chisq_norsd0, chisq_rsd0 = chisqrlts[0]
    for iomw in range(len(chisqrlts)):
        [om,w], chisq_norsd, chisq_rsd = chisqrlts[iomw]
        sig1 = np.sqrt(abs(chisq_norsd-chisq_norsd0)) * np.sign(chisq_norsd-chisq_norsd0)
        sig2 = np.sqrt(abs(chisq_rsd-chisq_rsd0)) * np.sign(chisq_rsd-chisq_rsd0)
        print addstr,'om/w = %.3f'%om,'/ %.3f'%w,':   No RSD/RSD chisq = %8.3f'%chisq_norsd+\
            ' (%6.1f'%sig1+'sigma)  /%8.3f'%chisq_rsd+\
            ' (%6.1f'%sig2+'sigma)'

def chisq_of_intxiatdiffr(curves_at_diffr, ers_at_diffr, mini_mu = 0, maxi_mu = 1000000, chisqmethod = 'minus'):
    numr = len(curves_at_diffr)
    nummu = len(curves_at_diffr[0])
    totchisq = 0.0
    numcurve = len(curves_at_diffr)
    for i_mu in range(mini_mu,min(nummu,maxi_mu)):
        #xiavg = sum([curves_at_diffr[ir][i_mu] for ir in range(numr)]) / float(numr)
        for ir in range(numr):
            if chisqmethod == 'div':
                A = curves_at_diffr[ir][i_mu]
                sigA = ers_at_diffr[ir][i_mu]
                B = curves_at_diffr[numcurve-1][i_mu]
                sigB = ers_at_diffr[numcurve-1][i_mu]
                val = A/B
                er  = val * np.sqrt((sigA/A)**2.0 + (sigB/B)**2.0 )
                totchisq += ( (val-1)/er )**2.0
            elif chisqmethod == 'minus':
                totchisq += ( (curves_at_diffr[ir][i_mu] - curves_at_diffr[numcurve-1][i_mu]) / \
                             np.sqrt(ers_at_diffr[ir][i_mu]**2.0+ers_at_diffr[numcurve-1][i_mu]**2.0)) **2.0            
    return totchisq
    
def chisq_of_intxiatdiffr_cov(binnedxis,printinfo=False, diffchisqmethod='use_lastrbin_as_ref'):
    nummock = len(binnedxis); numr = len(binnedxis[0]); nummu = len(binnedxis[0][0])
    totchisq, totlike = 0, 1
    covmats = []
    if diffchisqmethod == 'use_lastrbin_as_ref':
        ### Values of intxi(this bin) - intxi(last bin)
        diffintxis = [[[binnedxis[imock][ir][i_mu]-binnedxis[imock][numr-1][i_mu] 
                      for imock in range(nummock)] for i_mu in range(nummu)] for ir in range(numr-1)]
        for ir in range(numr-1):
            if printinfo:
                print 'rbin = ', ir
            if printinfo and False:
                print '########################\nDistance bin ', ir
                print 'diff(intxi) at different mocks:'
                for imock in range(nummock):
                    nowstr = ''
                    for i_mu in range(nummu):
                        nowstr += ('%.4f'%(diffintxis[ir][i_mu][imock])+' ')
                    print '\t',  nowstr
            ### Result of different mu. at different mocks; imock is the second index; so estimate covmat between mu;
            Xs = diffintxis[ir]
            covmat = get_covmat(Xs)
            covmats.append(copy.deepcopy(covmat))
            if printinfo:
                print '  covmat: '
                for row1 in range(nummu):
                    nowstr = ''
                    for row2 in range(nummu):
                        nowstr += ('%.6e'%(covmat[row1][row2])+' ')
                        #if row1 != row2:
                            #covmat[row1][row2] = 0
                    print '\t',  nowstr
                    
            ### mean value of diff(intxi); ...
            diffX = [ sum(diffintxis[ir][i_mu])/float(len(diffintxis[ir][i_mu])) for i_mu in range(nummu)]
            nowstr = ''
            for i_mu in range(nummu):
                nowstr += ('%.4f'%(diffX[i_mu])+' ')
            chisq, like = chisq_like_cov_xbar(diffX, covmat)
            if printinfo and False:
                print 'diffX:\n\t', nowstr
                print 'chisq, like = ', chisq, like
            totchisq += chisq; totlike *= like
    elif diffchisqmethod == 'use_weightedavg_as_ref':
	### in this case no need to do differential
        diffintxis = [[[binnedxis[imock][ir][i_mu]
                      for imock in range(nummock)] for i_mu in range(nummu)] for ir in range(numr)]
        for ir in range(numr):
                if printinfo:
                    print 'rbin = ', ir
                if printinfo and False:
                    print '########################\nDistance bin ', ir
                    print 'diff(intxi) at different mocks:'
                    for imock in range(nummock):
                        nowstr = ''
                        for i_mu in range(nummu):
                            nowstr += ('%.4f'%(diffintxis[ir][i_mu][imock])+' ')
                        print '\t',  nowstr
                ### Result of different mu. at different mocks; imock is the second index; so estimate covmat between mu;
                Xs = diffintxis[ir]
                covmat = get_covmat(Xs)
                covmats.append(copy.deepcopy(covmat))
                if printinfo:
                    print '  covmat: '
                    for row1 in range(nummu):
                        nowstr = ''
                        for row2 in range(nummu):
                            nowstr += ('%.6e'%(covmat[row1][row2])+' ')
                            #if row1 != row2:
                                #covmat[row1][row2] = 0
                        print '\t',  nowstr
                        
                ### mean value of diff(intxi); ...
                diffX = [ sum(diffintxis[ir][i_mu])/float(len(diffintxis[ir][i_mu])) for i_mu in range(nummu)]
                nowstr = ''
                for i_mu in range(nummu):
                    nowstr += ('%.4f'%(diffX[i_mu])+' ')
                
                chisq, like = chisq_like_cov_xbar(diffX, covmat)
                if printinfo and False:
                    print 'diffX:\n\t', nowstr
                    print 'chisq, like = ', chisq, like
                totchisq += chisq; totlike *= like
    return totchisq, totlike, covmats

imocks = range(1,9)
import numpy as np


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

execfile(pythonlibPATH+'/Tpcftools_plot.py')
