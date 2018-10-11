#!/home/xiaodongli/software/anaconda/bin/python
# Filename: Tpcftools.py 

import commands
pythonlibPATH = commands.getoutput('echo $pythonlibPATH')
pyfile=pythonlibPATH+'/'+'Tpcftools.py'

### Settings for smu

def smusettings_get_rmax_nbins_mubins(smusettings):
	return int(smusettings['smax']+0.5), int(smusettings['numsbin']+0.5), int(smusettings['nummubin']+0.5)

smusettings_sample = {
            'smin':0.0, 'smax':150.0, 'numsbin':150,
            'mumin':0.0, 'mumax':1.0, 'nummubin':120,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }

smusettings_smax60 = {
            'smin':0.0, 'smax':60.0, 'numsbin':60,
            'mumin':0.0, 'mumax':1.0, 'nummubin':120,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }

smusettings_smax51 = {
            'smin':0.0, 'smax':51.0, 'numsbin':51,
            'mumin':0.0, 'mumax':1.0, 'nummubin':120,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }

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

def smu__loadin(smufile, smusettings, icol='', skiprows=1):
    """		Load in 2pcf file in s, mu space
			if icol=='', load in all columns; else only load in icol column
    """
    data = np.loadtxt(smufile, skiprows=skiprows)
    #print smufile, len(data)

    if icol == '':
	Tpcfrlt = [[data[smu__rowofsmu(i_s,i_mu,smusettings)] for i_mu in range(smusettings['nummubin'])] for i_s in range(smusettings['numsbin'])]
    else:
	Tpcfrlt = [[data[smu__rowofsmu(i_s,i_mu,smusettings)][icol] for i_mu in range(smusettings['nummubin'])] for i_s in range(smusettings['numsbin'])]
    #Tpcfrlt = np.asarray(Tpcfrlt)
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

def xi_natural(dd, dr, rr):
    return dd/rr - 1.0
def xi_davis(dd, dr, rr):
    return dd/dr - 1.0
def xi_LS(dd, dr, rr):
    return (dd-2*dr)/rr + 1.0

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
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range_minimalmuge1:
	name =  'int--xi_normed_muge1'
    elif xifunction == int_s_squre_xi_normed__by_full_mu_range:
	name =  'int--s-squre-xi_normed'
    elif xifunction == int_s_xi_normed__by_full_mu_range:
	name =  'int--s-xi_normed'
    elif xifunction ==  intxi_FractionalS:
	name = 'int--xi'
    elif xifunction ==  int_Numcounts:
	name = 'int--Numcounts'
    else:
	name = 'int--xi'
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

# ismin_float=0.3, ismax_float=105, i_mumin=1, i_mumax=2)
## integral of xi, allowing use of a fraction of s
def int_Numcounts(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax):
    sumdd = 0; sumrr = 0
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
        rr = packed_count(RRlist, nowis, nowis, i_mumin, i_mumax)
        #print 'nowis = ', nowis
        sumdd += dd
        sumrr += rr
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, i_mumin, i_mumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, i_mumin, i_mumax)
    #print 'nowis = ', ismin-1
    sumdd += (dd * smin_tail)
    sumrr += (rr * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, i_mumin, i_mumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, i_mumin, i_mumax)
        #print 'nowis = ', ismax+1
        sumdd += (dd * smax_tail)
        sumrr += (rr * smax_tail)
    return sumdd - sumrr


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

def intxi_FractionalS_normed__by_full_mu_range_minimalmuge1(DDlist, DRlist, RRlist, ismin_float, ismax_float, i_mumin, i_mumax, 
	minimali_mumin=1,maximali_mumax=119):
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
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range_minimalmuge1:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range.minimalmuge1'
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
    elif xifunction == intxi_FractionalS_normed__by_full_mu_range_minimalmuge1:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range.minimalmuge1'
    elif xifunction == int_s_xi_normed__by_full_mu_range:
	nowstr += '.int_s_xi_normed__by_full_mu_range'
    elif xifunction == int_Numcounts:
	nowstr += '.int_Numcounts'
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



### something dealing with   number counts   and   merging of two 2pcf files...

def weighted_sum_of_counts(a, b, weiA, weiB):
    #if random.rand() > 0.999:
    #    print a, b, '\t', a*weiA, b*weiB, '\t', (a*weiA + b*weiB) / (weiA + weiB)
    return (a*weiA + b*weiB) / (weiA + weiB)

def merge_2pcf_files(gal_sumweiA, ran_sumweiA,
                     gal_sumweiB, ran_sumweiB,
                     #catnameA, galfileA, ranfileA, ## Info of the first 2pCF file
                     #catnameB, galfileB, ranfileB, ## Info of the second 2pCF file
                     TpcffileA, ## Name of 1st
                     TpcffileB, ## Name of 2rd
                     TpcffileC  ## Name of 3rd (the merged result)
                     ):
    dataA = np.loadtxt(TpcffileA)
    dataB = np.loadtxt(TpcffileB)
    
    nowf = open(TpcffileC, 'w')
    nowf.write(open(TpcffileA,'r').readline())
    for irow in range(len(dataA)):
        cosmin, cosmax, rmin, rmax, DDnormA, DRnormA, RRnormA, xinaturalA, xidavisA, xiLSA = dataA[irow]
        cosmin, cosmax, rmin, rmax, DDnormB, DRnormB, RRnormB, xinaturalB, xidavisB, xiLSB = dataB[irow]
        DDnormC = weighted_sum_of_counts(DDnormA, DDnormB, gal_sumweiA**2.0, gal_sumweiB**2.0)
        DRnormC = weighted_sum_of_counts(DRnormA, DRnormB, gal_sumweiA*ran_sumweiA, gal_sumweiB*ran_sumweiB)
        RRnormC = weighted_sum_of_counts(RRnormA, RRnormB, ran_sumweiA**2.0, ran_sumweiB**2.0)
        xinautral  = xi_natural(DDnormC, DRnormC, RRnormC) 
        xidavis    =   xi_davis(DDnormC, DRnormC, RRnormC) 
        xiLS       =      xi_LS(DDnormC, DRnormC, RRnormC) 
        if False:
            print 'cosmin, cosmax, rmin, rmax, DDnormA, DRnormA, RRnormA, xinaturalA, xidavisA, xiLSA = ',
            print cosmin, cosmax, rmin, rmax, DDnormA, DRnormA, RRnormA, xinaturalA, xidavisA, xiLSA
            print 'cosmin, cosmax, rmin, rmax, DDnormB, DRnormB, RRnormB, xinaturalB, xidavisB, xiLSB = ',
            print cosmin, cosmax, rmin, rmax, DDnormB, DRnormB, RRnormB, xinaturalB, xidavisB, xiLSB
            print 'DDnormA, DDnormB = ', DDnormA, DDnormB
            print 'DRnormA, DRnormB = ', DRnormA, DRnormB
            print 'RRnormA, RRnormB = ', RRnormA, RRnormB
            print 'gal_sumweiA, gal_sumweiB, ran_sumweiA, ran_sumweiB = ', gal_sumweiA,gal_sumweiB,ran_sumweiA,ran_sumweiB
            print cosmin, cosmax, rmin, rmax, '\t', xiLSA, xiLSB, '\t', xiLS
            print 'DDnormA, DDnormB, DDnormC, DRnormA, DRnormB, DRnormC, RRnormA, RRnormB, RRnormC = ',
            print DDnormA, DDnormB, DDnormC, DRnormA, DRnormB, DRnormC, RRnormA, RRnormB, RRnormC
            return 
        nowf.write(array_to_str_fmtted([cosmin, cosmax, rmin, rmax, DDnormC, DRnormC, RRnormC, 
                xinautral, xidavis, xiLS], fmt='%14.7e')+'\n')
    return TpcffileC


imocks = range(1,9)
import numpy as np
### Useful Shell Commands
execfile(pythonlibPATH+'/Tpcftools_cmds.py')

### Subroutins for smu intxi

execfile(pythonlibPATH+'/Tpcftools_smuintxi.py')

### Subroutins for plot

execfile(pythonlibPATH+'/Tpcftools_plot.py')

#py_Plot is the command to plot
def smu_smids(s1,s2):
 X = range(s1, s2+1)
 X = [((X[row]**3.0+X[row+1]**3.0)/2.0)**(1.0/3.0) for row in range(len(X)-1)]
 return X


def smu_xis(filename, 
	given_ipt_data = None,
	outputtofile = False, outputfilename=None,
	smax=51, nummubin=120,
	imumin = 0,
	ismax = None,
	is_sig_pi=False,
	sfact = 1, 
	make_plot = True, savefig=False, figname=None):
  smusettings = {
            'smin':0.0, 'smax':smax, 'numsbin':smax,
#            'mumin':0.0, 'mumax':1.0, 'nummubin':50,
            'mumin':0.0, 'mumax':1.0, 'nummubin':nummubin,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }
  smu__initsmusettings(smusettings)

  imumax = nummubin;
  if outputfilename == None:  outputfilename = filename+'.xi_s'
  if True:
		if given_ipt_data == None:
			data = smu__loadin(filename, smusettings)
		else:
			data = given_ipt_data
		DDlist, DRlist, RRlist = Xsfrom2ddata(data, [4,5,6])
		if ismax == None: ismax = smax
		if make_plot: fig, ax1 = figax()

                ### packed count of xi as a function of s
                now_s = 0;
                sasx = []; packedxiasy = [];
		sedges = range(smax+1)

                for now_s in range(ismax):
                    nowx=((sedges[now_s]**3.0+sedges[now_s+1]**3.0)/2.0)**(1.0/3.0);   #nowx = now_s+0.5
		    sasx.append(nowx)
		    #print now_s, nowx
		    if is_sig_pi:
			    	Y = []
			    	for nowxig in range(ismax-1):
				 nowpi = int(np.sqrt(nowx**2.0 - nowxig**2.0)) 
				 nowpi = max(nowpi,0)
				 Y.append(packedxi(DDlist,DRlist,RRlist,nowxig,nowxig+1,nowpi,nowpi+1))
			    	nowy = nowx**sfact*sum(Y) / (len(Y)+0.0)
				#nowy = np.log(np.abs(sum(Y)/len(Y)+0.0)) / np.log(10.0)
			    	packedxiasy.append(nowy)
		    else:
                                packedxiasy.append(nowx**sfact*packedxi(DDlist, DRlist, RRlist, now_s, now_s, imumin, imumax))
                    now_s += 1;
                #print sasx
  	        #print packedxiasy
		if outputtofile: np.savetxt(filename+'.xis', packedxiasy)
                if make_plot:
                        ax1.plot(sasx, packedxiasy, marker='o', markersize=1)
                        ax1.set_xlabel('$s\ [\\rm Mpc/h]$', fontsize=25)
                        ax1.set_ylabel('$s^{'+str(sfact)+'}\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                        ax1.set_xlim(0,ismax)
			ax1.set_title(separate_path_file(filename)[1])
		        if savefig: 
				if figname == None: figname = filename+'.png'
				fig.savefig(figname, format = 'png')
		        #plt.show()

  return packedxiasy

def smu_xis_ilus(sfact = 1,
                 s1=6,s2=50,
                 binbins=2,
                 ipt_omwlist = [(0.31,-1.0), (0.06,-1.0), (0.66,-1.0), (0.31,-2.0), (0.31,-0.3)],
                 catname2s = ['HR3', 'J08', 'data'],
		 catnames = ['DR12v4-LOWZ', 'DR12v4-CMASS'],
                 bin_iplot12=False,
                 savefig=True,
                 imumin =0,
                 fs=17):
    X = range(s1, s2+1)
    X = [((X[row]**3.0+X[row+1]**3.0)/2.0)**(1.0/3.0) for row in range(len(X)-1)]
    smids = [x for x in X]
    Xnew = [avgarray(X[row*binbins:(row+1)*(binbins)]) for row in range(len(X)/binbins)]; X=Xnew
    if True:
     print
     print '################## sfact = ', sfact, '##########################'
     print
     for bin_iplot12 in [False]:
      for catname2 in catname2s:
       if catname2 == 'data': 
             omwlist = ipt_omwlist
       else: 
             omwlist = [(0.26,-1.0)]
       for omw in omwlist:
        om, w = omw
        fig = plt.figure(figsize=(16,10));
        ax1 = fig.add_subplot(221);ax2 = fig.add_subplot(222);
        ax3 = fig.add_subplot(223);ax4 = fig.add_subplot(224);
        iplot = 0;
        for catname in catnames:
            for ibin in range(3):
                if catname2 == 'data':
                    filename = Tpcfrltfilename(cosmoconvertedfilename(
                        binsplittedfilename(datafile(catname), ibin+1), om, w), 51, 51, 120)
                    Y = smu_xis(filename, sfact=sfact, imumin=imumin, make_plot=False); Y = Y[s1:s2]
                else:
                    Ys = []
                    for imock in range(4):
                        filename = Tpcfrltfilename(binsplittedfilename(mockfile(catname, catname2, imock, 'RSD'), ibin+1), 150, 150, 120)
                        Y = smu_xis(filename, sfact=sfact, imumin=imumin, make_plot=False); Y = Y[s1:s2]
                        Ys.append([y for y in Y])
                    Y = avgarray_2d(Ys)
                binnumstr = str(ibin+1)
                Ynew = [avgarray(Y[row*binbins:(row+1)*(binbins)]) for row in range(len(Y)/binbins)]; Y=Ynew
                nowc = PLOT_COLOR_ARRAY[iplot]
                Ynorm = normto1(Y);
                if bin_iplot12:
                    if iplot == 1:
                        Ysave = [y for y in Y];
                        Ynormsave = [y for y in Ynorm]
                        iplot += 1
                        continue
                    elif iplot == 2:
                        Y = [(Y[row]+Ysave[row])*0.5 for row in range(len(Y))]
                        Ynorm = [(Ynorm[row]+Ynormsave[row])*0.5 for row in range(len(Ynorm))]
                        binnumstr = '2,3'
                legname = catname+'-bin'+binnumstr
                ax1.plot(X, Y, label = legname, lw=2, c=nowc)
                ax2.plot(X, Ynorm, label = legname, lw=2, c=nowc)
                if iplot == 0: 
                    refY = [y for y in Y]; refYnorm = [y for y in Ynorm]
                else:
                    difY = [Y[row]-refY[row] for row in range(len(Y))]
                    difYnorm = [Ynorm[row]-refYnorm[row] for row in range(len(Ynorm))]
                    ax3.plot(X, difY, label = legname + 'wrt bin1', lw=2, c=nowc)
                    ax4.plot(X, difYnorm, label = legname + 'wrt bin1', lw=2, c=nowc)
                iplot +=1;
        for ax in [ax1,ax2,ax3,ax4]:
                ax.grid();ax.set_xlabel('$s$ [Mpc/h]', fontsize=fs); 
                ax.legend(loc='best', frameon=False)
        ax1.set_ylabel('$s^{'+str(sfact)+'} \\xi(s)$', fontsize=fs);
        ax2.set_ylabel('$\\hat s^{'+str(sfact)+'} \\xi(s)$', fontsize=fs);
        ax3.set_ylabel('$\\delta s^{'+str(sfact)+'} \\xi(s)$', fontsize=fs);
        ax4.set_ylabel('$\\delta \\hat s^{'+str(sfact)+'} \\xi(s)$', fontsize=fs);
        fig.suptitle(catname2+'; Cosmology = '+str(om)+' '+str(w)+'; s = '+str(s1)+'-'+str(s2)+' Mpc/h; imumin='+str(imumin), fontsize=13)
        ax4.set_ylim(-0.15,0.15)
        if savefig:
            fig.savefig('3.xis_'+catname2+'.'+omwstr(om,w)+'.sfact'+str(sfact)+'.imumin'+str(imumin)+'.binbins'+str(binbins)+'.png', format='png')
        plt.show()

def smu_xis_loadmockrlt(imumin, 
		print_nummock=False,
		just_check_correctness=False):
		''' Load in the xis computed from all mocks; returning the dictonary.
			format: catname, catname2, RSDstr, ibin, imock'''
	 	### 1. Load in xis from mock 
		xis_mock = {}
		#xisfile = smu_xis_xisdir + 'MockResult.imumin%03i'%imumin+'.xis'
		xisfile = smu_xis_xisdir + 'MockResult.MDincluded.imumin%03i'%imumin+'.xis'
		for catname in catnamelist:
		  xis_mock[catname] = {}
		  for catname2 in catname2list:
			xis_mock[catname][catname2] = {}
			for RSDstr in ['noRSD', 'RSD']:
				xis_mock[catname][catname2][RSDstr] = {}
		   		for ibin in range(totbin):
				    xis_mock[catname][catname2][RSDstr][ibin]= []
		nowf = open(xisfile, 'r')
		print '  Open and read: ', xisfile
		while True:
			nowstr = nowf.readline()
			if nowstr == '': break
			nowstr = nowstr.split()
			if nowstr[0][0] == '#': continue
			catname, catname2, RSDstr,  imock, ibin = \
				nowstr[0] , nowstr[1], nowstr[2], int(nowstr[3]), int(nowstr[4])
			xis_mock[catname][catname2][RSDstr][ibin].append([float(nowstr[row]) for row in range(5,len(nowstr))])
			### Checking the correctness of xi loaded from file...
			if just_check_correctness and random.random() > 0.99: 
				X1 = [float(nowstr[row]) for row in range(5,len(nowstr))]
                                galfile = binsplittedfilename(mockfile(catname, catname2, imock, RSDstr, ), ibin+1,totbin)
                             	Tpcfrltfile = Tpcfrltfilename(galfile, smax_mock, smax_mock, nummubins)
				X2 = smu_xis(Tpcfrltfile, smax=smax_mock, sfact=0, make_plot=False, imumin=imumin);
				print len(X1), len(X2)
				sys.exit()
		nowf.close()
		if print_nummock:
		 print '  Load in Done. # of xis loaded: '
		 for catname in catnamelist:
		  print '#################'
                  for catname2 in catname2list:
			print '######'
                        for RSDstr in ['noRSD', 'RSD']:
				print 
                                for ibin in range(totbin):
					print '    ', catname, catname2, RSDsr, ibin,':', \
						 len(xis_mock[catname][catname2][RSDstr][ibin])
		return xis_mock


def smu_xis_loaddatarlt(xis_data_file, totbin =3):
                print '  Open and read:', xis_data_file
		xis_data = {}
                for catname in catnamelist:
                        xis_data[catname] =  {}
                        for ibin in range(totbin):
                                xis_data[catname][ibin] = {}
                for nowstr in open(xis_data_file, 'r').readlines():
                        nowstr = nowstr.split()
                        if nowstr[0][0] == '#': continue
                        catname, ibin, nowomwstr = nowstr[0:3]; ibin = int(ibin)-1  # DR12v4-LOWZ   1  om0.0600_w-2.5000;
                        xis =[float(nowstr[row]) for row in range(3,len(nowstr))]
                        xis_data[catname][ibin][nowomwstr] = [x for x in xis]
		return xis_data

#def xismu(data, sbin=150, mubin=120, s1=6, s2=40, DDicol=4, DRicol=5, RRicol=6, DDnorm=1, DRnorm=1, RRnorm=1, nmubin = None):
#    xis = [0 for row in range(sbin)]
#    ximu = [0 for row in range(mubin)]
#    xismu = [[0 for rowsub in range(mubin)] for rows in range(sbin)]
    
#    DDs, DRs, RRs = [0 for row in range(sbin)], [0 for row in range(sbin)], [0 for row in range(sbin)]
#    DDmu, DRmu, RRmu = [0 for row in range(mubin)], [0 for row in range(mubin)], [0 for row in range(mubin)]
#    for rows in range(sbin):
#        for rowmu in range(mubin):
#            irow = rows*mubin+rowmu
#            DD, DR, RR = data[irow][DDicol]/DDnorm, data[irow][DRicol]/DRnorm, data[irow][RRicol]/RRnorm
#            DDs[rows] += DD; DRs[rows] += DR; RRs[rows] += RR
#            DDmu[rowmu] += DD; DRmu[rowmu] += DR; RRmu[rowmu] += RR
#            xismu[rows][rowmu] = xi_LS(DD, DR, RR)
#    xis = [xi_LS(DDs[row], DRs[row], RRs[row]) for row in range(sbin)]
#    ximu = [xi_LS(DDmu[row], DRmu[row], RRmu[row]) for row in range(mubin)]
#    intximu = [ sum([xismu[rows][rowmu] for rows in range(s1,s2)]) for rowmu in range(mubin)]
#    return xis, ximu, intximu

def xismu(data, sbin=150, mubin=120, s1=6, s2=40, DDicol=4, DRicol=5, RRicol=6, DDnorm=1, DRnorm=1, RRnorm=1, 
          nmubin = None, imu1=None, rs=None, mus=None):
    ''' xi(s), xi(mu), AND \int xi(s,mu) ds as a function of mu
	Example:
		data = np.loadtxt(nowfile)
		nmubin = 15
		imu1 = 1
		rs, mus = [], []
		xis, ximu, intximu = xismu(data, 150, 120, 6, 40, 
	        DDicol=4,DRicol=5,RRicol=6, DDnorm=1, DRnorm=1, RRnorm=1, nmubin=nmubin, imu1=imu1, rs=rs, mus=mus,);
		
		fig = plt.figure(figsize=(14,6));
		ax1, ax2, ax3 = fig.add_subplot(131), fig.add_subplot(132), fig.add_subplot(133), 
		ax1.plot(rs, xis)
		ax2.plot(mus, ximu)
		ax3.plot(mus, intximu)'''
    if nmubin == None:
        xis = [0 for row in range(sbin)]
        ximu = [0 for row in range(mubin)]
        xismu = [[0 for rowsub in range(mubin)] for rows in range(sbin)]
        
        DDs, DRs, RRs = [0 for row in range(sbin)], [0 for row in range(sbin)], [0 for row in range(sbin)]
        DDmu, DRmu, RRmu = [0 for row in range(mubin)], [0 for row in range(mubin)], [0 for row in range(mubin)]
        for rows in range(sbin):
            for rowmu in range(mubin):
                irow = rows*mubin+rowmu
                DD, DR, RR = data[irow][DDicol]/DDnorm, data[irow][DRicol]/DRnorm, data[irow][RRicol]/RRnorm
                DDs[rows] += DD; DRs[rows] += DR; RRs[rows] += RR
                DDmu[rowmu] += DD; DRmu[rowmu] += DR; RRmu[rowmu] += RR
                xismu[rows][rowmu] = xi_LS(DD, DR, RR)
        xis = [xi_LS(DDs[row], DRs[row], RRs[row]) for row in range(sbin)]
        ximu = [xi_LS(DDmu[row], DRmu[row], RRmu[row]) for row in range(mubin)]
        intximu = [ sum([xismu[rows][rowmu] for rows in range(s1,s2)]) for rowmu in range(mubin)]
    else:
        xis = [0 for row in range(sbin)]
        ximu = [0 for row in range(nmubin)]
        xismu = [[0 for rowsub in range(nmubin)] for rows in range(sbin)]

        DDs, DRs, RRs = [0 for row in range(sbin)], [0 for row in range(sbin)], [0 for row in range(sbin)]
        DDmu, DRmu, RRmu = [0 for row in range(nmubin)], [0 for row in range(nmubin)], [0 for row in range(nmubin)]
        for rows in range(sbin):
            rowmu=0; irow1 = rows*mubin+rowmu
            rowmu=mubin-1; irow2 = rows*mubin+rowmu
            #print irow1, irow2
            DD, DR, RR = array_fracbin([data[row][DDicol]/DDnorm for row in range(irow1,irow2+1)], nmubin, imu1),\
                         array_fracbin([data[row][DRicol]/DRnorm for row in range(irow1,irow2+1)], nmubin, imu1),\
                         array_fracbin([data[row][RRicol]/RRnorm for row in range(irow1,irow2+1)], nmubin, imu1),
            #print len(DD)
            #return
            
            DDs[rows] += sum(DD); DRs[rows] += sum(DR); RRs[rows] += sum(RR);
            DDmu  = XplusY(DDmu, DD);
            DRmu  = XplusY(DRmu, DR);
            RRmu  = XplusY(RRmu, RR);
            #DRmu = DR; RRmu[rowmu] += RR
            for rowmu in range(nmubin):
                xismu[rows][rowmu] = xi_LS(DD[rowmu], DR[rowmu], RR[rowmu])
        xis = [xi_LS(DDs[row], DRs[row], RRs[row]) for row in range(sbin)]
        ximu = [xi_LS(DDmu[row], DRmu[row], RRmu[row]) for row in range(nmubin)]
        intximu = [ sum([xismu[rows][rowmu] for rows in range(s1,s2)]) for rowmu in range(nmubin)]
    nowrs, nowmus = get_mid_array1d(range(sbin+1)), get_mid_array1d(np.linspace(0,1,nmubin+1))
    if rs != None:
	try:
          if len(rs) == len(nowrs):
	     for row in range(len(nowrs)):
		rs[row] = nowrs[row]
          elif len(rs) == 0:
	     for row in range(len(nowrs)):
		rs.append(nowrs[row])
        except: pass
    if mus != None:
	try:
          if len(mus) == len(nowmus):
	     for row in range(len(nowmus)):
		mus[row] = nowmus[row]
          elif len(mus) == 0:
	     for row in range(len(nowmus)):
		mus.append(nowmus[row])
        except: pass
    return xis, ximu, intximu




execfile(pythonlibPATH+'/Tpcftools_smuxis.py')
execfile(pythonlibPATH+'/Tpcftools_smuximu.py')
execfile(pythonlibPATH+'/Tpcftools_smuximu_calcximus.py')
execfile(pythonlibPATH+'/Tpcftools_ChisqContour.py')

