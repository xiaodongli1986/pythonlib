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
	

def intxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax):
    sumxi = 0
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    if ismin > ismax or imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        sumxi += calc_xils(dd, dr, rr)
        #print s
    return sumxi

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def intximu_over_intxi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        #print 'nowis = ', nowis
        sumxi += calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    #print 'nowis = ', ismin-1
    sumxi += (calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        #print 'nowis = ', ismax+1
        sumxi += (calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        #print 'nowis = ', nowis
        sumxi += calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    #print 'nowis = ', ismin-1
    sumxi += (calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        #print 'nowis = ', ismax+1
        sumxi += (calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 



#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        sumxi += now_s_square*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s_square*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s_square*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi    

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_s3_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        sumxi += now_s_square**1.5 *calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s_square**1.5 *calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s_square**1.5 *calc_xils(dd, dr, rr) * smax_tail)
    return sumxi    

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_weightedmu_s2xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square*nowmu
		sumwei += nowxi*now_s_square
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square*nowmu * smin_tail
		sumwei += nowxi*now_s_square  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square*nowmu * smax_tail
	        sumwei += nowxi*now_s_square  * smax_tail
	        
    return sumxi / sumwei


#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_weightedmu_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*nowmu
		sumwei += nowxi
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*nowmu * smin_tail
		sumwei += nowxi  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*nowmu * smax_tail
	        sumwei += nowxi  * smax_tail
	        
    return sumxi / sumwei

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_weighteds_s2xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0; sumwei = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square**1.5
		sumwei += nowxi*now_s_square
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
    for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square**1.5 * smin_tail
		sumwei += nowxi*now_s_square  * smin_tail
    if not (ismax+1 >= numsbin):
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s_square=(s0*s0+s1*s1)/2.0;
        #print 'nowis = ', ismax+1
        for nowimu in range(imumin,imumax+1):
        	nowmu = mumin + (nowimu+0.5)*deltamu
		nowxi = calc_xils(DDlist[nowis][nowimu], DRlist[nowis][nowimu], RRlist[nowis][nowimu])
	        sumxi += nowxi*now_s_square**1.5 * smax_tail
	        sumwei += nowxi*now_s_square  * smax_tail
	        
    return sumxi / sumwei

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        sumxi += now_s*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (now_s*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (now_s*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi  

#intxi_FractionalS(DDlist1s[0][0][0], DRlist1s[0][0][0], RRlist1s[0][0][0], 
# ismin_float=0.3, ismax_float=105, imumin=1, imumax=2)
## integral of xi, allowing use of a fraction of s
def int_1overs_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax):
    sumxi = 0
    ## boundary as integer ...
    ismin = int(ismin_float+1.0); ismax = int(ismax_float-1.0);
    smin_tail = ismin - ismin_float; smax_tail = ismax_float - 1.0 - ismax  ## two tails (at minimal and maximal s)
    imumin = max(imumin,0); imumax = min(imumax, nummubin-1)
    ismin = max(ismin,0); ismax = min(ismax, numsbin-1)
    #print 'ismin, ismax, ismin_float, ismax_float = ',ismin,ismax,ismin_float,ismax_float,\
    #    '; smin_tail, smax_tail=',smin_tail,smax_tail
    if imumin > imumax:
        return 0
    for nowis in range(ismin, ismax+1):
        dd = packed_count(DDlist, nowis, nowis, imumin, imumax)
        dr = packed_count(DRlist, nowis, nowis, imumin, imumax)
        rr = packed_count(RRlist, nowis, nowis, imumin, imumax)
        s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        sumxi += 1.0/now_s*calc_xils(dd, dr, rr)
    ## considering small s tail
    dd = packed_count(DDlist, ismin-1, ismin-1, imumin, imumax)
    dr = packed_count(DRlist, ismin-1, ismin-1, imumin, imumax)
    rr = packed_count(RRlist, ismin-1, ismin-1, imumin, imumax)
    nowis=ismin-1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
    #print 'nowis = ', ismin-1
    sumxi += (1.0/now_s*calc_xils(dd, dr, rr) * smin_tail)
    ## considering large s tail
    if not (ismax+1 >= numsbin):
        dd = packed_count(DDlist, ismax+1, ismax+1, imumin, imumax)
        dr = packed_count(DRlist, ismax+1, ismax+1, imumin, imumax)
        rr = packed_count(RRlist, ismax+1, ismax+1, imumin, imumax)
        nowis=ismax+1; s0=smin+nowis*deltas;s1=s0+deltas;now_s=(s0+s1)/2.0;
        #print 'nowis = ', ismax+1
        sumxi += (1.0/now_s*calc_xils(dd, dr, rr) * smax_tail)
    return sumxi 

def intxi_rat(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax): 
	ismid = (ismin_float+ismax_float)/2.0
	return intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismid,  imumin, imumax) / intxi_FractionalS(DDlist, DRlist, RRlist, ismid,ismax_float, imumin, imumax)
	
def int_s_square_xi_rat(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax): 
	ismid = (ismin_float+ismax_float)/2.0
	return int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismid,  imumin, imumax) / int_s_square_xi(DDlist, DRlist, RRlist, ismid,ismax_float, imumin, imumax)

def intxi_FractionalS_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax, 
	minimalimumin=0,maximalimumax=119):
	return intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax) /\
		intxi_FractionalS(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimalimumin, maximalimumax)

def int_s_squre_xi_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax, 
	minimalimumin=0,maximalimumax=119):
	return int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax) /\
		int_s_square_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimalimumin, maximalimumax)

def int_s_xi_normed__by_full_mu_range(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax, 
	minimalimumin=0,maximalimumax=119):
	return int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, imumin, imumax) /\
		int_s_xi(DDlist, DRlist, RRlist, ismin_float, ismax_float, minimalimumin, maximalimumax)


def settingstr_2pcf(nowsmin=5, nowsmax=100, nownummubin=10, nowminmucut=0.1, RSDcor=True, flexiblesmin=False, flexiblesmax=False, 	
	diffchisqmethod = 'use_lastrbin_as_ref', calcquan = 0):
    nowstr = 's'+str(nowsmin)+'to'+str(nowsmax)+'--'+str(nownummubin)+'mubins--muge%.3f'%nowminmucut+'--RSDCor'+str(RSDcor)
    if flexiblesmin:
        nowstr += '--flexiblesmin'
    if flexiblesmax:
        nowstr += '--flexiblesmax'
    if diffchisqmethod == 'use_weightedavg_as_ref':
        nowstr += '.use_weightedavg_as_ref'
    if calcquan == int_s_square_xi:
        nowstr += '.calcquan--s-square-xi'
    elif calcquan == int_s_xi:
        nowstr += '.calcquan--s-xi'
    elif calcquan == int_s3_xi:
        nowstr += '.calcquan--s3-xi'
    elif calcquan == int_1overs_xi:
        nowstr += '.calcquan--1overs-xi'
    elif calcquan == intxi_rat:
        nowstr += '.calcquan--intxi_rat'
    elif calcquan == int_s_square_xi_rat:
        nowstr += '.calcquan--int_s_square_xi_rat'
    elif calcquan == int_weightedmu_s2xi:
        nowstr += '.int_weightedmu_s2xi'
    elif calcquan == int_weighteds_s2xi:
        nowstr += '.int_weighteds_s2xi'
    elif calcquan == int_weightedmu_xi:
        nowstr += '.int_weightedmu_xi'
    elif calcquan == intxi_FractionalS_normed__by_full_mu_range:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range'
    elif calcquan == int_s_squre_xi_normed__by_full_mu_range:
	nowstr += '.int_s_squre_xi_normed__by_full_mu_range'
    elif calcquan == int_s_xi_normed__by_full_mu_range:
	nowstr += '.int_s_xi_normed__by_full_mu_range'
    return nowstr

#def filename_2pcfcovmat(nowsmin=5, nowsmax=100, nownummubin=2, nowminmucut=0.1, rsdstr = 'norsd', nowdir='./2pcfdata/covmatchisqs-Oct31/', 
def filename_2pcfcovmat(nowsmin=5, nowsmax=100, nownummubin=2, nowminmucut=0.1, rsdstr = 'norsd', nowdir='./2pcfdata/covmatchisqs/', 
	extraprefix='', diffchisqmethod = 'use_lastrbin_as_ref', calcquan = 0):
    nowstr = nowdir+'/'+extraprefix+'covmats--s'+str(nowsmin)+'to'+str(nowsmax)+'--'+str(nownummubin)\
           +'mubins--muge%.3f'%nowminmucut+'--'+rsdstr+'.covmat'
    if diffchisqmethod == 'use_weightedavg_as_ref':
        nowstr += '.use_weightedavg_as_ref'
    if calcquan == int_s_square_xi:
        nowstr += '.calcquan--s-square-xi'
    elif calcquan == int_s_xi:
        nowstr += '.calcquan--s-xi'
    elif calcquan == int_s3_xi:
        nowstr += '.calcquan--s3-xi'
    elif calcquan == int_1overs_xi:
        nowstr += '.calcquan--1overs-xi'
    elif calcquan == intxi_rat:
        nowstr += '.calcquan--intxi_rat'
    elif calcquan == int_s_square_xi_rat:
        nowstr += '.calcquan--int_s_square_xi_rat'
    elif calcquan == int_weightedmu_s2xi:
        nowstr += '.int_weightedmu_s2xi'
    elif calcquan == int_weighteds_s2xi:
        nowstr += '.int_weighteds_s2xi'
    elif calcquan == int_weightedmu_xi:
        nowstr += '.int_weightedmu_xi'
    elif calcquan == intxi_FractionalS_normed__by_full_mu_range:
	nowstr += '.intxi_FractionalS_normed__by_full_mu_range'
    elif calcquan == int_s_squre_xi_normed__by_full_mu_range:
	nowstr += '.int_s_squre_xi_normed__by_full_mu_range'
    elif calcquan == int_s_xi_normed__by_full_mu_range:
	nowstr += '.int_s_xi_normed__by_full_mu_range'
    elif calcquan != intxi and calcquan != intxi_FractionalS:
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

def chisq_of_intxiatdiffr(curves_at_diffr, ers_at_diffr, minimu = 0, maximu = 1000000, chisqmethod = 'minus'):
    numr = len(curves_at_diffr)
    nummu = len(curves_at_diffr[0])
    totchisq = 0.0
    numcurve = len(curves_at_diffr)
    for imu in range(minimu,min(nummu,maximu)):
        #xiavg = sum([curves_at_diffr[ir][imu] for ir in range(numr)]) / float(numr)
        for ir in range(numr):
            if chisqmethod == 'div':
                A = curves_at_diffr[ir][imu]
                sigA = ers_at_diffr[ir][imu]
                B = curves_at_diffr[numcurve-1][imu]
                sigB = ers_at_diffr[numcurve-1][imu]
                val = A/B
                er  = val * np.sqrt((sigA/A)**2.0 + (sigB/B)**2.0 )
                totchisq += ( (val-1)/er )**2.0
            elif chisqmethod == 'minus':
                totchisq += ( (curves_at_diffr[ir][imu] - curves_at_diffr[numcurve-1][imu]) / \
                             np.sqrt(ers_at_diffr[ir][imu]**2.0+ers_at_diffr[numcurve-1][imu]**2.0)) **2.0            
    return totchisq
    
def chisq_of_intxiatdiffr_cov(binnedxis,printinfo=False, diffchisqmethod='use_lastrbin_as_ref'):
    nummock = len(binnedxis); numr = len(binnedxis[0]); nummu = len(binnedxis[0][0])
    totchisq, totlike = 0, 1
    covmats = []
    if diffchisqmethod == 'use_lastrbin_as_ref':
        ### Values of intxi(this bin) - intxi(last bin)
        diffintxis = [[[binnedxis[imock][ir][imu]-binnedxis[imock][numr-1][imu] 
                      for imock in range(nummock)] for imu in range(nummu)] for ir in range(numr-1)]
        for ir in range(numr-1):
            if printinfo:
                print 'rbin = ', ir
            if printinfo and False:
                print '########################\nDistance bin ', ir
                print 'diff(intxi) at different mocks:'
                for imock in range(nummock):
                    nowstr = ''
                    for imu in range(nummu):
                        nowstr += ('%.4f'%(diffintxis[ir][imu][imock])+' ')
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
            diffX = [ sum(diffintxis[ir][imu])/float(len(diffintxis[ir][imu])) for imu in range(nummu)]
            nowstr = ''
            for imu in range(nummu):
                nowstr += ('%.4f'%(diffX[imu])+' ')
            chisq, like = chisq_like_cov_xbar(diffX, covmat)
            if printinfo and False:
                print 'diffX:\n\t', nowstr
                print 'chisq, like = ', chisq, like
            totchisq += chisq; totlike *= like
    elif diffchisqmethod == 'use_weightedavg_as_ref':
	### in this case no need to do differential
        diffintxis = [[[binnedxis[imock][ir][imu]
                      for imock in range(nummock)] for imu in range(nummu)] for ir in range(numr)]
        for ir in range(numr):
                if printinfo:
                    print 'rbin = ', ir
                if printinfo and False:
                    print '########################\nDistance bin ', ir
                    print 'diff(intxi) at different mocks:'
                    for imock in range(nummock):
                        nowstr = ''
                        for imu in range(nummu):
                            nowstr += ('%.4f'%(diffintxis[ir][imu][imock])+' ')
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
                diffX = [ sum(diffintxis[ir][imu])/float(len(diffintxis[ir][imu])) for imu in range(nummu)]
                nowstr = ''
                for imu in range(nummu):
                    nowstr += ('%.4f'%(diffX[imu])+' ')
                
                chisq, like = chisq_like_cov_xbar(diffX, covmat)
                if printinfo and False:
                    print 'diffX:\n\t', nowstr
                    print 'chisq, like = ', chisq, like
                totchisq += chisq; totlike *= like
    return totchisq, totlike, covmats

def make_2pcf_plottings(nowsmins=[5],nowsmaxs=[100],iomws=[0,1,2,3,4],irbinlist=[0,1,2,3,4], noplot=False, plotRSDeff=False,
                   nowxrange=[],nowyrange=[],plotnoRSD=True,plotRSD=True,printsetttings=False, calcquan=intxi,
		   chisqmethod = 'minus',
                   iimocklist=range(len(imocks)),muedges=linspace(0,1,11), minmucut = -1, RSDcor=True, RSDcortozero=False,
		   consistsrange=False,
                   maxmucut=1000000,useourerrorbar=False,savfig=False,shifttozero=False,divto1=False, usecovmat=False,
                   reprocessed_shifttozero_covmat=False,reprocessed_divto1_covmat=False,
                   set_last_mubindiff_tozero = False, separate_RSDCor=False,
                   diffchisqmethod='use_lastrbin_as_ref', weight_for_avgdiff=[],
                   ConsiderVarVar=False, CorrectedCovmat=False, usegivencovmats=[], chisqprintinfo=True,
		   notitle = False):

    ### Range of s
    # Minimal/Maximal indice of s
    if calcquan == intxi or calcquan == packed_count:
	    nowismins = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmins]
	    nowismaxs = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmaxs]
    else:
	    nowismins = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmins]
	    nowismaxs = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmaxs]
#    else:
#    	print 'ERROR (make_2pcf_plottings)! Unknown calcquan!'
#    	return
    	
    ### Binning of mu
    imuedges = [[min(max(0,  int((muedges[row]-mumin)/deltamu)  ),nummubin), 
                 min(max(0,  int((muedges[row+1]-mumin)/deltamu)-1  ),nummubin)]
                 for row in range(len(muedges)-1)]
    # Minimal/Maximal cut of mu
    if minmucut > 0:
        minimu = get_imu(minmucut)
        imuedges2 = copy.deepcopy(imuedges); imuedges=[]
        for imuedge in imuedges2:
            if imuedge[1] < minimu:
                continue
            else:
                imuedges.append([max(imuedge[0],minimu),imuedge[1]])
        if printsetttings:
            print 'Re-define imuedge due to minimal mucut ', minmucut, '.\n\torig:', imuedges2, '\n\tnow:', imuedges
    if maxmucut < 1:
        maximu = get_imu(maxmucut)
        imuedges2 = copy.deepcopy(imuedges); imuedges=[]
        for imuedge in imuedges2:
            if imuedge[0] > maximu:
                continue
            else:
                imuedges.append([imuedge[0],min(maximu,imuedge[1])])
        if printsetttings:
            print 'Re-define imuedge due to maximal mucut ', maxmucut, '.\n\torig:', imuedges2, '\n\tnow:', imuedges
    # indice of mu for each binning; middle value at each bin
    imuranges = [range(imuedge[0],imuedge[1]+1) for imuedge in imuedges]
    muedgemids = [(get_mid_smu(0,imuedge[0])[1]+get_mid_smu(0,imuedge[1])[1])*0.5 for imuedge in imuedges]
    
    ### print settings of the program
    if printsetttings:
        print 'Settings:'
        print '\tlist of irbins: \n\t\t',irbinlist,'\n\t\t', [rbinstrs[irbin] for irbin in irbinlist]
        print '\tlist of mocks: \n\t\t',iimocklist,'\n\t\t', [imocks[iimock] for iimock in iimocklist]
        print '\tlist of muedges:  \n\t\t',muedges, '\n\t\t', imuedges
        print '\t\tBinning of imu:  \n\t\t\t', imuranges
        print '\t\t# of imu in each bin:  \n\t\t\t', [len(imuranges[row]) for row in range(len(imuranges))]
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
                    if calcquan == intxi or calcquan == packed_count:
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
                        [ calcquan(DD1, DR1, RR1, nowismin, nowismax, muedge[0], muedge[1]) for muedge in imuedges]
                    binnedxis2[i_of_mock][i_of_bin] = \
                        [ calcquan(DD2, DR2, RR2, nowismin, nowismax, muedge[0], muedge[1]) for muedge in imuedges]

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

### This plot 2pcf for boss data; a simple plot (no RSD correction, no covmat load in, no chisq computation, ...)
def make_2pcf_plottings_bossdata_simpleplot(nowsmins=[5],nowsmaxs=[150],iomws=[0],irbinlist=[0], noplot=False, plotRSDeff=False,
                   nowxrange=[],nowyrange=[],plotnoRSD=True,plotRSD=True,printsetttings=False, calcquan=intxi,
		   chisqmethod = 'minus',
                   iimocklist=range(len(imocks)),muedges=linspace(0,1,11), minmucut = -1, RSDcor=True, RSDcortozero=False,
		   consistsrange=False,
                   maxmucut=1000000,useourerrorbar=False,savfig=False,shifttozero=False,divto1=False, usecovmat=False,
                   reprocessed_shifttozero_covmat=False,reprocessed_divto1_covmat=False,
                   set_last_mubindiff_tozero = False, separate_RSDCor=False,
                   diffchisqmethod='use_lastrbin_as_ref', weight_for_avgdiff=[],
                   ConsiderVarVar=False, CorrectedCovmat=False, usegivencovmats=[], chisqprintinfo=True,
		   notitle = False, nomustd=False):

    ### Range of s
    # Minimal/Maximal indice of s
    #if calcquan == intxi or calcquan == packed_count:
#	    nowismins = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmins]#
#	    nowismaxs = [min(max(0,int((nows-smin)/deltas)),numsbin) for nows in nowsmaxs]
 #   else:
    nowismins = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmins]
    nowismaxs = [min(max(0,(nows-smin)/deltas),numsbin) for nows in nowsmaxs]
#    else:
#    	print 'ERROR (make_2pcf_plottings)! Unknown calcquan!'
#    	return
    	
    ### Binning of mu
    imuedges = [[min(max(0,  int((muedges[row]-mumin)/deltamu)  ),nummubin), 
                 min(max(0,  int((muedges[row+1]-mumin)/deltamu)-1  ),nummubin)]
                 for row in range(len(muedges)-1)]
    # Minimal/Maximal cut of mu
    if minmucut > 0:
        minimu = get_imu(minmucut)
        imuedges2 = copy.deepcopy(imuedges); imuedges=[]
        for imuedge in imuedges2:
            if imuedge[1] < minimu:
                continue
            else:
                imuedges.append([max(imuedge[0],minimu),imuedge[1]])
        if printsetttings:
            print 'Re-define imuedge due to minimal mucut ', minmucut, '.\n\torig:', imuedges2, '\n\tnow:', imuedges
    if maxmucut < 1:
        maximu = get_imu(maxmucut)
        imuedges2 = copy.deepcopy(imuedges); imuedges=[]
        for imuedge in imuedges2:
            if imuedge[0] > maximu:
                continue
            else:
                imuedges.append([imuedge[0],min(maximu,imuedge[1])])
        if printsetttings:
            print 'Re-define imuedge due to maximal mucut ', maxmucut, '.\n\torig:', imuedges2, '\n\tnow:', imuedges
    # indice of mu for each binning; middle value at each bin
    imuranges = [range(imuedge[0],imuedge[1]+1) for imuedge in imuedges]
    muedgemids = [(get_mid_smu(0,imuedge[0])[1]+get_mid_smu(0,imuedge[1])[1])*0.5 for imuedge in imuedges]
    
    ### print settings of the program
    if printsetttings:
        print 'Settings:'
        print '\tlist of irbins: \n\t\t',irbinlist,'\n\t\t', [rbinstrs[irbin] for irbin in irbinlist]
        print '\tlist of mocks: \n\t\t',iimocklist,'\n\t\t', [imocks[iimock] for iimock in iimocklist]
        print '\tlist of muedges:  \n\t\t',muedges, '\n\t\t', imuedges
        print '\t\tBinning of imu:  \n\t\t\t', imuranges
        print '\t\t# of imu in each bin:  \n\t\t\t', [len(imuranges[row]) for row in range(len(imuranges))]
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
            binnedxis = [[0 for irbin in irbinlist] for iimock in iimocklist]
            intxiavgs = [[0 for irbin in irbinlist] for iimock in iimocklist]
            
            ### for each cosmology, loop of different rbin (redshift)
            i_of_bin = -1; 
            xi_diffrs, xistd_diffrs = [],[]
            
            for irbin in irbinlist:
                i_of_bin +=1;
                zetalist,zeta_mu_arrays,mucompactarrays = [],[],[]
                ### loop of indice of mock
                i_of_mock = -1;
                if consistsrange:
                    nowsmin = nowsmins[i_of_smin]*rescalfac[iomw][irbin]
                    nowsmax = nowsmaxs[i_of_smin]*rescalfac[iomw][irbin]
                    if calcquan == intxi or calcquan == packed_count:
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
                        nowfile = get_2pcffile(imocks[iimock], '', rbinstrs[irbin], om, w)
                        zetalist = loadin_smu(nowfile);
                        DD, DR, RR = loadin_counts(nowfile)
                    else:
                        #zetalist1 = zetalist1s[iomw][irbin][iimock]
                        #zetalist2 = zetalist2s[iomw][irbin][iimock]
                        DD, DR, RR = DDlist1s[iomw][irbin][iimock], DRlist1s[iomw][irbin][iimock], RRlist1s[iomw][irbin][iimock]
                        
                    ### calculate integral xi,  save result to binnedxis
                    binnedxis[i_of_mock][i_of_bin] = \
                        [ calcquan(DD, DR, RR, nowismin, nowismax, muedge[0], muedge[1]) for muedge in imuedges]

            ### Normalization & Plotting & Chisq calculation 
            # Figure
            if not noplot:
                fig=plt.figure(figsize=(18,6));ax=fig.add_subplot(111); 
                ax.set_xlabel('$\\mu$',fontsize=18); ax.set_ylabel('$\\sum \\xi$',fontsize=18)
            # binned xi and its STD at different r
            mucompactarray_diffrs, mucompactarraystd_diffrs = [], []
            for i_of_bin in range(len(irbinlist)):
                mucompactarrays = []
                for i_of_mock in range(len(iimocklist)):
                    if divto1:
                        intxiavgs[i_of_mock][i_of_bin], binnedxis[i_of_mock][i_of_bin] = \
                        	normto1(binnedxis[i_of_mock][i_of_bin],returnavg=True)
                        
                    if shifttozero:
                        intxiavgs[i_of_mock][i_of_bin], binnedxis[i_of_mock][i_of_bin] = \
                        	shiftozero(binnedxis[i_of_mock][i_of_bin],returnavg=True)
                        
                    intxi = binnedxis[i_of_mock][i_of_bin]
                    mucompactarrays.append(copy.deepcopy(intxi));
		if not nomustd and len(iimocklist)>1: 
	                mucompactarray, mucompactarraystd = get_avgstd_array(mucompactarrays,ConsiderVarVar=ConsiderVarVar)
		else:
			mucompactarray = mucompactarrays[0]
			mucompactarraystd = [0 for row in range(len(mucompactarray))]
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
                mucompactarray_diffrs.append(copy.deepcopy(mucompactarray)); 
                mucompactarraystd_diffrs.append(copy.deepcopy(mucompactarraystd));
                
                if noplot:
                    continue
                    
                ### label of the noRSD/RSD curve
                labelstr = '${\\rm om/w=%.2f'%om+'/%.2f'%w+';\ \ '+rbinstrs[irbin]+\
                    '}$ BOSS Data: $\\bar\\sum_\\xi=%.2f'%(avgarray(mucompactarray))+'$'

                ### plot the curve
                if True: #BOSSLi: For BOSS data there is no "noRSD" or "RSD"; so we just keep one set of plot
                    if useourerrorbar:
                        ourerrorbar(ax, muedgemids, mucompactarray, mucompactarraystd, lw=lw, c=lclist[mod(ilc,len(lclist))],
                                ls=lslist[mod(ils,len(lslist))], capsize=10, label=labelstr, 
                                marker=markerlist[mod(iplot,len(markerlist))], ms=10, not_dof_div=False, polyorder=1)
			plt.show()
                    else:
                        ax.errorbar(muedgemids,mucompactarray,mucompactarraystd,label=labelstr,
                            lw=lw,ls=lslist[mod(ils,len(lslist))],capsize=10,c=lclist[mod(ilc,len(lclist))],
                            marker=markerlist[mod(iplot,len(markerlist))],markersize=10)
                iplot+=1;ilc+=1

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
            continue ### So far we just plot; no RSDCor, no covmat, no chisq calculation, ...
	
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
		    print 'hahaha'
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
#    if not usecovmat:
#        return chisqrlts, mucompactarray1_diffrs,mucompactarray2_diffrs
#    else:
#        return chisqrlts, covmats_norsd, covmats_rsd, mucompactarray1_diffrs,mucompactarray2_diffrs
