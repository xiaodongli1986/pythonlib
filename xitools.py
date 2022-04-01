
import numpy as np
import os, sys
import matplotlib.pyplot as plt


def Hz(omegam, w, h, z):
    return 100*h*np.sqrt(omegam*(1.0+z)**3.0 + (1.0-omegam)*(1.0+z)**(3.0*(1+w)))

def omegamz(om0, z):
    '''compute ratio of matter component at redshift z; assume lambda cdm
	H = Hz(om0, -1, 0.7, z);
	om = om0 * (1.0+z)**3
	return om/H**2.0*10000*0.7**2'''
    H = Hz(om0, -1, 0.7, z);
    om = om0 * (1.0+z)**3
    return om/H**2.0*10000*0.7**2

def Hz_Mpc_to_h(omegam, w, h, z):
    return 100*np.sqrt(omegam*(1.0+z)**3.0 + (1.0-omegam)*(1.0+z)**(3.0*(1+w)))

def comov_r(omegam, w, h, z):
    CONST_C = 3e8
    #x, y = scipy.integrate.quad(lambda x: 1.0/Hz(omegam, w, h, x), 0, z)
    x = np.linspace(0,z, 1000); y = [1.0/Hz(omegam, w, h, xx) for xx in x]
    return CONST_C * np.trapz(y, x, ) * h

def DA(omegam, w, h, z):
	return comov_r(omegam, w, h, z) / (1.0 + z)

def packarray2d(A, rat):
    As = [A[:,row::rat] for row in range(rat)]
    #print As
    nowlen =  min([len(X[0]) for X in As])
    B = As[0][:,:nowlen]
    for row in range(1,len(As)):
        B = B + As[row][:,:nowlen]
    return B

def meannorm(Y):
    if np.mean(Y) > 0:
        return Y / np.mean(Y)
    else:
        return Y / abs(np.mean(Y)) + 2

def smu__CosmoConvert(s,mu,DA1,DA2,H1,H2):
    ''' s1: angular direction; s2: LOS direction '''
    s2 = s*mu;
    s1 = np.sqrt(s*s - s2*s2)
    alpha1 = DA2 / DA1
    alpha2 = H1 / H2
    s_prime  =  np.sqrt((alpha1*s1)**2 + (alpha2*s2)**2)
    mu_prime =  alpha2*s2 / s_prime 
    return s_prime, mu_prime

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

def mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd, DAstd, DAnew, Hstd, Hnew,
        deltas1=0.2, deltamu1=1.0/600.0, deltas2=1.0, deltamu2=1.0/120.0, smin_mapping=1,smax_mapping=51,
        compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[], method='simple_bin'):
        #compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[9]):
    ''' method can be simple_bin or divided_pixel'''
    nummu1 = int(1.0/deltamu1 + 0.5)
    sbound1, mubound1 = smu__CosmoConvert(smax_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smax_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = max(sbound1, sbound2); nums1 = int(sbound / deltas1 + 0.5)
    smutab1 = np.array(smutabstd)

    sbound1, mubound1 = smu__CosmoConvert(smin_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smin_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = min(sbound1, sbound2); mins1 = int(sbound / deltas1)

    nums2, nummu2 = int(smax_mapping / deltas2 + 0.5), int(1.0/deltamu2 + 0.5)
    numrow3=len(smutabstd[0][0])
    smutab2 = np.array([[[0 for row3 in range(numrow3+1)] for row2 in range(nummu2)] for row1 in range(nums2)])

    #numtests = 0###DEBUG
    #print mins1, nums1, nummu1, nums2, nummu2
    if method == 'divided_pixel':
        smutab1_centers = np.array([[0 for row2 in range(nummu1)] for row1 in range(nums1)])
    for is1 in range(mins1, nums1):
        for iangle1 in range(nummu1):
            scenter, anglecenter = (is1+0.5)*deltas1, (iangle1+0.5)*deltamu1
            mucenter = 1.0 - anglecenter
            scenter2, mucenter2 = smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,)
            anglecenter2 = 1.0 - mucenter2
            is2 = int(scenter2  / deltas2  )
            iangle2 = int(anglecenter2 / deltamu2)
            if method == 'simple_bin':
                if is2 < nums2 and iangle2 < nummu2:
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutab1[is1][iangle1][row3]
                    if save_counts_row != None:
                        smutab2[is2][iangle2][save_counts_row] += 1
            elif method == 'divided_pixel':
                smutab1_centers[is1][iangle1] = [scenter2, anglecenter2, is2, iangle2]

    if method == 'divided_pixel':
        for is1 in range(mins1, nums1):
            for iangle1 in range(nummu1):
                scenter2, anglecenter2, is2, iangle2 = smutab1_centers[is1][iangle1]
                if not (is2 < nums2 and iangle2 < nummu2):
                    continue

                ### firstly, check boundary:
                sboundflag, muboundflag = False, False
                for is1_nearby in [max(is1-1,mins1), min(is1+1,nums1-1)]:
                    is2_b, iangle2_b = smutab1_centers[is1_nearby][iangle1][2], smutab1_centers[is1_nearby][iangle1][3]
                    if is2_b != is2:
                        sboundflag=True; is1_bound = is1_nearby; is2_bound = is2_b;
                        scenter2_bound = smutab1_centers[is1_nearby][iangle1][0]### Xiao-DOng : I think it should be 0??
                        #scenter2_bound = (is2_bound +0.5)*deltas2
                for iangle1_nearby in [max(iangle1-1,0), min(iangle1+1,nummu1-1)]:
                    is2_b, iangle2_b = smutab1_centers[is1][iangle1_nearby][2], smutab1_centers[is1][iangle1_nearby][3]
                    if iangle2_b != iangle2:
                        muboundflag=True; iangle1_bound = iangle1_nearby; iangle2_bound = iangle2_b
                        anglecenter2_bound = smutab1_centers[is1][iangle1_nearby][1]### Xiao-dong: I think it should be 1?
                        #anglecenter2_bound = (iangle2_bound+0.5)*deltamu2

                ### Then, treat them case by case...
                ## s, mu are all not near the boundary of tab 2
                if ((not sboundflag)and(not muboundflag)):
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutab1[is1][iangle1][row3]
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += 1

                ## s is near the boundary of tab2
                if sboundflag and (not muboundflag):
                    s = (is2 + is2_bound +1) * 0.5 * deltas2
                    if False:
                        rat = (s-scenter2) / (scenter2_bound-scenter2)
                    else:
                        #print 'scenter2, scenter2_bound = ', scenter2, scenter2_bound
                        #numtests += 1
                        #if numtests > 10:
                        #    sys.exit(1)
                        scenter3=(scenter2+scenter2_bound) * 0.5
                        ds = scenter3-scenter2
                        d1 = s-scenter2+ds
                        rat = d1 / (2*ds)
                        rat = min(rat, 1)

                    rat_bound = 1-rat
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutab1[is1][iangle1][row3]*rat
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
                    if is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2][row3] += smutab1[is1][iangle1][row3]*rat_bound
                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat_bound

                ## mu is near the boundary of tab2
                if muboundflag and (not sboundflag):
                    angle = (iangle2 + iangle2_bound +1) * 0.5 * deltamu2
                    if False:
                        rat = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
                    else:
                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
                        dangle = anglecenter3-anglecenter2
                        d1 = angle-anglecenter2+dangle
                        rat = d1 / (2*dangle)
                        rat = min(rat, 1)

                    rat_bound = 1-rat
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutab1[is1][iangle1][row3]*rat
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
                    if iangle2_bound < nummu2:
                      for row3 in compute_rows:
                        smutab2[is2][iangle2_bound][row3] += smutab1[is1][iangle1][row3]*rat_bound
                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat_bound

                ## both s, mu are near the boundy...
                if muboundflag and sboundflag:
                    s = (is2 + is2_bound+1) * 0.5 * deltas2
                    angle = (iangle2 + iangle2_bound+1) * 0.5 * deltamu2
                    if False:
                        rats = (s-scenter2) / (scenter2_bound-scenter2)
                        ratangle = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
                    else:
                        scenter3=(scenter2+scenter2_bound) * 0.5
                        ds = scenter3-scenter2
                        d1 = s-scenter2+ds
                        rats = d1 / (2*ds)
                        rats = min(rats, 1)

                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
                        dangle = anglecenter3-anglecenter2
                        d1 = angle-anglecenter2+dangle
                        ratangle = d1 / (2*dangle)
                        ratangle = min(ratangle, 1)
                    # original pixel
                    rat1 = rats*ratangle
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutab1[is1][iangle1][row3]*rat1
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat1
                    # diff s
                    rat2 = (1-rats)*ratangle
                    if is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2][row3] += smutab1[is1][iangle1][row3]*rat2
                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat2
                    # diff angle
                    rat3 = rats*(1-ratangle)
                    if iangle2_bound < nummu2:
                      for row3 in compute_rows:
                        smutab2[is2][iangle2_bound][row3] += smutab1[is1][iangle1][row3]*rat3
                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat3
                    # diff s and diff angle
                    rat4 = (1-rats)*(1-ratangle)
                    if iangle2_bound < nummu2 and is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2_bound][row3] += smutab1[is1][iangle1][row3]*rat4
                      if save_counts_row != None: smutab2[is2_bound][iangle2_bound][save_counts_row] += rat4

    if div_counts_rows != [] and save_counts_row != None:
        for is2 in range(nums2):
            for iangle2 in range(nummu2):
                for row3 in div_counts_rows:
                    smutab2[is2][iangle2][row3] /= smutab2[is2][iangle2][save_counts_row]
    return smutab2
def mapping_smudata_to_another_cosmology_simple(smutabstd, DAstd, DAnew, Hstd, Hnew, deltamu=1.0/120.0, simple_replacement=False, max_mubin=120,
                                         smin_mapping=1,smax_mapping=51):
    ''' mapping: so far only mapping the last coordinate!! '''
    import copy
    smutab2 = copy.deepcopy(smutabstd)

    for row1 in range(smin_mapping,smax_mapping):
      for row2 in range(len(smutabstd[row1])):
            smax, anglemax = row1, row2
            mumin = 1 - anglemax * deltamu
            if True:
            #smid = smax - 0.5
            #mumid = mumin + deltamu * 0.5
                smax2, mumin2 = smu__CosmoConvert(smax,mumin,DAnew,DAstd,Hnew,Hstd)
            else: # coordinate transformation using "middle" rather than "boundary" values of s and mu: not much improvement...
                smid = smax - 0.5
                mumid = mumin + deltamu * 0.5
                smid2, mumid2 = smu__CosmoConvert(smid,mumid,DAnew,DAstd,Hnew,Hstd)
                smax2, mumin2 = smid2+0.5, mumid2-deltamu*0.5
                #print smid, mumid, smid2, mumid2, smax2, mumin2
            anglemax2 = (1 - mumin2)/deltamu
            if not simple_replacement:
             x1, y1 = int(smax2), int(anglemax2)
             #if x1 <=5: x1 +=1
             if y1 == max_mubin-1: y1=y1-1
             x2, y2 = x1 + 1, y1 + 1
             #x3, y3 = x2 + 1, y2 + 1
             #print smax, anglemax, smax2, anglemax2
             smutab2[row1][row2] = \
                    LinearInterpolation_2d(x1,y1,x2,y2,
                                   smutabstd[x1][y1],smutabstd[x1][y2],smutabstd[x2][y1],\
                                   smutabstd[x2][y2],smax2, anglemax2)
            else:
             x1, y1 = int(smax2+0.5), int(anglemax2+0.5)
             smutab2[row1][row2] = smutabstd[x1][y1]
    return smutab2

class xismu:
    def __init__(self, filename=None, smax=150, sbin=150, mubin=120,
                DDnorm=None, DRnorm=None, RRnorm=None, DD=None, DR=None, RR=None, ):
        self.smax, self.sbin, self.mubin =smax, sbin, mubin
        if filename != None:
            self.DDnorm, self.DRnorm, self.RRnorm = [float(xx) for xx in open(filename, 'r').readline().split()[1:4]]
            self.data = np.loadtxt(filename, )
            self.DD, self.DR, self.RR = [self.data[:sbin*mubin,row].reshape(sbin,mubin) for row in [3,4,6]]
        else:
            self.DDnorm, self.DRnorm, self.RRnorm, self.DD, self.DR, self.RR = DDnorm, DRnorm, RRnorm, DD, DR, RR

        self.DD /= self.DDnorm; self.DR /= self.DRnorm; self.RR /= self.RRnorm

    def s_xis(self, ):
        X = np.linspace(0,self.smax,self.sbin+1)
        X = (X[1:] + X[:len(X)-1]) / 2.
        self.xis =  (self.DD[:,:].sum(1) - 2*self.DR[:,:].sum(1) + self.RR[:,:].sum(1)) / self.RR[:,:].sum(1)
        return X, self.xis

    def intximu(self, mumax=0.97, intxismin=6, intxismax=40, mubin_pack_rat=1):
        imumax = int(self.mubin*mumax)
        DD, DR, RR = self.DD[:,:imumax], self.DR[:,:imumax], self.RR[:,:imumax]
        DD, DR, RR = packarray2d(DD, mubin_pack_rat), packarray2d(DR, mubin_pack_rat), packarray2d(RR, mubin_pack_rat)
        xi = np.divide(DD-2*DR+RR,RR)

        ds = (self.smax-0)/float(self.sbin)
        is1, is2 = int((intxismin-0)/ds), int((intxismax-0)/ds)
        intximu = (xi[is1:is2+1,:].sum(0)*ds)[::-1]

        X = np.linspace(1-mumax,1,len(intximu)+1)
        X = (X[1:] + X[:len(X)-1]) / 2.
        return X, intximu
    def cosmo_conv(self, omstd= 0.3071, wstd = -1, omwrong = 0.5, wwrong = -1, redshift=None, smax_mapping=100, ):
            h = 0.6777
            redshift = max(redshift , 0.0001)
            DAstd, DAnew = DA(omstd, wstd, h, redshift), DA(omwrong, wwrong, h, redshift)
            Hzstd, Hznew = Hz(omstd, wstd, h, redshift), Hz(omwrong, wwrong, h, redshift)
            data = self.data.reshape(self.sbin, self.mubin, -1); data = data[:,::-1,:]
            data = mapping_smudata_to_another_cosmology_simple(data, DAstd, DAnew, Hzstd, Hznew, deltamu=1.0/self.mubin,
                                                       simple_replacement=False,smax_mapping=smax_mapping)
            data = data[:,::-1,:]
            sbin2=int(self.sbin*smax_mapping/self.smax)
            DDs ,DRs ,RRs = [data[:sbin2,:,row] for row in [3,4,6]]; #DDs = DDs.reshape(-1); DRs = DRs.reshape(-1); RRs = RRs.reshape(-1)
            return xismu(  filename=None, smax=smax_mapping, sbin=sbin2, mubin=self.mubin,
                    DDnorm=1, DRnorm=1, RRnorm=1, DD=DDs, DR=DRs, RR=RRs)
    def cosmo_conv_DenseToSparse(self, omstd=0.3071, wstd=-1, omwrong=0.5, wwrong=-1, redshift=None, 
        sbin2=150, mubin2=120, smin_mapping=1, smax_mapping=60, method='divided_pixel'):
            h = 0.6777
            redshift = max(redshift , 0.0001)
            DAstd, DAnew = DA(omstd, wstd, h, redshift), DA(omwrong, wwrong, h, redshift)
            Hzstd, Hznew = Hz(omstd, wstd, h, redshift), Hz(omwrong, wwrong, h, redshift)
            data = self.data.reshape(self.sbin, self.mubin, -1); data = data[:,::-1,:]
        
            data = mapping_smudata_to_another_cosmology_DenseToSparse(data, DAstd, DAnew, Hzstd, Hznew,
                deltas1=1.0/self.sbin, deltamu1=1.0/self.mubin, deltas2=1.0/sbin2, deltamu2=1.0/mubin2, 
                smin_mapping=smin_mapping,smax_mapping=smax_mapping,
                compute_rows=[3,4,6], save_counts_row=0, method=method)
            DDs,DRs ,RRs = [data[:,:,row] for row in [3,4,6]]; #DDs = DDs.reshape(-1); DRs = DRs.reshape(-1); RRs = RRs.reshape(-1)
            return xismu(  filename=None, smax=smax_mapping, sbin=sbin2, mubin=mubin2,
                    DDnorm=1, DRnorm=1, RRnorm=1, DD=DDs, DR=DRs, RR=RRs)

