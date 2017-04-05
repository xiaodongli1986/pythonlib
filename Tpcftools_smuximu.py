

def smu__CosmoConvert(s,mu,DA1,DA2,H1,H2):
    ''' s1: angular direction; s2: LOS direction '''
    s2 = s*mu;
    s1 = np.sqrt(s*s - s2*s2)
    alpha1 = DA2 / DA1
    alpha2 = H1 / H2
    s_prime  =  np.sqrt((alpha1*s1)**2 + (alpha2*s2)**2)
    mu_prime =  alpha2*s2 / s_prime 
    return s_prime, mu_prime

#DA1, H1 = DA(0.26, -1, 0.7, 0.5), Hz(0.26, -1, 0.7, 0.5); DA2, H2 = DA(0.21, -1, 0.7, 0.5), Hz(0.21, -1, 0.7, 0.5); 
#print smu__CosmoConvert(20, 0.5, DA1,DA2,H1,H2)

def mapping_smudata_to_another_cosmology(smutabstd, DAstd, DAnew, Hstd, Hnew, deltamu=1.0/120.0, simple_replacement=False, max_mubin=120,
                                         smin_mapping=1,smax_mapping=51):
    ''' mapping: so far only mapping the last coordinate!! '''
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
             for row3 in range(len(smutabstd[row1][row2])):            
	     #if True:
	     #	row3 = 9
                smutab2[row1][row2][row3] = \
                    LinearInterpolation_2d(x1,y1,x2,y2, 
                                   smutabstd[x1][y1][row3],smutabstd[x1][y2][row3],smutabstd[x2][y1][row3],\
                                   smutabstd[x2][y2][row3],smax2, anglemax2)
	    else:
             x1, y1 = int(smax2+0.5), int(anglemax2+0.5)
             for row3 in range(len(smutabstd[row1][row2])):            
              smutab2[row1][row2][row3] = smutabstd[x1][y1][row3]
    return smutab2

def mapping_smudata_to_another_cosmology_DenseToSparse_old(smutabstd, DAstd, DAnew, Hstd, Hnew, 
        deltas1=0.2, deltamu1=1.0/600.0, deltas2=1.0, deltamu2=1.0/120.0, smin_mapping=1,smax_mapping=51, 
        compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[], method='simple_bin'):
        #compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[9]):
    ''' method can be simple_bin or divided_pixel'''
    nummu1 = int(1.0/deltamu1 + 0.5)
    sbound1, mubound1 = smu__CosmoConvert(smax_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smax_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = max(sbound1, sbound2); nums1 = int(sbound / deltas1 + 0.5)
    
    sbound1, mubound1 = smu__CosmoConvert(smin_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smin_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = min(sbound1, sbound2); mins1 = int(sbound / deltas1)
    
    nums2, nummu2 = int(smax_mapping / deltas2 + 0.5), int(1.0/deltamu2 + 0.5)
    numrow3=len(smutabstd[0][0])
    smutab2 = [[[0 for row3 in range(numrow3+1)] for row2 in range(nummu2)] for row1 in range(nums2)]
    
    #print mins1, nums1, nummu1, nums2, nummu2
    if method == 'divided_pixel':
        smutabstd_centers = [[0 for row2 in range(nummu1)] for row1 in range(nums1)]
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
                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
                    if save_counts_row != None:
                        smutab2[is2][iangle2][save_counts_row] += 1
            elif method == 'divided_pixel':
                smutabstd_centers[is1][iangle1] = [scenter2, anglecenter2, is2, iangle2]
                
    if method == 'divided_pixel':
        for is1 in range(mins1, nums1):
            for iangle1 in range(nummu1):
                scenter2, anglecenter2, is2, iangle2 = smutabstd_centers[is1][iangle1]
                if not (is2 < nums2 and iangle2 < nummu2):
                    continue
                    
                ### firstly, check boundary:
                sboundflag, muboundflag = False, False
                for is1_nearby in [max(is1-1,mins1), min(is1+1,nums1-1)]:
                    is2_b, iangle2_b = smutabstd_centers[is1_nearby][iangle1][2], smutabstd_centers[is1_nearby][iangle1][3]
                    if is2_b != is2:
                        sboundflag=True; is1_bound = is1_nearby; is2_bound = is2_b; 
                        scenter2_bound = smutabstd_centers[is1_nearby][iangle1][1]
                for iangle1_nearby in [max(iangle1-1,0), min(iangle1+1,nummu1-1)]:
                    is2_b, iangle2_b = smutabstd_centers[is1][iangle1_nearby][2], smutabstd_centers[is1][iangle1_nearby][3]
                    if iangle2_b != iangle2:
                        muboundflag=True; iangle1_bound = iangle1_nearby; iangle2_bound = iangle2_b
                        anglecenter2_bound = smutabstd_centers[is1][iangle1_nearby][2]
                
                ### Then, treat them case by case...
                ## s, mu are all not near the boundary of tab 2
                if ((not sboundflag)and(not muboundflag)):
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += 1
                            
                ## s is near the boundary of tab2
                if sboundflag and (not muboundflag):
                    s = (is2 + is2_bound) * 0.5 * deltas2
                    if False:
                        rat = (s-scenter2) / (scenter2_bound-scenter2)
                    else:
                        scenter3=(scenter2+scenter2_bound) * 0.5
                        ds = scenter3-scenter2
                        d1 = s-scenter2+ds
                        rat = d1 / (2*ds)
                        rat = min(rat, 1)
                        
                    rat_bound = 1-rat
                    for row3 in compute_rows:
                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
                    if is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat_bound
                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat_bound
                            
                ## mu is near the boundary of tab2
                if muboundflag and (not sboundflag):
                    angle = (iangle2 + iangle2_bound) * 0.5 * deltamu2
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
                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
                    if iangle2_bound < nummu2:
                      for row3 in compute_rows:
                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat_bound
                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat_bound
                
                ## both s, mu are near the boundy...
                if muboundflag and sboundflag:
                    s = (is2 + is2_bound) * 0.5 * deltas2
                    angle = (iangle2 + iangle2_bound) * 0.5 * deltamu2
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
                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat1
                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat1
                    # diff s
                    rat2 = (1-rats)*ratangle
                    if is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat2
                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat2
                    # diff angle
                    rat3 = rats*(1-ratangle)
                    if iangle2_bound < nummu2:
                      for row3 in compute_rows:
                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat3
                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat3
                    # diff s and diff angle
                    rat4 = (1-rats)*(1-ratangle)
                    if iangle2_bound < nummu2 and is2_bound < nums2:
                      for row3 in compute_rows:
                        smutab2[is2_bound][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat4
                      if save_counts_row != None: smutab2[is2_bound][iangle2_bound][save_counts_row] += rat4
                            
    if div_counts_rows != [] and save_counts_row != None:
        for is2 in range(nums2):
            for iangle2 in range(nummu2):
                for row3 in div_counts_rows:
                    smutab2[is2][iangle2][row3] /= smutab2[is2][iangle2][save_counts_row]
    return smutab2

def mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd, DAstd, DAnew, Hstd, Hnew, 
        deltas1=0.2, deltamu1=1.0/600.0, deltas2=1.0, deltamu2=1.0/120.0, smin_mapping=1,smax_mapping=51, 
        compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[], method='simple_bin'):
        #compute_rows=[4,5,6,9], save_counts_row=0, div_counts_rows=[9]):
    ''' method can be simple_bin or divided_pixel'''
    nummu1 = int(1.0/deltamu1 + 0.5)
    sbound1, mubound1 = smu__CosmoConvert(smax_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smax_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = max(sbound1, sbound2); nums1 = int(sbound / deltas1 + 0.5)
    
    sbound1, mubound1 = smu__CosmoConvert(smin_mapping,0,DAnew,DAstd,Hnew,Hstd)
    sbound2, mubound2 = smu__CosmoConvert(smin_mapping,1,DAnew,DAstd,Hnew,Hstd)
    sbound = min(sbound1, sbound2); mins1 = int(sbound / deltas1)
    
    nums2, nummu2 = int(smax_mapping / deltas2 + 0.5), int(1.0/deltamu2 + 0.5)
    numrow3=len(smutabstd[0][0])
    smutab2 = [[[0 for row3 in range(numrow3+1)] for row2 in range(nummu2)] for row1 in range(nums2)]
    
    #numtests = 0###DEBUG
    #print mins1, nums1, nummu1, nums2, nummu2
    if method == 'divided_pixel':
        smutab1_centers = [[0 for row2 in range(nummu1)] for row1 in range(nums1)]
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


def smu_ximu_calcchisqs(
	omlist, wlist,### 1. List of omegam, w
	baseoutputfile, ### 2. Basic name of outputfile
	mumins = [0.01, 0.02, 0.03], ### 3. Range of mumin
	#smu__intxi__settings = smu__intxi__settings_std,
	smusettings_data = smusettings_smax51, 
	smusettings_mock = { 'smin':0.0, 'smax':150.0, 'numsbin':150, 'mumin':0.0, 'mumax':1.0, 'nummubin':120, 'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[] } , 
	s1=6.0, s2=40, NumMubin=21,   Smax = 51, Smax_covmock = 150, Smax_sysmock=150, ### 4. Things about s
	cov_catname2='HR3', cov_RSDstr = 'noRSD', cov_nummock=72, ### 5. settings for covmat and systematic
	xifunction = intxi_FractionalS,
	sys_catname2= 'J08', sys_RSDstr ='RSD', sys_catnamesuffix = '', syscor_imocks = range(4), polyfitdeg=3,
	make_plot=True, showfig=False, savefig=True, 
	totbin = 3, delta_time = 100, nummubins = 120, 	normfun=norm_OneSkip, ### Other miscellous
	outputdir = '/home/xiaodongli/software/', suffix = '',
        use_omw_as_sys = False,
	cosmo_rescaled_s = False, z_of_each_bin = [0.2154098242, 0.3156545036, 0.3869159269, 0.4785988874, 0.5408209467, 0.6186561000],
	consider_tilt_in_chisq = False,
	covmat_nondiag_rescale = False,
        usingmapping_for_nonstd_omw = False, basepar_for_smu_mapping = [0.26,-1],
	rebinxi=False, # In early settings we re-bin the values of DD/DR/RR
	polyfit_dxi = None,
	simple_replacement=False,
	use_DenseMapping=False, DM_nbins=750, DM_nummubin=600, DM_smax=150, method2str='divided_pixel',
		):
	ommin, ommax, wmin, wmax = min(omlist), max(omlist), min(wlist), max(wlist)
	#omwlist = sumlist([[[om,w] for om in omlist] for w in wlist])
	omwlist = [[0.31,-1.5]]+sumlist([[[om,w] for om in omlist] for w in wlist]) #DEBUG
	#sarray = smu_smids(s1, s2)
	#sarray = [ x**sfact for x in sarray ]
	for mumin in mumins:
		smu__intxi__settings = {
            	'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,# In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
         	      }
		smu__intxi__settings_orig = {
            	'xifunction': xifunction, 'smin':s1, 'smax':s2, 'mumin':mumin, 'mumax':1.0, 'nummubin':NumMubin,## In fact this name shall be 'nummuedge'; the actual number of bins is 1 smaller
         	      }
		nowchisqstr = baseoutputfile+'_'+smu__intxi_str(smu__intxi__settings)

		nowchisqstr += ('.CovMock_'+str(cov_nummock)+'_'+cov_catname2+'_'+cov_RSDstr)
		if cosmo_rescaled_s != False:
			CRbaseom, CRbasew = cosmo_rescaled_s
			nowchisqstr += ('.CosmoRescaledS_'+omwstr(CRbaseom,CRbasew))
		if use_omw_as_sys !=False:
			omsys, wsys = use_omw_as_sys
			nowchisqstr += '.force_bestfit_as_'+omwstr(omsys,wsys)
		if sys_catname2 != 'J08' and use_omw_as_sys == False:
			nowchisqstr += ('.SysCor' + sys_catname2)
		if sys_catnamesuffix != '':
			if sys_catnamesuffix == '-N':
				nowchisqstr += '.SysCor_NGC'
			elif sys_catnamesuffix == '-S':
				nowchisqstr += '.SysCor_SGC'
			else: 
				nowchisqstr += ('.SysCor_'+sys_catnamesuffix)
		if consider_tilt_in_chisq != False:
			nowchisqstr += '.Chisq_ConsiderTilt'
		if polyfitdeg != 3:
			nowchisqstr += ('.SysCor_poly'+str(polyfitdeg))
		if covmat_nondiag_rescale != False:
			nowchisqstr += ('.CovNDRescal%.3f'%covmat_nondiag_rescale)
			def covmat_adjust(covmat):
				for rowa in range(len(covmat)):
					for rowb in range(len(covmat)):
						if rowa != rowb:
							covmat[rowa][rowb] *= covmat_nondiag_rescale
				return covmat
		if usingmapping_for_nonstd_omw == True:
			nowchisqstr += '.SMU_Mapping'
			if simple_replacement:
				nowchisqstr += '.simple_replacement'
			if use_DenseMapping: 
				nowchisqstr += ('.DenseMapping.'+method2str)
				DM_smusettings = {'deltamu': 0,'deltas': 0,'mulist': [],'mumax': 1.0,'mumin': 0.0,\
					'nummubin': DM_nummubin,'numsbin': DM_nbins,'slist': [], 'smax': DM_smax,'smin': 0.0};
				smu__initsmusettings(DM_smusettings)
				smu__initsmusettings(smusettings_data)
				#print DM_smusettings, smusettings_data
		if rebinxi == True:
			nowchisqstr += '.rebinxi'
		if polyfit_dxi!= None:
			nowchisqstr += '.polyfit_dxi_deg'+str(polyfit_dxi)

	        ximudir = smu_ximu_covchisqdir
		#nowfile=outputdir+'/'+baseoutputfile+'_'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
		nowfile=outputdir+'/'+nowchisqstr+'.txt'; nowf = open(nowfile,'w')
		print '    Output chisq to file: ', nowfile
		nowf.write('### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor\n')
		### 6. Compute
		#### 6.1 Load in 
		##### 6.1.1 Mock
		#xis_mock = smu_xis_loadmockrlt(mumin, 
		#				      print_nummock=False)  #fmt: catname, catname2, RSDstr, ibin
		##### 6.1.2 Data
		#xis_data_file = smu_xis_ximudir+'/'+baseoutputfile+'.smax'+str(smax)+'.mumin%03i'%mumin+'.xis'
		#xis_data = smu_xis_loaddatarlt(xis_data_file) # fmt: catname, ibin, omwstr
		#print xis_data['DR12v4-CMASS'][1][omwstr(0.31,-1)]
		
		#### 6.2 Compuate chisqs for all ...
		chisqs_nosyscor, chisqs_syscor = {}, {} # fmt: omwstr, ibin
		iomw = 0
		covmats = []
		xisys_list = []
		dxisys_list = []

		t0 = time.clock(); t1 = t0; 
		if usingmapping_for_nonstd_omw:
			omstd, wstd = basepar_for_smu_mapping
		smutabstd_list = []
		for omw in omwlist:
			#print omw ## BOSSLi
			t2 = time.clock();
			if t2 - t1 > 10:
				print t2-t0, 'sec ellapses.   iomw/numomw = ', iomw,'/', len(omwlist), \
					'     rat = %.4f'%(float(iomw)/float(len(omwlist)))
				t1 = t2
			om,w = omw; nowomwstr=omwstr(om,w); chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr] = [], []
			i_redshiftbin = 0
			xihatdata_now = []
			
			for catname in ['DR12v4-LOWZ', 'DR12v4-CMASS', ]:
				for ibin in range(totbin):
					if usingmapping_for_nonstd_omw:
				          if catname == 'DR12v4-LOWZ':
					    zmedian = SixBinRedshift_NS[ibin]
				          elif  catname == 'DR12v4-CMASS':
					    zmedian = SixBinRedshift_NS[ibin+3]
					  #print 'ibin / zmedian =   ', ibin, zmedian
    				          DAstd, Hstd = DA(omstd, wstd, 0.7, zmedian), Hz(omstd, wstd, 0.7, zmedian)
    				          DAnew, Hnew = DA(om,    w,    0.7, zmedian), Hz(om,    w,    0.7, zmedian)
					if iomw == 0 and usingmapping_for_nonstd_omw:
						if not use_DenseMapping:
							smufile = Tpcfrltfilename(cosmoconvertedfilename(\
	 						 	binsplittedfilename(datafile(catname), ibin+1, totbin),omstd,wstd), \
							 	mubins=nummubins,nbins=Smax,rmax=Smax );
							smutabstd_list.append(smu__loadin(smufile,smusettings_data ))
						else:
							smufile = Tpcfrltfilename(cosmoconvertedfilename(\
	 						 	binsplittedfilename(datafile(catname), ibin+1, totbin),omstd,wstd), \
							 	mubins=DM_nummubin,nbins=DM_nbins,rmax=DM_smax );
							smutabstd_list.append(smu__loadin(smufile,DM_smusettings ))
							
						#print 'smax of the smutabstd_list: ', len(smu__loadin(smufile,smusettings_data ))
						#print smusettings_data
					if cosmo_rescaled_s != False:
						nowz = z_of_each_bin[i_redshiftbin]
						VolBase = DA(CRbaseom,CRbasew,0.7,nowz)**2.0 / Hz(CRbaseom,CRbasew,0.7,nowz)
						VolNow  = DA(om,w,0.7,nowz)**2.0 / Hz(om,w,0.7,nowz)
						ReScale = (VolNow/VolBase)**(1.0/3.0)
						CRs1, CRs2 = s1*ReScale, s2*ReScale
						smu__intxi__settings['smin'] = CRs1
						smu__intxi__settings['smax'] = CRs2
						#print om, w, nowz, ': ', CRs1, CRs2
					if i_redshiftbin in [0]:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
						#print smu__intxi_file
						if isfile(smu__intxi_file):
							xidata_base = smu__intxi_quickload(smu__intxi_file)
						else:
						  if not usingmapping_for_nonstd_omw:
							xidata_base = smu__intxi_calcwrite(smufile, smusettings_data,\
								writetofile=False, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
						  elif not use_DenseMapping:
							smudata = mapping_smudata_to_another_cosmology(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								deltamu=1.0/120.0, smin_mapping=1, smax_mapping=51,\
								simple_replacement=simple_replacement )         
							xidata_base = smu__intxi_calcwrite(smufile, smusettings_data, \
							writetofile=False, smudata=smudata, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
						  else:
							smudata = mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								method=method2str,\
								smin_mapping=1, smax_mapping=51,\
                    						deltas1=DM_smusettings['deltas'],  deltamu1=DM_smusettings['deltamu'],\
					                    	deltas2=smusettings_data['deltas'],deltamu2=smusettings_data['deltamu'])
							xidata_base = smu__intxi_calcwrite(smufile, smusettings_data, \
							writetofile=False, smudata=smudata, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
						## B. normalize the amplitude
						if polyfit_dxi== None:
						  xihatdata_base =   normfun(xidata_base)
						else:
						  xihatdata_base =   normfun(polyfitY(range(len(xidata_base)),xidata_base,polyfit_dxi))
							
						
						xihatdata_now.append([x for x in xihatdata_base])
						if iomw == 0:
							xicov_base = []
							for imock in range(cov_nummock):
								smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_covmock,rmax=Smax_covmock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
								if isfile(smu__intxi_file):
									xicov_base.append(smu__intxi_quickload(smu__intxi_file))
								else:
									xicov_base.append(smu__intxi_calcwrite(smufile, smusettings_mock,\
								writetofile=False, smu__intxi__settings=smu__intxi__settings_orig,rebinxi=rebinxi)[2])
								xihatcov_base  = [ normfun(X) for X in xicov_base]
								#if polyfit_dxi== None:
								# xihatcov_base  = [ normfun(X) for X in xicov_base]
								#else:
								# xihatcov_base  = [ normfun(polyfitY(range(len(X)),X,polyfit_dxi)) for X in xicov_base]
							covmats.append([])
							xisys_base = []
							if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_sysmock,rmax=Smax_sysmock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
									#print 'smu__intxi_file', smu__intxi_file # HB
									if isfile(smu__intxi_file):
										xisys_base.append(smu__intxi_quickload(smu__intxi_file))
									else:
										xisys_base.append(smu__intxi_calcwrite(smufile, smusettings_mock,\
										  writetofile=False, smu__intxi__settings=smu__intxi__settings_orig,rebinxi=rebinxi)[2])
							else:
								omsys, wsys = use_omw_as_sys
								smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                                                		if isfile(smu__intxi_file):
		                                                        xisys_base.append(smu__intxi_quickload(smu__intxi_file))
		                                                else:
                		                                        xisys_base.append(smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2])
							xisys_list.append(get_avg_array([normfun(X) for X in xisys_base]))
							dxisys_list.append([])
					else:
						## A. picku up the xis
						smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),om,w), mubins=nummubins,nbins=Smax,rmax=Smax );	smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
						#print smu__intxi_file; sys.exit()
						if isfile(smu__intxi_file):
							#print 'load in smu__intxi_file: ', smu__intxi_file
							xidata = smu__intxi_quickload(smu__intxi_file)
						else:
                                                  if not usingmapping_for_nonstd_omw:
                                                        xidata = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
                                                  elif not use_DenseMapping:
                                                        smudata2 = mapping_smudata_to_another_cosmology(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
                                                                deltamu=1.0/120.0, smin_mapping=1, smax_mapping=51,\
								simple_replacement=simple_replacement )
							#print  'smutabstd_list[i_redshiftbin][10][10][9], smudata[10][10][9] = ', smutabstd_list[i_redshiftbin][10][10][9], smudata[10][10][9]
                                                        xidata = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, \
                                                                smudata=smudata2, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
						  else: 
                                                        smudata2 = mapping_smudata_to_another_cosmology_DenseToSparse(smutabstd_list[i_redshiftbin], \
								DAstd, DAnew, Hstd, Hnew, \
								method=method2str,\
                                                                smin_mapping=1, smax_mapping=51,\
                    						deltas1=DM_smusettings['deltas'],  deltamu1=DM_smusettings['deltamu'],\
					                    	deltas2=smusettings_data['deltas'],deltamu2=smusettings_data['deltamu'])
                                                        xidata = smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, \
                                                                smudata=smudata2, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2]
						#print 'om, w, i_redshiftbin, normto1(xidata) = ', om, w, i_redshiftbin, normto1(xidata) #DEBUG
						## B. normalize the amplitude
						#xihatdata =   normfun(xidata)
						if polyfit_dxi== None:
						 xihatdata =   normfun(xidata)
						else:
						 xihatdata =   normfun(polyfitY(range(len(xidata)),xidata,polyfit_dxi))
						xihatdata_now.append([x for x in xihatdata])
						dxidata = get_diffarray(xihatdata_base, xihatdata)

						#if om == 0.26 and w == -1.0:
						#	print om,w, ':', ibin, xidata, xihatdata
						if iomw == 0:
							xisys = []
	                                                if not use_omw_as_sys:
								for imock in syscor_imocks:
									smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname+sys_catnamesuffix, sys_catname2, imock, sys_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_sysmock,rmax=Smax_sysmock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
									#print 'smu__intxi_file', smu__intxi_file # HB
									if isfile(smu__intxi_file):
										xisys.append(smu__intxi_quickload(smu__intxi_file))
									else:
										xisys.append(smu__intxi_calcwrite(smufile, smusettings_mock, writetofile=False, smu__intxi__settings=smu__intxi__settings_orig,rebinxi=rebinxi)[2])
	                                                else:
        	                                                 omsys, wsys = use_omw_as_sys
                	                                         smufile = Tpcfrltfilename(cosmoconvertedfilename(binsplittedfilename(datafile(catname), ibin+1, totbin),omsys,wsys), mubins=nummubins,nbins=Smax,rmax=Smax ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings)
                        	                                 if isfile(smu__intxi_file):
                                	                                xisys.append(smu__intxi_quickload(smu__intxi_file))
                                        	                 else:
                                                	                xisys.append(smu__intxi_calcwrite(smufile, smusettings_data, writetofile=False, smu__intxi__settings=smu__intxi__settings,rebinxi=rebinxi)[2])

							polyfitdeg = polyfitdeg
							xihatsys_base  = [ normfun(X) for X in xisys_base]
							xihatsys  = [ normfun(X) for X in xisys]
							#print get_avg_array(xisys)
							#print 'i_redshiftbin, xihatsys, xihatdata', i_redshiftbin, get_avg_array(xihatsys), xihatdata  # HB
							dxisys = get_diffarray(get_avg_array(xihatsys_base),get_avg_array(xihatsys))
							dxisys = polyfitY(range(len(dxisys)), dxisys, deg=polyfitdeg)

                                                        xisys_list.append(get_avg_array(xihatsys))
							dxisys_list.append([x for x in dxisys])

                                                        xicov = []
                                                        for imock in range(cov_nummock):
                                                                smufile = Tpcfrltfilename(binsplittedfilename(mockfile(catname, cov_catname2, imock, cov_RSDstr), ibin+1, totbin), mubins=nummubins,nbins=Smax_covmock,rmax=Smax_covmock ); smu__intxi_file = smu__intxi_filename(smufile, smu__intxi__settings=smu__intxi__settings_orig)
								if isfile(smu__intxi_file):
                                                                	xicov.append(smu__intxi_quickload(smu__intxi_file))
								else:
									xicov.append(smu__intxi_calcwrite(smufile, smusettings_mock, writetofile=False, smu__intxi__settings=smu__intxi__settings_orig,rebinxi=rebinxi)[2])

                                                                xihatcov  = [ normfun(X) for X in xicov]
								#if polyfit_dxi== None:
                                                                # xihatcov  = [ normfun(X) for X in xicov]
								#else:
                                                                # xihatcov  = [ normfun(polyfitY(range(len(X)),X,polyfit_dxi)) for X in xicov]
							X0=range(len(xihatcov[0]))
							dxicov  = [ get_diffarray(xihatcov_base[row], xihatcov[row]) \
								for row in range(len(xihatcov))]
							if covmat_nondiag_rescale == False:
								covmats.append(get_covmat(transpose(dxicov)))
							else:
								covmats.append(covmat_adjust(get_covmat(transpose(dxicov))))
							#print '\n\n ibin, covmat = ', ibin, get_covmat(transpose(dxicov))
						covmat = covmats[i_redshiftbin]
						#if iomw == 0: print np.mat(covmat).I # BOSSLi
						chisq_nosyscor, like = chisq_like_cov_xbar(dxidata, covmat)
						chisq_syscor, like   = chisq_like_cov_xbar(XplusY(dxidata,dxisys_list[i_redshiftbin],b=-1), covmat)
						 
						#print nowomwstr, '/ redshiftbin=', i_redshiftbin, '/ dxi = ', dxidata, chisq_nosyscor, chisq_syscor
						if consider_tilt_in_chisq:
							 xisys_avg = get_avg_array(xisys_list)
							 X = range(len(xisys_avg))
							 slop_ref = polyfit(X,xisys_avg,1)[0]

							 xidata_avg = get_avg_array(xihatdata_now)
							 slop_now = polyfit(X,xidata_avg,1)[0]
							 #print nowomwstr, slop_ref, slop_now

							 chisq_rescale_fact = (slop_now / slop_ref)**2.0
							 chisq_nosyscor /= chisq_rescale_fact
							 chisq_syscor /= chisq_rescale_fact
						#if i_redshiftbin in [5]:
						#	chisq_nosyscor, chisq_syscor = 0, 0
						chisqs_nosyscor[nowomwstr].append(chisq_nosyscor)  
						chisqs_syscor[nowomwstr].append(chisq_syscor) 
					i_redshiftbin += 1
			X1, X2 = chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
			str1, str2, strA, strB = \
				array_to_str(X1), array_to_str(X2), str(sum(X1)), str(sum(X2))
			#print str(mumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2
			nowf.write(str(mumin)+' '+nowomwstr+'   '+strA+' '+strB+'  '+ str1+'   '+str2+'\n')
			iomw += 1	
			#if iomw == 2:
			#	sys.exit()
			#print om, w, chisqs_nosyscor[nowomwstr], chisqs_syscor[nowomwstr]
		#print chisqs_nosyscor.keys()
		#sys.exit(' ')
		chisqlist_nosyscor = [ [  sum(chisqs_nosyscor[omwstr(om,w)])    for om in omlist] for w in wlist ]
		chisqlist_syscor =   [ [  sum(chisqs_syscor[omwstr(om,w)])      for om in omlist] for w in wlist ]
		nowf.close()
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
			fig.suptitle(nowchisqstr+'.mumin'+str(mumin))
			if savefig: 
				figname = outputdir+'/'+nowchisqstr+'.mumin'+str(mumin)+'.png';
				print '   Figure saved: ', figname
				plt.savefig(figname, fmt='png')
			if showfig: plt.show()

