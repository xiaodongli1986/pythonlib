import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as Axes3D
import commands, os, random, scipy.integrate, copy
import scipy.ndimage as ndimage

########################################
### Basic definition
########################################

### dataname bm (binned mu): information of binned mu
class binnedmuinfo:
    om  = 100.0; w = 100.0; 
    num_bin = 0; nbinlist = []; # how many bins; list of #-bin;
    num_inputed = [0 for row in range(num_bin)]; ## 
    # # of read in; may read in many results using same #-bins, take the average; ...
    mult_redshiftlist = [[] for row in range(num_bin)]; ## list of redshift bins
    mult_rlist = [[] for row in range(num_bin)]; ## list of r bins
    mult_muavlist= [[] for row in range(num_bin)]; ## list of muav bins
    mult_muerlist = [[] for row in range(num_bin)]; ## list of muaver bins
    def __init__(self, num_bin=14): ## initilize: set everything as empty list;
        self.num_bin = num_bin; self.om=100.0; self.w=100.0;
        self.num_inputed = [0 for row in range(self.num_bin)];
        self.mult_redshiftlist = [[] for row in range(self.num_bin)]; 
        self.mult_rlist = [[] for row in range(self.num_bin)];
        self.mult_muavlist= [[] for row in range(self.num_bin)];
        self.mult_muerlist = [[] for row in range(self.num_bin)];
        return
    def update(self, ibin, ipt_redshiftlist, ipt_rlist, ipt_muavlist, ipt_muerlist): ## update ibin-th result;
        if ibin >= self.num_bin or ibin < 0:
            print 'ERROR (binnedmuinfo.update): wrong ibin: ', ibin, self.num_bin
            return
        for combin in [(self.mult_redshiftlist, ipt_redshiftlist), 
                        (self.mult_rlist,ipt_rlist), 
                        (self.mult_muavlist,ipt_muavlist), 
                        (self.mult_muerlist,ipt_muerlist)]:
            mult_list, ipt_list = combin
            if mult_list[ibin] == []:
                mult_list[ibin] = copy.deepcopy(ipt_list)
            else:
                if len(ipt_list) != len(mult_list[ibin]):
                    print 'ERROR (binnedmuinfo.update): mismatch of list length! ibin, len1, len2 = ', \
                        ibin, len(ipt_list), len(mult_list[ibin])
                for row in range(len(mult_list[ibin])):
                    mult_list[ibin][row] = (mult_list[ibin][row]*self.num_inputed[ibin] + ipt_list[row]) \
                        / float(self.num_inputed[ibin]+1)
        self.num_inputed[ibin] += 1
        self.nbinlist = [len(self.mult_muavlist[row]) for row in range(len(self.mult_muavlist))]

########################################
### Loading in
########################################

### read in bm from a mufile
def get_bm_from_mufile(mudatafile, om0, w0, num_nbin, num_read, action_is_update = False):
    if not action_is_update:
        bm = binnedmuinfo(num_nbin)
    f = open(mudatafile,'r')
    for i_read in range(num_read):
        floatlist = str_to_numbers(f.readline());
        nowiom, nowiw, om, w = floatlist ### fmt of first line: iom, iw, om, w
        if abs(om-om0)>0.0001 or abs(w-w0)>0.0001:
            print 'ERROR (get_chisq)!: Mismatched omegams /ws: ', om,om0, w, w0
            return
        bm.om = om; bm.w = w;
        for ibin in range(num_nbin):
            floatlist=str_to_numbers(f.readline());
            nownbin = floatlist[0] ## nbin
            bm.update(ibin,str_to_numbers(f.readline()),str_to_numbers(f.readline()),
                      str_to_numbers(f.readline()),str_to_numbers(f.readline())) ### list of redshift, r, muav, muaver
    f.close()
    return bm

def loadin_bmlist(datanamestr, gfstr, mutypestr, omlist=[], wlist=[], idrop=1,num_read=1,imocklist=[0],ipatchlist=[1]):
    bmlists = [[[ 0 for iom in range(numom)] for iw in range(numw)] for idata in range(numdata) ] 
    idata = -1
    for imock in imocklist:
        for ipatch in ipatchlist:
            idata += 1
            for iom in range(len(omlist)):
                for iw in range(len(wlist)):
                    om = omlist[iom]; w=wlist[iw]
                    filename = get_mudatafile(dirstr, datanamestr, gfstr, mutypestr, imock, ipatch, 
                                                              om, w, idrop=idrop)
                    bmlists[idata][iw][iom] = get_bm_from_mufile(filename, om, w, num_ibin, num_read=num_read)
    return bmlists

########################################
### About muav list
########################################

### Get the difference between muav for bm1 & bm2 (bm2-bm1)
def get_multmuavlist_diff(bm1, bm2):
    num_bin = bm1.num_bin
    diffmuavlist = [];
    for ibin in range(num_bin):
        muavlist1 = bm1.mult_muavlist[ibin]
        muavlist2 = bm2.mult_muavlist[ibin]
        diffmuavlist.append([muavlist2[row]-muavlist1[row] for row in range(len(muavlist2))])
    return diffmuavlist

## Getting muav/muerlist from a list of bm;
# error estimated from variance
def get_muaverlists_from_bms(bmlist, calcer=True):
    numbm = len(bmlist)
    num_nbin = len(bmlist[0].mult_muavlist)
    mult_muavlists, mult_muerlists = [0 for row in range(num_nbin)], [0 for row in range(num_nbin)]
    for ibin in range(num_nbin):
        if calcer:
            mult_muavlists[ibin], mult_muerlists[ibin] = get_avgstd_array([bmlist[ibm].mult_muavlist[ibin] for ibm in range(numbm)])
        else:
            mult_muavlists[ibin] = get_avg_array([bmlist[ibm].mult_muavlist[ibin] for ibm in range(numbm)])
            mult_muerlists[ibin] = get_avg_array([bmlist[ibm].mult_muerlist[ibin] for ibm in range(numbm)])
    return mult_muavlists, mult_muerlists

def get_mult_chisqs_from_bmlist(bmlists, fixidata=False,fixedidata=0, calcer = True, diffchisqmethod='use_weightedavg_as_ref'):
    mult_chisqs = [[0 for iom in range(numom)] for iw in range(numw)]
    for iw in range(numw):
        for iom in range(numom):
            mult_muavlist, mult_muerlist = get_muaverlists_from_bms([bmlists[idata][iw][iom] for idata in range(numdata)], calcer)
            if fixidata:
                mult_muavlist = bmlists[fixedidata][iw][iom].mult_muavlist
            mult_chisqs[iw][iom] = get_chisqs_from_yyerlists(mult_muavlist, mult_muerlist,diffchisqmethod=diffchisqmethod)
    return mult_chisqs

########################################
### Calculating chisq
########################################

def calc_multchisq_April23(bm,ipt_ibinlist=[],ipt_diffmuavlist=[],ax='',extraleg='',plotibin = 13,use_polyfit=False):
        chisqlist = []; nbinlist = [];
        if ipt_ibinlist == []:
            ibinlist = range(bm.num_bin)
        else:
            ibinlist = ipt_ibinlist
        for ibin in ibinlist:
            muavlist = copy.deepcopy(bm.mult_muavlist[ibin])
            if ipt_diffmuavlist!= []:
                for row in range(len(muavlist)):
                    muavlist[row] -= ipt_diffmuavlist[ibin][row]
            muerlist = bm.mult_muerlist[ibin]
            nownbin = len(muavlist); summu = 0; sumwei = 0
            for j in range(nownbin):
                summu  += muavlist[j]/(muerlist[j])**2.0
                sumwei += + 1.0/(muerlist[j])**2.0

            muav = summu / sumwei
            if use_polyfit:
                if ibin >= 7:
                    pforder = 3
                elif ibin >= 3:
                    pforder = 2
                else:
                    pforder = 1
                coefficients = np.polyfit(redshiftlist, muavlist, pforder );
                polynomial = np.poly1d(coefficients);
                for j in range(len(muavlist)):
                    muavlist[j] = polynomial(redshiftlist[j])

            Y = [((muavlist[row]-muav)/muerlist[row])**2.0 for row in range(len(muavlist))]
            chisqlist.append(sum(Y))
            nbinlist.append(nownbin)
            if ax != '' and ibin == plotibin:
                ax.errorbar(self.mult_redshiftlist[ibin], muavlist, muerlist, label = extraleg+'$\\chi^2$: %.2f'%sum(Y))
        return chisqlist#sum(chisqlist)/float(len(chisqlist))#, chisqlist
    
def chisq_from_chisqlist(ipt_nbinlist, ipt_chisqlist, ipt_ibinlist=[]):
    if ipt_ibinlist == []:
        nbinlist = ipt_nbinlist
        chisqlist = ipt_chisqlist
    else:
        nbinlist = [ipt_nbinlist[ipt_ibinlist[row]] for row in range(len(ipt_ibinlist))]
        chisqlist = [ipt_chisqlist[ipt_ibinlist[row]] for row in range(len(ipt_ibinlist))]
    #weilist = [1.0/sqrt(nbinlist[row]) for row in range(len(nbinlist))]
    weilist = [1.0 for row in range(len(nbinlist))]
    sumchisq = 0; sumwei = 0

    for i in range(len(chisqlist)):
        sumchisq += chisqlist[i]*weilist[i]
        sumwei += weilist[i]
    chisq = sumchisq / sumwei
    return chisq

########################################
### Plotting
########################################

def ourerrorbar(ax, X, Y, YErr, lw=1, c='k', ls='-', capsize=5, label='', marker='o', ms=8, not_dof_div=False, polyorder=1):
    chisq = get_chisq_from_yyerlist(Y, YErr)
    if not not_dof_div:
        ax.errorbar(X[0],Y[0],YErr[0],lw=lw,c=c,ls=ls,capsize=capsize,marker=marker,ms=ms,
                label=label+': $\\chi^2/{\\rm dof}=%.2f'%(chisq/float(len(Y)))+'$')
    else:
        ax.errorbar(X[0],Y[0],YErr[0],lw=lw,c=c,ls=ls,capsize=capsize,marker=marker,ms=ms,
                label=label+': $\\chi^2=%.2f'%(chisq)+'$')
    ax.errorbar(X,Y,YErr, lw=lw,c=c,ls='.',capsize=capsize,marker=marker,ms=ms,)
    coef = np.polyfit(X, Y, polyorder)
    xl = X[0] - (X[1]-X[0]); xr = X[len(X)-1] + (X[len(X)-1]-X[len(X)-2])
    NewX = [xl]+X+[xr]
    ax.plot(NewX, [np.polyval(coef, x) for x in NewX], lw=lw, c=c, ls=ls)
    return xl, xr
    
def do_plot(dirstr_lab_list, datanamestr_lab_list, gfstr_lab_list, omwlist, imocklist, ipatchlist, add_extrafigname=True,
            figname='',figtitlestr='',figxsize=14,figysize=5,oneplot=False, quantypestr = 'mu',idrop=1, polyorder=1,
            calcer=False, not_dof_div=False, plot_sing_point=True, musub=0.0, ipt_doavg=True, nolegend=False,Plot_NoSub=False,
            MinorMutypestr='.noabs',MajorMutypestr='',printfilename=False,outputformat='jpg',  
            compctXlim=True,  compctXlim_released=True, legfs=12, tiltfs=15, labelfs=15, figtop=0.9,
            xlabel='$redshift$', ylabel='$\\bar\\mu$',
            Sub_FirstTwoMuDiff=False, Sub_FirstMU=False, NormTo0_MU=False ):
    if Sub_FirstTwoMuDiff and Sub_FirstMU:
        print 'ERROR! Sub_FirstTwoMuDiff, Sub_FirstMU shall not be both True!'
        return
    idata=-1; iplot=-1; i_plotdone = 0;
    lclist = ['r', 'b', 'g', 'k', 'y', 'c', 'grey', 'm']
    lslist = ['-', '--', '--','-.', ':']
    markerlist = ['o', '*', 'D', 'p', '+', 's']
    fig = plt.figure(figsize=(figxsize,figysize)); figtitle = figtitlestr
    if oneplot: # one plot: abs mu (you will see the magnitude of line is around 0.5)
        axA = fig.add_subplot(111); axmutypelist = [[axA, MajorMutypestr]]
    else: # two plot: abs mu and mu
        axA = fig.add_subplot(121); axB = fig.add_subplot(122); axmutypelist = [[axA, MinorMutypestr], [axB, MajorMutypestr]]
    Xls, Xrs = [], [];
    # loop: dirstr -> datanamestr -> gfstr -> omwstr -> (axmutype) -> imock -> ipatch
    for dirstr_lab in dirstr_lab_list:
        dirstr, dirlab = dirstr_lab
        for datanamestr_lab in datanamestr_lab_list:
            datanamestr, dtlab = datanamestr_lab; idata +=1; igf=-1;
            # Whether do average (over different mocks & patches)
            if datanamestr[0:4] == 'Data': # now multi relizations only for mocks; hope data having soon..
                doavg = False
            else:
                doavg = True and ipt_doavg and (len(imocklist)*len(ipatchlist)>1)
            for gfstr_lab in gfstr_lab_list:
                gfstr, gflab = gfstr_lab; igf+=1; iomw = -1; 
                for omw in omwlist:
                    om, w = omw; omwstr = ('om%.3f'%om+', w%.2f'%w); labelstr=''; iomw+=1; lw=3;                     
                    lc = lclist[np.mod(iplot,len(lclist))]; # color of the line
                    ls = lslist[np.mod(idata,len(lslist))]; # style of the line
                    marker = markerlist[np.mod(iplot,len(markerlist))] # style of marker
                    labelstr = ''
                    if dirlab != '':
                        labelstr += (dirlab+', ')
                    if dtlab != '':
                        labelstr += (dtlab+', ')
                    if gflab != '':
                        labelstr += (gflab+', ')
                    if len(omwlist) == 1:
                        if add_extrafigname:
                            figtitle += omwstr
                    else:
                        labelstr += omwstr
                    allfound = True
                    for axmutype in axmutypelist: 
                        ax, mutypestr = axmutype
                        if doavg:
                            Zs, MUs, MUERs = [], [], []
                            for imock in imocklist:
                                if not allfound: # skip if not all found!!!
                                    break
                                for ipatch in ipatchlist:
                                    filename = get_mudatafile(dirstr=dirstr,datanamestr=datanamestr, gfstr=gfstr, mutypestr=mutypestr, 
                                                    quantypestr=quantypestr, imock=imock, ipatch=ipatch, idrop=idrop, om=om,w=w)
                                    if printfilename:
                                        print filename
                                    if os.path.isfile(filename) == False:
                                        allfound = False;  break
                                    bm = get_bm_from_mufile(filename, om, w, num_nbin, num_read)
                                    X1 = bm.mult_redshiftlist[ibin]; Zs.append([X1[row] for row in range(len(X1))])
                                    X2 = bm.mult_muavlist[ibin]; MUs.append([X2[row] for row in range(len(X2))])
                                    X3 = bm.mult_muerlist[ibin]; MUERs.append([X3[row] for row in range(len(X3))])
                                    if plot_sing_point:
                                        ax.plot(X1, X2, ls=ls, c=lc, lw=0.1) # also show each single plot even if do_avg
                            if not calcer: # adopt the error bars available in the mu file
                                Z = get_avg_array(Zs)
                                MU = get_avg_array(MUs)
                                MUER = get_avg_array(MUERs)
                            else: # calculate error bar using the variance of results from different realizations
                                Z, ZEr = get_avgstd_array(Zs)
                                MU, MUER = get_avgstd_array(MUs)
                        else:
                            filename = get_mudatafile(dirstr=dirstr,datanamestr=datanamestr, gfstr=gfstr, mutypestr=mutypestr, 
                                                        quantypestr=quantypestr, imock=0, ipatch=ipatchlist[0], idrop=idrop, om=om,w=w)
                            if printfilename:
                                print filename
                            if os.path.isfile(filename) == False:
                                allfound = False
                            else:
                                bm = get_bm_from_mufile(filename, om, w, num_nbin, num_read)
                                Z = bm.mult_redshiftlist[ibin]; MU = bm.mult_muavlist[ibin]; MUER = bm.mult_muerlist[ibin]
                        if allfound: # make plot only if allfound
                            if iplot != 0:
                                MU = [MU[row]-musub for row in range(len(MU))]
                            if Sub_FirstMU:
                                if i_plotdone == 0:
                                    MUSUB = [MU[row] for row in range(len(MU))]
                                    MU = [0 for row in range(len(MU))]
                                else:
                                    SubedMU = [MU[row]-MUSUB[row] for row in range(len(MU))];
                                    if NormTo0_MU:
                                        MUAV = sum(SubedMU) / float(len(SubedMU))
                                        SubedMU = [SubedMU[row]-MUAV for row in range(len(SubedMU))]
                                    lc = lclist[np.mod(iplot,len(lclist))]; ls = lslist[np.mod(idata,len(lslist))];
                                    xl, xr = ourerrorbar(ax,Z,SubedMU,
                                        MUER,lw=lw,c=lc,ls=ls,label=labelstr+' (Curve1 subtracted)',
                                        marker=marker,not_dof_div=not_dof_div,polyorder=polyorder);
                                    Xls.append(xl);Xrs.append(xr); i_plotdone += 1; iplot+=1;
                            if Sub_FirstTwoMuDiff:
                                if i_plotdone == 0:
                                    MU1 = [MU[row] for row in range(len(MU))]; MUSUB = [0 for row in range(len(MU))]
                                elif i_plotdone == 1:
                                    MU2 = [MU[row] for row in range(len(MU))]; MUSUB = [MU2[row]-MU1[row] for row in range(len(MU))]
                                else:
                                    SubedMU = [MU[row]-MUSUB[row] for row in range(len(MU))];
                                    if NormTo0_MU:
                                        MUAV = sum(SubedMU) / float(len(SubedMU))
                                        SubedMU = [SubedMU[row]-MUAV for row in range(len(SubedMU))]
                                    lc = lclist[np.mod(iplot,len(lclist))]; ls = lslist[np.mod(idata,len(lslist))];
                                    xl, xr = ourerrorbar(ax,Z,SubedMU,
                                        MUER,lw=lw,c=lc,ls=ls,label=labelstr+' (Curve2-Curve1 subtracted)',
                                        marker=marker,not_dof_div=not_dof_div,polyorder=polyorder);
                                    Xls.append(xl);Xrs.append(xr); i_plotdone += 1;  iplot+=1;
                            if NormTo0_MU:
                                MUAV = sum(MU) / float(len(MU))
                                MU = [MU[row]-MUAV for row in range(len(MU))]
                            if ((not Sub_FirstMU) and (not Sub_FirstTwoMuDiff)) or Plot_NoSub \
                                or (Sub_FirstMU and i_plotdone==0) or (Sub_FirstTwoMuDiff and i_plotdone<=1):
                                lc = lclist[np.mod(iplot,len(lclist))]; ls = lslist[np.mod(idata,len(lslist))];
                                xl, xr = ourerrorbar(ax,Z,MU,MUER,lw=lw,c=lc,ls=ls,label=labelstr,
                                                     marker=marker,not_dof_div=not_dof_div,polyorder=polyorder);
                                Xls.append(xl);Xrs.append(xr); i_plotdone += 1; iplot += 1;
                            
    for axmutype in axmutypelist:
        ax, mutypestr = axmutype
        ax.grid(ls='-',c='g')
        ax.set_xlabel(xlabel,fontsize=labelfs);ax.set_ylabel(ylabel,fontsize=labelfs);
        if not nolegend:
            ax.legend(loc='best', prop={'size':legfs, 'family':"Times New Roman"},frameon=False)
        if Xls!=[] and Xrs!=[] and compctXlim:
            if compctXlim_released:
                ax.set_xlim(min(Xls),max(Xrs))
            else:
                ax.set_xlim(max(Xls),min(Xrs))
    fig.tight_layout()
    if figtitle != '':
        fig.suptitle(figtitle,fontsize=tiltfs); fig.subplots_adjust(top=figtop)
    if figname!='':
        print 'figure saved: ', figname
        fig.savefig(figname,format=outputformat)
    return fig, ax

####
## Plot the curve of mu for given set of data, cosmological parameters...   

### Contour plot
def plot_contour(ax, omlist, wlist, chisqlist, label='NO RSD',
                    ommin = 0.01, ommax = 0.6, wmin = -2.0, wmax = -0.0,  do_smooth=True, smsigma=0.5, 
                    extratitle = '', titleftsize=15, notitle = False,
                    sigA = 0.683, sigB = 0.954, sigC = 0.997,  nolegend = False, nolabel = False, legftsize=15,
                    noxticks = False, noyticks = False, showgrid = False, use_ratCL = True, plotformat = 1):
    if True:
        smsigma = smsigma
        if do_smooth:
            Z = ndimage.gaussian_filter(chisqlist, sigma=smsigma, order=0)
        else:
            Z = chisqlist
        if use_ratCL:
            chisqA, chisqB, chisqC = list_find_CL_chisqcut(chisqlist, [sigA, sigB, sigC])
            #chisqA, likeratioA = find_CL_chisqcut(Z, sigA)
            #chisqB, likeratioB = find_CL_chisqcut(Z, sigB)
            #chisqC, likeratioC = find_CL_chisqcut(Z, sigC)
            #print chisqA, likeratioA, chisqB, likeratioB, chisqC, likeratioC
        else:
            chisqA = 2.3; chisqB = 6.17; chisqC = 11.83;

        X = [-100,-99]; Y=[0,0]
        if plotformat == 1:
            ax.contourf(omlist, wlist, Z, [0, chisqC], colors='0.75')
            ax.contourf(omlist, wlist, Z, [0, chisqB], colors='0.65')
            ax.contourf(omlist, wlist, Z, [0, chisqA], colors='0.55' )
            ax.plot(X,Y,c='0.55',lw=10,label=label)
        elif plotformat == 2:
            ax.contour(omlist, wlist, Z, [chisqA, chisqB, chisqC], colors='r', linewidths = 2, linestyles = 'dashed')
            ax.plot(X,Y,c='r',lw=3,ls='--',label=label)
        elif plotformat == 3:
            ax.contour(omlist, wlist, Z, [chisqA, chisqB, chisqC], colors='b', linewidths = 2)
            ax.plot(X,Y,c='b',lw=3,ls='-',label=label)
            
        ax.scatter([0.26], [-1], marker = '+', c = 'g', s = 200, lw=2)        
        
        ax.set_xlim(ommin,ommax);   ax.set_ylim(wmin, wmax)
        if showgrid:
            grid(b=True, which='major', color='g', linestyle='-');
            grid(b=True, which='minor', color='g', linestyle='--');
        if not notitle:
            ax.set_title(extratitle, fontsize=titleftsize)
        if not nolegend:
            ax.legend(loc='upper left',prop={'size':legftsize},frameon=False)
        if not nolabel and not noxticks:
            ax.set_xlabel('$\Omega_m$', fontsize=26)
        if not nolabel and not noyticks:
            ax.set_ylabel('$w$', fontsize=26)
        if noyticks:
            for ylabel_i in ax.get_yticklabels():
                ylabel_i.set_fontsize(0.0)
                ylabel_i.set_visible(False)
        if noxticks:
            for xlabel_i in ax.get_xticklabels():
                xlabel_i.set_fontsize(0.0)
                xlabel_i.set_visible(False)
