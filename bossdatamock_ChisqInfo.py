
### paths for chisqs of our method and cosmomc chains
bossdatamock_bosschisqfiledir = '/home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/covmats/'
bossdatamock_cosmomcchaindir = '/home/xiaodongli/software/cosmomcs/12Oct_generic/chains/boss2pcfAP/'

### Chisqfiles for Dense1 & Dense2

Dense1chisqfile = '/home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/covmats/'+\
    'Dense1sub--25-51-Scan--NSGC-6bins-J08-RSD-refibin0--int--xi.40mu0.15to1.00.s6.00to40.00.1275chisqs.381files.Corrected'

Dense2chisqfile = '/home/xiaodongli/SparseFilaments/data/local_input/boss2pcf/data/covmats/'+\
    'Dense2sub--71-45-Scan--NSGC-6bins-J08-RSD-refibin0--int--xi.40mu0.15to1.00.s6.00to40.00.3195chisqs.381files.9skiped.Corrected'

###
### 2d interpolated grid of chisqs for our method
def Dense12_2d_Interpolation(numomw=100, 
                             output_MCMC_file=False, 
                             conrange =  [0.1, 0.55, -1.9, -0.6], 
                             smsigma=0.0,
			     cosmomcchaindir = bossdatamock_cosmomcchaindir,
			     Dense1chisqfile = Dense1chisqfile, 
			     Desne2chisqfile = Dense2chisqfile,
			     doplot=True,
			     Set_Dense2_as_Dense1=True,
			     suffix = ''):
    '''
        Many quantities will be returned (1 means Dense1; 2 means Dense2; 3 is 2d interpolation of Dense2)
    
    omlist1, wlist1, numom1, numw1, ommin1, ommax1, wmin1, wmax1, chisqlist1, \
            omlist2, wlist2, numom2, numw2, ommin2, ommax2, wmin2, wmax2, chisqlist2, \
            omlist3, wlist3, numom3, numw3, ommin3, ommax3, wmin3, wmax3, chisqlist3 
    '''
    
    conommin, conommax, conwmin, conwmax = conrange
    
    ### Dense1 Info
    omlist1, wlist1 = Dense1subscan_omlist, Dense1subscan_wlist
    numom1, numw1 = len(omlist1), len(wlist1)
    ommin1, ommax1, wmin1, wmax1 = min(omlist1), max(omlist1), min(wlist1), max(wlist1)
    chisqlist1 = get_2darray_from_1d(Xfromdata(loadtxt(Dense1chisqfile),2), numw1, numom1)

    if not Set_Dense2_as_Dense1:
	    ### Dense2 Info
	    omlist2, wlist2 = Dense2subscan_omlist, Dense2subscan_wlist
	    numom2, numw2 = len(omlist2), len(wlist2)
	    ommin2, ommax2, wmin2, wmax2 = min(omlist2), max(omlist2), min(wlist2), max(wlist2)
	    chisqlist2 = get_2darray_from_1d(Xfromdata(loadtxt(Dense2chisqfile),2), numw2, numom2)
    else:
	    omlist2, wlist2, numom2, numw2, ommin2, ommax2, wmin2, wmax2, chisqlist2 = \
			   omlist1, wlist1, numom1, numw1, ommin1, ommax1, wmin1, wmax1, chisqlist1

    ### Plot of the original contour from Dense1 and Dense2
    if False:
        fig, ax = figax()
        if False:
            plot_contour(ax, omlist1, wlist1,  chisqlist1, ommin=conommin, ommax=conommax, wmin=conwmin, wmax=conwmax, 
                         smsigma=0.0,
                         sigs=[sig1,sig2], plotformat=2,scatter_WMAP5=False,  )
        if False:
            plot_contour(ax, omlist2, wlist2,  chisqlist2, ommin=conommin, ommax=conommax, wmin=conwmin, wmax=conwmax, 
                         smsigma=0.0,
                         sigs=[sig1,sig2], plotformat=1,scatter_WMAP5=False,  )
        ax.grid()
        plt.show()

    #######################################################
    ### Interpolated Contour

    for numomw in [numomw]:
    #for nuomw in [80, 100, 150, 200, 250, 300]:
        numom3 = numw3 = numomw #numom2, numw2 #

        ommin3, ommax3, wmin3, wmax3 = ommin2, ommax2, wmin2, wmax2
        omlist3, wlist3 = np.linspace(ommin3, ommax3, numom3), np.linspace(wmin3, wmax3, numw3)
        omwlist3 = []
        chisqlist3 = []
        for w in wlist3:
            omwlist3 += [[om,w] for om in omlist3]

        X=omlist2;Y=wlist2;chisqlist=chisqlist2;
        X3,Y3=XYfromdata(omwlist3);
        chisqlist3_1d = LinearInterpolation_2d_UniformGrid(min(X), max(X), len(X), min(Y),max(Y),len(Y), 
                chisqlist, X3, Y3)
        chisqlist3 = get_2darray_from_1d(chisqlist3_1d, numw3, numom3)

        if output_MCMC_file:
            ### Output chisqs to file
            MCMCfile = bossdatamock_cosmomcchaindir+'Dense2_Interpolated_'+str(numom3)+suffix
            nowf = open(MCMCfile, 'w')
            for row in range(len(X3)):
                nowf.write(str(X3[row])+'\t'+str(Y3[row])+'\t'+str(chisqlist3_1d[row])+'\n')
            nowf.close()
            MCMC_cosmomc_fmt_convert(MCMCfile, '_1.txt', )
            print '### outputing chisqs to file: ', MCMCfile
            

    ### Plot of the contour of interpolated result
    if doplot:
	    fig, ax = figax(figxsize=14,figysize=8)
	    plot_contour(ax, omlist3, wlist3, chisqlist3, 
			 ommin=conommin, ommax=conommax, wmin=conwmin, wmax=conwmax, 
		         smsigma=smsigma,
		         sigs=[sig1,sig2], plotformat=1,scatter_WMAP5=False,  )
	    ax.grid()
	    plt.show()
    return  omlist1, wlist1, numom1, numw1, ommin1, ommax1, wmin1, wmax1, chisqlist1, \
            omlist2, wlist2, numom2, numw2, ommin2, ommax2, wmin2, wmax2, chisqlist2, \
            omlist3, wlist3, numom3, numw3, ommin3, ommax3, wmin3, wmax3, chisqlist3 

