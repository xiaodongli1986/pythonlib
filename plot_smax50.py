
#py_Plot is the command to plot

execfile('/home/xiaodongli/software/pythonlib/Tpcftools.py')


smusettings = {
            'smin':0.0, 'smax':50.0, 'numsbin':50,
            'mumin':0.0, 'mumax':1.0, 'nummubin':50,
#            'mumin':0.0, 'mumax':1.0, 'nummubin':120,
            'deltas':0, 'deltamu':0, 'slist':[], 'mulist':[]
               }

is_sig_pi = True
#is_sig_pi = False

#filename =  "data.xyzw.finesplitted.rmax50.50rbins.120mubins.2pcf_20160329"
#filename =  "data.xyzw.finesplitted.rmax50.50rbins.120mubins.2pcf_20160330"
#filename = "2pcf_aniso_smu_rmax50_50rbins_120mubins"
filename = "data.xyzw.finesplitted.rmax50.50rbins.sigpi.2pcf_20160330"


smu__initsmusettings(smusettings)

if True:

		DDlist, DRlist, RRlist = Xsfrom2ddata(smu__loadin(filename, smusettings), [4,5,6])
		ismax = 50
		deltais = 2
		no_s_sq_in_y = False
		fig, ax1 = figax()

                ### packed count of xi as a function of s
                ismin = 0;
                imumin = 0; imumax = nummubin;
                sasx = []; packedxiasy = [];

                while ismin < 50:
                    ismax = min(ismin+deltais,numsbin-1);
                    nows=(slist[ismin]+slist[ismax])/2.0;   sasx.append(nows)
                    if not no_s_sq_in_y:
			    if is_sig_pi:
			    	Y = []
			    	for nowsig in range(ismax-1):
				 nowpi = int(np.sqrt(nows**2.0 - nowsig**2.0)) 
				 nowpi = max(nowpi,0)
				 Y.append(packedxi(DDlist,DRlist,RRlist,nowsig,nowsig+1,nowpi,nowpi+1))
			    	packedxiasy.append(nows*nows*sum(Y) / (len(Y)+0.0))
			    else:
                                packedxiasy.append(nows*nows*packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
                    else:
                            packedxiasy.append(packedxi(DDlist, DRlist, RRlist, ismin, ismax, imumin, imumax))
                    ismin += deltais;
                if True:
                        ax1.plot(sasx, packedxiasy, marker='o', markersize=1)
                        ax1.set_xlabel('$s\ [\\rm Mpc/h]$', fontsize=25)
                        if not no_s_sq_in_y:
                                ax1.set_ylabel('$s^2\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                        else:
                                ax1.set_ylabel('$\\xi\ [\\rm Mpc/h]^2$', fontsize=25)
                        #ax1.set_xlim(min(sasx),max(sasx))
                        ax1.set_xlim(0,50)
		fig.savefig(filename+'.png', format = 'png')
