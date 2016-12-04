
#EXE=/home/xiaodongli/software/kstat/bin/3pcf_eq
execfile('/home/xiaodongli/software/pythonlib/stdA.py')

### Definitions
def TriCF(ddd, ddr, drr, rrr):
	return (ddd-3*ddr+3*drr-rrr)/rrr


### Read in data
def TriCFeq_readin(outfile, nmu=12, nbins=40, printinfo=False):
 data = [[[ 0 for imu2 in range(nmu)] for imu1 in range(nmu)] for ir in range(nbins)]
 nowf = open(outfile, 'r')
 numline = 0
 for ir in range(nbins):
	for imu1 in range(nmu):
		for imu2 in range(nmu):
			nowstrs = nowf.readline().split()
			ddd, ddr, drr, rrr, zeta1, zeta2 = [float(xx) for xx in nowstrs[6:12]]
			#print ir, imu1, imu2
			data[ir][imu1][imu2] = [ddd,ddr,drr,rrr, zeta1, zeta2]
			#print data[ir][imu1][imu2]
			numline += 1
 if printinfo: print 'TriCFeq_readin: Finishing read in ', numline, 'lines'
 nowf.close()
 return data

def TriCFeq_plot_funr(outfile, savefig=True, figname = None):
 data = TriCFeq_readin(outfile)

 Ys = [0 for row in range(nbins)]
 for ir in range(nbins):
	ddd = ddr = drr = rrr = 0
	for imu1 in range(nmu):
		for imu2 in range(nmu):
			ddd0, ddr0, drr0, rrr0 = data[ir][imu1][imu2]
			ddd+=ddd0; ddr+=ddr0; drr+=drr0; rrr+=rrr0
			Ys[ir] = [ddd,ddr,drr,rrr]
 DDD, DDR, DRR, RRR = Xsfromdata(Ys, [0,1,2,3])
 X = [ir+0.5 for ir in range(nbins)]
 fig, ax = figax()
 ax.plot(X,DDD, label='DDD')
 ax.plot(X,DDR, label='DDR')
 ax.plot(X,DRR, label='DRR')
 ax.plot(X,RRR, label='RRR')
 ax.legend(frameon=False)
 ax.grid()
 ax.set_xlabel('r [Mpc/h]', fontsize=18)
 ax.set_ylabel('3-point CF', fontsize=18)
 fig.show()
 if savefig: 
	if figname == None:
		figname = outfile+'.Fun_r.png'
	fig.savefig(figname, format='png')

