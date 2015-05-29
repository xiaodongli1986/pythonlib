execfile('/home/xiaodongli/software/pythonlib/stdA.py');

def BSK_plotinfo(ax, infofile, nbin=20, Xrange=[0,15], lw=2, c='k', ls='-'):
    print 'Loading ', infofile, '...'
    data = np.loadtxt(infofile);
    Len, Mu = XYfromdata(data, 2,3)
    Mu = get_absarray(Mu)
    Lenavs, Lenavers, Muavs, Muavers = binned_quan_meaner(Len, Mu, nbin=nbin, Xrange=Xrange)
    ax.errorbar(Lenavs, Muavs, Muavers, label=infofile, c=c, lw=lw, ls=ls);
    ax.legend(frameon=False)    
def BSK_plotnumberinfo(ax, infofile, nbin=20, Xrange=[0,15], lw=2, c='k', ls='-'):
    print 'Loading ', infofile, '...'
    data = np.loadtxt(infofile);
    Len, Mu = XYfromdata(data, 2,3)
    rlt = ax.hist(Len, bins=nbin, range=Xrange,label=infofile, color=c, linewidth=lw, linestyle=ls);
    ax.legend(frameon=False)    
    return rlt 
