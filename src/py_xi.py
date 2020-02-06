
import numpy as np
import os, sys
import matplotlib.pyplot as plt

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

printstr = 'Usage: py_xi filename [-smin 0] [-smax 150] [-sbin 150] [-mubin 120] [-intxismin 6] [-intxismax 40] [-savefig T/F] [-showfig T/F]'+\
        '\n\t[-mubin_pack_rat 1] reduces the bins by pack-up small bins; should be an integer 1,2,3,...; '+\
        '\n\t[-mumax 1] maximal value of mu in plot; can impose a cut (e.g. set -mumax 0.97) to remove FOG'+\
        '\n\t[-figname YOUR_FIGURE_FILE_NAME] '+\
        '\n\t[-figtitle YOUR_FIGURE_NAME] '+\
        '\n\t[-labs hand-given labels of curves, separted using comma; e.g. \'curve1,curve2,curve3\']'


print printstr+'\n'
args = sys.argv

if len(args) <=1:
    print printstr
    sys.exit()

filenames = os.popen('ls '+args[1]).read().split()

def tf_readin(str2):
     if str2[0] in ['t', 'T']:
                return True
     elif str2[0] in ['f', 'F']:
                return  False
     else:
                print 'Error (py_xi): intput option not recognizable! ', str2
                print printstr
                sys.exit()

smin=0; smax=150; sbin=150; mubin=120; intxismin=6; intxismax=40; mubin_pack_rat=1; fs=14; savefig = True; showfig = True; mumax = 1.0; 
figtitle=None; figname=None; labs=[];

for icmd in range(2,len(args)):
    if icmd%2 == 1: 
        continue
    else:
        str1, str2 = args[icmd], args[icmd+1]
        if str1 == '-smin':
            smin = float(str2)
        elif str1 == '-smax':
            smax = float(str2)
        elif str1 == '-sbin':
            sbin = int(str2)
        elif str1 == '-mubin':
            mubin = int(str2)
        elif str1 == '-intxismin':
            intxismin = int(str2)
        elif str1 == '-intxismax':
            intxismax = int(str2)
        elif str1 == '-mubin_pack_rat':
            mubin_pack_rat = int(str2)
        elif str1 == '-savefig':
            savefig = tf_readin(str2)
        elif str1 == '-mumax':
            mumax = float(str2)
        elif str1 == '-figtitle':
            figtitle = str2
        elif str1 == '-figname':
            figname = str2
        elif str1 == '-showfig':
            showfig = tf_readin(str2)
        elif str1 == '-labs':
            labs = str2.split(',')
        else:
            print 'Error (py_xi): cmd not recognizable! ', str1, str2
            print printstr
            sys.exit()

print 'Plot CF for ', len(filenames), 'files...'
for filename in filenames:
    print '\t',filename
print ' * smin/smax/sbin = ', smin, smax, sbin
print ' * mubin/intxismin/intxisma = ', mubin, intxismin, intxismax
print ' * mubin_pack_rat/mumax = ', mubin_pack_rat, mumax

fig, axs = plt.subplots(3,2,figsize=(18,15))
axs = axs.reshape(-1)

used_labs = []

for ifile, filename in enumerate(filenames):
    DDnorm, DRnorm, RRnorm = [float(xx) for xx in open(filename, 'r').readline().split()[1:4]]
    data = np.loadtxt(filename)
    DD, DR, RR = [data[:,row].reshape(sbin,mubin) for row in [3,4,6]]
    DD /= DDnorm; DR /= DRnorm; RR /= RRnorm

    DDs = DD[:,:].sum(1); 
    DRs = DR[:,:].sum(1); 
    RRs = RR[:,:].sum(1); 
    xis = (DDs - 2*DRs + RRs) / RRs
    
    try:
        lab = labs[ifile]
    except:
        lab = filename
    used_labs.append(lab)
    #print('lab = ', lab)

    X = np.linspace(smin,smax,sbin+1)
    X = (X[1:] + X[:len(X)-1]) / 2.
    axs[0].plot(X,xis,label=lab); axs[0].set_xlabel('s',fontsize=fs); axs[0].set_ylabel('$\\xi$',fontsize=fs); axs[0].set_yscale('symlog')
    axs[1].plot(X,meannorm(xis)); axs[1].set_xlabel('s',fontsize=fs); axs[1].set_ylabel('$\\xi$, normed',fontsize=fs); axs[1].set_yscale('symlog')
    axs[2].plot(X,X**2*xis); axs[2].set_xlabel('s',fontsize=fs); axs[2].set_ylabel('$s^2\\ \\xi$',fontsize=fs) 
    axs[3].plot(X,meannorm(X**2*xis)); axs[3].set_xlabel('s',fontsize=fs); axs[3].set_ylabel('$s^2\\ \\xi$, normed',fontsize=fs) 

    imumax = int(mubin*mumax) 
    DD, DR, RR = DD[:,:imumax], DR[:,:imumax], RR[:,:imumax] 
    DD, DR, RR = packarray2d(DD, mubin_pack_rat), packarray2d(DR, mubin_pack_rat), packarray2d(RR, mubin_pack_rat) 
    xi = np.divide(DD-2*DR+RR,RR)
    
    ds = (smax-smin)/float(sbin)
    is1, is2 = int((intxismin-smin)/ds), int((intxismax-smin)/ds)
    intximu = (xi[is1:is2+1,:].sum(0)*ds)[::-1]

    X = np.linspace(1-mumax,1,len(intximu)+1)
    X = (X[1:] + X[:len(X)-1]) / 2.
    axs[4].plot(X,intximu); axs[4].set_xlabel('$1-\\mu$',fontsize=fs); 
    axs[4].set_ylabel('$\\xi_{\\Delta s}, s_1='+str(intxismin)+',s_2='+str(intxismax)+'$',fontsize=fs)
    axs[5].plot(X,meannorm(intximu)); axs[5].set_xlabel('$1-\\mu$',fontsize=fs);    
    axs[5].set_ylabel('$\\hat\\xi_{\\Delta s}, s_1='+str(intxismin)+',s_2='+str(intxismax)+'$',fontsize=fs)

#print('used_labs = ', used_labs)
legfs = (60./  max([len(lab) for lab in used_labs]))**0.5 * 10; legfs = min(legfs, 16)
for ax in axs:
   ax.grid()
   ax.legend(frameon=False,fontsize=legfs)

if figtitle==None:
    fig.suptitle(figname,fontsize=16)
else:
    fig.suptitle(figtitle,fontsize=16)

if figname==None:
    figname = filenames[0]+'.'+str(len(filenames))+'files.plotxi.png'
else:
    figname = figname+'.png'

print 'figure saved: ', figname
fig.savefig(figname, format='png')
if showfig:
    plt.show()


