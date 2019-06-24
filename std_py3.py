
import numpy as np

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def mean_erbar_str(A, fmt='%.3f'):
    A0, A1 = mean(A), np.sqrt(A.var())
    return '$'+fmt%A0+'\\pm'+fmt%A1+'$'

def get_covmat(Xs, weis=[]):
        n = len(Xs)
        Xibars = [0 for row in range(n)]
        XiXjbars = [ [0 for row1 in range(n)] for row2 in range(n)]
        covmat = [ [0 for row1 in range(n)] for row2 in range(n)]
        ndat = len(Xs[0])
        if weis == []: weis = [1. for row in range(ndat)]
        sumwei = sum(weis)
        sumweisq = sum([wei**2.0 for wei in weis])
        for i in range(n):
                ## \bar Xi
                for idat in range(ndat):
                        Xibars[i] += Xs[i][idat]*weis[idat]
                Xibars[i] /= sumwei
                ## \bar (XiXj)
                for j in range(i,n):
                        for idat in range(ndat):
                                XiXjbars[i][j] += (Xs[i][idat] * Xs[j][idat]) * weis[idat]
                        XiXjbars[i][j] /= sumwei
        for i in range(n):
                for j in range(i,n):
                        XiXjbars[j][i] = XiXjbars[i][j]
        #print 'XiXjbars: '
        #for i in range(n):
        #       print XiXjbars[i]
        for i in range(n):
                for j in range(n):
                        covmat[i][j] = (XiXjbars[i][j] - Xibars[i]*Xibars[j]) * ndat / float(ndat-1.0)
#                               = E ((Xi -Xibar)(Xj-Xjbar) ) = E(XiXj) - E(Xi)*Xjbar - E(Xj) * Xibar + Xibar * Xjbar
        return covmat

