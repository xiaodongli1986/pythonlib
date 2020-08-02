import numpy as np
import matplotlib.pyplot as plt
from camb import model,initialpower
import camb,os,sys



maxkh = 4
npoints = 600
path = 'data/'



def linPk(Omega_m=0.307115,w=-1):
    print('Start producing power spectrum')
    h = 0.6777
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h*100, ombh2=0.048206*h*h, omch2=(Omega_m-0.048206)*h*h)
    pars.set_dark_energy(w=w)
    def checks8(As=2e-9):
        pars.InitPower.set_params(ns=0.96,As=As)
        pars.set_matter_power(redshifts=[0], kmax=100.,)
        pars.NonLinear = model.NonLinear_none
        results = camb.get_results(pars)
        s8 = results.get_sigma8()
        return(s8)
    As = 2e-9;Astmp=2.01e-9
    s8 = checks8(As);tmp=checks8(Astmp)
    target = 0.8228
    for i in range(10):
        if(abs(s8-target)<1e-5):break
        diff = (tmp-s8)/(Astmp-As)
        Astmp = As;tmp = s8
        As = As - (s8-target)/diff
        s8 = checks8(As)
        print('Iteration of Sigma8 ',s8)
    pars.InitPower.set_params(ns=0.96,As=As)
    pars.set_matter_power(redshifts=[0], kmax=100.,)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=maxkh, npoints = npoints,)
    return(As,s8,kh,pk)


def Pksave(Om,w,As,s8,kh,pk):
    As,s8 = str(round(As[0],13)),str(round(s8[0],6))
    filename = path+'Om'+str(Om)+'w'+str(w)+'As'+As+'s8'+s8+'.dat'
    data = np.vstack([kh,pk[0]]).T
    np.savetxt(filename,data)
    return 0 
    

omList = 0.2671+0.02*np.arange(5)
wList = np.array([-1.4,-1.2,-0.8,-0.6])
#omList,wList





As,s8,kh,pk = linPk(omList[1],wList[1])
Pksave(omList[1],wList[1],As,s8,kh,pk)




As,s8,kh,pk = linPk(omList[1],wList[2])
Pksave(omList[1],wList[2],As,s8,kh,pk)





As,s8,kh,pk = linPk(omList[3],wList[1])
Pksave(omList[3],wList[1],As,s8,kh,pk)





As,s8,kh,pk = linPk(omList[3],wList[2])
Pksave(omList[3],wList[2],As,s8,kh,pk)





for Omega in omList:
    As,s8,kh,pk = linPk(Omega,-1)
    Pksave(Omega,-1,As,s8,kh,pk)





for w in wList:
    As,s8,kh,pk = linPk(0.3071,w)
    Pksave(0.3071,w,As,s8,kh,pk)







