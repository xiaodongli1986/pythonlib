def Hz_CPL(omegam, w0, wa, h, z):
    return 100*h*np.sqrt(omegam*(1.0+z)**3.0 + (1.0-omegam)*(1.0+z)**(3.0*(1.0+w0+wa)) * exp(-3.0*z*wa/(1.0+z)))

def comov_r_CPL(omegam, w0,wa, h, z):
     x, y = scipy.integrate.quad(lambda x: 1.0/Hz_CPL(omegam, w0, wa, h, x), 0, z)
     return CONST_C* x * h 

def lcdm_qz(omegam, z=0):
    Hz = np.sqrt(omegam*(1.0+z)**3.0 + 1-omegam)
    dHdz = 3*omegam*(1.0+z)**2.0 / (2*Hz)
    q = dHdz / Hz * (1.0+z) - 1.0
    return q

def lcdm_dqdz(omegam, z=0, deltaz=0.01):
    dqdz = (lcdm_qz(omegam, z+deltaz) - lcdm_qz(omegam, z-deltaz)) / (2*deltaz)
    return dqdz
