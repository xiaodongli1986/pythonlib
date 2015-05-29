def lcdm_qz(omegam, z=0):
    Hz = np.sqrt(omegam*(1.0+z)**3.0 + 1-omegam)
    dHdz = 3*omegam*(1.0+z)**2.0 / (2*Hz)
    q = dHdz / Hz * (1.0+z) - 1.0
    return q

def lcdm_dqdz(omegam, z=0, deltaz=0.01):
    dqdz = (lcdm_qz(omegam, z+deltaz) - lcdm_qz(omegam, z-deltaz)) / (2*deltaz)
    return dqdz
