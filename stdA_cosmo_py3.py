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


class zerror:
	def __init__(self, sig=0.002, z_slope=0.):
		''' 
	redshift error, based on model of 
		z_error = sig*(1+z) * (z+z_slope*z)
	sig  is the base scattering; z_slope describes linear variation of z_error '''
		self.sig = sig; self.z_slope=z_slope
	def z_error(self, z):
		return self.sig*(1+z)*(1+self.z_slope*z)
	def vLOS_effective(self, z):
		''' dz = vLOS * (1+z) / c '''
		return self.z_error(z) * 3e5 / (1.+z)
	def r_error(self, z, om=0.3071, w=-1):
		dz = self.z_error(z)
		return comov_r( om, w, 0.7, z+dz) - comov_r( om, w, 0.7, z) 
