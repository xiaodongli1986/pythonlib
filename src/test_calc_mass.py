import numpy as np
Omega = 0.307115
HubbleParam = 0.692885
boxsize = 1024
nparticle = 512
Omega, HubbleParam = 0.21, 0.6786682


### particle mass
G = 6.67428*(10**-8) #cm3*g-1*s-2
M_sun = 1.98892*(10**33) #g
Mpc = 30.85677*(10**23) #cm
h = float(HubbleParam); Omega_m = float(Omega); boxsize=float(boxsize); nparticle=float(nparticle)

rho = (3*((100*h)**2)/(8*np.pi*G)) * ((10**10)*Mpc/M_sun) #M_sun*h2/Mpc-3
rho_h = rho/(h**2)
m_p1 = (((boxsize**3)*3*Omega_m*((100)**2))/(nparticle**3*8*np.pi*G)) * ((10**10) *Mpc /M_sun) #M_sun/h
m_p2 = (rho_h*Omega_m*(boxsize**3))/(nparticle**3)
print ('m_p=',m_p2)

