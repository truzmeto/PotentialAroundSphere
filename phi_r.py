import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#currently not needed
def deltaF(r0, r):
    if r*r == r0*r0:
        return 1.0
    else:
        return  0.0


def solve_phi(y, r, eps_in, eps_out, sigma, sigma_s, k0, q):

    rsq = r*r 
    unit_factor = 6987.0 #this factor ensures that potential

    eps_r = eps_out + (eps_in - eps_out)*np.exp(-rsq/sigma**2)
    de_dr = -2.0*(r)*(eps_in - eps_out)*np.exp(-rsq/sigma**2)
    lambda_r = 1.0 - np.exp(-rsq / (sigma + sigma_s)**2)
    #delta_r = deltaF(r0, r)

    phi, dphi_dr = y
    phi_new = [dphi_dr, #-(4*np.pi*unit_factor *q/(eps_r))*delta_r
               -k0*k0*lambda_r*phi/eps_r -
               dphi_dr*(de_dr/eps_r + 2/r)]
            
    return phi_new


unit_factor = 6987.0 #this factor ensures that potential
                     #unit is in (k_B*T / e), where T = 300 K room temp.


sigma = 10.0                          #atom size AA
q = 10.0                              #1.6*10^(-19)

#initilize
r = np.linspace(sigma, sigma+20, 1001)      #r range
p0 = unit_factor*q/(4* np.pi*sigma)   #initial value for phi in (kB T/e)
phi0 = [p0, 0.0]                      #initial values

eps_in = 4.0
eps_out = 79.0
sigma_s = 1.0                         #stern layer width
k0 = 0.5                              # 

sol = odeint(solve_phi, phi0, r, args=(eps_in, eps_out, sigma, sigma_s, k0, q))

plt.plot(r, sol[:, 0], 'b', label=r'$\phi(r)$')
plt.legend(loc='best')
plt.xlabel('r')
plt.ylabel(r'$\phi$')
plt.grid()
plt.show()
