import numpy as np
import pandas as pd
from scipy.integrate import quad

# **************************
# This program returns the delta_m_k value
# **************************
# (5/ln10)*[delz(1/H*(1/int_0^z'delz/H ) + 1/(1+z))] !!NOT USED HERE!!
# w_a = 0, w_o = -1 , omega_m = 0.3
# H(z) = sqrt(omega_m*(1+z)^3 + (1-omega_m))

# D_L = (1+z) integration(dz/H(z))LIMIT:[0,z]

const = 5/np.log(10)


def H_z(z):
    omega_m = 0.3
    omega_de = 1-omega_m
    hz = omega_m*(1+z)**3+omega_de
    return np.sqrt(hz)


def H_z_inverse(z):
    return 1/(H_z(z))

#d0 =0.   ; d1 =0. ;
#d0 =0.01 ; d1 =0. ;
#d0 =0.   ; d1 =0.01 ;
#print 'd0 = %s, d1 = %s '%(d0,d1)


def del_z(d0, d1, z):
    return (d0+d1*z)


def integration_delz_h(z):
    I = quad(H_z_inverse, 0, z)[0]
    return I


def D_l(z):
    return (1+z)*integration_delz_h(z)


def del_m(z, d0, d1):
    x = const*np.log(D_l(z+del_z(d0, d1, z))/D_l(z))
    return x


for i in np.linspace(0.05, 1.15, 12):
    del_m(i, 0., 0.01)
