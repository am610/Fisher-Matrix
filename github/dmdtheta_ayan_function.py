import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
import pandas as pd

Om = 0.3  # Omega Matter
Omatter = (1-Om)
w0 = -1
wa = 0
#w0 = -1.2;wa = -0.9;
#w0 = -0.8;wa = 0.9;
#a = (1+z)**(-1)
# 1
#print('w0, wa used : ', w0, wa)
# Not Required


def X(z):
    a = (1+z)**(-1)
    return Om*(a**-3)+((1-Om)*a**(-3)*(a**(-3*(w0+wa))*np.exp(3*wa*(a-1))))


def dxdm(z, wa, w0):
    a = (1+z)**(-1)
    zt = (1+z)
    w = (wa+w0)
    X3 = np.exp(3*wa*(a-1))
    X2 = zt**(3*w)
    X1 = X2*X3
    X4 = 1.-X1
    X = Om*(a**-3)+((1-Om)*a**(-3)*(a**(-3*(w0+wa))*np.exp(3*wa*(a-1))))
    A = (X**(-1.5))*X4*(zt**3)
    #A  = [1-zt**(-3*w)*np.exp(-3*wa*z/zt)]*zt**3
    return A

# 2


def dxdw0(z, wa, w0):
    a = (1+z)**(-1)
    zt = (1+z)
    ztt = z/zt  # z/1+z
    a_a = a**(3*wa)
    a_0 = a**(3*w0)
    e_a = np.exp(3*wa*(-ztt))
    A3 = (a_a**-1)*(-3*np.log(a)/a_0)
    A4 = e_a*A3
    X = Om*(a**-3)+((1-Om)*a**(-3)*(a**(-3*(w0+wa))*np.exp(3*wa*(a-1))))
    A = (X**(-1.5))*Omatter*zt**3*(A4)
    return A


# 3
def dxdwa(z, wa, w0):
    a = (1+z)**(-1)
    a2 = a-1
    zt = (1+z)
    ztt = z/zt
    a_a = a**(3*wa)
    a_0 = a**(3*w0)
    e_a = np.exp(3*wa*(-ztt))
    A3 = e_a*(-3*np.log(a)/a_a)
    A2 = 3*(a-1)*e_a/(a_a)
    A4 = A2+A3
    X = Om*(a**-3)+((1-Om)*a**(-3)*(a**(-3*(w0+wa))*np.exp(3*wa*(a-1))))
    A = (X**(-1.5))*(Omatter)*(zt**3)*(a_0**-1)*(A4)
    return A


# PRINTING CHECKS
#print 'Redshift  dx/dm      dx/dwo      dx/dwa'
#print '---------------------------------------'
# for i in (np.linspace(0.05,1.65,17)):
#          print("%.4f    %.4f    %.4f     %.4f " %(i,dxdm(i,0,-1),dxdw0(i,0,-1),dxdwa(i,0,-1)) )
# ************************************************************
# ************************************************************


# Integration Part : STARTS FROmatter BELOW

def H_0D_L(z):  # THIS WILL BE THEN INTEGRATED IN LINE NO. 106
    a = (1+z)**(-1)
    sec = (1-Om)*(a**-3)*(a**(-3*(w0+wa))*np.exp(3*(wa*(a-1))))
    fir = (Om*(a**-3))
    return((fir+sec)**(-0.5))


# Z BINS
# z = np.linspace(0.05,1.65,17) # bin means
# z = np.linspace(0.05,1.15,12) # LSST
z = np.array([0.3, 0.6, 1.0, 1.4, 1.8, 2.0, 2.6, 3.0, 3.3, 4.0])


# ADDITIONAL FUNCTION DECLARATION
def fu1(z):
    return dxdm(z, wa, w0)


def fu2(z):
    return dxdw0(z, wa, w0)


def fu3(z):
    return dxdwa(z, wa, w0)


# VARIABLE DECLARATION
I_om = []
I_wa = []
I_wo = []
hdl = []
prefactor1 = []
dx_domega = []
dx_dwa = []
dx_dw0 = []

# INTEGRATIONS RUNNING FROM LIMIT [0,Z_BIN[i]]
# THE fu FUNCTION CALLS THE OTHER FUNCTION : DX/DTHETA
# QUAD IS A PYTHON INTEGRATION FUNCTION WITH 3 ARGUMENTS :
# FUNCTION TO INTEGRATE, LOWER LIMIT, UPPER LIMIT

for ii in range(len(z[:])):
    #print quad(fu, 0,0.1)
    #       print 'z bin ',z_bin[ii]
    #       print 'z mean',z[ii]
    #       print 'fu1[z]',fu1(z[ii])
    #       print 'quad(fu1,0,%s)'%z[ii]
    #       print 'I_Om :',(quad(fu1,0,z[ii]))
    I_om.append(quad(fu1, 0, z[ii])[0])
    I_wo.append(quad(fu2, 0, z[ii])[0])
    I_wa.append(quad(fu3, 0, z[ii])[0])
    hdl.append((1+z[ii])*quad(H_0D_L, 0, z[ii])[0])  # H0_Dl INTEGRATION
    # PREFACTOR WITHOUT THE H0_Dl TERM. THE 1* STANDS FOR 'c'
    prefactor1.append((-5/np.log(10))*(1*(1+z[ii]))/2)
#       print 'PREFACTOR :',prefactor1[ii]
#print '*******************************************************'

# STORING THE VALUES IN ARRAYS
for i in range(len(I_om)):
    dx_domega.append((prefactor1[i]/hdl[i])*I_om[i])
    dx_dw0.append((prefactor1[i]/hdl[i])*I_wo[i])
    dx_dwa.append((prefactor1[i]/hdl[i])*I_wa[i])
    #print 'DX/Domega_m :',I_om[i]
    dx_dM = np.ones(17)
# WRITING TO A FILE
DXDT = pd.DataFrame(zip(z, dx_domega, dx_dw0, dx_dwa, dx_dM))
DXDT.columns = ['z', 'dmdom', 'dmdwo', 'dmdwa', 'dmdM']
DXDT.to_csv('Fisher_Table_SS.csv', index=False)


#print DXDT

# PRINT OUTPUT TO SCREEN
# print 'Redshift H0_Dl  dm/domegam      dm/dwo      dm/dwa'#    (3*dm/dwa)/(dm/dwo)'
#print '-------------------------------------------------------------------------'
# for i in range(len(I_om)):
#          print("%.4f    %.4f    %.4f    %.4f     %.4f " %(z[i],hdl[i],dx_domega[i],dx_dw0[i],dx_dwa[i]))#,3*dx_dwa[i]/dx_dw0[i]))

def func():
    return DXDT
