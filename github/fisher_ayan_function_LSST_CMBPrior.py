import pandas as pd
import numpy as np
np.set_printoptions(precision=4)

# Produced by dmdtheta_ayan_function.py
P = pd.read_csv('Fisher_Table_LSST.csv')
dmdom = np.array(P['dmdom'])
dmdwo = np.array(P['dmdwo'])
dmdwa = np.array(P['dmdwa'])
dmdM = np.array(P['dmdM'])
z = np.array(P['z'])
dmdtheta = np.array(zip(dmdom, dmdwa, dmdwo, dmdM))
N_parms = np.shape(dmdtheta)[1]

# ****LSST*****
s = 300*np.ones(10)
s2 = np.array([150, 150])
s2 = np.ndarray.tolist(s2)
s = np.ndarray.tolist(s)
N = np.array(s+s2)
sigma_int_sq = 0.15*0.15
sigma_sys_sq = (0.02*(1+z)/2.)*(0.02*(1+z)/2.)
sigma_sq = sigma_int_sq/N + sigma_sys_sq
# ****LSST*****

F = np.zeros([N_parms, N_parms])
# since its scalar product, so order of the terms doesn't matter
dmdM_dmdM = (np.sum(np.multiply(np.multiply(dmdM, dmdM), (1/sigma_sq))))
dmdom_dmdom = (np.sum(np.multiply(np.multiply(dmdom, dmdom), (1/sigma_sq))))
dmdwa_dmdwa = (np.sum(np.multiply(np.multiply(dmdwa, dmdwa), (1/sigma_sq))))
dmdwo_dmdwo = (np.sum(np.multiply(np.multiply(dmdwo, dmdwo), (1/sigma_sq))))
dmdM_dmdom = (np.sum(np.multiply(np.multiply(dmdM, dmdom), (1/sigma_sq))))
dmdM_dmdwo = (np.sum(np.multiply(np.multiply(dmdM, dmdwo), (1/sigma_sq))))
dmdM_dmdwa = (np.sum(np.multiply(np.multiply(dmdM, dmdwa), (1/sigma_sq))))
dmdom_dmdwa = (np.sum(np.multiply(np.multiply(dmdom, dmdwa), (1/sigma_sq))))
dmdom_dmdwo = (np.sum(np.multiply(np.multiply(dmdom, dmdwo), (1/sigma_sq))))
dmdwa_dmdwo = (np.sum(np.multiply(np.multiply(dmdwa, dmdwo), (1/sigma_sq))))

#dtheta_dtheta = pd.DataFrame(dmdM_dmdM,dmdomegam_dmdomegam)

F_I = np.array([dmdM_dmdM, dmdM_dmdom, dmdM_dmdwo, dmdM_dmdwa, dmdM_dmdom, dmdom_dmdom, dmdom_dmdwo, dmdom_dmdwa,
                dmdM_dmdwo, dmdom_dmdwo, dmdwo_dmdwo, dmdwa_dmdwo, dmdM_dmdwa, dmdom_dmdwa, dmdwa_dmdwo, dmdwa_dmdwa])

Fisher_Matrix = F_I.reshape(4, 4)
Cov_Matrix = np.linalg.pinv(Fisher_Matrix)
C_ij = np.sqrt(np.diag(Cov_Matrix))
# Original prints
#print np.matrix(F_I)
#print '\n\nFisher Matrix \n',Fisher_Matrix
#print '\n\n*****************************\n\n'
#print 'Covariance Matrix \n',Cov_Matrix
#print '\n\n*****************************\n\n'
#print 'sigma_ii \n', C_ij

# ***************************
# CMB PRIOR
F_OmOm = 26282.777992634321
F_Omw0 = -7211.5106699248900
F_Omwa = -1974.3936330616734
F_w0w0 = 1978.7058337978979
F_w0wa = 541.73728345787083
F_wawa = 148.31880478414203
F1 = [F_OmOm, F_Omw0, F_Omwa, F_Omw0, F_w0w0, F_w0wa, F_Omwa, F_w0wa, F_wawa]
F1 = np.reshape(F1, (3, 3))
A1 = np.zeros(16)
A1 = np.reshape(A1, (4, 4))
for row in range(1, np.shape(A1)[0]):
    for col in range(1, np.shape(A1)[1]):
        A1[row][col] += F1[row-1][col-1]
F_Ext = A1


Fisher_Matrix_PRIOR = Fisher_Matrix + F_Ext
Cov_Matrix_PRIOR = np.linalg.pinv(Fisher_Matrix_PRIOR)
C_ij_PRIOR = np.sqrt(np.abs(np.diag(Cov_Matrix_PRIOR)))
#print np.matrix(F_I)

print '\n\n The ORDER IS : dM_dM, domega_m_domega_m, dw0_dw0, dwa_dwa]\n\n'
print '\n\nFisher Matrix with prior \n', Fisher_Matrix_PRIOR
#print '\n\n*****************************\n\n'
print 'Covariance Matrix with prior\n', Cov_Matrix_PRIOR
#print '\n\n*****************************\n\n'
print 'sigma_ii with prior\n', C_ij_PRIOR


def fisher_matrix():
    return Fisher_Matrix_PRIOR


def covariance_matrix():
    return Cov_Matrix_PRIOR
