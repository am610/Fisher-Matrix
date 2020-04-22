import pandas as pd
import numpy as np 
np.set_printoptions(precision=4)

P=  pd.read_csv('Fisher_Table.csv')
dmdom = np.array(P['dmdom'])
dmdwo = np.array(P['dmdwo'])
dmdwa = np.array(P['dmdwa'])
dmdM = np.array(P['dmdM'])
z = np.array(P['z'])
dmdtheta = np.array(zip(dmdom,dmdwa,dmdwo,dmdM))
N_parms = np.shape(dmdtheta)[1]
N=np.array([300.0, 35.0, 64.0, 95.0, 124.0, 150.0, 171.0, 183.0, 179.0, 170.0, 155.0, 142.0, 130.0, 119.0,107.0, 94.0, 80.0]) 
sigma_stat_sq = 0.12*0.12
sigma_sys_sq  = (0.02*(1+z)/2.7)*(0.02*(1+z)/2.7)
sigma_sq      = sigma_stat_sq/N + sigma_sys_sq
F = np.zeros([N_parms,N_parms])

#dmdM_dmdM = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdM,dmdM),(1/sigma_sq))),name='dmdM_dmdM'))
#dmdom_dmdom = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdom,dmdom),(1/sigma_sq))),name='dmdomegam_dmdomegam'))
#dmdwa_dmdwa =pd.DataFrame(pd.Series (np.sum(np.multiply(np.multiply(dmdwa,dmdwa),(1/sigma_sq))),name='dmdwa_dmdwa'))
#dmdwo_dmdwo = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdwo,dmdwo),(1/sigma_sq))),name='dmdwo_dmwo'))
#dmdM_dmdom = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdM,dmdom),(1/sigma_sq))),name='dmdM_dmdom'))
#dmdM_dmdwo = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdM,dmdwo),(1/sigma_sq))),name='dmdM_dmdwo'))
#dmdM_dmdwa = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdM,dmdwa),(1/sigma_sq))),name='dmdM_dmdwa'))
#dmdom_dmdwa = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdom,dmdwa),(1/sigma_sq))),name='dmdom_dmdwa'))
#dmdom_dmdwo = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdom,dmdwo),(1/sigma_sq))),name='dmdom_dmdwo'))
#dmdwa_dmdwo = pd.DataFrame(pd.Series(np.sum(np.multiply(np.multiply(dmdwa,dmdwo),(1/sigma_sq))),name='dmdwa_dmdwo'))




dmdM_dmdM = (np.sum(np.multiply(np.multiply(dmdM,dmdM),(1/sigma_sq)))) # since its scalar product, so order of the terms doesn't matter
dmdom_dmdom =  (np.sum(np.multiply(np.multiply(dmdom,dmdom),(1/sigma_sq))))
dmdwa_dmdwa =  (np.sum(np.multiply(np.multiply(dmdwa,dmdwa),(1/sigma_sq))))
dmdwo_dmdwo =  (np.sum(np.multiply(np.multiply(dmdwo,dmdwo),(1/sigma_sq))))
dmdM_dmdom =  (np.sum(np.multiply(np.multiply(dmdM,dmdom),(1/sigma_sq))))
dmdM_dmdwo =  (np.sum(np.multiply(np.multiply(dmdM,dmdwo),(1/sigma_sq))))
dmdM_dmdwa =  (np.sum(np.multiply(np.multiply(dmdM,dmdwa),(1/sigma_sq))))
dmdom_dmdwa =  (np.sum(np.multiply(np.multiply(dmdom,dmdwa),(1/sigma_sq))))
dmdom_dmdwo =  (np.sum(np.multiply(np.multiply(dmdom,dmdwo),(1/sigma_sq))))
dmdwa_dmdwo =  (np.sum(np.multiply(np.multiply(dmdwa,dmdwo),(1/sigma_sq))))             

#dtheta_dtheta = pd.DataFrame(dmdM_dmdM,dmdomegam_dmdomegam)

F_I = np.array([dmdM_dmdM,dmdM_dmdom,dmdM_dmdwo,dmdM_dmdwa,dmdM_dmdom,dmdom_dmdom,dmdom_dmdwo,dmdom_dmdwa,dmdM_dmdwo,dmdom_dmdwo,dmdwo_dmdwo,dmdwa_dmdwo,dmdM_dmdwa,dmdom_dmdwa,dmdwa_dmdwo,dmdwa_dmdwa])

Fisher_Matrix = F_I.reshape(4,4)
Cov_Matrix    = np.linalg.pinv(Fisher_Matrix)
C_ij  = np.sqrt(np.diag(Cov_Matrix))
#print np.matrix(F_I)
print '\n\nFisher Matrix \n',Fisher_Matrix
print '\n\n*****************************\n\n'
print 'Covariance Matrix \n',Cov_Matrix
print '\n\n*****************************\n\n'
print 'sigma_ii \n', C_ij
#***************************
# F_Ext
# Adding a prior to Omega_matter
#F_[dMdom_dMdom] = [1,1]
# PRIOR OF 0.03
F_Ext = 1/(0.03)**2
v = Fisher_Matrix[1,1] + F_Ext
Fisher_Matrix_PRIOR = Fisher_Matrix
Fisher_Matrix_PRIOR[1,1] =v
Cov_Matrix_PRIOR    = np.linalg.pinv(Fisher_Matrix_PRIOR)
C_ij_PRIOR  = np.sqrt(np.abs(np.diag(Cov_Matrix_PRIOR)))
#print np.matrix(F_I)
print '\n\n The ORDER IS : dM_dM, domega_m_domega_m, dw0_dw0, dwa_dwa]\n\n'
print '\n\nFisher Matrix with prior \n',Fisher_Matrix_PRIOR
print '\n\n*****************************\n\n'
print 'Covariance Matrix with prior\n',Cov_Matrix_PRIOR
print '\n\n*****************************\n\n'
print 'sigma_ii with prior\n', C_ij_PRIOR                               

def fisher_matrix():
    return Fisher_Matrix_PRIOR

def covariance_matrix():
    return cov_Matrix_PRIOR
