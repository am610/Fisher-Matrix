#!/usr/bin/python
import del_m_fisher_bias_full as m  # CHECK WHICH D0,D1 IS USED
import numpy as np
#import fisher_ayan_function as f
import dmdtheta_ayan_function as dm
import fisher_ayan_function_LSST as f
#import fisher_ayan_function_LSST_CMBPrior as f
import pandas as pd
import sys
#print sys.argv[1]
# apply = sys.argv[2] # which bin to apply with dz and omit the rest
z = np.linspace(0.05, 1.15, 12)
dm_dp = dm.func()
sum = []
#d0 = [0.,0.01,0.]
#d1 = [0.,0.,0.01]
do_d1 = [[0., 0.], [0.01, 0.], [0., 0.01]]
cs1 = []
cs2 = []
cs3 = []

#print 'Case 1  : d0, d1  %s'%(do_d1[0])
for i in (z):
    #print i
    x1 = m.del_m(i, do_d1[0][0], do_d1[0][1])  # + 1.4/(1+i)*(m.del_z(z,d0,d1))
#    print x1
    cs1.append(x1)
#print 'Case 2  : d0, d1  %s'%(do_d1[1])
for i in (z):
    x2 = m.del_m(i, do_d1[1][0], do_d1[1][1])
#    print x2
    cs2.append(x2)
#print 'Case 3  : d0, d1  %s'%(do_d1[2])
for i in (z):
    x3 = m.del_m(i, do_d1[2][0], do_d1[2][1])
#    print x3
    cs3.append(x3)
print '***********************\nFull expressions : No Approxiamtion'
print '************************\n\ndel_m with 3 different sets of d0,d1\n'

print "0,0       0.01,0    0,0.01\n"
for i in range(len(cs2)):
    print("%.4f    %.4f    %.4f  ") % (cs1[i], cs2[i], cs3[i])
print '***********************'
s = 300*np.ones(10)
s2 = np.array([150, 150])
s2 = np.ndarray.tolist(s2)
s = np.ndarray.tolist(s)
N = np.array(s+s2)
sigma_int_sq = 0.15*0.15
sigma_sys_sq = (0.02*(1+z)/2.)*(0.02*(1+z)/2.)
sigma_sq = sigma_int_sq/N + sigma_sys_sq

print 'N :  ', N  # OK
print 'sigma_square : ', sigma_sq
if sys.argv[1] == '1':
    flag = cs1
elif sys.argv[1] == '2':
    flag = cs2
elif sys.argv[1] == '3':
    flag = cs3
# flag = cs3 # **** FLAG SET
del_m_sigma_sq = flag/sigma_sq  # cs1 is used
cov = f.covariance_matrix()
XX_om = 0.
XX_wa = 0.
XX_w0 = 0.0
XX_M = 0.
for n in range(len(dm_dp)):
    XX_om += dm_dp['dmdom'][n]*del_m_sigma_sq[n]
    XX_w0 += dm_dp['dmdwo'][n]*del_m_sigma_sq[n]
    XX_wa += dm_dp['dmdwa'][n]*del_m_sigma_sq[n]
    XX_M += dm_dp['dmdM'][n] * del_m_sigma_sq[n]
arr = [XX_M, XX_om, XX_w0, XX_wa]
if flag == cs1:
    print '\n********\nUsing d0 , d1 :0.0, 0\n********\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[0][0], do_d1[0][1])
    ide = str(do_d1[0][0])+'_'+str(do_d1[0][1])
elif flag == cs2:
    print '\n********\nUsing d0 , d1 :0.01, 0\n********\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[1][0], do_d1[1][1])
    ide = str(do_d1[1][0])+'_'+str(do_d1[1][1])
elif flag == cs3:
    print '\n********\nUsing d0 , d1 :0., 0.01\n********\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[2][0], do_d1[2][1])
    ide = str(do_d1[2][0])+'_'+str(do_d1[2][1])

print '\ncov : \n', cov
print '\ndel_m / sigma_sq\n', del_m_sigma_sq
print 'array of 4 vectors\n %s' % arr
del_p2 = cov.dot(arr)
print '\ndel_p2 : \n', (del_p2)
name = 'del_p_LSST_%s.txt' % ide
np.savetxt(name, np.matrix(del_p2), fmt='%.4f')
