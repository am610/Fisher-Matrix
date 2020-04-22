#!/usr/bin/python

# ADDITION OF STRETCH FACTOR IN BIAS TERM. LINE 26 ONWARDS
# TWO ARGUMENTS PROVIDED :
# 1 = OPTIONS : [1-7]
# WHICH BINS TO APPLY THE del_z AND WHICH COMBINATION
# OF d0,d1 TO APPLY
# THE FIRST 3 BASES APPLY D0,D1 = 0,(0.01,0),(0,0.01),
# UNIFORMLY ON ALL THE 12 z BINS
# THE NEXT 4 CASES IN ARGV[1] SELECTIVELY APPLIES THE TWO
# NON-ZERO d0,d1 CASES TO THE 12 BINS

# THE SECOND ARGUMENT : OPTIONS [1-12]
# WHICH BIN TO SELECT OUT OF 12 TO APPLY
# del_z OR WHICH BIN TO EXCLUDE FOR APPLYING
# THE del_z


import del_m_fisher_bias_full as m  # CHECK WHICH D0,D1 IS USED
import numpy as np
#import fisher_ayan_function as f
import dmdtheta_ayan_function as dm
#import fisher_ayan_function_LSST as f
import fisher_ayan_function_LSST_CMBPrior as f
import pandas as pd
import sys


#print sys.argv[1]
apply = int(sys.argv[2])-1  # which bin to apply with dz and omit the rest
print '****** z - Bin to work on : ', apply+1
counter = 0
z = np.linspace(0.05, 1.15, 12)
dm_dp = dm.func()
sum = []
#d0 = [0.,0.01,0.]
#d1 = [0.,0.,0.01]
do_d1 = [[0., 0.], [0.01, 0.], [0., 0.01]]
cs1 = []
cs2 = []
cs3 = []
cs4 = []
cs5 = []
cs6 = []
cs7 = []
cols = ['d0', 'd1', 'bin_used', 'how']
dat1 = pd.DataFrame(columns=cols)
dat2 = pd.DataFrame(columns=cols)
dat3 = pd.DataFrame(columns=cols)
dat4 = pd.DataFrame(columns=cols)
dat5 = pd.DataFrame(columns=cols)
dat6 = pd.DataFrame(columns=cols)
dat7 = pd.DataFrame(columns=cols)


#print 'Case 1  : d0, d1  %s'%(do_d1[0])
for i in (z):
    #print i
    x1 = m.del_m(i, do_d1[0][0], do_d1[0][1]) + \
        (1.4/(1+i)*(m.del_z(do_d1[0][0], do_d1[0][1], i)))
#    print x1
    cs1.append(x1)
dat1 = dat1.append({'d0': str(do_d1[2][0]), 'd1': str(
    do_d1[2][1]), 'bin_used': str(0), 'how': str(0)}, ignore_index=True)
#print 'Case 2  : d0, d1  %s'%(do_d1[1])
for i in (z):
    x2 = m.del_m(i, do_d1[1][0], do_d1[1][1]) + \
        (1.4/(1+i)*(m.del_z(do_d1[1][0], do_d1[1][1], i)))
   # print'Stretch factor : ', (1.4/(1+i)*(m.del_z(do_d1[1][0],do_d1[1][1],i)))
    cs2.append(x2)
dat2 = dat2.append({'d0': str(do_d1[1][0]), 'd1': str(
    do_d1[1][1]), 'bin_used': str(0), 'how': str(0)}, ignore_index=True)
#print 'Case 3  : d0, d1  %s'%(do_d1[2])
for i in (z):
    x3 = m.del_m(i, do_d1[2][0], do_d1[2][1]) + \
        (1.4/(1+i)*(m.del_z(do_d1[2][0], do_d1[2][1], i)))
    #print'Stretch factor : ', (1.4/(1+i)*(m.del_z(do_d1[2][0],do_d1[2][1],i)))
    #print'Difference : ',x3-(1.4/(1+i)*(m.del_z(do_d1[2][0],do_d1[2][1],i)))
    cs3.append(x3)
dat3 = dat3.append({'d0': str(do_d1[2][0]), 'd1': str(
    do_d1[2][1]), 'bin_used': str(0), 'how': str(0)}, ignore_index=True)
for i in (z):  # APPLY ONLY ON ONE BIN
    if apply == counter:
        #        print 'Condition match'
        #        print 'counter :',counter
        #        print 'apply :',apply
        x4 = m.del_m(i, do_d1[1][0], do_d1[1][1]) + (1.4/(1+i)
                                                     * (m.del_z(do_d1[1][0], do_d1[1][1], i)))  # 0.01,0
    else:
        x4 = m.del_m(i, do_d1[0][0], do_d1[0][1])
    counter += 1
    cs4.append(x4)
dat4 = dat4.append({'d0': str(do_d1[1][0]), 'd1': str(
    do_d1[1][1]), 'bin_used': str(apply+1), 'how': str(1)}, ignore_index=True)
counter = 0
for i in (z):  # APPLY ONLY ON ONE  BIN
    if apply == counter:
        #        print 'Condition match'
        #        print 'counter :',counter
        #        print 'apply :',apply
        x5 = m.del_m(i, do_d1[2][0], do_d1[2][1]) + (1.4/(1+i)
                                                     * (m.del_z(do_d1[2][0], do_d1[2][1], i)))  # 0.0,0.01
    else:
        x5 = m.del_m(i, do_d1[0][0], do_d1[0][1])
    counter += 1
    cs5.append(x5)
dat5 = dat5.append({'d0': str(do_d1[2][0]), 'd1': str(
    do_d1[2][1]), 'bin_used': str(apply+1), 'how': str(1)}, ignore_index=True)
counter = 0
for i in (z):  # APPLY ON ALL BIN EXCEPT ONE
    #    print 'All BIn Except One : 6'
    if apply != counter:
        x6 = m.del_m(i, do_d1[1][0], do_d1[1][1]) + (1.4/(1+i)
                                                     * (m.del_z(do_d1[1][0], do_d1[1][1], i)))  # 0.01,0
    else:
        x6 = m.del_m(i, do_d1[0][0], do_d1[0][1])
    counter += 1
    cs6.append(x6)
dat6 = dat6.append({'d0': str(do_d1[1][0]), 'd1': str(
    do_d1[1][1]), 'bin_used': str(apply+1), 'how': str(11)}, ignore_index=True)
counter = 0
for i in (z):  # APPLY ON ALL BIN EXCEPT ONE
    #    print 'All BIn Except One : 7'
    if apply != counter:
        x7 = m.del_m(i, do_d1[2][0], do_d1[2][1]) + (1.4/(1+i)
                                                     * (m.del_z(do_d1[2][0], do_d1[2][1], i)))  # 0.0,0.01
    else:
        x7 = m.del_m(i, do_d1[0][0], do_d1[0][1])
    counter += 1
    cs7.append(x7)
dat7 = dat7.append({'d0': (do_d1[2][0]), 'd1': str(
    do_d1[2][1]), 'bin_used': str(apply+1), 'how': str(11)}, ignore_index=True)


#print '***********************\nFull expressions : No Approxiamtion'
#print '************************\n\ndel_m with 3 different sets of d0,d1\n'

#print "0,0       0.01,0    0,0.01\n"
# for i in range(len(cs2)):
#    print("%.4f    %.4f    %.4f  ")%(cs1[i],cs2[i],cs3[i])
#print '***********************'
s = 300*np.ones(10)
s2 = np.array([150, 150])
s2 = np.ndarray.tolist(s2)
s = np.ndarray.tolist(s)
N = np.array(s+s2)
sigma_int_sq = 0.15*0.15
sigma_sys_sq = (0.02*(1+z)/2.)*(0.02*(1+z)/2.)
sigma_sq = sigma_int_sq/N + sigma_sys_sq

# print 'N :  ',N # OK
#print 'sigma_square : ',sigma_sq
if sys.argv[1] == '1':
    flag = cs1
    print_d = do_d1[0]
    dat = dat1
elif sys.argv[1] == '2':
    flag = cs2
    print_d = do_d1[1]
    dat = dat2
elif sys.argv[1] == '3':
    flag = cs3
    print_d = do_d1[2]
    dat = dat3
elif sys.argv[1] == '4':  # only one bin [0.01,0]
    flag = cs4
    print_d = do_d1[1]
    dat = dat4
elif sys.argv[1] == '5':  # only one bin [0,0.01]
    flag = cs5
    print_d = do_d1[2]
    dat = dat5
elif sys.argv[1] == '6':  # except one bin [0.01,0]
    flag = cs6
    print_d = do_d1[1]
    dat = dat6
elif sys.argv[1] == '7':  # except one bin [0,0.01]
    flag = cs7
    print_d = do_d1[2]
    dat = dat7

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
elif flag == cs4:
    print '\n********\nUsing d0 , d1 :0.01, 0\n****1 bin del_z****\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[1][0], do_d1[1][1])
    ide = str(do_d1[1][0])+'_'+str(do_d1[1][1])+'_1_bin_'+str(apply+1)
elif flag == cs5:
    print '\n********\nUsing d0 , d1 :0., 0.01\n****1 bin del_z****\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[2][0], do_d1[2][1])
    ide = str(do_d1[2][0])+'_'+str(do_d1[2][1])+'_1_bin_'+str(apply+1)
elif flag == cs6:
    print '\n********\nUsing d0 , d1 :0.01, 0\n****11 bins del_z****\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[1][0], do_d1[1][1])
    ide = str(do_d1[1][0])+'_'+str(do_d1[1][1])+'_11_bin'+'_'+str(apply+1)
elif flag == cs7:
    print '\n********\nUsing d0 , d1 :0., 0.01\n****11 bins del_z****\n'  # Change L 40 also
    comment = 'd0 = %s, d1 = %s' % (do_d1[2][0], do_d1[2][1])
    ide = str(do_d1[2][0])+'_'+str(do_d1[2][1])+'_11_bin'+'_'+str(apply+1)

#print '\ncov : \n',cov
#print '\ndel_m / sigma_sq\n',del_m_sigma_sq
#print 'array of 4 vectors\n %s'%arr
del_p2 = cov.dot(arr)
columns = ['M', 'omegam', 'w0', 'wa']
delt = pd.DataFrame(columns=columns)
delt = delt.append({'M': del_p2[0], 'omegam':  del_p2[1],
                    'w0':  del_p2[2], 'wa':  del_p2[3]}, ignore_index=True)
print '\ndel_p2 (using d0,d1 = %s ): %s\n' % (print_d, del_p2)
data2 = pd.concat([delt, dat], axis=1)
name = 'del_p_LSST_CMBPrior_Stretch_%s.txt' % ide
name2 = 'del_p_LSST_CMBPrior_Stretch_LONG_%s.txt' % ide
np.savetxt(name, np.matrix(del_p2), fmt='%.4f')
data2.to_csv(name2, index=False)
