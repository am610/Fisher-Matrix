from matplotlib.patches import Ellipse
import matplotlib as mpl
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import Fisher_Ellipse_LSST_CMBPrior as F
alpha = 1.52
# angle = F.tan_two_theta() # In radian
params = F.par(F.F[2], F.F[3], F.F[4])
# xy   = centers = w_o,w_a = -1,0
#width = sigma_w_oo
#height= sigma_w_aa
# angle = sigma_w_ao/sqrt(sigma_w_oo*sigma_w_aa)  = 0.5716
#print '5 vectors', F.F
print 'Width : ', params[1]
print 'Height :', params[0]
print 'Angle :', -math.degrees(params[2])

ell = mpl.patches.Ellipse(xy=[-1, 0], width=params[1], height=params[0], angle=-
                          math.degrees(params[2]), facecolor='none', edgecolor='b', label='1-sigma')
#ell = mpl.patches.Ellipse(xy=[-1,0], width=2*1.52*np.sqrt(F.F[2]), height=2*1.52*np.sqrt(F.F[3]), angle = math.degrees(F.par(F.F[2],F.F[3],F.F[4])[2]),facecolor='none',edgecolor='b',label='1-sigma')
#ell = mpl.patches.Ellipse(xy=[-1,0], width=4.28, height=209.2915, angle = 0.5716*180/np.pi,facecolor='darkgray',edgecolor='b',label='1-sigma')


fig, ax = plt.subplots()
ax.add_patch(ell)
ax.set_aspect('auto')
ax.set_facecolor('none')
ax.autoscale()  # ''tight
ax.scatter(F.F[1], F.F[0], c='black', marker='o', label='w0-wa')
# ax.axvline(x=-1,alpha=0.5,c='black',ls='--')
# ax.axhline(y=0,alpha=0.5,c='black',ls='--')
ax.set_xlabel('w_0')
ax.set_ylabel('w_a')
# ax.set_xlim(-1.7,-0.7)
ax.set_xlim(-1.35, -0.75)
ax.set_ylim(-0.9, 0.9)
# ax.set_ylim(-3,3.0)
#ax.axvline(x=-0.871,ls = '--',c='green');ax.axvline(x=-1.129,c='green',ls='--',label='del w_0 ~ %.2f'%(2*alpha*np.sqrt(F.F[2])));
#ax.axhline(y=0.9,ls = '--',c='cyan');ax.axhline(y=-0.9,ls='--',c='cyan',label='del w_a ~ %.2f'%(2*alpha*np.sqrt(F.F[3])))
ax.tick_params(labeltop=True, labelright=True)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xticks(np.arange(-1.35, -0.75, 0.03))
plt.yticks(np.arange(-0.9, 0.9, 0.1))


# x_slope =[4.35,-6.3,-1]   #[-9.6,7.7,-1]# [4.35,-6.3,-1]
# y_slope = [-12.3,12.3,0]#[7.7,-7.3,0]#[-12.3,12.3,0]
# plt.plot(x_slope,y_slope)


r = pd.read_csv('del_p2_Long.dat')  # READING in the d_w0,d_wa data
del_w0_1 = r.iloc[12:24]['del_w0']  # d0,d1 = (0,0.01) , the x axis
del_wa_1 = r.iloc[12:24]['del_wa']  # d0,d1 = (0,0.01) , the y axis

del_w0_1_11 = r.iloc[3]['del_w0']
del_wa_1_11 = r.iloc[3]['del_wa']

del_w0_2 = r.iloc[36:]['del_w0']  # d0,d1 = (0.01,0)
del_wa_2 = r.iloc[36:]['del_wa']  # d0,d1 = (0.01,0)

del_w0_2_11 = r.iloc[27]['del_w0']
del_wa_2_11 = r.iloc[27]['del_wa']

for i in range(12):
    a1 = plt.arrow(-1, 0, del_w0_1.iloc[i], del_wa_1.iloc[i],
                   color='red', label='d0,d1=[0,0.01]')
    # a2=plt.arrow(-1,0,del_w0_2.iloc[i],del_wa_2.iloc[i],color='green',label='d0,d1=[0.01,0]')
    #plt.legend([a1,a2,], ['d0,d1=[0,0.01]','d0,d1=[0.0,0.0]',])
a_1_11 = plt.arrow(-1, 0, del_w0_1_11, del_wa_1_11, color='magenta')
# a_2_11=plt.arrow(-1,0,del_w0_2_11,del_wa_2_11,color='greenyellow')

plt.scatter(-1, 0, c='red', marker='.', label='d0,d1=[0.0,0.01]')
# plt.scatter(-1,0,c='green',marker='.',label='d0,d1=[0.01,0.0]')

plt.legend()
plt.savefig('Ellipse_LSST_CMBPrior_1.png')
plt.show()
