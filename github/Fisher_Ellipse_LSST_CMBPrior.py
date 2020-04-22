import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# The Marginalized [2X2] subset of COVARIANCE Matrix.
# Here Marginalized over M, omega_m
# Obtained from : fisher_ayan_function_LSST_CMBPrior.py
F = [0, -1, 2.8902e-01, 1.6886e-02, -6.8702e-02]
# F = [0, -1,  3.5020e-01,  7.3036e-03, -3.1243e-02]  # Used
# F = [wa,wo,wawa,wowo,wowa]  wa, w0, sigma_aa, sigma_00, sigma_a0
alpha = 1.52    # arxiv : 0906.4123. 1-s CL
print 'The Height and Width roughly is %.2f  %.2f' % (
    2*alpha*np.sqrt(F[3]), 2*alpha*np.sqrt(F[2]))


def a_square(sigma_x, sigma_y, sigma_xy):
    F = sigma_x + sigma_y
    G = (sigma_x - sigma_y)**2
    H = (G/4) + sigma_xy**2
    A = F/2 + np.sqrt(H)
    #print 'A^2 :',A
    return A


def b_square(sigma_x, sigma_y, sigma_xy):
    F = sigma_x + sigma_y
    G = (sigma_x - sigma_y)**2
    H = (G/4) + sigma_xy**2
    B = F/2 - np.sqrt(H)
    #print 'F, G, H',F,G,H
    #print np.sqrt(H)
    #print 'B^2 :',B
    return B


def tan_two_theta(sigma_x, sigma_y, sigma_xy):
    print 'theta :', -math.degrees(0.5*math.atan(2*sigma_xy/(sigma_x-sigma_y)))
    return(0.5*math.atan(2*sigma_xy/(sigma_x-sigma_y)))


def area(sigma_x, sigma_y, sigma_xy):
    return np.pi*(np.sqrt(a_square(sigma_x, sigma_y, sigma_xy))*alpha)*(np.sqrt(b_square(sigma_x, sigma_y, sigma_xy))*alpha)


def par(sigma_x, sigma_y, sigma_xy):
    a_star = np.sqrt(np.abs(a_square(sigma_x, sigma_y, sigma_xy)))*2*1.52
    b_star = np.sqrt(np.abs(b_square(sigma_x, sigma_y, sigma_xy)))*2*1.52
    two_theta = tan_two_theta(sigma_x, sigma_y, sigma_xy)
    param = [a_star, b_star, two_theta]
    #param = [2*alpha*np.sqrt(F[3]),2*alpha*np.sqrt(F[2]),two_theta]

    #param = [a_star/alpha,b_star/alpha,two_theta]
    #print 'a_star, b_star, 1/2*atan(2sigma_xy/sigma_x^2-sigma_y^2) : \n',param[0],param[1],param[2]
    return param

#tan_two_theta(7.3036e-03,3.5020e-01, -3.1243e-02)
#par(7.3036e-03,3.5020e-01, -3.1243e-02)
