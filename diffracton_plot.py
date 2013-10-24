"""Generates a pcolor+quiver plot showing the structure of a single
   diffracton.  This is the first figure in the paper."""

import matplotlib.pyplot as plt
import numpy as np
path = './scaling/'

def load_data(i):
    K = np.loadtxt(path+'K.txt')
    rho = np.loadtxt(path+'rho.txt')
    eps  = np.loadtxt(path+'sw'+str(i)+'_q0.txt')
    xmom = np.loadtxt(path+'sw'+str(i)+'_q1.txt')
    ymom = np.loadtxt(path+'sw'+str(i)+'_q2.txt')

    sigma = np.exp(K*eps[-1])-1  #Stress
    u = xmom/rho
    v = ymom/rho

    return sigma,eps,u,v
        
i_diffracton = 4 # Select which diffracton to plot
sigma,eps,u,v = load_data(i_diffracton)
x = np.linspace(0,300,9600)
y = np.linspace(0,1,128)
X,Y = np.meshgrid(x,y,indexing='ij')
x_skip = 8
y_skip = 4
plt.clf()
plt.pcolormesh(x,y,eps.T,cmap='Blues')
plt.quiver(X[::x_skip,::y_skip],Y[::x_skip,::y_skip],u[::x_skip,::y_skip],v[::x_skip,::y_skip],pivot='center',units='width',headaxislength=5,scale=7)
plt.xlim((153,159.5))
