# Generate Figure 1 of the manuscript, showing the structure of a diffracton.
import matplotlib.pyplot as plt
import numpy as np

def load_data():
    K = np.loadtxt('K.txt')
    rho = np.loadtxt('rho.txt')
    eps = []
    sig = []
    u = []
    v = []
    emax = []
    smax = []
    umax = []
    vmax = []
    for i in range(1,8):
        print i
        eps.append(np.loadtxt('sw'+str(i)+'_q0.txt'))
        u.append(np.loadtxt('sw'+str(i)+'_q1.txt'))
        v.append(np.loadtxt('sw'+str(i)+'_q2.txt'))
        sig.append(np.exp(K*eps[-1])-1)

        emax.append(np.max(eps[-1],0))
        smax.append(np.max(sig[-1],0))
        umax.append(np.max(np.abs(u[-1]),0))
        vmax.append(np.max(np.abs(v[-1]),0))
    return sig,eps,u,v
        
sig,eps,u,v = load_data()
x = np.linspace(0,300,9600)
y = np.linspace(0,1,128)
X,Y = np.meshgrid(x,y,indexing='ij')
i_diffracton = 3
x_skip = 8; y_skip = 4
plt.clf()
plt.pcolormesh(x,y,eps[i_diffracton].T,cmap='Blues')
plt.quiver(X[::x_skip,::y_skip],Y[::x_skip,::y_skip],u[i_diffracton][::x_skip,::y_skip],v[i_diffracton][::x_skip,::y_skip],pivot='center',units='width',headaxislength=5,scale=7)
plt.xlim((153,159.5))
