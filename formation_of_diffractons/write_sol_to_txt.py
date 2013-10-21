from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

# material parameters
E1=5./8;   p1=8./5
E2=5./2;   p2=2./5
# interface parameters
alphax=0.5; deltax=1000.0
alphay=0.5; deltay=1.0

def write_sol(frame):
    sol=Solution(frame,file_format='petsc',read_aux=False)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; my=len(y)
    q=sol.state.q
    np.savetxt('surface_plots/q0_frame'+str(frame)+'.txt',q[0,:,:])
    np.savetxt('surface_plots/q1_frame'+str(frame)+'.txt',q[1,:,:])
    np.savetxt('surface_plots/q2_frame'+str(frame)+'.txt',q[2,:,:])

def write_aux(frame):
    sol=Solution(frame,file_format='petsc',read_aux=False)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers;
    aux=setaux(x,y)
    np.savetxt('surface_plots/rho.txt',aux[0,:,:])
    np.savetxt('surface_plots/K.txt',aux[1,:,:])

def setaux(x,y):
    r"""Creates a matrix representing every grid cell in the domain, 
    whose size is len(x),len(y)
    Each entry of the matrix contains a vector of size 3 with:
    The material density p
    The young modulus E
    A flag indicating which material the grid is made of.
    The domain pattern is a checkerboard."""

    aux = np.empty((2,len(x),len(y)), order='F')
    # xfrac and yfrac are x and y relative to deltax and deltay resp.
    xfrac=x-np.floor(x/deltax)*deltax
    yfrac=y-np.floor(y/deltay)*deltay
    # create a meshgrid out of xfrac and yfrac
    [yyfrac,xxfrac]=np.meshgrid(yfrac,xfrac)
    # density 
    aux[0,:,:]=p1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)\
        +p1*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay)\
        +p2*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay)\
        +p2*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay)
    #Young modulus
    aux[1,:,:]=E1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)\
        +E1*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay)\
        +E2*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay)\
        +E2*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay)
    
    return aux

if __name__== "__main__":
    
    write_sol(400)
    write_aux(400)

