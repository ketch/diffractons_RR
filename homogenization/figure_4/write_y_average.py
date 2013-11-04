from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

def write_y_average(frame,file_prefix,path,name):
    sol=Solution(frame,file_format='petsc',path=path,read_aux=False,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; my=len(y)
    q=sol.state.q
    
    f=open(name+'.txt','w')
    f.writelines(str(xc)+" "+str(sum(q[0,i,:])/my)+"\n" for i,xc in enumerate(x))
    f.close()

if __name__== "__main__":

    write_y_average(120,'claw_p','./_output/_p/','FV_stress')
