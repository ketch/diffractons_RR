from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import matplotlib.cm as cm
import numpy as np
import os

def plot_p(frame,skip=1,file_prefix='claw_p',path='./_output/_p',slices_limits=None):
    sol_clawpack_finer=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix='clawpack_finer')
    sol_sharpclaw_finer=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix='sharpclaw_finer')
    sol_sharpclaw_coarser=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix='sharpclaw_coarser')

    x_finer=sol_sharpclaw_finer.state.grid.x.centers; y_finer=sol_sharpclaw_finer.state.grid.y.centers; my_finer=len(y_finer)
    x_coarser=sol_sharpclaw_coarser.state.grid.x.centers; y_coarser=sol_sharpclaw_coarser.state.grid.y.centers; my_coarser=len(y_coarser)
    
    if frame < 10:
        str_frame = "00"+str(frame)
    elif frame < 100:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    p_clawpack_finer=sol_clawpack_finer.state.q[0,:,:]
    p_sharpclaw_finer=sol_sharpclaw_finer.state.q[0,:,:]
    p_sharpclaw_coarser=sol_sharpclaw_coarser.state.q[0,:,:]

    pl.figure(figsize=(8,3))
    p1,=pl.plot(x_finer,p_clawpack_finer[:,my_finer/4.],'-b',linewidth=2) # clawpack finer
    p2,=pl.plot(x_finer,p_sharpclaw_finer[:,my_finer/4.],'--r',linewidth=2) # sharpclaw finer
    p3,=pl.plot(x_coarser,p_sharpclaw_coarser[:,my_coarser/4.],'-.g',linewidth=2) # coarser finer
    pl.legend([p1,p2,p3],['Clawpack with Nx=32, Ny=128','SharpClaw with Nx=32, Ny=128','SharpClaw with Nx=16, Ny=64'],loc='upper left')
    
    pl.title("t= "+str(sol_clawpack_finer.state.t),fontsize=20)
    #pl.xlabel('x',fontsize=20)
    #pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    if slices_limits is not None:
        pl.axis([slices_limits[0],slices_limits[1],slices_limits[2],slices_limits[3]])
    pl.savefig('./_plots_to_paper/counter-interaction_015_comparison.eps')
    pl.close()
   
if __name__== "__main__":
    if not os.path.exists('./_plots_to_paper'): os.mkdir('_plots_to_paper')
    
    print('**********************')
    print('**********************')
    print('Plotting p ...')
    plot_p(frame=15,slices_limits=[22.5,42.5,-0.01,0.1])
        
