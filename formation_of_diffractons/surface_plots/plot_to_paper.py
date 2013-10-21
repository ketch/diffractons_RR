from clawpack.petclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D as pl3
from matplotlib import rc
import matplotlib.cm as cm
import numpy as np
import os

def plot_p(frame,zlimits=None,xshift=0.):
    if not os.path.exists('./_plots_paper'): os.mkdir('./_plots_paper')
    sol=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers; x=x+xshift
    mx=len(x); my=len(y)    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    else:
        str_frame = "0"+str(frame)

    p=sol.state.q[0,:,:]

    y_index_max = p.max(0).argmax()
    x_index_max = p[:,y_index_max].argmax()
    xm = x[x_index_max]

    # plot slices
    pl.figure(figsize=(8,3))
    #pl.gcf().subplots_adjust(left=0.15)
    pl.plot(x,p[:,3*my/4.],'-r',linewidth=2)
    pl.plot(x,p[:,my/4.],'-b',linewidth=2)
    pl.plot(-x,p[:,3*my/4.],'-r',linewidth=2)
    pl.plot(-x,p[:,my/4.],'-b',linewidth=2)
    pl.title("t= "+str(sol.state.t),fontsize=20)
    pl.xlabel('x',fontsize=20)
    #pl.ylabel('Stress',fontsize=20)
    pl.xticks(size=20); pl.yticks(size=20)
    pl.xlim([xm-60, xm+15])
    if zlimits is not None:
        pl.ylim(zlimits[0]-0.05*zlimits[1],zlimits[1]*1.05)
    pl.savefig('./_plots_paper/sw_formation_frame'+str_frame+'_slices.eps')
    pl.close()

def get_extremum_p(frame):
    sol_p=Solution(frame,file_format='petsc',read_aux=False,path='./_output/_p/',file_prefix='claw_p')
    sig=sol_p.state.q[0,:]
    return [np.min(sig), np.max(sig)]

if __name__== "__main__":
    if not os.path.exists('./_plots_paper'): os.mkdir('./_plots_paper')

    frames = [0, 10, 20, 100, 200, 400]
    xshift = [0, 0, 0, 0, 200, 400]
    zlim_p = get_extremum_p(0)

    print('**********************')
    print('Plotting ...')
    index = 0
    for i in frames:
        plot_p(i,zlimits=zlim_p,xshift=xshift[index])
        index+=1
        print ('frame '+str(i)+' plotted')

